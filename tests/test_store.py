"""
Unit tests for the chunked storage layer and quantization.
No GWAS-VCF files required.
"""

import tempfile
from pathlib import Path

import numpy as np
import pytest

from pleiodb.store import ChunkedMatrix
from pleiodb.quantize import (
    encode_z, decode_z, encode_neff, decode_neff,
    encode_eaf, decode_eaf, reconstruct_beta_se,
    Z_NA, NEFF_NA,
)


# ---------------------------------------------------------------------------
# Quantization round-trips
# ---------------------------------------------------------------------------

class TestZQuant:
    def test_roundtrip_random(self):
        z = np.random.normal(0, 3, 10000).astype(np.float32)
        rec = decode_z(encode_z(z))
        # max quantization error = 0.5 / Z_SCALE = 0.005
        np.testing.assert_allclose(rec, z, atol=0.005)

    def test_nan_sentinel(self):
        z = np.array([np.nan, 1.5, np.nan], dtype=np.float32)
        enc = encode_z(z)
        assert enc[0] == Z_NA
        assert enc[2] == Z_NA
        dec = decode_z(enc)
        assert np.isnan(dec[0])
        assert abs(dec[1] - 1.5) < 0.01

    def test_clip(self):
        z = np.array([1000.0, -1000.0], dtype=np.float32)
        enc = encode_z(z)
        assert enc[0] == np.iinfo(np.int16).max
        assert enc[1] == np.iinfo(np.int16).min + 1

    def test_sign_preserved(self):
        z = np.array([-3.14, 3.14], dtype=np.float32)
        dec = decode_z(encode_z(z))
        assert dec[0] < 0
        assert dec[1] > 0


class TestNeffQuant:
    def test_roundtrip(self):
        neff = np.array([1000.0, 50000.0, 500000.0, 1e6], dtype=np.float32)
        rec = decode_neff(encode_neff(neff))
        np.testing.assert_allclose(rec, neff, rtol=0.001)  # < 0.1% error

    def test_nan_sentinel(self):
        neff = np.array([np.nan, 10000.0], dtype=np.float32)
        enc = encode_neff(neff)
        assert enc[0] == NEFF_NA
        dec = decode_neff(enc)
        assert np.isnan(dec[0])
        assert abs(dec[1] - 10000.0) / 10000.0 < 0.001

    def test_nonpositive(self):
        neff = np.array([0.0, -5.0], dtype=np.float32)
        enc = encode_neff(neff)
        assert np.all(enc == NEFF_NA)


class TestReconstructBetaSE:
    def test_reconstruction(self):
        beta_true = np.array([0.05, -0.1], dtype=np.float32)
        se_true = np.array([0.01, 0.02], dtype=np.float32)
        eaf = np.array([0.4, 0.3], dtype=np.float32)
        # Neff = 1 / (2 * p * (1-p) * se^2)
        neff = 1.0 / (2.0 * eaf * (1.0 - eaf) * se_true**2)
        z = beta_true / se_true

        beta_r, se_r = reconstruct_beta_se(
            z.reshape(2, 1),
            neff.reshape(2, 1),
            eaf,
        )
        np.testing.assert_allclose(se_r[:, 0], se_true, rtol=0.001)
        np.testing.assert_allclose(beta_r[:, 0], beta_true, rtol=0.01)


# ---------------------------------------------------------------------------
# ChunkedMatrix: write → read round-trip
# ---------------------------------------------------------------------------

class TestChunkedMatrix:
    def _make(self, tmpdir, shape=(1000, 200), chunk=(100, 50)):
        base = Path(tmpdir) / "test"
        rng = np.random.default_rng(42)
        data = rng.integers(-100, 100, size=shape, dtype=np.int16)
        mat = ChunkedMatrix(base, shape, np.int16, chunk)
        mat.open_write()
        CV, CT = chunk
        V, T = shape
        for vi in range((V + CV - 1) // CV):
            for ti in range((T + CT - 1) // CT):
                v0, v1 = vi * CV, min((vi + 1) * CV, V)
                t0, t1 = ti * CT, min((ti + 1) * CT, T)
                mat.write_chunk(vi, ti, data[v0:v1, t0:t1])
        mat.close_write()
        return mat, data

    def test_full_block(self, tmp_path):
        mat, data = self._make(tmp_path)
        result = mat.get_block(0, data.shape[0], 0, data.shape[1])
        np.testing.assert_array_equal(result, data)

    def test_partial_block(self, tmp_path):
        mat, data = self._make(tmp_path)
        result = mat.get_block(50, 150, 30, 80)
        np.testing.assert_array_equal(result, data[50:150, 30:80])

    def test_single_row(self, tmp_path):
        mat, data = self._make(tmp_path)
        result = mat.get_block(7, 8, 0, data.shape[1])
        np.testing.assert_array_equal(result[0], data[7])

    def test_single_col(self, tmp_path):
        mat, data = self._make(tmp_path)
        result = mat.get_block(0, data.shape[0], 13, 14)
        np.testing.assert_array_equal(result[:, 0], data[:, 13])

    def test_cross_chunk_boundary(self, tmp_path):
        mat, data = self._make(tmp_path, shape=(1000, 200), chunk=(100, 50))
        # Spans two v-chunks and two t-chunks
        result = mat.get_block(80, 130, 40, 70)
        np.testing.assert_array_equal(result, data[80:130, 40:70])
