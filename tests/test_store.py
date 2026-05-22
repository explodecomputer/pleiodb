"""
Unit tests for the chunked storage layer, quantization, and ALID utilities.
No GWAS-VCF files required.
"""

import tempfile
import textwrap
from pathlib import Path

import numpy as np
import pytest

from pleiodb.alid import compress_allele, is_compressed, canonical_alid
from pleiodb.build import TraitInfo, load_trait_list
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


# ---------------------------------------------------------------------------
# ALID utilities (issue #8)
# ---------------------------------------------------------------------------

class TestCompressAllele:
    def test_short_allele_unchanged(self):
        """Alleles ≤ 20 chars must not be modified."""
        short = "ACGTACGTACGTACGTACGT"  # exactly 20 chars
        assert compress_allele(short) == short

    def test_single_char_unchanged(self):
        assert compress_allele("A") == "A"

    def test_long_allele_compressed(self):
        """Allele > 20 chars → {allele[:8]}~{sha256[:4]}, total 13 chars."""
        long_allele = "A" * 21
        result = compress_allele(long_allele)
        assert len(result) == 13
        assert result[:8] == "AAAAAAAA"
        assert result[8] == "~"
        assert len(result[9:]) == 4

    def test_compressed_format_hex(self):
        """The 4-char suffix must be valid lowercase hex."""
        import re
        long_allele = "GCTAGCTAGCTAGCTAGCTAGCTA"  # 24 chars
        result = compress_allele(long_allele)
        suffix = result[9:]
        assert re.match(r"^[0-9a-f]{4}$", suffix), f"Not hex: {suffix!r}"

    def test_two_distinct_long_alleles_differ(self):
        """Two different long alleles at the same position must compress differently."""
        a = "A" * 25
        b = "T" * 25
        assert compress_allele(a) != compress_allele(b)

    def test_compression_deterministic(self):
        """Same allele always produces the same compressed form."""
        allele = "GCTAGCTAGCTAGCTAGCTAGCTA"
        assert compress_allele(allele) == compress_allele(allele)

    def test_boundary_21_chars(self):
        """21-char allele is the first that gets compressed."""
        assert len(compress_allele("A" * 20)) == 20  # unchanged
        assert len(compress_allele("A" * 21)) == 13  # compressed


class TestIsCompressed:
    def test_plain_allele_not_compressed(self):
        assert not is_compressed("ACGT")

    def test_compressed_allele_detected(self):
        allele = "A" * 21
        compressed = compress_allele(allele)
        assert is_compressed(compressed)

    def test_tilde_in_plain_allele(self):
        """Tilde is the marker; it cannot appear in a real genomic allele."""
        assert is_compressed("ACGT~1234")


class TestCanonicalAlid:
    def test_canonical_order_unchanged(self):
        """When a1 ≤ a2 alphabetically, no flip."""
        alid, flipped = canonical_alid("1", 100, "A", "T")
        assert alid == "1:100_A_T"
        assert not flipped

    def test_non_canonical_order_flipped(self):
        """When a1 > a2 alphabetically, alleles are swapped."""
        alid, flipped = canonical_alid("1", 100, "T", "A")
        assert alid == "1:100_A_T"
        assert flipped

    def test_long_alleles_compressed_in_alid(self):
        """Long alleles are compressed in the resulting ALID string."""
        long = "A" * 25
        short = "G"
        alid, flipped = canonical_alid("3", 500, short, long)
        compressed = compress_allele(long)
        # short "G" > compress("AAAA...") lexicographically? Let's just verify
        # the ALID contains the compressed form and has correct structure.
        assert compressed in alid or short in alid
        # Both parts should be in the alid
        c_long = compress_allele(long)
        parts = alid.split("_", 2)  # chrom:pos, a1, a2
        assert c_long in (parts[1], parts[2])

    def test_flip_with_long_alleles(self):
        """EAF-flip signal is correct even when alleles are compressed."""
        long_a = "AAAAAAAAAAAAAAAAAAAAAAA"  # 23 chars
        long_b = "TTTTTTTTTTTTTTTTTTTTTTT"  # 23 chars  → T... > A... so flip expected
        alid1, f1 = canonical_alid("1", 1, long_a, long_b)
        alid2, f2 = canonical_alid("1", 1, long_b, long_a)
        # Swapping input order should produce the same ALID but opposite flip
        assert alid1 == alid2
        assert f1 != f2


# ---------------------------------------------------------------------------
# Traits input file loader (issue #9)
# ---------------------------------------------------------------------------

def _write_traits_tsv(tmp_path, content: str) -> Path:
    """Write a traits TSV fixture and return its Path."""
    p = tmp_path / "traits.tsv"
    p.write_text(textwrap.dedent(content))
    return p


class TestLoadTraitList:
    """Tests for load_trait_list() — the public traits-file parser."""

    # ------------------------------------------------------------------
    # Happy-path
    # ------------------------------------------------------------------

    def test_parses_header_and_all_columns(self, tmp_path):
        """A well-formed file returns correct TraitInfo objects."""
        tsv = _write_traits_tsv(tmp_path, """\
            trait_id\ttrait_name\tN\tK\tvcf_path\tbuild
            t1\tHeight\t337000\t\t/data/t1.vcf.gz\thg19
            t2\tCAD\t60801\t0.34\t/data/t2.vcf.gz\thg19
        """)
        traits = load_trait_list(tsv)
        assert len(traits) == 2

        t1 = traits[0]
        assert t1.trait_id == "t1"
        assert t1.trait_name == "Height"
        assert t1.N == 337000
        assert t1.K is None          # empty K → continuous
        assert t1.vcf_path == "/data/t1.vcf.gz"
        assert t1.vcf_build == "hg19"

        t2 = traits[1]
        assert t2.trait_id == "t2"
        assert t2.N == 60801
        assert abs(t2.K - 0.34) < 1e-9   # binary

    def test_build_column_optional(self, tmp_path):
        """build column may be absent; vcf_build should be None."""
        tsv = _write_traits_tsv(tmp_path, """\
            trait_id\ttrait_name\tN\tK\tvcf_path
            t1\tBMI\t500000\t\t/data/t1.vcf.gz
        """)
        traits = load_trait_list(tsv)
        assert traits[0].vcf_build is None

    def test_k_absent_column_all_continuous(self, tmp_path):
        """No K column at all → every trait is continuous (K=None)."""
        tsv = _write_traits_tsv(tmp_path, """\
            trait_id\ttrait_name\tN\tvcf_path
            t1\tBMI\t500000\t/data/t1.vcf.gz
        """)
        traits = load_trait_list(tsv)
        assert traits[0].K is None

    def test_k_valid_boundary_values(self, tmp_path):
        """K strictly inside (0, 1) is accepted."""
        for k in ("0.01", "0.5", "0.99"):
            tsv = _write_traits_tsv(tmp_path, f"""\
                trait_id\ttrait_name\tN\tK\tvcf_path
                t1\tTrait\t10000\t{k}\t/data/t.vcf.gz
            """)
            traits = load_trait_list(tsv)
            assert traits[0].K == pytest.approx(float(k))

    def test_comment_and_blank_lines_skipped(self, tmp_path):
        """Lines starting with # and blank lines are ignored."""
        tsv = _write_traits_tsv(tmp_path, """\
            # comment line
            trait_id\ttrait_name\tN\tK\tvcf_path

            t1\tBMI\t500000\t\t/data/t1.vcf.gz
        """)
        traits = load_trait_list(tsv)
        assert len(traits) == 1

    # ------------------------------------------------------------------
    # Error cases
    # ------------------------------------------------------------------

    def test_missing_n_column_raises(self, tmp_path):
        """No N column → ValueError at parse time."""
        tsv = _write_traits_tsv(tmp_path, """\
            trait_id\ttrait_name\tvcf_path
            t1\tBMI\t/data/t1.vcf.gz
        """)
        with pytest.raises(ValueError, match="[Nn]"):
            load_trait_list(tsv)

    def test_empty_n_value_raises(self, tmp_path):
        """Empty N cell → ValueError."""
        tsv = _write_traits_tsv(tmp_path, """\
            trait_id\ttrait_name\tN\tvcf_path
            t1\tBMI\t\t/data/t1.vcf.gz
        """)
        with pytest.raises(ValueError):
            load_trait_list(tsv)

    def test_k_equals_zero_raises(self, tmp_path):
        """K = 0 is not a valid case fraction."""
        tsv = _write_traits_tsv(tmp_path, """\
            trait_id\ttrait_name\tN\tK\tvcf_path
            t1\tCAD\t10000\t0\t/data/t.vcf.gz
        """)
        with pytest.raises(ValueError, match="[Kk]"):
            load_trait_list(tsv)

    def test_k_equals_one_raises(self, tmp_path):
        """K = 1 is not a valid case fraction."""
        tsv = _write_traits_tsv(tmp_path, """\
            trait_id\ttrait_name\tN\tK\tvcf_path
            t1\tCAD\t10000\t1\t/data/t.vcf.gz
        """)
        with pytest.raises(ValueError, match="[Kk]"):
            load_trait_list(tsv)

    def test_k_out_of_range_raises(self, tmp_path):
        """K > 1 raises ValueError."""
        tsv = _write_traits_tsv(tmp_path, """\
            trait_id\ttrait_name\tN\tK\tvcf_path
            t1\tCAD\t10000\t1.5\t/data/t.vcf.gz
        """)
        with pytest.raises(ValueError, match="[Kk]"):
            load_trait_list(tsv)
