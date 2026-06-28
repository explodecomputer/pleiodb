"""
Unit tests for the chunked storage layer, quantization, and ALID utilities.
No GWAS-VCF files required.
"""

import tempfile
import textwrap
from pathlib import Path

import numpy as np
import pytest

from pleiodb.alid import compress_allele, is_compressed, canonical_alid, parse_alid
from pleiodb.build import TraitInfo, load_trait_list, derive_neff
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
        for ti in range((T + CT - 1) // CT):
            for vi in range((V + CV - 1) // CV):
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


# ---------------------------------------------------------------------------
# ALID parser (issue #13 — needed by db.py to read variants.tsv)
# ---------------------------------------------------------------------------

class TestParseAlid:
    """parse_alid returns raw (chrom, pos, a1, a2) without canonicalisation."""

    def test_standard_snp(self):
        chrom, pos, a1, a2 = parse_alid("1:100_A_T")
        assert chrom == "1"
        assert pos == 100
        assert a1 == "A"
        assert a2 == "T"

    def test_chr_prefix_preserved(self):
        chrom, pos, a1, a2 = parse_alid("chr22:1234567_C_G")
        assert chrom == "chr22"
        assert pos == 1234567
        assert a1 == "C"
        assert a2 == "G"

    def test_compressed_allele_round_trips(self):
        """A canonical ALID that contains a compressed allele must parse back."""
        long_allele = "A" * 25
        compressed = compress_allele(long_allele)   # e.g. "AAAAAAAA~xxxx"
        alid = f"3:500_G_{compressed}"
        chrom, pos, a1, a2 = parse_alid(alid)
        assert chrom == "3"
        assert pos == 500
        assert a1 == "G"
        assert a2 == compressed

    def test_indel(self):
        chrom, pos, a1, a2 = parse_alid("10:1206798_C_CAAT")
        assert a1 == "C" and a2 == "CAAT"

    def test_invalid_format_raises(self):
        with pytest.raises(ValueError):
            parse_alid("missing_colon")

    def test_roundtrip_with_canonical_alid(self):
        """canonical_alid → parse_alid must recover the same alleles."""
        alid_str, _ = canonical_alid("7", 12345, "G", "ACGT")
        _, _, a1, a2 = parse_alid(alid_str)
        # canonical_alid compresses if needed — alleles should match compressed forms
        assert a1 == compress_allele("ACGT") or a2 == compress_allele("ACGT")


# ---------------------------------------------------------------------------
# var_y estimation (issue #10)
# ---------------------------------------------------------------------------

from pleiodb.build import compute_neff_study, estimate_var_y


class TestComputeNeffStudy:
    def test_continuous_trait(self):
        assert compute_neff_study(100_000, None) == 100_000

    def test_continuous_large(self):
        assert compute_neff_study(461_460, None) == 461_460

    def test_binary_trait(self):
        # 4 * N * K * (1 - K)
        N, K = 60_801, 0.34
        expected = 4 * N * K * (1 - K)
        assert abs(compute_neff_study(N, K) - expected) < 1e-6

    def test_binary_k_half_gives_n(self):
        """K = 0.5 → Neff = 4 * N * 0.25 = N."""
        assert abs(compute_neff_study(10_000, 0.5) - 10_000) < 1e-6


class TestEstimateVarY:
    def _make_arrays(self, n=500, var_y_true=1.0, neff_study=50_000, seed=42):
        """Generate synthetic SE and EAF consistent with a given var_y."""
        rng = np.random.default_rng(seed)
        eaf = rng.uniform(0.05, 0.95, n).astype(np.float32)
        # SE = sqrt(var_y / (2 * eaf * (1-eaf) * neff_study))
        se = np.sqrt(var_y_true / (2 * eaf * (1 - eaf) * neff_study)).astype(np.float32)
        return se, eaf

    def test_returns_tuple(self):
        se, eaf = self._make_arrays()
        result = estimate_var_y(se, eaf, 50_000)
        assert isinstance(result, tuple) and len(result) == 2

    def test_continuous_trait_recovers_var_y(self):
        """For normalised beta (var_y = 1), estimate ≈ 1.0."""
        se, eaf = self._make_arrays(n=1000, var_y_true=1.0, neff_study=50_000)
        var_y_est, n_used = estimate_var_y(se, eaf, 50_000)
        assert n_used > 900
        assert abs(var_y_est - 1.0) < 0.05, f"var_y = {var_y_est:.4f}, expected ≈ 1.0"

    def test_binary_trait_recovers_var_y(self):
        """For binary trait, var_y ≈ K*(1-K); Neff_study = 4*N*K*(1-K) handles it."""
        K = 0.34
        var_y_true = K * (1 - K)           # ≈ 0.2244
        N = 60_000
        neff_study = 4 * N * K * (1 - K)   # = 4*N*var_y_true
        se, eaf = self._make_arrays(n=1000, var_y_true=var_y_true, neff_study=neff_study)
        var_y_est, _ = estimate_var_y(se, eaf, neff_study)
        assert abs(var_y_est - var_y_true) < 0.02, (
            f"var_y = {var_y_est:.4f}, expected ≈ {var_y_true:.4f}"
        )

    def test_maf_filter_low_eaf(self):
        """Variants with EAF < 0.01 must be excluded."""
        rng = np.random.default_rng(0)
        eaf = np.array([0.005, 0.5, 0.5], dtype=np.float32)   # first is too rare
        se = np.full(3, 0.01, dtype=np.float32)
        _, n_used = estimate_var_y(se, eaf, 10_000)
        assert n_used == 2

    def test_maf_filter_high_eaf(self):
        """Variants with EAF > 0.99 must be excluded."""
        eaf = np.array([0.5, 0.995, 0.5], dtype=np.float32)
        se = np.full(3, 0.01, dtype=np.float32)
        _, n_used = estimate_var_y(se, eaf, 10_000)
        assert n_used == 2

    def test_nan_se_excluded(self):
        """Variants with NaN SE must not contribute."""
        eaf = np.array([0.5, 0.5, 0.5], dtype=np.float32)
        se = np.array([np.nan, 0.01, 0.01], dtype=np.float32)
        _, n_used = estimate_var_y(se, eaf, 10_000)
        assert n_used == 2

    def test_zero_qualifying_variants_raises(self):
        """All variants filtered → ValueError."""
        eaf = np.array([0.001, 0.999], dtype=np.float32)   # both outside MAF window
        se = np.full(2, 0.01, dtype=np.float32)
        with pytest.raises(ValueError, match="[Vv]ariant"):
            estimate_var_y(se, eaf, 10_000)

    def test_n_variants_var_y_count(self):
        """n_used must equal the number of contributing variants."""
        eaf = np.array([0.1, 0.2, 0.3, 0.4, 0.005], dtype=np.float32)
        se = np.full(5, 0.01, dtype=np.float32)
        _, n_used = estimate_var_y(se, eaf, 10_000)
        assert n_used == 4   # last one (EAF=0.005) excluded


# ---------------------------------------------------------------------------
# derive_neff (issue #11) — Neff from SE, not VCF SS
# ---------------------------------------------------------------------------

class TestDeriveNeff:
    """Unit tests for derive_neff: Neff[v] = var_y / (SE[v]² × 2·EAF[v]·(1−EAF[v]))."""

    def test_formula_correct(self):
        """Spot-check the formula against hand-computed values."""
        var_y = 1.0
        se = np.array([0.01, 0.02], dtype=np.float32)
        eaf = np.array([0.5, 0.3], dtype=np.float32)
        result = derive_neff(var_y, se, eaf)
        expected_0 = 1.0 / (0.01**2 * 2 * 0.5 * 0.5)   # = 10 000
        expected_1 = 1.0 / (0.02**2 * 2 * 0.3 * 0.7)   # ≈ 5952
        assert abs(float(result[0]) - expected_0) / expected_0 < 1e-4
        assert abs(float(result[1]) - expected_1) / expected_1 < 1e-4

    def test_nan_se_propagates_nan(self):
        """NaN SE → NaN Neff; other entries are unaffected."""
        se = np.array([np.nan, 0.01], dtype=np.float32)
        eaf = np.array([0.5, 0.5], dtype=np.float32)
        result = derive_neff(1.0, se, eaf)
        assert np.isnan(result[0])
        assert np.isfinite(result[1])

    def test_zero_eaf_produces_nan(self):
        """EAF = 0 → denominator = 0 → NaN (no divide-by-zero exception)."""
        se = np.array([0.01], dtype=np.float32)
        eaf = np.array([0.0], dtype=np.float32)
        result = derive_neff(1.0, se, eaf)
        assert np.isnan(result[0])

    def test_one_eaf_produces_nan(self):
        """EAF = 1 → denominator = 0 → NaN."""
        se = np.array([0.01], dtype=np.float32)
        eaf = np.array([1.0], dtype=np.float32)
        result = derive_neff(1.0, se, eaf)
        assert np.isnan(result[0])

    def test_roundtrip_with_known_neff(self):
        """SE derived from Neff → derive_neff must recover that Neff."""
        var_y = 1.5
        neff_true = 50_000.0
        eaf = np.array([0.2, 0.4, 0.6], dtype=np.float32)
        se = np.sqrt(var_y / (2.0 * eaf * (1.0 - eaf) * neff_true)).astype(np.float32)
        result = derive_neff(var_y, se, eaf)
        for v in result:
            assert abs(float(v) - neff_true) / neff_true < 1e-3

    def test_output_dtype_float32(self):
        """Return array must be float32."""
        se = np.array([0.01], dtype=np.float32)
        eaf = np.array([0.5], dtype=np.float32)
        result = derive_neff(1.0, se, eaf)
        assert result.dtype == np.float32

    def test_all_finite_values_positive(self):
        """All finite Neff values must be positive."""
        se = np.full(50, 0.01, dtype=np.float32)
        eaf = np.linspace(0.05, 0.95, 50).astype(np.float32)
        result = derive_neff(1.0, se, eaf)
        assert np.all(result[np.isfinite(result)] > 0)
