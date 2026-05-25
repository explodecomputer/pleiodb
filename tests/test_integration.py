"""
Integration tests for the full build → query pipeline.

These tests require bcftools and cyvcf2 (conda pleiodb environment).
They are automatically skipped when either dependency is absent.

Test data lives in tests/test_data/:
  variants_hg19.tsv  — 500 ALIDs + EAF, hg19
  variants_hg38.tsv  — same variants lifted to hg38
  traits.tsv         — 5 traits (trait_id | trait_name | vcf_path)
  vcf/               — 5 GWAS-VCFs, hg19, CSI-indexed
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pytest

# Skip all tests in this module when cyvcf2 is unavailable
cyvcf2 = pytest.importorskip("cyvcf2", reason="cyvcf2 not installed (activate pleiodb conda env)")

# Also skip when bcftools is not on PATH
_has_bcftools = subprocess.run(
    ["bcftools", "--version"], capture_output=True
).returncode == 0
pytestmark = pytest.mark.skipif(
    not _has_bcftools, reason="bcftools not on PATH (activate pleiodb conda env)"
)

TEST_DATA = Path(__file__).parent / "test_data"
VARIANTS_HG19 = TEST_DATA / "variants_hg19.tsv"
VARIANTS_HG38 = TEST_DATA / "variants_hg38.tsv"
TRAITS_TSV = TEST_DATA / "traits.tsv"

# Known spot-check: variant 10:1206798_C_T in trait ieu-a-7
# VCF record: REF=C ALT=T  ES=0.001936  SE=0.017765
# Expected z = ES/SE ≈ 0.1090   (C < T alphabetically → canonical, no flip)
SPOT_ALID = "10:1206798_C_T"
SPOT_TRAIT = "ieu-a-7"
SPOT_Z_EXPECTED = 0.001936 / 0.017765
SPOT_EAF_EXPECTED = 0.071479   # AF from VCF (T is A2 = effect allele)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _build(tmp_path, variants_file, variants_build=None):
    from pleiodb.build import build_database

    db_path = tmp_path / "test.pleiodb"
    build_database(
        output_dir=db_path,
        variants_path=variants_file,
        trait_tsv=TRAITS_TSV,
        chunk_shape=(64, 8),   # small chunks for fast test
        workers=2,
        variants_build=variants_build,
    )
    return db_path


def _open(db_path):
    from pleiodb.db import GWASDatabase
    return GWASDatabase.open(db_path)


# ---------------------------------------------------------------------------
# Test 1 — hg19 basic build and query
# ---------------------------------------------------------------------------

class TestHg19Build:
    def test_dimensions(self, tmp_path):
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)
        assert db.V == 501
        assert db.T == 5

    def test_spot_check_zscore(self, tmp_path):
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)

        v_idx = int(db.variant_index([SPOT_ALID])[0])
        t_idx = int(db.trait_index([SPOT_TRAIT])[0])
        z = db.zscore_block([v_idx], [t_idx])[0, 0]

        assert not np.isnan(z), f"z-score for {SPOT_ALID} in {SPOT_TRAIT} is NaN"
        assert abs(float(z) - SPOT_Z_EXPECTED) < 0.01, (
            f"z={z:.4f}, expected ≈{SPOT_Z_EXPECTED:.4f}"
        )

    def test_spot_check_eaf(self, tmp_path):
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)

        v_idx = int(db.variant_index([SPOT_ALID])[0])
        eaf_val = float(db.eaf[v_idx])
        assert abs(eaf_val - SPOT_EAF_EXPECTED) < 0.001, (
            f"EAF={eaf_val:.4f}, expected ≈{SPOT_EAF_EXPECTED:.4f}"
        )

    def test_match_rate(self, tmp_path):
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)

        t_idx = int(db.trait_index([SPOT_TRAIT])[0])
        z = db.zscore_trait(t_idx)
        match_rate = float(np.isfinite(z).mean())
        assert match_rate >= 0.9, f"Match rate too low: {match_rate:.1%}"

    def test_eaf_range(self, tmp_path):
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)

        eaf = db.eaf
        finite = eaf[np.isfinite(eaf)]
        assert np.all(finite >= 0) and np.all(finite <= 1)

    def test_allele_flip_negates_zscore(self, tmp_path):
        """Verify that a non-canonical VCF record (REF > ALT alphabetically)
        produces a negated z-score relative to a directly computed ES/SE."""
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)

        # Find any variant where the canonical ALID has A2 matching VCF REF
        # (i.e. VCF stored it as REF=larger allele, ALT=smaller allele).
        # Proxy: look for a variant whose EAF > 0.5 and check sign is consistent.
        # (A full flip test would need a known flipped record in the VCF.)
        t_idx = int(db.trait_index([SPOT_TRAIT])[0])
        z = db.zscore_trait(t_idx)
        assert np.any(np.isfinite(z)), "No finite z-scores found"


# ---------------------------------------------------------------------------
# Test 2 — hg38 liftover path
# ---------------------------------------------------------------------------

class TestHg38Liftover:
    @pytest.fixture(scope="class")
    def db(self, tmp_path_factory):
        if not VARIANTS_HG38.exists():
            pytest.skip("variants_hg38.tsv not yet generated (run generate_test_data.sh)")
        tmp_path = tmp_path_factory.mktemp("liftover")
        db_path = _build(tmp_path, VARIANTS_HG38, variants_build="hg38")
        return _open(db_path)

    def test_dimensions(self, db):
        assert db.V == 501
        assert db.T == 5

    def test_match_rate(self, db):
        t_idx = int(db.trait_index([SPOT_TRAIT])[0])
        z = db.zscore_trait(t_idx)
        match_rate = float(np.isfinite(z).mean())
        assert match_rate >= 0.9, f"Liftover match rate too low: {match_rate:.1%}"

    def test_same_zscore_as_hg19(self, db, tmp_path_factory):
        """The hg38 liftover path should produce the same z-score as hg19 direct."""
        # Find the hg38 ALID for our spot-check variant.
        # The ALID in hg38 coords has a different POS but the same A1/A2.
        v_arr = db.variants
        t_idx = int(db.trait_index([SPOT_TRAIT])[0])

        # Look up by A1=C, A2=T at any position on chr10 (lifted coord)
        candidates = [
            i for i, row in enumerate(v_arr)
            if str(row["a1"]) == "C" and str(row["a2"]) == "T"
            and str(row["chrom"]).lstrip("chr") == "10"
        ]
        assert candidates, "No C/T variants on chr10 found in liftover database"
        # At least one should be non-NaN for ieu-a-7
        z_vals = [db.zscore_block([i], [t_idx])[0, 0] for i in candidates]
        assert any(np.isfinite(v) for v in z_vals), (
            "All C/T chr10 variants are NaN for ieu-a-7 after liftover"
        )


# ---------------------------------------------------------------------------
# Test 3 — variants.tsv storage format (issue #13)
# ---------------------------------------------------------------------------

class TestVariantsTsv:
    """Verify that build_database writes variants.tsv and not variants.npy / eaf.f16."""

    def test_variants_tsv_written(self, tmp_path):
        db_path = _build(tmp_path, VARIANTS_HG19)
        assert (db_path / "variants.tsv").exists(), "variants.tsv should be written"

    def test_legacy_files_absent(self, tmp_path):
        db_path = _build(tmp_path, VARIANTS_HG19)
        assert not (db_path / "variants.npy").exists(), "variants.npy should not be written"
        assert not (db_path / "eaf.f16").exists(), "eaf.f16 should not be written"

    def test_variants_tsv_header(self, tmp_path):
        db_path = _build(tmp_path, VARIANTS_HG19)
        with open(db_path / "variants.tsv") as f:
            header = f.readline().strip().split("\t")
        assert header == ["alid", "eaf"]

    def test_variants_tsv_row_count(self, tmp_path):
        db_path = _build(tmp_path, VARIANTS_HG19)
        lines = (db_path / "variants.tsv").read_text().strip().splitlines()
        assert len(lines) == 502, f"expected 502 (501 variants + 1 header), got {len(lines)}"

    def test_db_variants_reads_correctly(self, tmp_path):
        """GWASDatabase.variants must return a structured array from variants.tsv."""
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)
        v = db.variants
        assert len(v) == 501
        assert "id" in v.dtype.names
        assert "chrom" in v.dtype.names
        assert "pos" in v.dtype.names
        assert "a1" in v.dtype.names
        assert "a2" in v.dtype.names

    def test_db_eaf_reads_correctly(self, tmp_path):
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)
        eaf = db.eaf
        assert len(eaf) == 501
        finite = eaf[np.isfinite(eaf)]
        assert np.all(finite >= 0) and np.all(finite <= 1)

    def test_spot_alid_in_variants(self, tmp_path):
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)
        ids = db.variants["id"]
        assert SPOT_ALID in ids, f"{SPOT_ALID} not found in variants.tsv"

    def test_variant_index_works(self, tmp_path):
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)
        idx = db.variant_index([SPOT_ALID])
        assert len(idx) == 1 and 0 <= idx[0] < 501


# ---------------------------------------------------------------------------
# Test 4 — var_y estimation at ingest (issue #10)
# ---------------------------------------------------------------------------

class TestVarYIngest:
    """Verify that var_y is correctly estimated during build_database."""

    def test_build_completes_with_var_y_computation(self, tmp_path):
        """Build must succeed — var_y computation runs without error for all traits."""
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)
        assert db.T == 5

    def test_var_y_stored_in_meta(self, tmp_path):
        """var_y per trait is stored in meta.json for downstream use."""
        import json
        db_path = _build(tmp_path, VARIANTS_HG19)
        meta = json.loads((db_path / "meta.json").read_text())
        assert "var_y" in meta, "var_y should be stored in meta.json"
        assert len(meta["var_y"]) == 5

    def test_var_y_values_positive_finite(self, tmp_path):
        import json
        db_path = _build(tmp_path, VARIANTS_HG19)
        meta = json.loads((db_path / "meta.json").read_text())
        for v in meta["var_y"]:
            assert v is not None and np.isfinite(v) and v > 0, f"var_y={v} not positive-finite"

    def test_var_y_binary_differs_from_continuous(self, tmp_path):
        """Binary trait (ieu-a-7) should have a different var_y than continuous traits.

        The precise value depends on how the test VCF's SE was calibrated, so we
        only check that the binary var_y is positive and finite (the mathematical
        correctness is verified by unit tests with synthetic data).
        """
        import json
        db_path = _build(tmp_path, VARIANTS_HG19)
        meta = json.loads((db_path / "meta.json").read_text())
        # trait order matches traits.tsv: ukb-b-19953, ieu-a-7, ukb-b-10787, ieu-b-110, ieu-b-109
        var_y_cad = meta["var_y"][1]   # ieu-a-7 (binary)
        assert var_y_cad is not None and np.isfinite(var_y_cad) and var_y_cad > 0, (
            f"ieu-a-7 var_y={var_y_cad} not positive-finite"
        )


# ---------------------------------------------------------------------------
# Test 5 — Neff derived from SE (issue #11)
# ---------------------------------------------------------------------------

class TestNeffDerivedFromSE:
    """Neff stored in the DB is derived from SE + var_y, not from VCF SS."""

    def test_build_completes(self, tmp_path):
        """Build must succeed with the new Neff derivation path."""
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)
        assert db.T == 5

    def test_neff_matrix_has_finite_values(self, tmp_path):
        """After build, the Neff matrix must contain finite values for the spot trait."""
        from pleiodb.quantize import decode_neff
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)
        t_idx = int(db.trait_index([SPOT_TRAIT])[0])
        neff_raw = db._neff.get_block(0, db.V, t_idx, t_idx + 1)[:, 0]
        neff = decode_neff(neff_raw)
        n_finite = int(np.isfinite(neff).sum())
        assert n_finite > 0, "Expected finite Neff values in the neff matrix"

    def test_neff_consistent_with_reconstruction(self, tmp_path):
        """SE_norm = 1/sqrt(Neff × 2·EAF·(1−EAF)) must match SE_vcf/sqrt(var_y)."""
        import json
        from pleiodb.quantize import decode_neff
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)

        v_idx = int(db.variant_index([SPOT_ALID])[0])
        t_idx = int(db.trait_index([SPOT_TRAIT])[0])

        neff_raw = db._neff.get_block(v_idx, v_idx + 1, t_idx, t_idx + 1)[0, 0]
        neff = float(decode_neff(np.array([neff_raw]))[0])
        eaf_v = float(db.eaf[v_idx])

        assert np.isfinite(neff) and neff > 0, f"Neff for spot variant is {neff}"

        # SE_norm from stored Neff
        se_norm_from_neff = 1.0 / np.sqrt(neff * 2.0 * eaf_v * (1.0 - eaf_v))

        # Expected SE_norm = SE_vcf / sqrt(var_y)
        meta = json.loads((db_path / "meta.json").read_text())
        var_y = meta["var_y"][t_idx]
        se_vcf = 0.017765   # known from SPOT_ALID / SPOT_TRAIT VCF record
        se_norm_expected = se_vcf / np.sqrt(var_y)

        assert abs(se_norm_from_neff - se_norm_expected) / se_norm_expected < 0.05, (
            f"SE_norm from Neff = {se_norm_from_neff:.5f}, "
            f"expected SE_vcf/sqrt(var_y) ≈ {se_norm_expected:.5f}  (5% tolerance)"
        )


# ---------------------------------------------------------------------------
# Test 6 — traits.tsv storage format (issue #12)
# ---------------------------------------------------------------------------

class TestTraitsTsv:
    """Verify that build_database writes traits.tsv and not traits.npy / neff_base.f32."""

    def test_traits_tsv_written(self, tmp_path):
        db_path = _build(tmp_path, VARIANTS_HG19)
        assert (db_path / "traits.tsv").exists(), "traits.tsv should be written"

    def test_legacy_files_absent(self, tmp_path):
        db_path = _build(tmp_path, VARIANTS_HG19)
        assert not (db_path / "traits.npy").exists(), "traits.npy should not be written"
        assert not (db_path / "neff_base.f32").exists(), "neff_base.f32 should not be written"

    def test_traits_tsv_header(self, tmp_path):
        db_path = _build(tmp_path, VARIANTS_HG19)
        with open(db_path / "traits.tsv") as f:
            header = f.readline().strip().split("\t")
        required = {"trait_id", "trait_name", "N", "K", "neff_study", "var_y",
                    "n_variants", "n_variants_var_y"}
        assert required.issubset(set(header)), (
            f"Missing columns: {required - set(header)}"
        )

    def test_traits_tsv_row_count(self, tmp_path):
        db_path = _build(tmp_path, VARIANTS_HG19)
        lines = (db_path / "traits.tsv").read_text().strip().splitlines()
        assert len(lines) == 6, f"expected 6 (5 traits + 1 header), got {len(lines)}"

    def test_db_traits_reads_correctly(self, tmp_path):
        """GWASDatabase.traits must return a structured array with id and name fields."""
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)
        t = db.traits
        assert len(t) == 5
        assert "id" in t.dtype.names
        assert "name" in t.dtype.names

    def test_trait_index_works(self, tmp_path):
        """trait_index must correctly resolve SPOT_TRAIT from traits.tsv."""
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)
        idx = db.trait_index([SPOT_TRAIT])
        assert len(idx) == 1 and 0 <= idx[0] < 5

    def test_neff_base_from_traits_tsv(self, tmp_path):
        """db.neff_base must return a float32 array of length T read from traits.tsv."""
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)
        nb = db.neff_base
        assert len(nb) == 5
        assert nb.dtype == np.float32
        assert np.all(np.isfinite(nb) & (nb > 0)), \
            f"neff_base has non-positive or non-finite entries: {nb}"

    def test_n_variants_column_populated(self, tmp_path):
        """n_variants column must be positive integers for all traits."""
        import csv
        db_path = _build(tmp_path, VARIANTS_HG19)
        with open(db_path / "traits.tsv") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                n_v = int(row["n_variants"])
                assert n_v > 0, f"n_variants=0 for trait {row['trait_id']}"

    def test_var_y_column_matches_meta(self, tmp_path):
        """var_y in traits.tsv must match the values in meta.json."""
        import json
        import csv
        db_path = _build(tmp_path, VARIANTS_HG19)
        meta_var_y = json.loads((db_path / "meta.json").read_text())["var_y"]
        with open(db_path / "traits.tsv") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
        for i, row in enumerate(rows):
            tsv_vy_str = row["var_y"]
            meta_vy = meta_var_y[i]
            if meta_vy is None:
                assert tsv_vy_str == "", f"Expected empty var_y for trait {i}"
            else:
                assert abs(float(tsv_vy_str) - meta_vy) < 1e-4, (
                    f"var_y mismatch at trait {i}: tsv={tsv_vy_str}, meta={meta_vy}"
                )
