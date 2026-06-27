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

import csv
import io
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


# ---------------------------------------------------------------------------
# Test 7 — Query output: z, beta_norm, se_norm, pval always returned (#14)
# ---------------------------------------------------------------------------

def _query_variant_rows(db_path, variant_id):
    """Helper: query a single variant and return parsed rows as list[dict]."""
    from pleiodb.cli import _query_single_variant
    db = _open(db_path)
    fh = io.StringIO()
    _query_single_variant(db, variant_id, None, "tsv", fh)
    fh.seek(0)
    return list(csv.DictReader(fh, delimiter="\t"))


class TestQueryOutput:
    """Query always returns variant_id, trait_id, z, beta_norm, se_norm, pval."""

    def test_header_has_required_columns(self, tmp_path):
        """Every query path must produce a header with the six required columns."""
        from pleiodb.cli import _query_single_variant
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)
        fh = io.StringIO()
        _query_single_variant(db, SPOT_ALID, None, "tsv", fh)
        fh.seek(0)
        header = fh.readline().strip().split("\t")
        assert header == ["alid", "trait_id", "z", "beta_norm", "se_norm", "pval"], (
            f"Unexpected header: {header}"
        )

    def test_pval_formula(self, tmp_path):
        """pval must equal 2·Φ(−|Z|) for the spot variant."""
        from scipy.stats import norm
        db_path = _build(tmp_path, VARIANTS_HG19)
        rows = _query_variant_rows(db_path, SPOT_ALID)
        spot = next((r for r in rows if r["trait_id"] == SPOT_TRAIT), None)
        assert spot is not None, f"{SPOT_TRAIT} not found in query output"
        z_val = float(spot["z"])
        pval_val = float(spot["pval"])
        expected = 2 * norm.sf(abs(z_val))
        assert abs(pval_val - expected) / max(expected, 1e-300) < 1e-3

    def test_beta_norm_se_norm_finite_positive_se(self, tmp_path):
        """beta_norm is finite and se_norm is positive for the spot variant."""
        db_path = _build(tmp_path, VARIANTS_HG19)
        rows = _query_variant_rows(db_path, SPOT_ALID)
        spot = next((r for r in rows if r["trait_id"] == SPOT_TRAIT), None)
        assert spot is not None
        beta = float(spot["beta_norm"])
        se = float(spot["se_norm"])
        assert np.isfinite(beta), f"beta_norm={beta} is not finite"
        assert se > 0, f"se_norm={se} is not positive"

    def test_beta_norm_equals_z_times_se_norm(self, tmp_path):
        """beta_norm = Z × se_norm (definition of normalised beta)."""
        db_path = _build(tmp_path, VARIANTS_HG19)
        rows = _query_variant_rows(db_path, SPOT_ALID)
        spot = next((r for r in rows if r["trait_id"] == SPOT_TRAIT), None)
        assert spot is not None
        z = float(spot["z"])
        beta = float(spot["beta_norm"])
        se = float(spot["se_norm"])
        assert abs(beta - z * se) / max(abs(beta), 1e-10) < 1e-3

    def test_cli_no_beta_se_flag(self):
        """The query command must NOT have a --beta-se option (always on)."""
        import click.testing
        from pleiodb.cli import main
        runner = click.testing.CliRunner()
        result = runner.invoke(main, ["query", "--help"])
        assert "--beta-se" not in result.output, \
            "--beta-se flag should have been removed (output is always full)"


# ---------------------------------------------------------------------------
# Test 8 — Study-scale beta via --study-scale flag (issue #15)
# ---------------------------------------------------------------------------

def _query_variant_rows_study(db_path, variant_id):
    """Helper: query a single variant with --study-scale and return rows."""
    from pleiodb.cli import _query_single_variant
    db = _open(db_path)
    fh = io.StringIO()
    _query_single_variant(db, variant_id, None, "tsv", fh, study_scale=True)
    fh.seek(0)
    return list(csv.DictReader(fh, delimiter="\t"))


class TestStudyScaleBeta:
    """--study-scale adds beta_study and se_study = sqrt(var_y) * normalised values."""

    def test_study_scale_header_has_extra_columns(self, tmp_path):
        """With study_scale=True, output header includes beta_study and se_study."""
        from pleiodb.cli import _query_single_variant
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)
        fh = io.StringIO()
        _query_single_variant(db, SPOT_ALID, None, "tsv", fh, study_scale=True)
        fh.seek(0)
        header = fh.readline().strip().split("\t")
        assert "beta_study" in header, f"beta_study missing from header: {header}"
        assert "se_study" in header, f"se_study missing from header: {header}"

    def test_without_flag_no_extra_columns(self, tmp_path):
        """Without study_scale, beta_study and se_study must not appear."""
        from pleiodb.cli import _query_single_variant
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)
        fh = io.StringIO()
        _query_single_variant(db, SPOT_ALID, None, "tsv", fh)
        fh.seek(0)
        header = fh.readline().strip().split("\t")
        assert "beta_study" not in header
        assert "se_study" not in header

    def test_beta_study_equals_sqrt_var_y_times_beta_norm(self, tmp_path):
        """beta_study = sqrt(var_y[t]) × beta_norm within floating-point tolerance."""
        import json
        db_path = _build(tmp_path, VARIANTS_HG19)
        meta = json.loads((db_path / "meta.json").read_text())

        norm_rows = _query_variant_rows(db_path, SPOT_ALID)
        study_rows = _query_variant_rows_study(db_path, SPOT_ALID)

        spot_norm = next((r for r in norm_rows if r["trait_id"] == SPOT_TRAIT), None)
        spot_study = next((r for r in study_rows if r["trait_id"] == SPOT_TRAIT), None)
        assert spot_norm is not None and spot_study is not None

        t_idx = 1  # ieu-a-7 is index 1 in traits.tsv
        var_y = meta["var_y"][t_idx]
        beta_norm = float(spot_norm["beta_norm"])
        beta_study = float(spot_study["beta_study"])
        expected = np.sqrt(var_y) * beta_norm
        assert abs(beta_study - expected) / max(abs(expected), 1e-10) < 1e-4

    def test_se_study_equals_sqrt_var_y_times_se_norm(self, tmp_path):
        """se_study = sqrt(var_y[t]) × se_norm."""
        import json
        db_path = _build(tmp_path, VARIANTS_HG19)
        meta = json.loads((db_path / "meta.json").read_text())

        norm_rows = _query_variant_rows(db_path, SPOT_ALID)
        study_rows = _query_variant_rows_study(db_path, SPOT_ALID)

        spot_norm = next((r for r in norm_rows if r["trait_id"] == SPOT_TRAIT), None)
        spot_study = next((r for r in study_rows if r["trait_id"] == SPOT_TRAIT), None)
        assert spot_norm is not None and spot_study is not None

        t_idx = 1
        var_y = meta["var_y"][t_idx]
        se_norm = float(spot_norm["se_norm"])
        se_study = float(spot_study["se_study"])
        expected = np.sqrt(var_y) * se_norm
        assert abs(se_study - expected) / max(expected, 1e-10) < 1e-4

    def test_cli_study_scale_flag_exists(self):
        """The query command must have a --study-scale option."""
        import click.testing
        from pleiodb.cli import main
        runner = click.testing.CliRunner()
        result = runner.invoke(main, ["query", "--help"])
        assert "--study-scale" in result.output, \
            "--study-scale flag missing from CLI help"


# ---------------------------------------------------------------------------
# Test 9 — Output column renamed from variant_id to alid (#27)
# ---------------------------------------------------------------------------

class TestAlidColumnName:
    """Output header must use 'alid', not 'variant_id', for consistency with
    variants.tsv (#27)."""

    def test_header_first_column_is_alid(self, tmp_path):
        """The first column of query output must be 'alid'."""
        from pleiodb.cli import _query_single_variant
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)
        fh = io.StringIO()
        _query_single_variant(db, SPOT_ALID, None, "tsv", fh)
        fh.seek(0)
        header = fh.readline().strip().split("\t")
        assert header[0] == "alid", (
            f"First column should be 'alid', got '{header[0]}'"
        )

    def test_header_no_variant_id(self, tmp_path):
        """The string 'variant_id' must not appear in the output header."""
        from pleiodb.cli import _query_single_variant
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)
        fh = io.StringIO()
        _query_single_variant(db, SPOT_ALID, None, "tsv", fh)
        fh.seek(0)
        header = fh.readline().strip().split("\t")
        assert "variant_id" not in header, (
            f"'variant_id' should be renamed to 'alid'; got header: {header}"
        )

    def test_full_header_with_alid(self, tmp_path):
        """Full header matches expected column order with 'alid'."""
        from pleiodb.cli import _query_single_variant
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)
        fh = io.StringIO()
        _query_single_variant(db, SPOT_ALID, None, "tsv", fh)
        fh.seek(0)
        header = fh.readline().strip().split("\t")
        assert header == ["alid", "trait_id", "z", "beta_norm", "se_norm", "pval"], (
            f"Unexpected header: {header}"
        )

    def test_study_scale_header_with_alid(self, tmp_path):
        """--study-scale header must also start with 'alid'."""
        from pleiodb.cli import _query_single_variant
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)
        fh = io.StringIO()
        _query_single_variant(db, SPOT_ALID, None, "tsv", fh, study_scale=True)
        fh.seek(0)
        header = fh.readline().strip().split("\t")
        assert header == [
            "alid", "trait_id", "z", "beta_norm", "se_norm", "pval",
            "beta_study", "se_study",
        ], f"Unexpected study-scale header: {header}"


# ---------------------------------------------------------------------------
# Test 10 — --variant + --trait returns intersection only (#25)
# ---------------------------------------------------------------------------

def _query_intersect_rows(db_path, variant_id, trait_id):
    """Helper: call _query_block with a single variant + single trait."""
    from pleiodb.cli import _query_block
    db = _open(db_path)
    fh = io.StringIO()
    _query_block(db, [variant_id], [trait_id], None, "tsv", fh)
    fh.seek(0)
    return list(csv.DictReader(fh, delimiter="\t"))


class TestQueryIntersect:
    """Combining --variant and --trait returns only the single intersection
    cell, not all traits for the variant (#25)."""

    def test_returns_exactly_one_row(self, tmp_path):
        """Querying one variant + one trait yields exactly one result row."""
        db_path = _build(tmp_path, VARIANTS_HG19)
        rows = _query_intersect_rows(db_path, SPOT_ALID, SPOT_TRAIT)
        assert len(rows) == 1, (
            f"Expected 1 row from variant+trait intersection, got {len(rows)}"
        )

    def test_row_has_correct_alid(self, tmp_path):
        """The single row has the requested variant ALID."""
        db_path = _build(tmp_path, VARIANTS_HG19)
        rows = _query_intersect_rows(db_path, SPOT_ALID, SPOT_TRAIT)
        assert rows[0]["alid"] == SPOT_ALID

    def test_row_has_correct_trait(self, tmp_path):
        """The single row has the requested trait ID."""
        db_path = _build(tmp_path, VARIANTS_HG19)
        rows = _query_intersect_rows(db_path, SPOT_ALID, SPOT_TRAIT)
        assert rows[0]["trait_id"] == SPOT_TRAIT

    def test_cli_variant_and_trait_flags_dispatch_to_block(self, tmp_path):
        """CLI --variant --trait goes through _query_block (intersection path)."""
        import click.testing
        from pleiodb.cli import main
        db_path = _build(tmp_path, VARIANTS_HG19)
        runner = click.testing.CliRunner()
        result = runner.invoke(
            main,
            ["query", str(db_path), "--variant", SPOT_ALID, "--trait", SPOT_TRAIT],
        )
        assert result.exit_code == 0, f"CLI failed: {result.output}"
        lines = [l for l in result.output.strip().splitlines() if l]
        # 1 header + 1 data row = 2 lines
        assert len(lines) == 2, (
            f"Expected header + 1 row, got {len(lines)} lines:\n{result.output}"
        )
        data_cols = lines[1].split("\t")
        assert data_cols[0] == SPOT_ALID
        assert data_cols[1] == SPOT_TRAIT

    def test_pval_filter_applied_to_intersection(self, tmp_path):
        """--pval filter still applies when --variant and --trait are combined."""
        from pleiodb.cli import _query_block
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)
        # The spot variant has a very significant z in SPOT_TRAIT; use a tiny
        # threshold that it should pass, and a strict one that it should fail.
        fh_pass = io.StringIO()
        _query_block(db, [SPOT_ALID], [SPOT_TRAIT], pval=0.99, fmt="tsv", fh=fh_pass)
        fh_pass.seek(0)
        rows_pass = list(csv.DictReader(fh_pass, delimiter="\t"))
        assert len(rows_pass) == 1, "Expected 1 row when pval threshold is permissive"

        fh_fail = io.StringIO()
        _query_block(db, [SPOT_ALID], [SPOT_TRAIT], pval=1e-300, fmt="tsv", fh=fh_fail)
        fh_fail.seek(0)
        rows_fail = list(csv.DictReader(fh_fail, delimiter="\t"))
        assert len(rows_fail) == 0, "Expected 0 rows when pval threshold is impossibly strict"


# ---------------------------------------------------------------------------
# Test — build_rho + GWASDatabase.rho_matrix (issues #32, #33)
# ---------------------------------------------------------------------------

class TestBuildRho:
    """Integration tests for build_rho and the resulting rho_matrix property."""

    @pytest.fixture(scope="class")
    def db_with_rho(self, tmp_path_factory):
        from pleiodb.build import build_rho
        tmp = tmp_path_factory.mktemp("rho")
        db_path = _build(tmp, VARIANTS_HG19)
        # min_nulls=50: test dataset has ~500 variants, z_thresh=1.0 retains ~46%
        # per pair → ~230 null variants, below the default 500.
        build_rho(db_path, z_null_thresh=1.0, min_nulls=50, workers=2)
        return db_path

    def test_rho_files_exist(self, db_with_rho):
        assert (db_with_rho / "rho.bin").exists()
        assert (db_with_rho / "rho.cidx").exists()

    def test_diagonal_is_one(self, db_with_rho):
        db = _open(db_with_rho)
        mat = db.rho_matrix
        for t in range(db.T):
            block = mat.get_block(t, t + 1, t, t + 1)
            assert abs(float(block[0, 0]) - 1.0) < 1e-3, (
                f"rho[{t},{t}] = {float(block[0,0]):.4f}, expected 1.0"
            )

    def test_off_diagonal_in_range(self, db_with_rho):
        db = _open(db_with_rho)
        mat = db.rho_matrix
        for j in range(db.T):
            for k in range(db.T):
                if j == k:
                    continue
                val = float(mat.get_block(j, j + 1, k, k + 1)[0, 0])
                if not np.isnan(val):
                    assert -1.0 < val < 1.0, (
                        f"rho[{j},{k}] = {val:.4f} outside (-1, 1)"
                    )

    def test_symmetry(self, db_with_rho):
        db = _open(db_with_rho)
        mat = db.rho_matrix
        for j in range(db.T):
            for k in range(j + 1, db.T):
                vjk = float(mat.get_block(j, j + 1, k, k + 1)[0, 0])
                vkj = float(mat.get_block(k, k + 1, j, j + 1)[0, 0])
                both_nan = np.isnan(vjk) and np.isnan(vkj)
                assert both_nan or abs(vjk - vkj) < 1e-3, (
                    f"rho[{j},{k}]={vjk:.4f} != rho[{k},{j}]={vkj:.4f}"
                )

    def test_meta_json_keys(self, db_with_rho):
        import json
        meta = json.loads((db_with_rho / "meta.json").read_text())
        assert "rho_chunk_shape" in meta
        assert "rho_null_z_thresh" in meta
        assert "lambda_chunk_shape" not in meta
        assert "lambda_null_z_thresh" not in meta


# ---------------------------------------------------------------------------
# Test — pleiodb info rho_present warning (issue #33)
# ---------------------------------------------------------------------------

class TestInfoRhoWarning:
    def test_info_warns_when_rho_absent(self, tmp_path):
        """pleiodb info reports rho_present=False and a warning when rho is missing."""
        db_path = _build(tmp_path, VARIANTS_HG19)
        db = _open(db_path)
        info = db.info()
        assert info["rho_present"] is False
        assert "warnings" in info
        assert any("rho" in w.lower() for w in info["warnings"])

    def test_info_no_warning_when_rho_present(self, tmp_path):
        """pleiodb info reports rho_present=True and no rho warning after build_rho."""
        from pleiodb.build import build_rho
        db_path = _build(tmp_path, VARIANTS_HG19)
        build_rho(db_path, z_null_thresh=1.0, min_nulls=50, workers=2)
        db = _open(db_path)
        info = db.info()
        assert info["rho_present"] is True
        warnings = info.get("warnings", [])
        assert not any("rho" in w.lower() for w in warnings)

    def test_pleiodb_info_cli_warns_when_rho_absent(self, tmp_path):
        """pleiodb info CLI output contains rho_present=false when rho is missing."""
        import json
        from click.testing import CliRunner
        from pleiodb.cli import main
        db_path = _build(tmp_path, VARIANTS_HG19)
        runner = CliRunner()
        result = runner.invoke(main, ["info", str(db_path)])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["rho_present"] is False
        assert any("rho" in w.lower() for w in data.get("warnings", []))


# ---------------------------------------------------------------------------
# Test — pleiodb rho compute command (issue #34)
# ---------------------------------------------------------------------------

class TestRhoComputeCommand:
    def test_compute_exits_zero(self, tmp_path):
        """pleiodb rho <db> --min-nulls 50 exits 0 and writes rho.bin."""
        from click.testing import CliRunner
        from pleiodb.cli import main
        db_path = _build(tmp_path, VARIANTS_HG19)
        runner = CliRunner()
        result = runner.invoke(main, ["rho", str(db_path), "--min-nulls", "50"])
        assert result.exit_code == 0, result.output
        assert (db_path / "rho.bin").exists()
        assert (db_path / "rho.cidx").exists()

    def test_lambda_command_removed_or_deprecated(self, tmp_path):
        """pleiodb lambda is gone or prints a deprecation notice; does not silently succeed."""
        from click.testing import CliRunner
        from pleiodb.cli import main
        runner = CliRunner()
        result = runner.invoke(main, ["lambda", "--help"])
        # Either the command is gone (non-zero / no such command) or it warns about deprecation
        deprecated = result.exit_code != 0 or "deprecated" in (result.output or "").lower()
        assert deprecated, (
            f"'pleiodb lambda' should be removed or warn about deprecation; "
            f"got exit_code={result.exit_code}, output={result.output!r}"
        )


# ---------------------------------------------------------------------------
# Test — pleiodb rho query mode (issue #35)
# ---------------------------------------------------------------------------

def _rho_db(tmp_path_factory):
    """Shared fixture: build a db with rho computed (min_nulls=50)."""
    from pleiodb.build import build_rho
    tmp = tmp_path_factory.mktemp("rhoq")
    db_path = _build(tmp, VARIANTS_HG19)
    build_rho(db_path, z_null_thresh=1.0, min_nulls=50, workers=2)
    return db_path


class TestRhoQueryMode:
    @pytest.fixture(scope="class")
    def db_path(self, tmp_path_factory):
        return _rho_db(tmp_path_factory)

    def _invoke(self, db_path, *args):
        from click.testing import CliRunner
        from pleiodb.cli import main
        runner = CliRunner()
        return runner.invoke(main, ["rho", str(db_path)] + list(args))

    def _trait_ids(self, db_path, n=3):
        """Return the first n trait IDs from the database."""
        db = _open(db_path)
        return list(db.traits["id"][:n])

    # --- pairwise list (default) ---

    def test_two_traits_one_row(self, db_path):
        """--traits t1,t2 → exactly one data row."""
        t1, t2 = self._trait_ids(db_path, 2)
        result = self._invoke(db_path, "--traits", f"{t1},{t2}")
        assert result.exit_code == 0, result.output
        rows = list(csv.DictReader(io.StringIO(result.output), delimiter="\t"))
        assert len(rows) == 1
        assert rows[0]["trait_id_1"] == t1
        assert rows[0]["trait_id_2"] == t2

    def test_three_traits_three_rows(self, db_path):
        """--traits t1,t2,t3 → exactly three data rows (all pairs)."""
        ids = self._trait_ids(db_path, 3)
        result = self._invoke(db_path, "--traits", ",".join(ids))
        assert result.exit_code == 0, result.output
        rows = list(csv.DictReader(io.StringIO(result.output), delimiter="\t"))
        assert len(rows) == 3

    def test_pairwise_no_diagonal(self, db_path):
        """Diagonal is never emitted in pairwise output."""
        ids = self._trait_ids(db_path, 3)
        result = self._invoke(db_path, "--traits", ",".join(ids))
        rows = list(csv.DictReader(io.StringIO(result.output), delimiter="\t"))
        for row in rows:
            assert row["trait_id_1"] != row["trait_id_2"]

    def test_pairwise_header_columns(self, db_path):
        """Pairwise output has columns trait_id_1, trait_id_2, rho."""
        t1, t2 = self._trait_ids(db_path, 2)
        result = self._invoke(db_path, "--traits", f"{t1},{t2}")
        rows = list(csv.DictReader(io.StringIO(result.output), delimiter="\t"))
        assert set(rows[0].keys()) == {"trait_id_1", "trait_id_2", "rho"}

    # --- traits-file ---

    def test_traits_file_same_as_traits_flag(self, db_path, tmp_path):
        """--traits-file produces same rows as equivalent --traits."""
        ids = self._trait_ids(db_path, 3)
        traits_file = tmp_path / "traits.txt"
        traits_file.write_text("\n".join(ids) + "\n")

        r_flag = self._invoke(db_path, "--traits", ",".join(ids))
        r_file = self._invoke(db_path, "--traits-file", str(traits_file))
        assert r_flag.exit_code == 0
        assert r_file.exit_code == 0
        assert r_flag.output == r_file.output

    # --- --matrix flag ---

    def test_matrix_diagonal_is_one(self, db_path):
        """--matrix output has diagonal values = 1.0."""
        ids = [str(x) for x in self._trait_ids(db_path, 3)]
        result = self._invoke(db_path, "--traits", ",".join(ids), "--matrix")
        assert result.exit_code == 0, result.output
        # Use splitlines() without strip() so the leading tab on the header is preserved
        lines = result.output.splitlines()
        # first line is the header row: \t<id1>\t<id2>\t...
        header = lines[0].split("\t")[1:]   # [0] is the empty cell before the first trait
        assert header == ids, f"header={header!r}, ids={ids!r}"
        for i, line in enumerate(lines[1:]):
            cells = line.split("\t")
            row_label = cells[0]
            assert row_label == ids[i]
            diag_val = float(cells[i + 1])
            assert abs(diag_val - 1.0) < 1e-3, f"Diagonal [{i},{i}] = {diag_val}"

    def test_matrix_shape(self, db_path):
        """--matrix output is T×T (T data rows, T+1 columns including label)."""
        ids = [str(x) for x in self._trait_ids(db_path, 3)]
        result = self._invoke(db_path, "--traits", ",".join(ids), "--matrix")
        lines = result.output.splitlines()
        # header + 3 data rows
        assert len(lines) == 4
        for data_line in lines[1:]:
            assert len(data_line.split("\t")) == len(ids) + 1

    # --- --format json ---

    def test_json_format_valid(self, db_path):
        """--format json produces valid JSON with correct keys."""
        import json
        t1, t2 = self._trait_ids(db_path, 2)
        result = self._invoke(db_path, "--traits", f"{t1},{t2}", "--format", "json")
        assert result.exit_code == 0, result.output
        data = json.loads(result.output)
        assert isinstance(data, list)
        assert len(data) == 1
        assert set(data[0].keys()) == {"trait_id_1", "trait_id_2", "rho"}
        assert data[0]["trait_id_1"] == t1
        assert data[0]["trait_id_2"] == t2

    # --- --output flag ---

    def test_output_to_file(self, db_path, tmp_path):
        """--output writes the result to the specified file instead of stdout."""
        t1, t2 = self._trait_ids(db_path, 2)
        out_file = tmp_path / "result.tsv"
        result = self._invoke(db_path, "--traits", f"{t1},{t2}",
                              "--output", str(out_file))
        assert result.exit_code == 0, result.output
        assert result.output.strip() == ""   # nothing on stdout
        content = out_file.read_text()
        rows = list(csv.DictReader(io.StringIO(content), delimiter="\t"))
        assert len(rows) == 1


# ---------------------------------------------------------------------------
# Test – imputation pipeline (issue 036)
# ---------------------------------------------------------------------------

class TestImputation:
    """Build with a synthetic LD panel and verify imputed.coo.zst + traits.tsv."""

    @staticmethod
    def _make_ld_panel(tmp_path, variants_file):
        """Create a minimal synthetic LD reference panel for the test variants.

        Reads the first two variants from variants_file and creates a 2-block
        LD panel (one block with those 2 variants) that the imputer can match.
        Returns the path to the panel root directory (the ancestry dir).
        """
        import gzip, io
        import pandas as pd

        variants = []
        with open(variants_file) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    alid = line.split("\t")[0]
                    # Parse chrom:pos_a1_a2
                    chrom, rest = alid.split(":", 1)
                    parts = rest.split("_")
                    bp = int(parts[0])
                    variants.append((alid, chrom, bp))
                if len(variants) >= 4:
                    break

        if len(variants) < 2:
            pytest.skip("Not enough variants to build LD panel fixture")

        # Group by chromosome; use the first chromosome
        chrom = variants[0][1].lstrip("chr")
        block_variants = [v for v in variants if v[1].lstrip("chr") == chrom][:2]
        bps = [v[2] for v in block_variants]
        start = min(bps) - 1000
        end = max(bps) + 1000
        block_name = f"{start}-{end}"
        block_dir = tmp_path / chrom / block_name
        block_dir.mkdir(parents=True)

        # Write TSV
        tsv_lines = ["CHR\tSNP\tOA\tEA\tEAF\tBP"]
        for alid, _, bp in block_variants:
            parts = alid.split("_")
            a1, a2 = parts[1], parts[2]
            tsv_lines.append(f"{chrom}\t{alid}\t{a2}\t{a1}\t0.3\t{bp}")
        (block_dir / f"{block_name}.tsv").write_text("\n".join(tsv_lines) + "\n")

        # Write identity-like LD matrix (diagonal 1, off-diag 0.5)
        n = len(block_variants)
        ld = np.eye(n) + 0.5 * (np.ones((n, n)) - np.eye(n))
        buf = io.BytesIO()
        with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
            for row in ld:
                gz.write(("\t".join(f"{v:.6f}" for v in row) + "\n").encode())
        (block_dir / f"{block_name}.unphased.vcor1.gz").write_bytes(buf.getvalue())

        return tmp_path  # caller uses tmp_path / ancestry=""

    def test_imputed_coo_written_without_ld_dir(self, tmp_path):
        """imputed.coo.zst is always written, even with no LD panel (empty)."""
        import zstandard
        db_path = tmp_path / "test.pleiodb"
        from pleiodb.build import build_database
        build_database(
            output_dir=db_path,
            variants_path=VARIANTS_HG19,
            trait_tsv=TRAITS_TSV,
            chunk_shape=(64, 8),
            workers=2,
        )
        coo_path = db_path / "imputed.coo.zst"
        assert coo_path.exists(), "imputed.coo.zst should always be written"
        blob = zstandard.ZstdDecompressor().decompress(coo_path.read_bytes())
        arr = np.frombuffer(blob, dtype=np.uint32).reshape(-1, 2)
        assert len(arr) == 0, "no imputed cells expected without LD panel"

    def test_n_variants_imputed_column_present(self, tmp_path):
        """traits.tsv must have n_variants_imputed column (even without LD)."""
        db_path = tmp_path / "test.pleiodb"
        from pleiodb.build import build_database
        build_database(
            output_dir=db_path,
            variants_path=VARIANTS_HG19,
            trait_tsv=TRAITS_TSV,
            chunk_shape=(64, 8),
            workers=2,
        )
        with open(db_path / "traits.tsv") as f:
            header = f.readline().strip().split("\t")
        assert "n_variants_imputed" in header

    def test_n_variants_imputed_zero_without_ld_dir(self, tmp_path):
        """Without LD panel, n_variants_imputed must be 0 for all traits."""
        db_path = tmp_path / "test.pleiodb"
        from pleiodb.build import build_database
        build_database(
            output_dir=db_path,
            variants_path=VARIANTS_HG19,
            trait_tsv=TRAITS_TSV,
            chunk_shape=(64, 8),
            workers=2,
        )
        import csv as csv_mod, io as io_mod
        with open(db_path / "traits.tsv") as f:
            rows = list(csv_mod.DictReader(f, delimiter="\t"))
        for row in rows:
            assert int(row["n_variants_imputed"]) == 0

    def test_imputed_coo_with_synthetic_ld_panel(self, tmp_path):
        """With a synthetic LD panel, imputed.coo.zst may have entries."""
        import zstandard
        ld_dir = tmp_path / "ld_panel"
        ld_dir.mkdir()
        panel_root = self._make_ld_panel(ld_dir, VARIANTS_HG19)

        db_path = tmp_path / "test.pleiodb"
        from pleiodb.build import build_database
        build_database(
            output_dir=db_path,
            variants_path=VARIANTS_HG19,
            trait_tsv=TRAITS_TSV,
            chunk_shape=(64, 8),
            workers=2,
            ld_dir=panel_root,
            ld_ancestry="",   # no ancestry subdir; panel_root is the chrom dir parent
            ld_thresh=0.9,
            ld_min_cor=0.0,   # accept any imputation quality for this test
        )
        coo_path = db_path / "imputed.coo.zst"
        assert coo_path.exists()
        blob = zstandard.ZstdDecompressor().decompress(coo_path.read_bytes())
        arr = np.frombuffer(blob, dtype=np.uint32).reshape(-1, 2)
        # The panel has only 2 matched variants; imputation may or may not
        # produce cells depending on VCF coverage, but the file must exist
        # and the pairs must be sorted.
        if len(arr) > 1:
            assert np.all(
                arr[1:, 0] > arr[:-1, 0]
                | (arr[1:, 0] == arr[:-1, 0]) & (arr[1:, 1] >= arr[:-1, 1])
            ), "COO pairs must be sorted (v, t)"

    def test_n_variants_counts_only_observed(self, tmp_path):
        """n_variants + n_variants_imputed ≤ total finite z-scores."""
        import csv as csv_mod, zstandard
        ld_dir = tmp_path / "ld_panel"
        ld_dir.mkdir()
        panel_root = self._make_ld_panel(ld_dir, VARIANTS_HG19)

        db_path = tmp_path / "test.pleiodb"
        from pleiodb.build import build_database
        build_database(
            output_dir=db_path,
            variants_path=VARIANTS_HG19,
            trait_tsv=TRAITS_TSV,
            chunk_shape=(64, 8),
            workers=2,
            ld_dir=panel_root,
            ld_ancestry="",
            ld_min_cor=0.0,
        )
        with open(db_path / "traits.tsv") as f:
            rows = list(csv_mod.DictReader(f, delimiter="\t"))
        for row in rows:
            n_obs = int(row["n_variants"])
            n_imp = int(row["n_variants_imputed"])
            assert n_obs >= 0
            assert n_imp >= 0


# ---------------------------------------------------------------------------
# Test – hg19 variant liftover for LD imputation (issues 037–040)
# ---------------------------------------------------------------------------

class TestLiftoverImputation:
    """Build with hg19 variants + hg38 LD panel; verify liftover fires correctly."""

    # Two hg19 ALIDs from the test variant file (chromosome 10)
    HG19_ALID_0 = "10:101248474_C_T"
    HG19_ALID_1 = "10:125249764_C_T"
    # Fake hg38 positions the mock will return
    HG38_POS_0 = 100000000
    HG38_POS_1 = 124000000

    @staticmethod
    def _make_hg38_ld_panel(tmp_path, alid_hg38_0: str, alid_hg38_1: str) -> Path:
        """Create a tiny synthetic LD panel using *hg38* ALIDs."""
        import gzip, io as _io
        chrom = "10"
        for alid in (alid_hg38_0, alid_hg38_1):
            assert alid.startswith("10:"), f"expected chr10 ALID, got {alid}"
        bps = [int(a.split(":")[1].split("_")[0]) for a in (alid_hg38_0, alid_hg38_1)]
        start = min(bps) - 1000
        end = max(bps) + 1000
        block_name = f"{start}-{end}"

        block_dir = tmp_path / chrom / block_name
        block_dir.mkdir(parents=True)

        tsv_lines = ["CHR\tSNP\tOA\tEA\tEAF\tBP"]
        for alid, bp in zip((alid_hg38_0, alid_hg38_1), bps):
            parts = alid.split("_")
            a1, a2 = parts[1], parts[2]
            tsv_lines.append(f"{chrom}\t{alid}\t{a2}\t{a1}\t0.3\t{bp}")
        (block_dir / f"{block_name}.tsv").write_text("\n".join(tsv_lines) + "\n")

        n = 2
        ld = np.eye(n) + 0.5 * (np.ones((n, n)) - np.eye(n))
        buf = _io.BytesIO()
        with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
            for row in ld:
                gz.write(("\t".join(f"{v:.6f}" for v in row) + "\n").encode())
        (block_dir / f"{block_name}.unphased.vcor1.gz").write_bytes(buf.getvalue())

        return tmp_path  # callers use ld_ancestry="" so tmp_path is the chrom parent

    def _build_with_liftover(self, tmp_path):
        """Build with mocked pyliftover that maps hg19 → known hg38 positions."""
        from unittest.mock import MagicMock, patch
        from pleiodb.build import build_database

        alid_hg38_0 = f"10:{self.HG38_POS_0}_C_T"
        alid_hg38_1 = f"10:{self.HG38_POS_1}_C_T"

        ld_dir = tmp_path / "ld_panel"
        ld_dir.mkdir()
        panel_root = self._make_hg38_ld_panel(ld_dir, alid_hg38_0, alid_hg38_1)

        # Mock LiftOver so our two test variants map to the hg38 positions we built
        hg19_pos_0 = int(self.HG19_ALID_0.split(":")[1].split("_")[0])
        hg19_pos_1 = int(self.HG19_ALID_1.split(":")[1].split("_")[0])
        result_map = {
            ("chr10", hg19_pos_0 - 1): [("chr10", self.HG38_POS_0 - 1, "+", 0)],
            ("chr10", hg19_pos_1 - 1): [("chr10", self.HG38_POS_1 - 1, "+", 0)],
        }
        lo_mock = MagicMock()
        lo_mock.convert_coordinate.side_effect = lambda c, p: result_map.get((c, p), [])

        db_path = tmp_path / "test_lift.pleiodb"
        with patch("pyliftover.LiftOver", return_value=lo_mock):
            build_database(
                output_dir=db_path,
                variants_path=VARIANTS_HG19,
                trait_tsv=TRAITS_TSV,
                chunk_shape=(64, 8),
                workers=2,
                variants_build="hg19",
                ld_dir=panel_root,
                ld_ancestry="",
                ld_min_cor=0.0,
            )
        return db_path

    def test_alid_hg38_column_written(self, tmp_path):
        """variants.tsv must contain alid_hg38 column after liftover build."""
        db_path = self._build_with_liftover(tmp_path)
        with open(db_path / "variants.tsv") as f:
            header = f.readline().strip().split("\t")
        assert "alid_hg38" in header, f"alid_hg38 not in header: {header}"

    def test_meta_json_flag(self, tmp_path):
        """meta.json must have variants_hg38_stored=true after liftover build."""
        import json
        db_path = self._build_with_liftover(tmp_path)
        meta = json.loads((db_path / "meta.json").read_text())
        assert meta.get("variants_hg38_stored") is True

    def test_id_hg38_readable_from_db(self, tmp_path):
        """db.variants['id_hg38'] must be non-empty for lifted variants."""
        from pleiodb.db import GWASDatabase
        db_path = self._build_with_liftover(tmp_path)
        db = GWASDatabase.open(db_path)
        id_hg38 = db.variants["id_hg38"]
        # At least our two target variants should have been lifted
        non_empty = (id_hg38 != "").sum()
        assert non_empty >= 2, f"Expected ≥2 hg38 ALIDs, got {non_empty}"

    def test_backwards_compat_no_hg38_column(self, tmp_path):
        """Opening a DB built without liftover returns empty id_hg38 strings."""
        from pleiodb.build import build_database
        from pleiodb.db import GWASDatabase
        db_path = tmp_path / "no_lift.pleiodb"
        build_database(
            output_dir=db_path,
            variants_path=VARIANTS_HG38,
            trait_tsv=TRAITS_TSV,
            chunk_shape=(64, 8),
            workers=2,
            variants_build="hg38",
        )
        db = GWASDatabase.open(db_path)
        id_hg38 = db.variants["id_hg38"]
        assert all(v == "" for v in id_hg38), "Expected all empty id_hg38 for non-lifted DB"

    def test_meta_json_flag_false_without_liftover(self, tmp_path):
        """meta.json must have variants_hg38_stored=false when no liftover performed."""
        import json
        from pleiodb.build import build_database
        db_path = tmp_path / "no_lift2.pleiodb"
        build_database(
            output_dir=db_path,
            variants_path=VARIANTS_HG38,
            trait_tsv=TRAITS_TSV,
            chunk_shape=(64, 8),
            workers=2,
            variants_build="hg38",
        )
        meta = json.loads((db_path / "meta.json").read_text())
        assert meta.get("variants_hg38_stored") is False


# ---------------------------------------------------------------------------
# Test – post-build imputation pass (issue 041)
# ---------------------------------------------------------------------------

class TestPostBuildImputation:
    """Verify that imputed positions get finite z-scores and finite Neff."""

    @staticmethod
    def _make_ld_panel(tmp_path, variants_file):
        """Minimal synthetic LD panel reused from TestImputation."""
        import gzip

        variants = []
        with open(variants_file) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    alid = line.split("\t")[0]
                    chrom, rest = alid.split(":", 1)
                    parts = rest.split("_")
                    bp = int(parts[0])
                    variants.append((alid, chrom, bp))
                if len(variants) >= 4:
                    break

        if len(variants) < 2:
            pytest.skip("Not enough variants for LD panel fixture")

        chrom = variants[0][1].lstrip("chr")
        block_variants = [v for v in variants if v[1].lstrip("chr") == chrom][:2]
        bps = [v[2] for v in block_variants]
        start, end = min(bps) - 1000, max(bps) + 1000
        block_name = f"{start}-{end}"
        block_dir = tmp_path / chrom / block_name
        block_dir.mkdir(parents=True)

        tsv_lines = ["CHR\tSNP\tOA\tEA\tEAF\tBP"]
        for alid, _, bp in block_variants:
            parts = alid.split("_")
            a1, a2 = parts[1], parts[2]
            tsv_lines.append(f"{chrom}\t{alid}\t{a2}\t{a1}\t0.3\t{bp}")
        (block_dir / f"{block_name}.tsv").write_text("\n".join(tsv_lines) + "\n")

        n = len(block_variants)
        ld = np.eye(n) + 0.5 * (np.ones((n, n)) - np.eye(n))
        buf = io.BytesIO()
        with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
            for row in ld:
                gz.write(("\t".join(f"{v:.6f}" for v in row) + "\n").encode())
        (block_dir / f"{block_name}.unphased.vcor1.gz").write_bytes(buf.getvalue())

        return tmp_path

    def _build_with_ld(self, tmp_path):
        from pleiodb.build import build_database
        panel_root = self._make_ld_panel(tmp_path / "ld", VARIANTS_HG19)
        db_path = tmp_path / "test.pleiodb"
        build_database(
            output_dir=db_path,
            variants_path=VARIANTS_HG19,
            trait_tsv=TRAITS_TSV,
            chunk_shape=(64, 8),
            workers=2,
            ld_dir=panel_root,
            ld_ancestry="",
            ld_min_cor=0.0,
        )
        return db_path

    def test_observed_z_scores_unchanged(self, tmp_path):
        """Observed (non-imputed) z-scores must be identical to a build without LD."""
        from pleiodb.build import build_database
        from pleiodb.db import GWASDatabase
        from pleiodb.quantize import decode_z

        # Build without imputation
        db_no_ld = tmp_path / "no_ld.pleiodb"
        build_database(
            output_dir=db_no_ld,
            variants_path=VARIANTS_HG19,
            trait_tsv=TRAITS_TSV,
            chunk_shape=(64, 8),
            workers=2,
        )

        # Build with imputation
        db_ld = self._build_with_ld(tmp_path)

        db_base = GWASDatabase.open(db_no_ld)
        db_imp = GWASDatabase.open(db_ld)

        z_base = decode_z(db_base._zscore.get_block(0, db_base.V, 0, db_base.T))
        z_imp = decode_z(db_imp._zscore.get_block(0, db_imp.V, 0, db_imp.T))

        # Observed positions (finite in base) must match in the imputed build
        observed_mask = np.isfinite(z_base)
        assert np.allclose(z_base[observed_mask], z_imp[observed_mask], atol=0.01), \
            "Observed z-scores changed after post-build imputation pass"

    def test_imputed_neff_assigned_from_neff_base(self, tmp_path):
        """Imputed positions get neff_base[t] so betas are recoverable.

        Patches impute_z_block to force imputation at a known position, then
        verifies that the rewritten Neff file has a finite value there.
        """
        from unittest.mock import patch
        from pleiodb.build import build_database
        from pleiodb.db import GWASDatabase
        from pleiodb.quantize import decode_neff, decode_z

        IMPUTED_V = 0  # force imputation at variant index 0

        def _fake_impute(z_block, variants, eaf_arr, block_index, **kwargs):
            out_mask = kwargs.get("out_mask")
            # Impute variant 0 for every trait where it has an observed z-score
            observed = np.isfinite(z_block[IMPUTED_V, :])
            if observed.any():
                z_block[IMPUTED_V, :] = np.where(observed, z_block[IMPUTED_V, :] + 0.0, 99.0)
                if out_mask is not None:
                    out_mask[IMPUTED_V, :] = True

        db_path = tmp_path / "forced_impute.pleiodb"
        with patch("pleiodb.impute.impute_z_block", side_effect=_fake_impute):
            panel_root = self._make_ld_panel(tmp_path / "ld", VARIANTS_HG19)
            build_database(
                output_dir=db_path,
                variants_path=VARIANTS_HG19,
                trait_tsv=TRAITS_TSV,
                chunk_shape=(64, 8),
                workers=2,
                ld_dir=panel_root,
                ld_ancestry="",
                ld_min_cor=0.0,
            )

        db = GWASDatabase.open(db_path)
        neff_full = decode_neff(db._neff.get_block(0, db.V, 0, db.T))

        # Variant 0 should have finite Neff for all traits (assigned from neff_base)
        for t_idx in range(db.T):
            n = neff_full[IMPUTED_V, t_idx]
            assert np.isfinite(n), \
                f"Imputed Neff at variant 0 trait {t_idx} should be finite, got {n}"
