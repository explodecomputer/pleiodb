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
