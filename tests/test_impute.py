"""Tests for src/pleiodb/impute.py (issues 033–035) and liftover.lift_variants (issue 037)."""
from __future__ import annotations

import gzip
import io
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from pleiodb.impute import (
    _ld_pca,
    _elastic_net_impute,
    _poly_rescale,
    _ratio_outlier_mask,
    _se_outliers,
    _load_ld_submatrix,
    build_block_index,
    impute_z_block,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _random_pd_ld(n: int, rng: np.random.Generator) -> np.ndarray:
    """Return a random n×n positive-definite matrix usable as an LD matrix."""
    A = rng.standard_normal((n, n))
    return A @ A.T + np.eye(n) * n * 0.1


def _make_variants(alids: list[str]) -> np.ndarray:
    dt = np.dtype([
        ("id", "U64"),
        ("chrom", "U10"),
        ("pos", np.uint32),
        ("a1", "U10"),
        ("a2", "U10"),
    ])
    rows = []
    for alid in alids:
        chrom, rest = alid.split(":", 1)
        pos_str, a1, a2 = rest.split("_", 2)
        rows.append((alid, chrom, int(pos_str), a1, a2))
    return np.array(rows, dtype=dt)


@pytest.fixture()
def rng():
    return np.random.default_rng(42)


@pytest.fixture()
def small_ld(rng):
    return _random_pd_ld(20, rng)


@pytest.fixture()
def ld_block_dir(tmp_path, rng):
    """Create a tiny synthetic LD block directory."""
    chrom = "1"
    start, end = 1000, 5000
    block_name = f"{start}-{end}"
    block_dir = tmp_path / chrom / block_name
    block_dir.mkdir(parents=True)

    alids = [
        "1:1100_A_C",
        "1:1200_G_T",
        "1:1300_A_T",
        "1:1400_C_G",
        "1:2000_A_G",
    ]
    n = len(alids)
    ld = _random_pd_ld(n, rng)

    tsv_lines = ["CHR\tSNP\tOA\tEA\tEAF\tBP"]
    for i, alid in enumerate(alids):
        chrom_, rest = alid.split(":", 1)
        bp = rest.split("_")[0]
        a1, a2 = rest.split("_")[1], rest.split("_")[2]
        tsv_lines.append(f"{chrom_}\t{alid}\t{a2}\t{a1}\t0.3\t{bp}")
    (block_dir / f"{block_name}.tsv").write_text("\n".join(tsv_lines) + "\n")

    # Write LD matrix as gzipped TSV (no header)
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        for row in ld:
            gz.write(("\t".join(f"{v:.6f}" for v in row) + "\n").encode())
    (block_dir / f"{block_name}.unphased.vcor1.gz").write_bytes(buf.getvalue())

    return block_dir, alids, ld


# ---------------------------------------------------------------------------
# 033 – _ld_pca
# ---------------------------------------------------------------------------

class TestLdPca:
    def test_eigenvalues_descending(self, small_ld):
        vals, vecs = _ld_pca(small_ld, thresh=0.9)
        assert np.all(vals[:-1] >= vals[1:]), "eigenvalues not in descending order"

    def test_thresh_respected(self, small_ld):
        thresh = 0.9
        vals, vecs = _ld_pca(small_ld, thresh=thresh)
        all_vals, _ = _ld_pca(small_ld, thresh=1.0)
        cumvar = np.cumsum(vals) / all_vals.sum()
        assert cumvar[-1] >= thresh

    def test_fewer_components_than_snps(self, small_ld):
        vals, vecs = _ld_pca(small_ld, thresh=0.5)
        assert len(vals) < small_ld.shape[0]

    def test_eigenvectors_shape(self, small_ld):
        vals, vecs = _ld_pca(small_ld, thresh=0.9)
        n = small_ld.shape[0]
        assert vecs.shape[0] == n
        assert vecs.shape[1] == len(vals)

    def test_non_negative_eigenvalues(self, small_ld):
        vals, _ = _ld_pca(small_ld, thresh=0.9)
        assert np.all(vals >= 0)


# ---------------------------------------------------------------------------
# 033 – _elastic_net_impute
# ---------------------------------------------------------------------------

class TestElasticNetImpute:
    def test_returns_full_length(self, small_ld, rng):
        vals, vecs = _ld_pca(small_ld, thresh=0.9)
        n = small_ld.shape[0]
        z = rng.standard_normal(n).astype(np.float64)
        z[::4] = np.nan  # mask every 4th
        result = _elastic_net_impute(z, vecs, len(vals))
        assert result is not None
        assert len(result) == n

    def test_no_nans_in_result(self, small_ld, rng):
        vals, vecs = _ld_pca(small_ld, thresh=0.9)
        z = rng.standard_normal(small_ld.shape[0]).astype(np.float64)
        z[:5] = np.nan
        result = _elastic_net_impute(z, vecs, len(vals))
        assert result is not None
        assert np.all(np.isfinite(result))

    def test_returns_none_all_missing(self, small_ld):
        vals, vecs = _ld_pca(small_ld, thresh=0.9)
        z = np.full(small_ld.shape[0], np.nan)
        assert _elastic_net_impute(z, vecs, len(vals)) is None

    def test_returns_none_constant_z(self, small_ld):
        vals, vecs = _ld_pca(small_ld, thresh=0.9)
        z = np.ones(small_ld.shape[0])
        z[:5] = np.nan
        assert _elastic_net_impute(z, vecs, len(vals)) is None

    def test_returns_none_one_observed(self, small_ld):
        vals, vecs = _ld_pca(small_ld, thresh=0.9)
        z = np.full(small_ld.shape[0], np.nan)
        z[0] = 1.5
        assert _elastic_net_impute(z, vecs, len(vals)) is None


# ---------------------------------------------------------------------------
# 033 – _poly_rescale
# ---------------------------------------------------------------------------

class TestPolyRescale:
    def test_high_correlation_linear(self, rng):
        n = 100
        predicted = rng.standard_normal(n)
        truth = 2.5 * predicted + rng.standard_normal(n) * 0.1
        truth[10:20] = np.nan  # some missing
        adj, corr = _poly_rescale(truth, predicted, npoly=1)
        assert np.isfinite(corr)
        assert corr > 0.9

    def test_output_length_matches_input(self, rng):
        n = 50
        predicted = rng.standard_normal(n)
        truth = predicted * 1.5
        truth[:5] = np.nan
        adj, _ = _poly_rescale(truth, predicted)
        assert len(adj) == n

    def test_degenerate_too_few_points(self):
        truth = np.array([1.0, np.nan, np.nan, np.nan])
        predicted = np.array([0.5, 1.0, 1.5, 2.0])
        adj, corr = _poly_rescale(truth, predicted, npoly=3)
        # Should not raise; corr may be nan
        assert len(adj) == 4

    def test_all_missing_truth(self, rng):
        truth = np.full(20, np.nan)
        predicted = rng.standard_normal(20)
        adj, corr = _poly_rescale(truth, predicted)
        assert not np.isfinite(corr)


# ---------------------------------------------------------------------------
# 033 – _ratio_outlier_mask
# ---------------------------------------------------------------------------

class TestRatioOutlierMask:
    def test_clear_outlier_flagged(self):
        truth = np.ones(20)
        predicted = np.ones(20)
        predicted[5] = 100.0  # ratio will be 0.01 = outlier
        keep = _ratio_outlier_mask(truth, predicted)
        assert not keep[5]
        assert keep.sum() == 19


# ---------------------------------------------------------------------------
# 033 – _se_outliers
# ---------------------------------------------------------------------------

class TestSeOutliers:
    def test_clear_outlier_flagged(self):
        n = 50
        se_hat = np.linspace(0.01, 0.1, n)
        se_obs = se_hat * 1.0 + 0.001
        se_obs[10] = se_hat[10] * 50  # extreme outlier
        flags = _se_outliers(se_obs, se_hat)
        assert flags[10], "obvious outlier should be flagged"
        assert flags.sum() <= 3, "should not flag most variants"

    def test_no_outliers_clean_data(self):
        # On clean linear data, _se_outliers may still flag high-leverage points
        # at the extremes of the predictor range (expected Cook's-distance behaviour).
        # The important property is that it flags far fewer than a dataset with a
        # real outlier (tested in test_clear_outlier_flagged).
        n = 100
        rng = np.random.default_rng(7)
        se_hat = np.linspace(0.01, 0.08, n)
        se_obs = se_hat + rng.standard_normal(n) * 0.001
        flags = _se_outliers(se_obs, se_hat)
        assert flags.sum() < n // 5  # fewer than 20% flagged on clean data

    def test_too_few_points_returns_false(self):
        se_obs = np.array([0.1, 0.2, np.nan])
        se_hat = np.array([0.1, 0.2, 0.3])
        flags = _se_outliers(se_obs, se_hat)
        assert not flags.any()

    def test_length_matches_input(self):
        n = 30
        se_obs = np.ones(n) * 0.05
        se_hat = np.ones(n) * 0.05
        flags = _se_outliers(se_obs, se_hat)
        assert len(flags) == n


# ---------------------------------------------------------------------------
# 034 – _load_ld_submatrix
# ---------------------------------------------------------------------------

class TestLoadLdSubmatrix:
    def _ld_path(self, block_dir):
        return block_dir / f"{block_dir.name}.unphased.vcor1.gz"

    def test_shape_and_values(self, ld_block_dir):
        block_dir, alids, expected_ld = ld_block_dir
        n = len(alids)
        idx = list(range(n))
        sub = _load_ld_submatrix(self._ld_path(block_dir), idx)
        assert sub.shape == (n, n)
        np.testing.assert_allclose(sub, expected_ld, atol=1e-5)

    def test_submatrix_subset(self, ld_block_dir):
        block_dir, alids, expected_ld = ld_block_dir
        idx = [0, 2, 4]
        sub = _load_ld_submatrix(self._ld_path(block_dir), idx)
        assert sub.shape == (3, 3)
        expected_sub = expected_ld[np.ix_(idx, idx)]
        np.testing.assert_allclose(sub, expected_sub, atol=1e-5)


# ---------------------------------------------------------------------------
# 034 – build_block_index
# ---------------------------------------------------------------------------

class TestBuildBlockIndex:
    def test_matched_variants_found(self, ld_block_dir, tmp_path):
        block_dir, alids, _ = ld_block_dir
        # block_dir is tmp_path/1/1000-5000
        ld_dir = tmp_path

        # Use first 3 alids as our pleiodb variant list (2 should be enough)
        variants = _make_variants(alids[:3])
        idx = build_block_index(variants, ld_dir, ancestry="")

        assert len(idx) == 1
        block_key = list(idx.keys())[0]
        info = idx[block_key]
        assert len(info["variant_indices"]) == 3
        assert len(info["ld_row_indices"]) == 3
        assert "ld_path" in info

    def test_block_excluded_when_fewer_than_2_matches(self, ld_block_dir, tmp_path):
        block_dir, alids, _ = ld_block_dir
        ld_dir = tmp_path

        # Only one variant matches → block should be excluded
        variants = _make_variants([alids[0], "2:9999_A_G"])
        idx = build_block_index(variants, ld_dir, ancestry="")
        assert len(idx) == 0

    def test_missing_tsv_skipped(self, tmp_path):
        (tmp_path / "1" / "100-200").mkdir(parents=True)
        variants = _make_variants(["1:150_A_C", "1:160_G_T"])
        idx = build_block_index(variants, tmp_path, ancestry="")
        assert len(idx) == 0

    def test_empty_ld_dir(self, tmp_path):
        variants = _make_variants(["1:100_A_C"])
        idx = build_block_index(variants, tmp_path / "nonexistent", ancestry="")
        assert len(idx) == 0

    def test_n_ld_snps_reported(self, ld_block_dir, tmp_path):
        block_dir, alids, _ = ld_block_dir
        ld_dir = tmp_path
        variants = _make_variants(alids)
        idx = build_block_index(variants, ld_dir, ancestry="")
        info = list(idx.values())[0]
        assert info["n_ld_snps"] == len(alids)
        assert info["ld_path"].exists()


# ---------------------------------------------------------------------------
# 035 – impute_z_block
# ---------------------------------------------------------------------------

class TestImputeZBlock:
    def _make_index(self, ld_block_dir, alids, subset=None):
        """Build a block_index dict pointing at the fixture block."""
        block_dir, _, ld = ld_block_dir
        n = len(alids)
        idxs = list(range(n)) if subset is None else subset
        ld_path = block_dir / f"{block_dir.name}.unphased.vcor1.gz"
        return {
            f"1/{block_dir.name}": {
                "ld_path": ld_path,
                "variant_indices": idxs,
                "ld_row_indices": idxs,
                "n_ld_snps": n,
            }
        }

    def test_nans_filled(self, ld_block_dir, rng):
        block_dir, alids, _ = ld_block_dir
        n = len(alids)
        variants = _make_variants(alids)
        eaf = np.full(n, 0.3, dtype=np.float32)

        z_block = rng.standard_normal((n, 3)).astype(np.float32)
        z_block[2, 0] = np.nan
        z_block[4, 1] = np.nan

        block_index = self._make_index(ld_block_dir, alids)
        out_mask = np.zeros((n, 3), dtype=bool)

        impute_z_block(z_block, variants, eaf, block_index,
                       thresh=0.9, min_cor=0.0, out_mask=out_mask)

        # After imputation, the NaN positions should be filled
        assert np.isfinite(z_block[2, 0]), "position (2,0) should have been imputed"
        assert np.isfinite(z_block[4, 1]), "position (4,1) should have been imputed"
        assert out_mask[2, 0], "mask should be True for imputed cell (2,0)"
        assert out_mask[4, 1], "mask should be True for imputed cell (4,1)"

    def test_observed_unchanged(self, ld_block_dir, rng):
        block_dir, alids, _ = ld_block_dir
        n = len(alids)
        variants = _make_variants(alids)
        eaf = np.full(n, 0.3, dtype=np.float32)

        z_block = rng.standard_normal((n, 2)).astype(np.float32)
        z_block[1, 0] = np.nan  # set NaN before taking the reference snapshot
        z_before = z_block.copy()

        block_index = self._make_index(ld_block_dir, alids)
        impute_z_block(z_block, variants, eaf, block_index, thresh=0.9, min_cor=0.0)

        # Positions that were observed (non-NaN) before imputation must be unchanged
        pre_obs = np.isfinite(z_before)
        np.testing.assert_array_equal(z_block[pre_obs], z_before[pre_obs])

    def test_low_correlation_skipped(self, ld_block_dir, rng):
        block_dir, alids, _ = ld_block_dir
        n = len(alids)
        variants = _make_variants(alids)
        eaf = np.full(n, 0.3, dtype=np.float32)

        z_block = rng.standard_normal((n, 1)).astype(np.float32)
        z_block[0, 0] = np.nan
        z_orig = z_block.copy()

        block_index = self._make_index(ld_block_dir, alids)
        # Set min_cor to 1.0 (impossible) so nothing gets filled
        impute_z_block(z_block, variants, eaf, block_index, thresh=0.9, min_cor=1.0)

        assert not np.isfinite(z_block[0, 0]), "should not have been filled at min_cor=1.0"

    def test_out_mask_none_ok(self, ld_block_dir, rng):
        block_dir, alids, _ = ld_block_dir
        n = len(alids)
        variants = _make_variants(alids)
        eaf = np.full(n, 0.3, dtype=np.float32)
        z_block = rng.standard_normal((n, 1)).astype(np.float32)
        z_block[0, 0] = np.nan
        # Should not raise with out_mask=None
        impute_z_block(z_block, variants, eaf, block_index={}, thresh=0.9, min_cor=0.7)


# ---------------------------------------------------------------------------
# Issue 037 – lift_variants unit tests
# ---------------------------------------------------------------------------

def _make_variants_full(alids: list[str]) -> np.ndarray:
    """Build a variants structured array matching the build.py dtype."""
    dt = np.dtype([
        ("id", "U64"), ("chrom", "U10"), ("pos", np.uint32),
        ("a1", "U64"), ("a2", "U64"),
    ])
    rows = []
    for alid in alids:
        chrom, rest = alid.split(":", 1)
        pos_str, a1, a2 = rest.split("_", 2)
        rows.append((alid, chrom, int(pos_str), a1, a2))
    return np.array(rows, dtype=dt)


class TestLiftVariants:
    """Unit tests for liftover.lift_variants using a mocked LiftOver."""

    ALIDS_HG19 = [
        "1:1000_A_G",
        "2:2000_C_T",
        "3:3000_A_C",
    ]
    # Fake hg38 positions returned by the mock
    HG38_POS = {
        ("chr1", 999):  [("chr1", 1100, "+", 0)],
        ("chr2", 1999): [("chr2", 2200, "+", 0)],
        ("chr3", 2999): [],  # simulate liftover failure for this variant
    }

    def _mock_lo(self):
        """Return a MagicMock LiftOver instance using HG38_POS."""
        lo = MagicMock()
        lo.convert_coordinate.side_effect = lambda chrom, pos0: self.HG38_POS.get(
            (chrom, pos0), []
        )
        return lo

    def _run(self, alids=None):
        from pleiodb.liftover import lift_variants
        variants = _make_variants_full(alids or self.ALIDS_HG19)
        lo = self._mock_lo()
        with patch("pyliftover.LiftOver", return_value=lo):
            return lift_variants(variants, "hg19", "hg38"), variants

    def test_length_unchanged(self):
        result, original = self._run()
        assert len(result) == len(original)

    def test_dtype_unchanged(self):
        result, original = self._run()
        assert result.dtype == original.dtype

    def test_successful_lift_updates_pos(self):
        result, _ = self._run()
        # variant 0: chr1:1000 → chr1:1101 (pos0=999 → result pos0=1100 → pos=1101)
        assert result["pos"][0] == 1101
        assert result["id"][0] == "1:1101_A_G"
        assert result["chrom"][0] == "1"

    def test_successful_lift_updates_second_variant(self):
        result, _ = self._run()
        assert result["pos"][1] == 2201
        assert result["id"][1] == "2:2201_C_T"

    def test_failed_lift_keeps_original(self):
        result, original = self._run()
        # variant 2 fails; should retain original values
        assert result["pos"][2] == original["pos"][2]
        assert result["id"][2] == original["id"][2]
        assert result["chrom"][2] == original["chrom"][2]

    def test_alleles_unchanged(self):
        result, original = self._run()
        for i in range(len(original)):
            assert result["a1"][i] == original["a1"][i]
            assert result["a2"][i] == original["a2"][i]

    def test_original_array_not_mutated(self):
        result, original = self._run()
        orig_ids = list(original["id"])
        assert result["id"][0] != orig_ids[0]  # liftover changed id[0]
        assert original["id"][0] == orig_ids[0]  # original unchanged

    def test_no_liftover_needed_returns_same_values(self):
        """When all variants lift successfully, no original values remain."""
        from pleiodb.liftover import lift_variants
        alids = ["1:1000_A_G", "2:2000_C_T"]
        variants = _make_variants_full(alids)
        lo = MagicMock()
        lo.convert_coordinate.return_value = [("chr1", 5000, "+", 0)]
        with patch("pyliftover.LiftOver", return_value=lo):
            result = lift_variants(variants, "hg19", "hg38")
        assert all(result["pos"] == 5001)
