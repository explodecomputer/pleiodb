"""
Unit tests for estimate_rho_cml (issue #31).

These tests do NOT require bcftools or cyvcf2 — they use synthetic numpy data only.
"""

from __future__ import annotations

import numpy as np
import pytest

from pleiodb import estimate_rho_cml


RNG = np.random.default_rng(42)

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _bivariate_null(rho: float, n: int = 5000, z_thresh: float = 1.0,
                    rng=None) -> tuple[np.ndarray, np.ndarray]:
    """Draw n variants from a bivariate normal with correlation rho,
    then keep only rows where |z_j| < z_thresh AND |z_k| < z_thresh."""
    if rng is None:
        rng = RNG
    # draw until we have enough null variants
    zs = []
    while sum(len(b) for b in zs) < n:
        batch = rng.multivariate_normal([0, 0], [[1, rho], [rho, 1]], size=n * 5)
        mask = (np.abs(batch[:, 0]) < z_thresh) & (np.abs(batch[:, 1]) < z_thresh)
        zs.append(batch[mask])
    combined = np.vstack(zs)[:n]
    return combined[:, 0].astype(np.float32), combined[:, 1].astype(np.float32)


# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------

class TestEstimateRhoCml:
    def test_nan_when_too_few_null_variants(self):
        """With fewer than min_nulls variants the function returns NaN."""
        z = np.zeros(10, dtype=np.float32)
        result = estimate_rho_cml(z, z, z_thresh=1.0, min_nulls=500)
        assert np.isnan(result)

    def test_min_nulls_override(self):
        """min_nulls can be lowered so small datasets return a finite estimate."""
        z_j, z_k = _bivariate_null(0.0, n=50)
        result = estimate_rho_cml(z_j, z_k, z_thresh=1.0, min_nulls=50)
        assert np.isfinite(result)

    def test_zero_correlation(self):
        """Independent normals → estimate ≈ 0."""
        z_j, z_k = _bivariate_null(0.0, n=5000)
        result = estimate_rho_cml(z_j, z_k, z_thresh=1.0)
        assert abs(float(result)) < 0.05, f"Expected ≈0, got {result:.4f}"

    def test_positive_correlation(self):
        """ρ = 0.5 is recovered within ±0.05."""
        z_j, z_k = _bivariate_null(0.5, n=5000)
        result = estimate_rho_cml(z_j, z_k, z_thresh=1.0)
        assert abs(float(result) - 0.5) < 0.05, f"Expected ≈0.5, got {result:.4f}"

    def test_negative_correlation(self):
        """ρ = -0.3 is recovered within ±0.05."""
        z_j, z_k = _bivariate_null(-0.3, n=5000)
        result = estimate_rho_cml(z_j, z_k, z_thresh=1.0)
        assert abs(float(result) - (-0.3)) < 0.05, f"Expected ≈-0.3, got {result:.4f}"

    def test_result_in_range(self):
        """Output is always in the open interval (-1, 1) for finite results."""
        z_j, z_k = _bivariate_null(0.3, n=5000)
        result = estimate_rho_cml(z_j, z_k, z_thresh=1.0)
        assert -1.0 < float(result) < 1.0

    def test_different_z_thresh_consistent(self):
        """Estimates at z_thresh=0.5 and z_thresh=1.5 are within ±0.1 of each other."""
        # Use the same underlying population (same rng seed) but different thresholds
        rng = np.random.default_rng(99)
        z_j_lo, z_k_lo = _bivariate_null(0.4, n=5000, z_thresh=0.5, rng=rng)
        rng2 = np.random.default_rng(99)
        z_j_hi, z_k_hi = _bivariate_null(0.4, n=5000, z_thresh=1.5, rng=rng2)
        r_lo = estimate_rho_cml(z_j_lo, z_k_lo, z_thresh=0.5)
        r_hi = estimate_rho_cml(z_j_hi, z_k_hi, z_thresh=1.5)
        assert abs(float(r_lo) - float(r_hi)) < 0.1, (
            f"thresh=0.5 → {r_lo:.4f}, thresh=1.5 → {r_hi:.4f}"
        )
