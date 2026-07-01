"""
Conditional maximum-likelihood estimator for the sample-overlap-weighted
phenotypic correlation (rho) between two GWAS traits.

Reference: Forde et al., mr.simss::est_lambda
  https://github.com/amandaforde/mr.simss/blob/main/R/est_lambda.R

Public API
----------
estimate_rho_cml(z_j, z_k, z_thresh=1.0, min_nulls=500) -> float
    Returns the CML estimate of rho, or NaN if fewer than min_nulls
    null variants are available.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import multivariate_normal

# ---------------------------------------------------------------------------
# Normalising-constant grid
# ---------------------------------------------------------------------------

_GRID_N = 2000  # number of ρ points in (−1, 1)
_norm_grid_cache: dict[float, tuple[np.ndarray, np.ndarray]] = {}


def _build_norm_grid(z_thresh: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Precompute P(|X| < z_thresh, |Y| < z_thresh ; ρ) for a fine grid of ρ.

    Uses the bivariate normal CDF (exact, via scipy) evaluated at the four
    corners of the truncation square.  Returns (rho_grid, log_prob_grid).
    """
    if z_thresh in _norm_grid_cache:
        return _norm_grid_cache[z_thresh]

    rho_grid = np.linspace(-1 + 1e-6, 1 - 1e-6, _GRID_N)
    log_prob = np.empty(_GRID_N)

    z = float(z_thresh)
    corners = np.array([
        [ z,  z],
        [ z, -z],
        [-z,  z],
        [-z, -z],
    ])

    for i, r in enumerate(rho_grid):
        cov = [[1.0, r], [r, 1.0]]
        # P(X ≤ z, Y ≤ z) - P(X ≤ z, Y ≤ -z) - P(X ≤ -z, Y ≤ z) + P(X ≤ -z, Y ≤ -z)
        # = P(|X| < z, |Y| < z)
        cdf = multivariate_normal.cdf(corners, mean=[0.0, 0.0], cov=cov)
        prob = cdf[0] - cdf[1] - cdf[2] + cdf[3]
        log_prob[i] = math.log(max(prob, 1e-300))

    _norm_grid_cache[z_thresh] = (rho_grid, log_prob)
    return rho_grid, log_prob


# ---------------------------------------------------------------------------
# Public function
# ---------------------------------------------------------------------------

def estimate_rho_cml(
    z_j: np.ndarray,
    z_k: np.ndarray,
    z_thresh: float = 1.0,
    min_nulls: int = 500,
) -> float:
    """
    Conditional maximum-likelihood estimate of rho for a single trait pair.

    Parameters
    ----------
    z_j, z_k : float32 arrays
        Z-scores already filtered to |z_j| < z_thresh AND |z_k| < z_thresh.
    z_thresh : float
        The threshold used for filtering (required to compute the normalising
        constant).  Default 1.0.
    min_nulls : int
        Minimum number of null variants required.  Returns NaN when
        len(z_j) < min_nulls.  Default 500; set lower for small test datasets.

    Returns
    -------
    float
        CML estimate of rho in (−1, 1), or NaN if n < min_nulls.
    """
    z_j = np.asarray(z_j, dtype=np.float64)
    z_k = np.asarray(z_k, dtype=np.float64)
    n = len(z_j)

    if n < min_nulls:
        return float("nan")

    # Sufficient statistics
    A = float(np.dot(z_j, z_j))
    B = float(np.dot(z_j, z_k))
    C = float(np.dot(z_k, z_k))

    # Precompute normalising-constant grid once for this z_thresh
    rho_grid, log_prob_grid = _build_norm_grid(z_thresh)

    def neg_log_lik(rho: float) -> float:
        rho = float(np.clip(rho, -1 + 1e-9, 1 - 1e-9))
        one_minus_r2 = 1.0 - rho * rho
        # bivariate normal log-density terms (up to additive constants):
        #   −(A − 2ρB + C) / 2(1−ρ²)  − n/2 · log(1−ρ²)
        ll = -(A - 2.0 * rho * B + C) / (2.0 * one_minus_r2)
        ll -= 0.5 * n * math.log(one_minus_r2)
        # subtract log normalising constant (conditional correction)
        lp = float(np.interp(rho, rho_grid, log_prob_grid))
        ll -= n * lp
        return -ll  # minimise negative log-likelihood

    result = minimize_scalar(neg_log_lik, bounds=(-1 + 1e-6, 1 - 1e-6),
                             method="bounded")
    return float(result.x)
