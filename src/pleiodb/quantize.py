"""
Quantization helpers for GWAS statistics.

Z-score  → int16   (scale=100, NA=INT16_MIN, clipped ±327)
Neff     → uint16  (log2-encoded, 11 fractional bits, NA=0xFFFF)
EAF      → float16 (direct cast; effect allele = A2, alphabetically second)
Lambda   → float16 (direct cast)
"""

from __future__ import annotations

import numpy as np

# ---- z-score ----
Z_SCALE: int = 100          # stored = round(z * Z_SCALE)
Z_NA: int = np.iinfo(np.int16).min  # -32768
Z_MAX: int = np.iinfo(np.int16).max  # 32767
Z_MIN_VALID: int = Z_NA + 1  # -32767

# ---- neff ----
NEFF_FRAC_BITS: int = 11           # log2(N) stored with 2^11 fractional resolution
NEFF_NA: int = 0xFFFF              # sentinel for missing
_NEFF_SCALE: float = 2.0**NEFF_FRAC_BITS


def encode_z(z: np.ndarray) -> np.ndarray:
    """float32 z-scores → int16.  NaN maps to Z_NA."""
    out = np.where(
        np.isnan(z),
        Z_NA,
        np.clip(np.round(z * Z_SCALE), Z_MIN_VALID, Z_MAX),
    ).astype(np.int16)
    return out


def decode_z(arr: np.ndarray) -> np.ndarray:
    """int16 → float32.  Z_NA maps to NaN."""
    out = arr.astype(np.float32)
    out[arr == Z_NA] = np.nan
    out[arr != Z_NA] /= Z_SCALE
    return out


def encode_neff(neff: np.ndarray) -> np.ndarray:
    """
    float32 Neff → uint16.  NaN / non-positive → NEFF_NA.
    Stored value = round(log2(neff) * 2^NEFF_FRAC_BITS), clamped to [0, NEFF_NA-1].
    Recoverable relative error ≈ 0.034 %.
    """
    valid = np.isfinite(neff) & (neff > 0)
    stored = np.where(
        valid,
        np.clip(np.round(np.log2(np.where(valid, neff, 1.0)) * _NEFF_SCALE), 0, NEFF_NA - 1),
        NEFF_NA,
    ).astype(np.uint16)
    return stored


def decode_neff(arr: np.ndarray) -> np.ndarray:
    """uint16 → float32.  NEFF_NA maps to NaN."""
    out = np.full(arr.shape, np.nan, dtype=np.float32)
    mask = arr != NEFF_NA
    out[mask] = np.power(2.0, arr[mask].astype(np.float64) / _NEFF_SCALE).astype(np.float32)
    return out


def encode_eaf(eaf: np.ndarray) -> np.ndarray:
    """float32 → float16 (direct cast; NaN preserved)."""
    return eaf.astype(np.float16)


def decode_eaf(arr: np.ndarray) -> np.ndarray:
    return arr.astype(np.float32)


def reconstruct_beta_se(
    z: np.ndarray, neff: np.ndarray, eaf: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """
    Reconstruct (beta, se) from stored quantities.

    SE  = 1 / sqrt(Neff * 2 * p * (1 - p))
    beta = z * SE

    z    : (..., T) float32
    neff : (..., T) float32
    eaf  : (...,)   float32  — effect allele freq (A2); broadcast over T when shape is (V,)
    """
    if eaf.ndim < z.ndim:
        eaf = eaf[..., np.newaxis]
    denom = neff * 2.0 * eaf * (1.0 - eaf)
    se = np.where(denom > 0, 1.0 / np.sqrt(denom), np.nan)
    beta = z * se
    return beta, se
