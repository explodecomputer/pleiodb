"""
GWAS-VCF reader.  Extracts z-scores and effective sample sizes for a
pre-specified variant list.

Follows the GWAS-VCF spec (Lyon et al. 2021):
  FORMAT fields: ES (effect size), SE (standard error), SS (sample size)
  Optional:       EZ (pre-computed z-score), LP (−log10 p)
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Sequence

import numpy as np

log = logging.getLogger(__name__)

_MISSING = np.nan


def _open_vcf(path: str):
    try:
        import cyvcf2  # type: ignore
        return cyvcf2.VCF(path)
    except ImportError:
        raise ImportError("cyvcf2 is required for reading GWAS-VCF files: pip install cyvcf2")


def _variant_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"{chrom}:{pos}:{ref}:{alt}"


def read_vcf(
    vcf_path: str | Path,
    variant_ids: Sequence[str],
    id_col: str = "ID",
) -> tuple[np.ndarray, np.ndarray]:
    """
    Parse a GWAS-VCF file and return (z_scores, neff) float32 arrays aligned
    to `variant_ids`.  Missing entries are NaN.

    Parameters
    ----------
    vcf_path    : path to (possibly bgzipped + tabix-indexed) GWAS-VCF
    variant_ids : ordered list of variant identifiers to extract
    id_col      : "ID" to match on VCF ID field, or "CHRPOSREFALT" to
                  match on chr:pos:ref:alt key
    """
    path = str(vcf_path)
    n = len(variant_ids)
    z_out = np.full(n, _MISSING, dtype=np.float32)
    neff_out = np.full(n, _MISSING, dtype=np.float32)

    lookup: dict[str, int] = {vid: i for i, vid in enumerate(variant_ids)}

    vcf = _open_vcf(path)
    seen = 0

    for rec in vcf:
        if id_col == "ID":
            key = rec.ID or _variant_key(rec.CHROM, rec.POS, rec.REF, rec.ALT[0])
        else:
            key = _variant_key(rec.CHROM, rec.POS, rec.REF, rec.ALT[0])

        idx = lookup.get(key)
        if idx is None:
            continue

        # --- z-score ---
        try:
            ez = rec.format("EZ")
            if ez is not None and not np.isnan(ez[0][0]):
                z_out[idx] = float(ez[0][0])
            else:
                es = rec.format("ES")
                se = rec.format("SE")
                if es is not None and se is not None:
                    es_v, se_v = float(es[0][0]), float(se[0][0])
                    if se_v > 0:
                        z_out[idx] = es_v / se_v
        except (TypeError, IndexError, ValueError):
            pass

        # --- effective sample size ---
        try:
            ss = rec.format("SS")
            if ss is not None:
                neff_out[idx] = float(ss[0][0])
            else:
                ns = rec.INFO.get("NS")
                if ns:
                    neff_out[idx] = float(ns)
        except (TypeError, AttributeError, ValueError):
            pass

        seen += 1
        if seen == n:
            break

    vcf.close()
    return z_out, neff_out


def read_vcf_region(
    vcf_path: str | Path,
    variant_ids: Sequence[str],
    region: str,
    id_col: str = "ID",
) -> tuple[np.ndarray, np.ndarray]:
    """Like read_vcf but restricts to a tabix region string 'chr:start-end'."""
    path = str(vcf_path)
    n = len(variant_ids)
    z_out = np.full(n, _MISSING, dtype=np.float32)
    neff_out = np.full(n, _MISSING, dtype=np.float32)
    lookup: dict[str, int] = {vid: i for i, vid in enumerate(variant_ids)}

    vcf = _open_vcf(path)
    for rec in vcf(region):
        key = rec.ID if id_col == "ID" else _variant_key(rec.CHROM, rec.POS, rec.REF, rec.ALT[0])
        idx = lookup.get(key)
        if idx is None:
            continue
        try:
            es = rec.format("ES")
            se = rec.format("SE")
            if es is not None and se is not None:
                es_v, se_v = float(es[0][0]), float(se[0][0])
                if se_v > 0:
                    z_out[idx] = es_v / se_v
        except (TypeError, IndexError, ValueError):
            pass
        try:
            ss = rec.format("SS")
            if ss is not None:
                neff_out[idx] = float(ss[0][0])
        except (TypeError, ValueError):
            pass
    vcf.close()
    return z_out, neff_out
