"""
GWAS-VCF reader.  Extracts z-scores and effective sample sizes for a
pre-specified variant list.

Follows the GWAS-VCF spec (Lyon et al. 2021):
  FORMAT fields: ES (effect size), SE (standard error), SS (sample size)
  Optional:       EZ (pre-computed z-score), LP (−log10 p)
"""

from __future__ import annotations

import logging
import os
import subprocess
import tempfile
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
    key_lookup: dict[str, int] | None = None,
    regions_file: str | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Parse a GWAS-VCF file and return (z_scores, neff) float32 arrays aligned
    to `variant_ids`.  Missing entries are NaN.

    Parameters
    ----------
    vcf_path     : path to (possibly bgzipped + tabix-indexed) GWAS-VCF
    variant_ids  : ordered list of variant identifiers to extract
    id_col       : "ID" to match on VCF ID field, or "CHRPOSREFALT" to
                   match on chr:pos:ref:alt key.  Ignored when key_lookup is
                   provided.
    key_lookup   : optional pre-built mapping of 'chrom:pos:ref:alt' (in the
                   VCF's own coordinate system) → index into variant_ids.
                   Supply this when the variant list and VCF are on different
                   genome builds (positions already lifted by the caller).
    regions_file : optional path to a CHROM-POS TSV for bcftools -R pre-filter.
                   Requires bcftools on PATH and a tabix-indexed VCF; falls back
                   to full-scan if bcftools fails.
    """
    path = str(vcf_path)
    n = len(variant_ids)
    z_out = np.full(n, _MISSING, dtype=np.float32)
    neff_out = np.full(n, _MISSING, dtype=np.float32)

    if key_lookup is not None:
        lookup = key_lookup
    else:
        lookup = {vid: i for i, vid in enumerate(variant_ids)}

    tmp_vcf: str | None = None
    if regions_file is not None:
        fd, tmp_vcf = tempfile.mkstemp(suffix=".vcf")
        os.close(fd)
        try:
            subprocess.run(
                ["bcftools", "view", "-R", regions_file, "-Ov", path, "-o", tmp_vcf],
                check=True,
                stderr=subprocess.DEVNULL,
            )
            vcf = _open_vcf(tmp_vcf)
        except (subprocess.CalledProcessError, FileNotFoundError) as exc:
            log.warning("bcftools pre-filter failed (%s), falling back to full VCF scan", exc)
            os.unlink(tmp_vcf)
            tmp_vcf = None
            vcf = _open_vcf(path)
    else:
        vcf = _open_vcf(path)

    seen = 0

    for rec in vcf:
        if key_lookup is not None:
            key = _variant_key(rec.CHROM, rec.POS, rec.REF, rec.ALT[0])
        elif id_col == "ID":
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
    if tmp_vcf is not None:
        try:
            os.unlink(tmp_vcf)
        except OSError:
            pass
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
