"""
GWAS-VCF reader.  Extracts z-scores and effective sample sizes for a
pre-specified variant list using positional matching (CHROM:POS → allele check).

Follows the GWAS-VCF spec (Lyon et al. 2021):
  FORMAT fields: ES (effect size), SE (standard error), SS (sample size)
  Optional:       EZ (pre-computed z-score)

Allele convention: effect allele = A2 (alphabetically second).  When the VCF
record has genome-ref REF > ALT alphabetically (non-canonical orientation),
the z-score is negated so it reflects the effect of A2.
"""

from __future__ import annotations

import logging
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from .alid import compress_allele

log = logging.getLogger(__name__)

_MISSING = np.nan


def _open_vcf(path: str):
    try:
        import cyvcf2  # type: ignore
        return cyvcf2.VCF(path)
    except ImportError:
        raise ImportError("cyvcf2 is required: install via conda (bioconda channel)")


def _has_index(vcf_path: str) -> bool:
    return Path(vcf_path + ".csi").exists() or Path(vcf_path + ".tbi").exists()


def _build_regions_file(pos_lookup: dict[str, list]) -> str:
    """Write a CHROM-POS TSV temp file from pos_lookup keys for bcftools -R."""
    seen: set[tuple[str, int]] = set()
    for key in pos_lookup:
        chrom, pos_str = key.split(":", 1)
        chrom_bare = chrom.lstrip("chr")
        seen.add((chrom_bare, int(pos_str)))

    def _sort_key(cp: tuple[str, int]):
        try:
            return (int(cp[0]), cp[1])
        except ValueError:
            return (999, cp[1])

    rows = sorted(seen, key=_sort_key)
    fd, path = tempfile.mkstemp(suffix=".tsv", prefix="pleiodb_regions_")
    with os.fdopen(fd, "w") as fh:
        for chrom, pos in rows:
            fh.write(f"{chrom}\t{pos}\n")
            fh.write(f"chr{chrom}\t{pos}\n")
    return path


def read_vcf(
    vcf_path: str | Path,
    pos_lookup: dict[str, list[tuple[str, str, int]]],
    regions_file: str | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Parse a GWAS-VCF and return (z_scores, se, neff) float32 arrays aligned to
    the variant list encoded in *pos_lookup*.

    Parameters
    ----------
    vcf_path    : path to (bgzipped + CSI/TBI-indexed) GWAS-VCF
    pos_lookup  : ``{chrom:pos → [(a1, a2, row_idx), ...]}`` — both bare and
                  chr-prefixed chromosome forms should be present.  Built from
                  the ALID variant list (direct coords) or from make_lifted_lookup
                  (liftover coords).
    regions_file : pre-built CHROM-POS TSV for bcftools -R.  If None, a temp
                  file is built from pos_lookup keys and used when the VCF has
                  a CSI or TBI index.  Falls back to full cyvcf2 scan on any
                  bcftools failure.

    Returns
    -------
    z_out   : float32(V) — z-scores (NaN = missing)
    se_out  : float32(V) — per-variant SE from VCF FORMAT/SE field (NaN = missing)
    neff_out: float32(V) — Neff from VCF SS/NS field (NaN = missing; deprecated,
              will be replaced in a later version)
    """
    path = str(vcf_path)
    n = max((idx for entries in pos_lookup.values() for _, _, idx in entries),
            default=-1) + 1
    if n == 0:
        empty = np.array([], dtype=np.float32)
        return empty, empty.copy(), empty.copy()

    z_out = np.full(n, _MISSING, dtype=np.float32)
    se_out = np.full(n, _MISSING, dtype=np.float32)
    neff_out = np.full(n, _MISSING, dtype=np.float32)

    own_regions_file = False
    tmp_vcf: str | None = None

    try:
        # --- bcftools pre-filter -------------------------------------------
        if regions_file is None and _has_index(path):
            regions_file = _build_regions_file(pos_lookup)
            own_regions_file = True

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
                log.debug("bcftools pre-filter failed (%s), falling back to full scan", exc)
                if tmp_vcf and os.path.exists(tmp_vcf):
                    os.unlink(tmp_vcf)
                tmp_vcf = None
                vcf = _open_vcf(path)
        else:
            vcf = _open_vcf(path)

        # --- iterate and match by position + alleles -----------------------
        for rec in vcf:
            chrom = rec.CHROM
            key = f"{chrom}:{rec.POS}"
            candidates = pos_lookup.get(key)
            if not candidates:
                continue

            vcf_ref = compress_allele(rec.REF)
            vcf_alt_raw = rec.ALT[0] if rec.ALT else None
            if vcf_alt_raw is None:
                continue
            vcf_alt = compress_allele(vcf_alt_raw)

            for a1, a2, idx in candidates:
                if vcf_ref == a1 and vcf_alt == a2:
                    flip = False
                elif vcf_ref == a2 and vcf_alt == a1:
                    flip = True
                else:
                    continue  # allele mismatch — silent NaN

                z = _extract_z(rec)
                if z is not None:
                    z_out[idx] = -z if flip else z

                se = _extract_se(rec)
                if se is not None:
                    se_out[idx] = se   # SE is always positive; no flip needed

                neff = _extract_neff(rec)
                if neff is not None:
                    neff_out[idx] = neff

        vcf.close()

    finally:
        if tmp_vcf and os.path.exists(tmp_vcf):
            try:
                os.unlink(tmp_vcf)
            except OSError:
                pass
        if own_regions_file and regions_file and os.path.exists(regions_file):
            try:
                os.unlink(regions_file)
            except OSError:
                pass

    return z_out, se_out, neff_out


def _extract_z(rec) -> float | None:
    try:
        ez = rec.format("EZ")
        if ez is not None:
            v = float(ez[0][0])
            if not np.isnan(v):
                return v
    except (TypeError, IndexError, ValueError):
        pass
    try:
        es = rec.format("ES")
        se = rec.format("SE")
        if es is not None and se is not None:
            es_v, se_v = float(es[0][0]), float(se[0][0])
            if se_v > 0:
                return es_v / se_v
    except (TypeError, IndexError, ValueError):
        pass
    return None


def _extract_se(rec) -> float | None:
    """Return SE from FORMAT/SE field; None if absent or zero."""
    try:
        se = rec.format("SE")
        if se is not None:
            v = float(se[0][0])
            if v > 0:
                return v
    except (TypeError, IndexError, ValueError):
        pass
    return None


def _extract_neff(rec) -> float | None:
    try:
        ss = rec.format("SS")
        if ss is not None:
            return float(ss[0][0])
    except (TypeError, ValueError):
        pass
    try:
        ns = rec.INFO.get("NS")
        if ns:
            return float(ns)
    except (TypeError, AttributeError, ValueError):
        pass
    return None
