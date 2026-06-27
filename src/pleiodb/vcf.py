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
) -> tuple[np.ndarray, np.ndarray]:
    """
    Parse a GWAS-VCF and return (z_scores, se) float32 arrays aligned to
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
    z_out  : float32(V) — z-scores (NaN = missing).
             Computed from FORMAT/EZ if present, else from FORMAT/ES ÷ FORMAT/SE.
    se_out : float32(V) — per-variant SE from VCF FORMAT/SE field (NaN = missing).
             Always positive; no allele-flip sign change needed.

    Notes
    -----
    The VCF ``SS``/``NS`` field (sample size) is **not** read.  Effective sample
    size per variant is derived from SE and var_y at build time (see
    :func:`pleiodb.build.derive_neff`).
    """
    path = str(vcf_path)
    n = max((idx for entries in pos_lookup.values() for _, _, idx in entries),
            default=-1) + 1
    if n == 0:
        empty = np.array([], dtype=np.float32)
        return empty, empty.copy()

    z_out = np.full(n, _MISSING, dtype=np.float32)
    se_out = np.full(n, _MISSING, dtype=np.float32)

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

        vcf.close()
        return z_out, se_out

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

def _lift_region(
    chrom: str,
    start: int,
    end: int,
    from_build: str,
    to_build: str,
) -> "tuple[str, int, int] | None":
    """Lift a genomic interval from *from_build* to *to_build*.

    Returns ``(new_bare_chrom, new_start, new_end)`` on success, or ``None``
    when either endpoint fails liftover (callers should return ``{}``).
    """
    try:
        from pyliftover import LiftOver  # type: ignore
    except ImportError:
        raise ImportError("pyliftover is required for cross-build region queries")

    from .liftover import normalise_build
    lo = LiftOver(normalise_build(from_build), normalise_build(to_build))

    bare = chrom.lstrip("chr")
    chrom_in = f"chr{bare}"

    r_start = lo.convert_coordinate(chrom_in, start - 1)
    r_end = lo.convert_coordinate(chrom_in, end - 1)

    if not r_start or not r_end:
        return None

    new_chrom = r_start[0][0].lstrip("chr")
    new_start = int(r_start[0][1]) + 1
    new_end = int(r_end[0][1]) + 1
    return new_chrom, min(new_start, new_end), max(new_start, new_end)


def read_vcf_region(
    vcf_path: str | Path,
    chrom: str,
    start: int,
    end: int,
    vcf_build: "str | None" = None,
) -> dict[str, float]:
    """Extract z-scores for all variants in a genomic region from a GWAS-VCF.

    Parameters
    ----------
    vcf_path  : path to a GWAS-VCF (bgzipped + indexed, or plain VCF)
    chrom     : chromosome in hg38 coordinates (bare, e.g. "1"; "chr1" ok)
    start     : region start (1-based, inclusive) in hg38 coordinates
    end       : region end   (1-based, inclusive) in hg38 coordinates
    vcf_build : genome build of the VCF (hg19/hg38/aliases).  When provided
                and different from hg38, the region is lifted from hg38 to
                *vcf_build* before issuing the bcftools query.  None / hg38
                → no liftover performed.

    Returns
    -------
    dict mapping ``"{bare_chrom}:{pos}_{ref}_{alt}"`` → z-score (effect of ALT).
    An empty dict is returned when no variants are found, liftover fails, or
    on any read error.

    The keys are compatible with ``build_block_allele_lookup``: callers should
    check ``is_flipped`` in the allele lookup and negate z when True.
    """
    from .liftover import normalise_build, builds_differ

    path = str(vcf_path)
    result: dict[str, float] = {}
    bare = chrom.lstrip("chr")

    # Liftover block coordinates when the VCF uses a different build from hg38
    if vcf_build is not None and builds_differ("hg38", vcf_build):
        lifted = _lift_region(bare, start, end, from_build="hg38", to_build=vcf_build)
        if lifted is None:
            log.warning(
                "read_vcf_region: liftover hg38→%s failed for %s:%d-%d – returning {}",
                vcf_build, chrom, start, end,
            )
            return result
        bare, start, end = lifted

    def _collect(recs) -> None:
        for rec in recs:
            if not rec.ALT:
                continue
            z = _extract_z(rec)
            if z is None:
                continue
            ref = compress_allele(rec.REF)
            alt = compress_allele(rec.ALT[0])
            rec_bare = rec.CHROM.lstrip("chr")
            result[f"{rec_bare}:{rec.POS}_{ref}_{alt}"] = z

    if not _has_index(path):
        vcf = _open_vcf(path)
        try:
            _collect(
                rec for rec in vcf
                if rec.CHROM.lstrip("chr") == bare and start <= rec.POS <= end
            )
        finally:
            vcf.close()
        return result

    # Indexed: try bcftools region filter with both bare and chr-prefixed chrom forms
    bcftools_ok = False
    for region_chrom in (bare, f"chr{bare}"):
        region = f"{region_chrom}:{start}-{end}"
        fd, tmp_vcf = tempfile.mkstemp(suffix=".vcf")
        os.close(fd)
        try:
            subprocess.run(
                ["bcftools", "view", "-r", region, "-Ov", path, "-o", tmp_vcf],
                stderr=subprocess.DEVNULL,
            )
            bcftools_ok = True
            vcf = _open_vcf(tmp_vcf)
            try:
                recs = list(vcf)
            finally:
                vcf.close()
            if recs:
                _collect(recs)
                return result
        except (OSError, FileNotFoundError) as exc:
            log.debug("bcftools region query failed (%s) – will fall back", exc)
            break
        finally:
            try:
                os.unlink(tmp_vcf)
            except OSError:
                pass

    if not bcftools_ok:
        # bcftools not available: full cyvcf2 scan filtered by position
        vcf = _open_vcf(path)
        try:
            _collect(
                rec for rec in vcf
                if rec.CHROM.lstrip("chr") == bare and start <= rec.POS <= end
            )
        finally:
            vcf.close()

    return result


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


