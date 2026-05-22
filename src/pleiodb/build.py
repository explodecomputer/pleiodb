"""
Build a GWASDatabase from a variant list and a TSV of trait → GWAS-VCF paths.

Variant input format (tab-separated, no header):
    ALID  EAF
    e.g.  10:101558746_G_T  0.56079

    ALID = CHROM:POS_A1_A2 where A1 ≤ A2 alphabetically (canonical order).
    EAF  = frequency of the effect allele (A2).  Second column is optional;
           omit or use NaN when frequencies are unavailable.
    Alleles longer than 20 characters are compressed to
    ``{allele[:8]}~{sha256(allele)[:4]}`` (see :mod:`pleiodb.alid`).

Trait TSV format (tab-separated, **with header row**):
    trait_id  trait_name  N  [K]  vcf_path  [build]

    N      = total sample size (integer, required).
    K      = case fraction in (0, 1), required for binary traits.
             Absent or empty → trait treated as continuous.
    build  = genome build of the VCF (hg19 / hg38); optional.

Processing strategy:
  Traits are processed in batches of T_BATCH.  Within each batch, VCF files
  are read in parallel (thread pool).  This keeps peak memory at:
    V × T_BATCH × (4 + 4) bytes  ≈  100k × 512 × 8 ≈ 400 MB
"""

from __future__ import annotations

import json
import logging
import os
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import zstandard as zstd

from .alid import canonical_alid, compress_allele
from .quantize import encode_z, encode_neff, encode_eaf, Z_NA
from .store import ChunkedMatrix
from .vcf import read_vcf
from .liftover import builds_differ, make_lifted_lookup, normalise_build

log = logging.getLogger(__name__)
_CCTX = zstd.ZstdCompressor(level=3, threads=-1)

DEFAULT_CHUNK = (512, 512)


# ---------------------------------------------------------------------------
# Trait metadata
# ---------------------------------------------------------------------------

@dataclass
class TraitInfo:
    """Per-trait metadata loaded from the traits input TSV."""
    trait_id: str
    trait_name: str
    N: int              # total sample size (required)
    K: float | None     # case fraction in (0, 1); None = continuous trait
    vcf_path: str
    vcf_build: str | None


def load_trait_list(trait_tsv: str | Path) -> list[TraitInfo]:
    """Parse the traits input TSV and return a list of :class:`TraitInfo`.

    File format (tab-separated, **with header row**):

    .. code-block:: text

        trait_id  trait_name  N  [K]  vcf_path  [build]

    Column rules
    ------------
    - Header is required; columns are identified by name, not position.
    - ``N`` (integer sample size) is required for every trait.
    - ``K`` (case fraction) is optional: if the column is absent or the cell
      is empty the trait is treated as continuous (``K=None``).
    - ``K`` must be strictly inside the open interval (0, 1).
    - ``build`` is optional.
    - Lines starting with ``#`` and blank lines are silently skipped.
    """
    path = Path(trait_tsv)
    header: list[str] | None = None
    traits: list[TraitInfo] = []

    with open(path) as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            if header is None:
                # First non-blank, non-comment line is the header
                header = [col.strip() for col in line.split("\t")]
                if "N" not in header:
                    raise ValueError(
                        f"traits TSV {path} is missing required 'N' column "
                        f"(found columns: {header})"
                    )
                continue

            cells = [c.strip() for c in line.split("\t")]
            row: dict[str, str] = dict(zip(header, cells))

            trait_id = row.get("trait_id", "")
            trait_name = row.get("trait_name", "")
            vcf_path = row.get("vcf_path", "")
            vcf_build_raw = row.get("build", "")
            vcf_build = vcf_build_raw if vcf_build_raw else None

            # --- N (required) ---
            n_raw = row.get("N", "").strip()
            if not n_raw:
                raise ValueError(
                    f"Trait '{trait_id}' has an empty N value in {path}"
                )
            try:
                N = int(n_raw)
            except ValueError:
                raise ValueError(
                    f"Trait '{trait_id}' has non-integer N value '{n_raw}' in {path}"
                )

            # --- K (optional) ---
            k_raw = row.get("K", "").strip() if "K" in header else ""
            K: float | None
            if k_raw:
                try:
                    K = float(k_raw)
                except ValueError:
                    raise ValueError(
                        f"Trait '{trait_id}' has non-numeric K value '{k_raw}' in {path}"
                    )
                if not (0.0 < K < 1.0):
                    raise ValueError(
                        f"Trait '{trait_id}' has K={K} outside (0, 1) in {path}"
                    )
            else:
                K = None

            traits.append(TraitInfo(
                trait_id=trait_id,
                trait_name=trait_name,
                N=N,
                K=K,
                vcf_path=vcf_path,
                vcf_build=vcf_build,
            ))

    if header is None:
        raise ValueError(f"traits TSV {path} appears to be empty")

    return traits
DEFAULT_T_BATCH = 512
DEFAULT_WORKERS = 8
DEFAULT_PVAL_THRESHOLDS = [5e-8, 1e-5]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_raw_alid(alid: str) -> tuple[str, int, str, str]:
    """Parse 'CHROM:POS_A1_A2' → (chrom, pos, raw_a1, raw_a2) without canonicalisation.

    The returned alleles are in the original input order; call
    ``canonical_alid`` to obtain the canonical form and flip flag.
    """
    colon_idx = alid.index(":")
    chrom = alid[:colon_idx]
    rest = alid[colon_idx + 1:]
    parts = rest.split("_", 2)
    if len(parts) != 3:
        raise ValueError(f"Cannot parse ALID: {alid!r}  (expected CHROM:POS_A1_A2)")
    pos = int(parts[0])
    return chrom, pos, parts[1], parts[2]


def _load_variants(
    variants_path: str | Path,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Read variant input file (ALID | EAF).

    Returns
    -------
    variants : structured array with fields (id, chrom, pos, a1, a2)
    eaf      : float32 array, NaN where frequency is absent
    """
    dt = np.dtype([
        ("id",    "U64"),
        ("chrom", "U10"),
        ("pos",   np.uint32),
        ("a1",    "U64"),   # compressed alleles are at most 13 chars; 64 is generous
        ("a2",    "U64"),
    ])
    rows = []
    eaf_list = []
    with open(variants_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            raw_alid = parts[0]
            raw_eaf = float(parts[1]) if len(parts) > 1 and parts[1].strip() else np.nan

            chrom, pos, raw_a1, raw_a2 = _parse_raw_alid(raw_alid)
            alid_str, was_flipped = canonical_alid(chrom, pos, raw_a1, raw_a2)
            # Extract compressed alleles from the canonical ALID string
            _, allele_part = alid_str.split(":", 1)
            _, ca1, ca2 = allele_part.split("_", 2)
            eaf = (1.0 - raw_eaf) if (was_flipped and not np.isnan(raw_eaf)) else raw_eaf

            rows.append((alid_str, chrom, pos, ca1, ca2))
            eaf_list.append(eaf)

    return np.array(rows, dtype=dt), np.array(eaf_list, dtype=np.float32)




def _build_pos_lookup(
    variants: np.ndarray,
) -> dict[str, list[tuple[str, str, int]]]:
    """Build ``{chrom:pos → [(a1, a2, row_idx), ...]}`` from the variant array.

    Both bare and chr-prefixed CHROM forms are inserted.
    """
    lookup: dict[str, list[tuple[str, str, int]]] = {}
    for i, row in enumerate(variants):
        chrom = str(row["chrom"])
        pos = str(int(row["pos"]))
        a1 = str(row["a1"])
        a2 = str(row["a2"])
        chrom_bare = chrom.lstrip("chr")
        for c in (chrom_bare, f"chr{chrom_bare}"):
            lookup.setdefault(f"{c}:{pos}", []).append((a1, a2, i))
    return lookup


def _make_regions_file(pos_lookup: dict[str, list]) -> str:
    """Write a sorted CHROM-POS TSV temp file from pos_lookup keys for bcftools -R."""
    seen: set[tuple[str, int]] = set()
    for key in pos_lookup:
        chrom, pos_str = key.split(":", 1)
        seen.add((chrom.lstrip("chr"), int(pos_str)))

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
    log.info("Regions file: %d positions → %s", len(rows), path)
    return path


def _fetch_trait(args: tuple) -> tuple[int, np.ndarray, np.ndarray]:
    """Worker: read one VCF and return (local_t_idx, z, neff)."""
    local_t_idx, trait_id, vcf_path, pos_lookup, regions_file, n_variants = args
    log.info("  reading %-30s  %s", trait_id, vcf_path)
    z, neff = read_vcf(vcf_path, pos_lookup, regions_file=regions_file)
    # Pad / trim to exactly n_variants (pos_lookup may cover more rows than V
    # if the lifted lookup has extra entries from chr/bare duplication)
    z_out = np.full(n_variants, np.nan, dtype=np.float32)
    neff_out = np.full(n_variants, np.nan, dtype=np.float32)
    length = min(len(z), n_variants)
    z_out[:length] = z[:length]
    neff_out[:length] = neff[:length]
    n_hit = int(np.isfinite(z_out).sum())
    log.info("  done    %-30s  %d/%d variants matched", trait_id, n_hit, n_variants)
    return local_t_idx, z_out, neff_out


# ---------------------------------------------------------------------------
# Main build routine
# ---------------------------------------------------------------------------

def build_database(
    output_dir: str | Path,
    variants_path: str | Path,
    trait_tsv: str | Path,
    chunk_shape: tuple[int, int] = DEFAULT_CHUNK,
    t_batch: int = DEFAULT_T_BATCH,
    workers: int = DEFAULT_WORKERS,
    pval_thresholds: Sequence[float] = DEFAULT_PVAL_THRESHOLDS,
    overwrite: bool = False,
    variants_build: str | None = None,
) -> None:
    out = Path(output_dir)
    if out.exists() and not overwrite:
        raise FileExistsError(f"{out} already exists; pass overwrite=True to replace")
    out.mkdir(parents=True, exist_ok=True)
    (out / "masks").mkdir(exist_ok=True)

    variants, eaf = _load_variants(variants_path)
    trait_pairs = load_trait_list(trait_tsv)
    V = len(variants)
    T = len(trait_pairs)
    CV, CT = chunk_shape

    canon_build: str | None = normalise_build(variants_build) if variants_build else None
    log.info(
        "Building pleiodb: V=%d  T=%d  chunks=%s  variants_build=%s",
        V, T, chunk_shape, canon_build or "unspecified",
    )

    # --- Build pos_lookup and regions file for same-build traits ------------
    direct_pos_lookup = _build_pos_lookup(variants)
    direct_regions_file = _make_regions_file(direct_pos_lookup)

    # --- Per-build liftover cache -------------------------------------------
    _lifted_cache: dict[str, dict] = {}
    _lifted_regions_cache: dict[str, str] = {}

    def _lookup_for(vcf_build: str | None):
        if canon_build is None or vcf_build is None:
            return direct_pos_lookup, direct_regions_file
        vcf_build_norm = normalise_build(vcf_build)
        if vcf_build_norm == canon_build:
            return direct_pos_lookup, direct_regions_file
        if vcf_build_norm not in _lifted_cache:
            _lifted_cache[vcf_build_norm] = make_lifted_lookup(
                variants, from_build=canon_build, to_build=vcf_build_norm
            )
            _lifted_regions_cache[vcf_build_norm] = _make_regions_file(
                _lifted_cache[vcf_build_norm]
            )
        return _lifted_cache[vcf_build_norm], _lifted_regions_cache[vcf_build_norm]

    # ---- Write variant / trait metadata ------------------------------------
    np.save(out / "variants.npy", variants)
    trait_dt = np.dtype([("id", "U64"), ("name", "U256")])
    trait_arr = np.array(
        [(t.trait_id, t.trait_name) for t in trait_pairs], dtype=trait_dt
    )
    np.save(out / "traits.npy", trait_arr)

    # ---- Initialise chunked matrices ---------------------------------------
    zscore_mat = ChunkedMatrix(out / "zscore", (V, T), np.int16, chunk_shape)
    neff_mat = ChunkedMatrix(out / "neff", (V, T), np.uint16, chunk_shape)
    zscore_mat.open_write()
    neff_mat.open_write()

    neff_sum = np.zeros(T, dtype=np.float64)
    neff_count = np.zeros(T, dtype=np.int64)

    # ---- Main ingestion loop -----------------------------------------------
    for t_block_start in range(0, T, CT):
        t_block_end = min(t_block_start + CT, T)
        batch_traits = trait_pairs[t_block_start:t_block_end]
        B = len(batch_traits)

        z_block = np.full((V, B), np.nan, dtype=np.float32)
        neff_block = np.full((V, B), np.nan, dtype=np.float32)

        args = [
            (j, t.trait_id, t.vcf_path, *_lookup_for(t.vcf_build), V)
            for j, t in enumerate(batch_traits)
        ]
        n_done = 0
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(_fetch_trait, a): a[0] for a in args}
            for fut in as_completed(futures):
                j, z, neff = fut.result()
                z_block[:, j] = z
                neff_block[:, j] = neff
                n_done += 1
                log.info("  progress: %d/%d traits complete", t_block_start + n_done, T)

        valid = np.isfinite(neff_block)
        neff_sum[t_block_start:t_block_end] += np.nansum(neff_block, axis=0)
        neff_count[t_block_start:t_block_end] += valid.sum(axis=0)

        for v_block_start in range(0, V, CV):
            v_block_end = min(v_block_start + CV, V)
            z_chunk = encode_z(z_block[v_block_start:v_block_end, :])
            neff_chunk = encode_neff(neff_block[v_block_start:v_block_end, :])
            zscore_mat.write_chunk(v_block_start // CV, t_block_start // CT, z_chunk)
            neff_mat.write_chunk(v_block_start // CV, t_block_start // CT, neff_chunk)

        log.info("  traits %d–%d done", t_block_start, t_block_end - 1)

    zscore_mat.close_write()
    neff_mat.close_write()

    # ---- Neff base ---------------------------------------------------------
    neff_base = np.where(
        neff_count > 0, neff_sum / neff_count, np.nan
    ).astype(np.float32)
    neff_base.tofile(out / "neff_base.f32")

    # ---- EAF ---------------------------------------------------------------
    encode_eaf(eaf).tofile(out / "eaf.f16")

    # ---- Significance masks ------------------------------------------------
    log.info("Building significance masks (%s)…", pval_thresholds)
    _build_masks(out, zscore_mat, V, T, pval_thresholds)

    # ---- Clean up temp files -----------------------------------------------
    for rf in [direct_regions_file] + list(_lifted_regions_cache.values()):
        try:
            os.unlink(rf)
        except OSError:
            pass

    # ---- Metadata ----------------------------------------------------------
    meta = {
        "V": V,
        "T": T,
        "chunk_shape": list(chunk_shape),
        "pval_thresholds": list(pval_thresholds),
        "z_scale": 100,
        "neff_encoding": "log2_u16_frac11",
        "format_version": 2,
        "variants_build": canon_build,
    }
    (out / "meta.json").write_text(json.dumps(meta, indent=2))
    log.info("Build complete → %s", out)


# ---------------------------------------------------------------------------
# Significance masks
# ---------------------------------------------------------------------------

def _build_masks(
    out: Path,
    zscore_mat: ChunkedMatrix,
    V: int,
    T: int,
    thresholds: Sequence[float],
) -> None:
    from .quantize import decode_z
    from scipy.stats import norm  # type: ignore

    z_thresholds = {p: abs(norm.ppf(p / 2)) for p in thresholds}
    hits: dict[float, list[tuple[np.ndarray, np.ndarray]]] = {p: [] for p in thresholds}
    CV, CT = zscore_mat.chunk_shape

    for vi in range(zscore_mat.n_v_chunks):
        for ti in range(zscore_mat.n_t_chunks):
            raw = zscore_mat.get_raw_chunk(vi, ti)
            z = decode_z(raw)
            v_off = vi * CV
            t_off = ti * CT
            for p, zt in z_thresholds.items():
                hv, ht = np.where(np.abs(z) >= zt)
                if len(hv):
                    hits[p].append((
                        (hv + v_off).astype(np.uint32),
                        (ht + t_off).astype(np.uint32),
                    ))
        log.info("mask scan: v-block %d/%d", vi + 1, zscore_mat.n_v_chunks)

    for p, batches in hits.items():
        if batches:
            v_all = np.concatenate([b[0] for b in batches])
            t_all = np.concatenate([b[1] for b in batches])
        else:
            v_all = np.array([], dtype=np.uint32)
            t_all = np.array([], dtype=np.uint32)

        order = np.lexsort((t_all, v_all))
        pairs = np.column_stack([v_all[order], t_all[order]])
        blob = _CCTX.compress(pairs.tobytes())
        mask_name = f"{p:.0e}".replace("+", "")
        (out / "masks" / f"{mask_name}.coo.zst").write_bytes(blob)
        log.info("  mask %s: %d hits", mask_name, len(v_all))


# ---------------------------------------------------------------------------
# Lambda (sample overlap) computation
# ---------------------------------------------------------------------------

def build_lambda(
    db_path: str | Path,
    n_null_per_trait: int = 5000,
    z_null_thresh: float = 3.0,
    chunk_shape: tuple[int, int] = (256, 256),
    workers: int = 8,
) -> None:
    from .db import GWASDatabase
    from .quantize import decode_z

    db = GWASDatabase.open(db_path)
    out = db.path
    T = db.T
    zscore = db._zscore

    lam_mat = ChunkedMatrix(out / "lambda", (T, T), np.float16, chunk_shape)
    lam_mat.open_write()

    CT = chunk_shape[1]
    for ti in range(lam_mat.n_v_chunks):
        t0i = ti * CT
        t1i = min(t0i + CT, T)
        col_z = decode_z(zscore.get_block(0, db.V, t0i, t1i))

        for tj in range(lam_mat.n_t_chunks):
            t0j = tj * CT
            t1j = min(t0j + CT, T)
            row_z = decode_z(zscore.get_block(0, db.V, t0j, t1j))

            null_i = np.abs(col_z) < z_null_thresh
            null_j = np.abs(row_z) < z_null_thresh

            block = np.full((t1i - t0i, t1j - t0j), np.nan, dtype=np.float32)
            for a in range(t1i - t0i):
                mask_a = null_i[:, a]
                for b in range(t1j - t0j):
                    mask_ab = mask_a & null_j[:, b]
                    n = mask_ab.sum()
                    if n >= 30:
                        za = col_z[mask_ab, a]
                        zb = row_z[mask_ab, b]
                        block[a, b] = float(np.corrcoef(za, zb)[0, 1])

            lam_mat.write_chunk(ti, tj, block.astype(np.float16))
        log.info("lambda: t-block %d/%d done", ti + 1, lam_mat.n_v_chunks)

    lam_mat.close_write()

    meta = json.loads((out / "meta.json").read_text())
    meta["lambda_chunk_shape"] = list(chunk_shape)
    meta["lambda_null_z_thresh"] = z_null_thresh
    (out / "meta.json").write_text(json.dumps(meta, indent=2))
    log.info("Lambda matrix written")
