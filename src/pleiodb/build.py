"""
Build a GWASDatabase from a variant list and a TSV of trait → GWAS-VCF paths.

Variant input format (tab-separated, no header):
    ALID  EAF
    e.g.  10:101558746_G_T  0.56079

    ALID = CHROM:POS_A1_A2 where A1 < A2 alphabetically.
    EAF  = frequency of the effect allele (A2).  Second column is optional;
           omit or use NaN when frequencies are unavailable.

Trait TSV format (tab-separated, no header):
    trait_id  trait_name  vcf_path  [build]

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
from pathlib import Path
from typing import Sequence

import numpy as np
import zstandard as zstd

from .quantize import encode_z, encode_neff, encode_eaf, Z_NA
from .store import ChunkedMatrix
from .vcf import read_vcf
from .liftover import builds_differ, make_lifted_lookup, normalise_build

log = logging.getLogger(__name__)
_CCTX = zstd.ZstdCompressor(level=3, threads=-1)

DEFAULT_CHUNK = (512, 512)
DEFAULT_T_BATCH = 512
DEFAULT_WORKERS = 8
DEFAULT_PVAL_THRESHOLDS = [5e-8, 1e-5]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_alid(alid: str) -> tuple[str, int, str, str]:
    """Parse 'CHROM:POS_A1_A2' → (chrom, pos, a1, a2).

    Normalises to canonical orientation (a1 < a2 alphabetically) on the fly.
    Returns alleles in canonical order; the caller is responsible for flipping
    EAF when normalisation was applied.
    """
    colon_idx = alid.index(":")
    chrom = alid[:colon_idx]
    rest = alid[colon_idx + 1:]
    parts = rest.split("_", 2)
    if len(parts) != 3:
        raise ValueError(f"Cannot parse ALID: {alid!r}  (expected CHROM:POS_A1_A2)")
    pos = int(parts[0])
    raw_a1, raw_a2 = parts[1], parts[2]
    if raw_a1 <= raw_a2:
        return chrom, pos, raw_a1, raw_a2
    else:
        return chrom, pos, raw_a2, raw_a1  # caller must flip EAF


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
        ("a1",    "U512"),
        ("a2",    "U512"),
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

            chrom, pos, a1, a2 = _parse_alid(raw_alid)
            canonical_alid = f"{chrom}:{pos}_{a1}_{a2}"
            # If alleles were flipped during normalisation, flip EAF too
            was_flipped = (raw_alid != canonical_alid and
                           raw_alid.split(":", 1)[1].split("_", 1)[1] != f"{a1}_{a2}")
            eaf = (1.0 - raw_eaf) if (was_flipped and not np.isnan(raw_eaf)) else raw_eaf

            rows.append((canonical_alid, chrom, pos, a1, a2))
            eaf_list.append(eaf)

    return np.array(rows, dtype=dt), np.array(eaf_list, dtype=np.float32)


def _load_trait_list(
    trait_tsv: str | Path,
) -> list[tuple[str, str, str, str | None]]:
    """
    Returns list of (trait_id, vcf_path, trait_name, vcf_build_or_None).

    TSV columns (tab-separated, no header):
        1  trait_id    required
        2  trait_name  optional (empty string if absent)
        3  vcf_path    required
        4  build       optional (hg19 / hg38 / GRCh37 / GRCh38)
    """
    pairs = []
    with open(trait_tsv) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            trait_name = parts[1].strip() if len(parts) > 1 else ""
            vcf_path = parts[2].strip() if len(parts) > 2 else ""
            vcf_build = parts[3].strip() if len(parts) > 3 and parts[3].strip() else None
            pairs.append((parts[0], vcf_path, trait_name, vcf_build))
    return pairs


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
    trait_pairs = _load_trait_list(trait_tsv)
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
        [(tid, name) for tid, _, name, _ in trait_pairs], dtype=trait_dt
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
            (j, trait_id, vcf_path, *_lookup_for(vcf_build), V)
            for j, (trait_id, vcf_path, _, vcf_build) in enumerate(batch_traits)
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
