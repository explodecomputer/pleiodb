"""
Command-line interface for pleiodb.

Commands:
  pleiodb build   — create database from GWAS-VCF files
  pleiodb lambda  — compute/add sample-overlap matrix
  pleiodb query   — query z-scores, beta, SE
  pleiodb info    — show database summary
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

import click
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _open_db(db_path: str):
    from .db import GWASDatabase
    return GWASDatabase.open(db_path)


def _print_table(header: list[str], rows: list[list], fmt: str, fh=None):
    if fh is None:
        fh = sys.stdout
    if fmt == "tsv":
        print("\t".join(header), file=fh)
        for row in rows:
            print("\t".join(str(x) for x in row), file=fh)
    elif fmt == "json":
        out = [dict(zip(header, row)) for row in rows]
        print(json.dumps(out, allow_nan=True), file=fh)
    else:
        raise ValueError(f"Unknown format: {fmt}")


# ---------------------------------------------------------------------------
# CLI root
# ---------------------------------------------------------------------------

@click.group()
def main():
    """pleiodb — compact GWAS summary-statistic storage and query tool."""


# ---------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------

@main.command()
@click.argument("output_dir")
@click.option("--variants", "-v", required=True,
              help="TSV: ALID  EAF  (ALID = CHROM:POS_A1_A2, alphabetical alleles)")
@click.option("--traits", "-t", required=True,
              help="TSV: trait_id  trait_name  vcf_path  [build]")
@click.option("--chunk-v", default=512, show_default=True,
              help="Chunk size along variant axis")
@click.option("--chunk-t", default=512, show_default=True,
              help="Chunk size along trait axis")
@click.option("--workers", "-j", default=8, show_default=True,
              help="Parallel VCF readers")
@click.option("--pval", multiple=True, default=["5e-8", "1e-5"],
              show_default=True, help="P-value threshold(s) for significance masks")
@click.option("--overwrite", is_flag=True, default=False)
@click.option(
    "--variants-build", default=None,
    help="Genome build of the variant list (hg19/hg38/GRCh37/GRCh38). "
         "When set, VCF files whose traits TSV specifies a different build "
         "will have their coordinates lifted over automatically.",
)
def build(output_dir, variants, traits, chunk_v, chunk_t, workers, pval, overwrite,
          variants_build):
    """Build a pleiodb database from GWAS-VCF files."""
    from .build import build_database
    thresholds = [float(p) for p in pval]
    build_database(
        output_dir=output_dir,
        variants_path=variants,
        trait_tsv=traits,
        chunk_shape=(chunk_v, chunk_t),
        workers=workers,
        pval_thresholds=thresholds,
        overwrite=overwrite,
        variants_build=variants_build,
    )
    click.echo(f"Database written to {output_dir}")


# ---------------------------------------------------------------------------
# lambda
# ---------------------------------------------------------------------------

@main.command("lambda")
@click.argument("db_path")
@click.option("--workers", "-j", default=8, show_default=True)
@click.option("--null-thresh", default=3.0, show_default=True,
              help="Max |z| to consider a variant 'null' for correlation")
def compute_lambda(db_path, workers, null_thresh):
    """Compute sample-overlap (lambda) matrix and add to database."""
    from .build import build_lambda
    build_lambda(db_path, z_null_thresh=null_thresh, workers=workers)
    click.echo("Lambda matrix complete.")


# ---------------------------------------------------------------------------
# query
# ---------------------------------------------------------------------------

@main.command()
@click.argument("db_path")
@click.option("--variant", "-s", default=None,
              help="Single ALID → all traits (or just --trait if also given)")
@click.option("--trait", "-t", default=None,
              help="Single trait ID → all variants (or just --variant if also given)")
@click.option("--region", "-r", default=None,
              help="Genomic region chr:start-end → all traits")
@click.option("--variants-file", default=None,
              help="File with variant IDs (one per line)")
@click.option("--traits-file", default=None,
              help="File with trait IDs (one per line)")
@click.option("--pval", default=None, type=float,
              help="Return only associations with p ≤ this threshold")
@click.option("--output", "-o", default=None,
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "json"]), show_default=True)
@click.option("--study-scale", is_flag=True, default=False,
              help="Add beta_study and se_study columns (original phenotype scale). "
                   "beta_study = sqrt(var_y) × beta_norm.  For binary traits, "
                   "this recovers the log-OR scale.")
def query(db_path, variant, trait, region, variants_file, traits_file, pval, output,
          fmt, study_scale):
    """Query z-scores, normalised beta/SE, and p-values from the database.

    Every row in the output contains: alid, trait_id, z, beta_norm,
    se_norm, pval.  beta_norm and se_norm are on the var(y)=1 scale;
    pval = 2·Φ(−|Z|) (two-sided).

    With --study-scale, two additional columns beta_study and se_study are
    appended, scaled by sqrt(var_y[t]) to recover original phenotype units.
    """
    db = _open_db(db_path)
    fh = open(output, "w") if output else sys.stdout

    try:
        if pval is not None and variant is None and trait is None and region is None \
                and variants_file is None and traits_file is None:
            _query_pval(db, pval, fmt, fh, study_scale=study_scale)

        elif variant is not None and trait is not None:
            # Both specified → return the single intersection cell
            _query_block(db, [variant], [trait], pval, fmt, fh, study_scale=study_scale)

        elif variant is not None:
            _query_single_variant(db, variant, pval, fmt, fh, study_scale=study_scale)

        elif trait is not None:
            _query_single_trait(db, trait, pval, fmt, fh, study_scale=study_scale)

        elif region is not None:
            _query_region(db, region, pval, fmt, fh, study_scale=study_scale)

        elif variants_file or traits_file:
            v_ids = _read_id_file(variants_file) if variants_file else None
            t_ids = _read_id_file(traits_file) if traits_file else None
            _query_block(db, v_ids, t_ids, pval, fmt, fh, study_scale=study_scale)

        else:
            raise click.UsageError(
                "Provide at least one of: --variant, --trait, --region, "
                "--variants-file/--traits-file, --pval"
            )
    finally:
        if output:
            fh.close()


def _read_id_file(path: str) -> list[str]:
    return [l.strip() for l in open(path) if l.strip() and not l.startswith("#")]


def _variant_label(db, v_idx: int) -> str:
    v = db.variants
    return str(v["id"][v_idx])


def _trait_label(db, t_idx: int) -> str:
    return str(db.traits["id"][t_idx])


_QUERY_HEADER = ["alid", "trait_id", "z", "beta_norm", "se_norm", "pval"]
_STUDY_EXTRA = ["beta_study", "se_study"]


def _compute_pval(z_val: float) -> float:
    """Two-sided p-value: 2·Φ(−|Z|)."""
    from scipy.stats import norm
    return float(2 * norm.sf(abs(z_val)))


def _header(study_scale: bool) -> list[str]:
    return _QUERY_HEADER + (_STUDY_EXTRA if study_scale else [])


def _make_row(vid, tid, z_val, beta_val, se_val,
              study_scale: bool = False, var_y_t: float = np.nan) -> list:
    row = [
        vid, tid,
        round(z_val, 4),
        round(beta_val, 6),
        round(se_val, 6),
        _compute_pval(z_val),
    ]
    if study_scale:
        sqrt_vy = float(np.sqrt(var_y_t)) if np.isfinite(var_y_t) else np.nan
        row += [
            round(sqrt_vy * beta_val, 6) if np.isfinite(sqrt_vy) else float("nan"),
            round(sqrt_vy * se_val, 6)   if np.isfinite(sqrt_vy) else float("nan"),
        ]
    return row


def _query_single_variant(db, vid, pval, fmt, fh, study_scale: bool = False):
    v_idx = int(db.variant_index([vid])[0])
    z = db.zscore_variant(v_idx)
    beta_arr, se_arr = db.beta_se_variant(v_idx)
    var_y = db.var_y if study_scale else None
    rows = []
    for t_idx in range(db.T):
        z_val = float(z[t_idx])
        if np.isnan(z_val):
            continue
        if pval is not None and _compute_pval(z_val) > pval:
            continue
        vy = float(var_y[t_idx]) if study_scale else np.nan
        rows.append(_make_row(
            vid, _trait_label(db, t_idx), z_val,
            float(beta_arr[t_idx]), float(se_arr[t_idx]),
            study_scale=study_scale, var_y_t=vy,
        ))
    _print_table(_header(study_scale), rows, fmt, fh)


def _query_single_trait(db, tid, pval, fmt, fh, study_scale: bool = False):
    t_idx = int(db.trait_index([tid])[0])
    z = db.zscore_trait(t_idx)
    vy = float(db.var_y[t_idx]) if study_scale else np.nan
    rows = []
    for v_idx in range(db.V):
        z_val = float(z[v_idx])
        if np.isnan(z_val):
            continue
        if pval is not None and _compute_pval(z_val) > pval:
            continue
        beta_v, se_v = db.beta_se_block([v_idx], [t_idx])
        rows.append(_make_row(
            _variant_label(db, v_idx), tid, z_val,
            float(beta_v[0, 0]), float(se_v[0, 0]),
            study_scale=study_scale, var_y_t=vy,
        ))
    _print_table(_header(study_scale), rows, fmt, fh)


def _query_region(db, region_str, pval, fmt, fh, study_scale: bool = False):
    chrom, rest = region_str.split(":")
    start, end = (int(x) for x in rest.split("-"))
    v_idx, z_mat = db.zscore_region(chrom, start, end)
    var_y = db.var_y if study_scale else None
    rows = []
    for i, vi in enumerate(v_idx):
        vid = _variant_label(db, int(vi))
        for t_idx in range(db.T):
            z_val = float(z_mat[i, t_idx])
            if np.isnan(z_val):
                continue
            if pval is not None and _compute_pval(z_val) > pval:
                continue
            beta_v, se_v = db.beta_se_block([int(vi)], [t_idx])
            vy = float(var_y[t_idx]) if study_scale else np.nan
            rows.append(_make_row(
                vid, _trait_label(db, t_idx), z_val,
                float(beta_v[0, 0]), float(se_v[0, 0]),
                study_scale=study_scale, var_y_t=vy,
            ))
    _print_table(_header(study_scale), rows, fmt, fh)


def _query_block(db, v_ids, t_ids, pval, fmt, fh, study_scale: bool = False):
    v_idx = np.arange(db.V, dtype=np.int64) if v_ids is None else db.variant_index(v_ids)
    t_idx = np.arange(db.T, dtype=np.int64) if t_ids is None else db.trait_index(t_ids)
    z_mat = db.zscore_block(v_idx, t_idx)
    var_y = db.var_y if study_scale else None
    rows = []
    for i, vi in enumerate(v_idx):
        vid = _variant_label(db, int(vi))
        for j, ti in enumerate(t_idx):
            z_val = float(z_mat[i, j])
            if np.isnan(z_val):
                continue
            if pval is not None and _compute_pval(z_val) > pval:
                continue
            beta_v, se_v = db.beta_se_block([int(vi)], [int(ti)])
            vy = float(var_y[ti]) if study_scale else np.nan
            rows.append(_make_row(
                vid, _trait_label(db, int(ti)), z_val,
                float(beta_v[0, 0]), float(se_v[0, 0]),
                study_scale=study_scale, var_y_t=vy,
            ))
    _print_table(_header(study_scale), rows, fmt, fh)


def _query_pval(db, pval, fmt, fh, study_scale: bool = False):
    v_idx, t_idx, z_vals = db.query_significant(pval=pval)
    var_y = db.var_y if study_scale else None
    rows = []
    for vi, ti, zv in zip(v_idx, t_idx, z_vals):
        z_val = float(zv)
        beta_v, se_v = db.beta_se_block([int(vi)], [int(ti)])
        vy = float(var_y[ti]) if study_scale else np.nan
        rows.append(_make_row(
            _variant_label(db, int(vi)), _trait_label(db, int(ti)), z_val,
            float(beta_v[0, 0]), float(se_v[0, 0]),
            study_scale=study_scale, var_y_t=vy,
        ))
    _print_table(_header(study_scale), rows, fmt, fh)


# ---------------------------------------------------------------------------
# info
# ---------------------------------------------------------------------------

@main.command()
@click.argument("db_path")
def info(db_path):
    """Display database metadata and storage summary."""
    db = _open_db(db_path)
    d = db.info()
    click.echo(json.dumps(d, indent=2))


if __name__ == "__main__":
    main()
