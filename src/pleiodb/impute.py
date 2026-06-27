"""LD-based z-score imputation using elastic net on LD eigenvectors.

Ported from genotype-phenotype-map/pipeline_steps/imputation_method.R.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import scipy.linalg
from sklearn.linear_model import ElasticNetCV

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Issue 033 – pure-function math kernel (no filesystem I/O)
# ---------------------------------------------------------------------------

def _ld_pca(ld_matrix: np.ndarray, thresh: float = 0.9) -> tuple[np.ndarray, np.ndarray]:
    """Eigendecomposition of an LD matrix, components in descending variance order.

    Returns (eigenvalues, eigenvectors) truncated to the minimum number of
    components whose cumulative variance fraction reaches *thresh*.
    """
    vals, vecs = scipy.linalg.eigh(ld_matrix)
    # eigh returns ascending order; flip to descending
    vals = vals[::-1]
    vecs = vecs[:, ::-1]
    # Force non-negative (numerical noise can produce tiny negatives)
    vals = np.maximum(vals, 0.0)
    cumvar = np.cumsum(vals) / (vals.sum() or 1.0)
    n_comp = int(np.searchsorted(cumvar, thresh)) + 1
    n_comp = min(n_comp, len(vals))
    return vals[:n_comp], vecs[:, :n_comp]


def _elastic_net_impute(
    z: np.ndarray,
    eigenvectors: np.ndarray,
    n_comp: int,
) -> np.ndarray | None:
    """Predict z-scores for all positions (observed + missing) via elastic net.

    Fits on observed (non-NaN) positions using the first *n_comp* eigenvector
    columns as features; predicts for every position.

    Returns the prediction array, or None if fitting is not possible.
    """
    mask_obs = np.isfinite(z)
    if mask_obs.sum() < 2:
        return None

    E = eigenvectors[:, :n_comp]
    E_obs = E[mask_obs]
    z_obs = z[mask_obs]

    if len(np.unique(z_obs)) <= 1:
        return None

    n_obs = int(mask_obs.sum())
    cv = min(5, max(2, n_obs - 1))
    try:
        model = ElasticNetCV(l1_ratio=0.5, fit_intercept=False, cv=cv, max_iter=2000)
        model.fit(E_obs, z_obs)
        return model.predict(E).astype(np.float64)
    except Exception as exc:  # noqa: BLE001
        log.debug("elastic net failed: %s", exc)
        return None


def _poly_rescale(
    truth: np.ndarray,
    predicted: np.ndarray,
    npoly: int = 3,
) -> tuple[np.ndarray, float]:
    """Rescale *predicted* to match the scale of *truth* via polynomial regression.

    Removes Cook's-distance outliers (two-pass SD threshold), fits a degree-*npoly*
    polynomial through the origin, evaluates it on all of *predicted*.

    Returns (rescaled_array, pearson_r_between_rescaled_and_truth).
    r is computed only on the non-outlier observed positions.
    """
    obs_mask = np.isfinite(truth)
    if obs_mask.sum() < max(npoly + 1, 2):
        return predicted.copy(), float("nan")

    t_obs = truth[obs_mask]
    p_obs = predicted[obs_mask]

    # Two-pass outlier removal on ratio truth/predicted
    keep = _ratio_outlier_mask(t_obs, p_obs)
    if keep.sum() < max(npoly + 1, 2):
        return predicted.copy(), float("nan")

    t_clean = t_obs[keep]
    p_clean = p_obs[keep]

    # Reduce polynomial degree if not enough clean points remain
    while npoly > 1 and len(t_clean) < npoly + 1:
        npoly -= 1
    if len(t_clean) < 2:
        return predicted.copy(), float("nan")

    try:
        coeffs = np.polyfit(p_clean, t_clean, npoly)
    except np.linalg.LinAlgError:
        return predicted.copy(), float("nan")

    adjusted = np.polyval(coeffs, predicted)

    # Correlation on non-outlier observed positions
    adj_obs = np.polyval(coeffs, p_obs[keep])
    if len(adj_obs) < 2 or np.std(adj_obs) == 0 or np.std(t_clean) == 0:
        corr = float("nan")
    else:
        corr = float(np.corrcoef(adj_obs, t_clean)[0, 1])

    return adjusted, corr


def _ratio_outlier_mask(truth: np.ndarray, predicted: np.ndarray) -> np.ndarray:
    """Boolean keep-mask: True = not an outlier.

    Outliers are positions where truth/predicted deviates > 3 SD from the
    median in two consecutive passes (matching the R implementation).
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(predicted != 0, truth / predicted, np.nan)

    keep = np.ones(len(ratio), dtype=bool)
    for _ in range(2):
        r = ratio[keep]
        if len(r) < 2:
            break
        m = np.nanmedian(r)
        s = np.nanstd(r)
        if s == 0:
            break
        keep[keep] = np.abs(r - m) <= 3 * s

    return keep


def _se_outliers(
    se_obs: np.ndarray,
    se_hat: np.ndarray,
    outthresh: float = 3.0,
) -> np.ndarray:
    """Return a boolean mask of SE outliers detected via Cook's-distance analogue.

    Fits OLS se_obs ~ se_hat, computes Cook's distances, then flags positions
    whose Cook's distance is > *outthresh* SDs above the median (two passes).
    Only positions where both se_obs and se_hat are finite are evaluated;
    others are always False.
    """
    valid = np.isfinite(se_obs) & np.isfinite(se_hat)
    outliers = np.zeros(len(se_obs), dtype=bool)

    if valid.sum() < 4:
        return outliers

    x = se_hat[valid]
    y = se_obs[valid]

    # OLS: y = a*x + b
    X = np.column_stack([x, np.ones_like(x)])
    try:
        coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    except np.linalg.LinAlgError:
        return outliers

    y_hat = X @ coeffs
    residuals = y - y_hat
    n, p = len(y), 2
    mse = np.sum(residuals ** 2) / max(n - p, 1)

    # Hat matrix diagonal (leverage)
    try:
        XtX_inv = np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError:
        return outliers
    h = np.einsum("ij,jk,ik->i", X, XtX_inv, X)

    denom = p * mse * (1 - h) ** 2
    cook_d = np.where(denom > 0, (residuals ** 2 * h) / denom, 0.0)

    # Two-pass SD threshold on Cook's distances
    keep = np.ones(n, dtype=bool)
    for _ in range(2):
        d = cook_d[keep]
        if len(d) < 2:
            break
        m = np.median(d)
        s = np.std(d)
        if s == 0:
            break
        keep[keep] = np.abs(d - m) <= outthresh * s

    flag = ~keep  # outlier positions among valid entries
    valid_indices = np.where(valid)[0]
    outliers[valid_indices[flag]] = True
    return outliers


# ---------------------------------------------------------------------------
# Issue 034 – LD block index (filesystem scanning + variant matching)
# ---------------------------------------------------------------------------

def build_block_index(
    variants: np.ndarray,
    ld_dir: Path,
    ancestry: str = "EUR",
) -> dict[str, dict]:
    """Map pleiodb variants to LD reference panel blocks.

    Supports two LD panel layouts:

    * **Flat** (production): ``ld_dir/ancestry/{chr}/{block}.tsv`` and
      ``{block}.unphased.vcor1.gz`` as siblings in the chromosome directory.
    * **Nested** (test fixtures): ``ld_dir/ancestry/{chr}/{block}/{block}.tsv``
      and ``{block}.unphased.vcor1.gz`` inside a per-block subdirectory.

    Returns a dict keyed by ``"{chrom}/{block_name}"`` string:
    ``{'ld_path': Path, 'variant_indices': list[int], 'ld_row_indices': list[int],
    'n_ld_snps': int}``

    Only blocks with ≥2 matched variants are included.
    """
    import pandas as pd  # local import – not always installed in minimal envs

    panel_dir = Path(ld_dir) / ancestry
    if not panel_dir.is_dir():
        log.warning("LD panel directory not found: %s", panel_dir)
        return {}

    # Build a lookup from ALID → pleiodb row index
    alid_to_idx: dict[str, int] = {
        str(variants["id"][i]): i for i in range(len(variants))
    }

    index: dict[str, dict] = {}

    for chrom_dir in sorted(panel_dir.iterdir()):
        if not chrom_dir.is_dir():
            continue
        chrom = chrom_dir.name

        # Collect (block_name, tsv_path, ld_path) pairs from both layouts.
        # Use a dict so that if both layouts produce the same block_name, the
        # flat layout takes precedence (it's the production layout).
        blocks: dict[str, tuple[Path, Path]] = {}

        # Nested layout: block subdirectories
        for item in chrom_dir.iterdir():
            if not item.is_dir():
                continue
            block_name = item.name
            tsv_p = item / f"{block_name}.tsv"
            ld_p = item / f"{block_name}.unphased.vcor1.gz"
            if tsv_p.exists() and ld_p.exists():
                blocks[block_name] = (tsv_p, ld_p)

        # Flat layout: TSV files directly in chrom_dir (overrides nested)
        for tsv_p in chrom_dir.glob("*.tsv"):
            block_name = tsv_p.stem
            ld_p = chrom_dir / f"{block_name}.unphased.vcor1.gz"
            if ld_p.exists():
                blocks[block_name] = (tsv_p, ld_p)

        for block_name, (tsv_path, ld_path) in sorted(blocks.items()):
            try:
                tsv = pd.read_csv(tsv_path, sep="\t", usecols=["SNP"])
            except Exception as exc:  # noqa: BLE001
                log.warning("Could not read %s: %s – skipping", tsv_path, exc)
                continue

            snp_ids: list[str] = tsv["SNP"].tolist()
            variant_indices: list[int] = []
            ld_row_indices: list[int] = []

            for row_i, snp_id in enumerate(snp_ids):
                # LD panel SNP IDs may have a "chr" prefix; strip for comparison.
                bare_snp = snp_id.lstrip("chr") if snp_id.startswith("chr") else snp_id
                if bare_snp in alid_to_idx:
                    variant_indices.append(alid_to_idx[bare_snp])
                    ld_row_indices.append(row_i)

            if len(variant_indices) >= 2:
                key = f"{chrom}/{block_name}"
                ldeig_path = tsv_path.with_suffix(".ldeig.rds")
                index[key] = {
                    "ld_path": ld_path,
                    "tsv_path": tsv_path,
                    "ldeig_path": ldeig_path if ldeig_path.exists() else None,
                    "variant_indices": variant_indices,
                    "ld_row_indices": ld_row_indices,
                    "n_ld_snps": len(snp_ids),
                }

    log.info(
        "LD block index built: %d blocks matched across %d variants",
        len(index),
        len(alid_to_idx),
    )
    return index


def build_block_allele_lookup(
    tsv_path: Path | str,
) -> dict[str, tuple[int, bool]]:
    """Build an allele-ID → (LD-row-index, is_flipped) lookup for a single LD block.

    Reads the block's TSV file (the same file used by build_block_index) and
    returns a dict with *both* allele orientations for every variant, so that
    z-scores from read_vcf_region can be matched regardless of REF/ALT order.

    Parameters
    ----------
    tsv_path : path to the block's ``.tsv`` file in the LD reference panel

    Returns
    -------
    dict keyed by ``"{bare_chrom}:{pos}_{a}_{b}"`` mapping to
    ``(ld_row_index, is_flipped)``.  Two entries are added per variant:
    - canonical orientation (A1 ≤ A2) → ``(i, False)``
    - reversed orientation               → ``(i, True)``

    When ``is_flipped=True``, callers must negate the z-score before using it
    as a training value (the VCF gave the effect for the non-canonical allele).
    """
    import pandas as pd
    from .alid import parse_alid

    tsv_path = Path(tsv_path)
    try:
        tsv = pd.read_csv(tsv_path, sep="\t", usecols=["SNP"])
    except Exception as exc:  # noqa: BLE001
        log.warning("Could not read block TSV %s: %s", tsv_path, exc)
        return {}

    lookup: dict[str, tuple[int, bool]] = {}
    for row_i, snp_id in enumerate(tsv["SNP"]):
        bare_snp = snp_id.lstrip("chr") if isinstance(snp_id, str) and snp_id.startswith("chr") else str(snp_id)
        try:
            chrom, pos, a1, a2 = parse_alid(bare_snp)
        except ValueError:
            continue
        # Canonical orientation (A1 ≤ A2 alphabetically in the SNP ID)
        lookup[bare_snp] = (row_i, False)
        # Reversed orientation
        lookup[f"{chrom}:{pos}_{a2}_{a1}"] = (row_i, True)

    return lookup


def _load_ld_submatrix(ld_path: Path, ld_row_indices: list[int]) -> np.ndarray:
    """Load the full LD matrix from *ld_path* and return the submatrix.

    *ld_path* points directly to the ``.unphased.vcor1.gz`` file.
    The submatrix rows and columns correspond to *ld_row_indices*.
    """
    import pandas as pd

    full = pd.read_csv(ld_path, header=None, delimiter="\t").values.astype(np.float64)
    idx = np.array(ld_row_indices)
    return full[np.ix_(idx, idx)]


def _load_ld_full(ld_path: Path) -> np.ndarray:
    """Load the complete N_ld_snps × N_ld_snps LD correlation matrix."""
    import pandas as pd
    return pd.read_csv(ld_path, header=None, delimiter="\t").values.astype(np.float64)


_LDEIG_MAX_K = 250  # eigenvectors stored in the npz cache (covers ≥97% variance)

_RSCRIPT_EXTRACT = """
args <- commandArgs(trailingOnly=TRUE)
rds_path <- args[1]; out_path <- args[2]; max_k <- as.integer(args[3])
x <- readRDS(rds_path)
k <- min(max_k, length(x$values))
vals <- x$values                              # all eigenvalues
vecs <- x$vectors[, seq_len(k), drop=FALSE]  # first k eigenvectors
con <- file(out_path, "wb")
writeBin(as.integer(length(vals)), con, size=4L)
writeBin(as.integer(k),            con, size=4L)
writeBin(as.double(vals),          con, size=8L)
writeBin(as.double(as.vector(vecs)), con, size=8L)
close(con)
"""


def _load_block_eigenvectors(
    ldeig_path: Path,
    thresh: float = 0.9,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (eigenvalues[:k], eigenvectors[:, :k]) from a ``.ldeig.rds`` file.

    Caches the result alongside the source file as a ``.ldeig.npz`` binary
    (written atomically; safe under concurrent workers).  The cache stores
    all eigenvalues and the first :data:`_LDEIG_MAX_K` eigenvectors.

    On first access the R script is called via Rscript (~5 s).  Cached reads
    take < 0.2 s regardless of block size.

    Falls back to loading the ``.unphased.vcor1.gz`` LD matrix and running PCA
    if ``Rscript`` is unavailable.
    """
    import subprocess
    import tempfile

    npz_path = ldeig_path.with_suffix(".npz")

    def _load_npz() -> tuple[np.ndarray, np.ndarray]:
        data = np.load(str(npz_path))
        vals: np.ndarray = data["values"].astype(np.float64)
        vecs: np.ndarray = data["vectors"].astype(np.float64)
        total = float(np.maximum(vals, 0).sum()) or 1.0
        cumvar = np.cumsum(np.maximum(vals, 0)) / total
        k = int(np.searchsorted(cumvar, thresh)) + 1
        k = min(k, vecs.shape[1])
        return vals[:k], vecs[:, :k]

    if npz_path.exists():
        try:
            return _load_npz()
        except Exception:
            npz_path.unlink(missing_ok=True)

    # Cache miss: call Rscript to extract, save, reload
    fd, tmp_rscript = tempfile.mkstemp(suffix=".R")
    fd2, tmp_bin = tempfile.mkstemp(suffix=".bin")
    try:
        import os
        os.write(fd, _RSCRIPT_EXTRACT.encode())
        os.close(fd)
        os.close(fd2)

        result = subprocess.run(
            ["Rscript", tmp_rscript, str(ldeig_path), tmp_bin, str(_LDEIG_MAX_K)],
            capture_output=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"Rscript failed: {result.stderr.decode()[:200]}")

        with open(tmp_bin, "rb") as fh:
            n = int(np.frombuffer(fh.read(4), dtype=np.int32)[0])
            k_stored = int(np.frombuffer(fh.read(4), dtype=np.int32)[0])
            vals = np.frombuffer(fh.read(n * 8), dtype=np.float64)
            vecs = np.frombuffer(fh.read(n * k_stored * 8), dtype=np.float64).reshape(
                n, k_stored, order="F"
            )

        # Atomic write: temp → rename so concurrent workers don't corrupt the cache.
        # np.savez_compressed appends .npz to the given path, so use a .tmp base.
        fd3, tmp_base = tempfile.mkstemp(suffix=".tmp", dir=npz_path.parent)
        os.close(fd3)
        tmp_npz = Path(tmp_base + ".npz")
        try:
            np.savez_compressed(tmp_base, values=vals, vectors=vecs)
            tmp_npz.replace(npz_path)
        except Exception:
            tmp_npz.unlink(missing_ok=True)
            raise
        finally:
            Path(tmp_base).unlink(missing_ok=True)

    except (FileNotFoundError, RuntimeError) as exc:
        log.warning("_load_block_eigenvectors: Rscript unavailable (%s); using LD matrix PCA", exc)
        return None, None  # caller falls back to ld_matrix PCA
    finally:
        Path(tmp_rscript).unlink(missing_ok=True)
        Path(tmp_bin).unlink(missing_ok=True)

    return _load_npz()


def _tsv_region(tsv_path: Path) -> tuple[str, int, int]:
    """Return (bare_chrom, start_bp, end_bp) from a block's LD panel TSV file."""
    import pandas as pd
    df = pd.read_csv(tsv_path, sep="\t", usecols=["CHR", "BP"])
    chrom = str(df["CHR"].iloc[0])
    return chrom, int(df["BP"].min()), int(df["BP"].max())


def _build_dense_z(
    vcf_result: dict[str, float],
    allele_lookup: dict[str, tuple[int, bool]],
    n_ld_snps: int,
) -> np.ndarray:
    """Map VCF z-scores to LD panel row indices using the allele lookup.

    Returns a float64 array of length *n_ld_snps*.  Positions not found in
    the VCF or not in the allele lookup remain NaN.  When ``is_flipped=True``
    the z-score sign is negated to align with the canonical ALID orientation
    (effect on A2).
    """
    z_dense = np.full(n_ld_snps, np.nan, dtype=np.float64)
    for allele_id, z in vcf_result.items():
        entry = allele_lookup.get(allele_id)
        if entry is None:
            continue
        row_idx, is_flipped = entry
        z_dense[row_idx] = -z if is_flipped else z
    return z_dense


# ---------------------------------------------------------------------------
# Issue 035 – public impute_z_block API
# ---------------------------------------------------------------------------

def _impute_block_process(args: tuple) -> tuple[list, int, int, int]:
    """Picklable worker for ProcessPoolExecutor.

    Receives a pre-extracted z_mini (n_block_variants × B) so the full V×T
    matrix never crosses the process boundary.  Returns fills as a list of
    (missing_bool_mask, j, fill_z_float32) where the mask indexes into z_mini
    rows (i.e. into the block's variant set).

    When *block_vcf_paths* is non-empty the worker switches to **dense mode**:
    it reads VCF regions for traits with missing data (using a per-block
    ThreadPoolExecutor), assembles a dense z-score vector covering all LD panel
    variants, and trains the elastic net on that full set.  This improves
    imputation quality because training uses ~7 k variants instead of the
    ~58 pleiodb variants that fall in a typical block.

    Sparse mode (block_vcf_paths empty or tsv_path unavailable) is identical
    to the original behaviour.
    """
    (
        block_key, ld_path_str, li, z_mini, block_eaf,
        thresh, min_cor, tsv_path_str, ldeig_path_str, block_vcf_paths, vcf_threads,
    ) = args
    ld_path = Path(ld_path_str)
    li_arr = np.array(li, dtype=np.intp)
    B = z_mini.shape[1]

    fills: list = []
    n_imputed = 0
    n_skipped_size = 0
    n_skipped_cor = 0

    if (~np.isfinite(block_eaf)).any():
        log.warning(
            "%s: %d variant(s) have NaN EAF",
            block_key, int((~np.isfinite(block_eaf)).sum()),
        )

    use_dense = bool(block_vcf_paths and tsv_path_str)

    # ---- Load eigenvectors -------------------------------------------------
    # Dense mode: load precomputed eigenvectors from .ldeig.rds (via cache).
    # Sparse mode: load LD submatrix and compute PCA in-process.
    eigenvalues: np.ndarray | None = None
    eigenvectors: np.ndarray | None = None
    n_ld_snps: int | None = None

    if use_dense and ldeig_path_str:
        try:
            eigenvalues, eigenvectors = _load_block_eigenvectors(
                Path(ldeig_path_str), thresh
            )
            if eigenvectors is not None:
                n_ld_snps = eigenvectors.shape[0]
        except Exception as exc:  # noqa: BLE001
            log.warning("%s: eigenvector load failed (%s) – falling back to LD matrix", block_key, exc)

    if eigenvectors is None:
        # Fall back: load LD matrix and compute PCA
        try:
            ld_mat = _load_ld_full(ld_path) if use_dense else _load_ld_submatrix(ld_path, li)
        except Exception as exc:  # noqa: BLE001
            log.warning("Could not load LD matrix for %s: %s – skipping", block_key, exc)
            return fills, 0, 1, 0
        try:
            eigenvalues, eigenvectors = _ld_pca(ld_mat, thresh)
        except Exception as exc:  # noqa: BLE001
            log.debug("PCA failed for %s: %s – skipping", block_key, exc)
            return fills, 0, 1, 0
        n_ld_snps = ld_mat.shape[0] if use_dense else None
        del ld_mat

    n_comp = len(eigenvalues)

    # ---- Dense mode: read VCF regions in parallel --------------------------
    vcf_results: dict[int, dict[str, float]] = {}
    allele_lookup: dict[str, tuple[int, bool]] = {}

    if use_dense:
        from concurrent.futures import ThreadPoolExecutor, as_completed as _as_completed
        from .vcf import read_vcf_region

        allele_lookup = build_block_allele_lookup(Path(tsv_path_str))
        try:
            chrom, region_start, region_end = _tsv_region(Path(tsv_path_str))
        except Exception as exc:  # noqa: BLE001
            log.warning("%s: could not read block region from TSV (%s) – falling back to sparse", block_key, exc)
            use_dense = False
        else:
            def _read_one(jvb: tuple) -> tuple[int, dict]:
                j, vcf_path, vcf_build = jvb
                try:
                    return j, read_vcf_region(vcf_path, chrom, region_start, region_end, vcf_build=vcf_build)
                except Exception as exc2:  # noqa: BLE001
                    log.debug("%s trait-col %d: VCF read failed (%s)", block_key, j, exc2)
                    return j, {}

            with ThreadPoolExecutor(max_workers=max(1, vcf_threads)) as pool:
                futs = {pool.submit(_read_one, jvb): jvb[0] for jvb in block_vcf_paths}
                for fut in _as_completed(futs):
                    j_res, vres = fut.result()
                    if vres:
                        vcf_results[j_res] = vres

    # ---- Per-trait imputation ----------------------------------------------
    for j in range(B):
        z_col = z_mini[:, j].astype(np.float64)
        missing = ~np.isfinite(z_col)

        if missing.sum() == 0:
            continue

        if use_dense and j in vcf_results:
            # Dense mode: build full-LD z vector, fit on it, extract at pleiodb positions
            z_dense = _build_dense_z(vcf_results[j], allele_lookup, n_ld_snps)
            z_pred_full = _elastic_net_impute(z_dense, eigenvectors, n_comp)
            if z_pred_full is None:
                n_skipped_size += 1
                continue
            z_pred = z_pred_full[li_arr]
        else:
            # Sparse mode: fit on pleiodb positions only
            if (~missing).sum() < 2:
                n_skipped_size += 1
                continue
            evecs = eigenvectors[li_arr, :] if use_dense else eigenvectors
            z_pred = _elastic_net_impute(z_col, evecs, n_comp)
            if z_pred is None:
                n_skipped_size += 1
                continue

        z_adj, corr = _poly_rescale(z_col, z_pred, npoly=3)

        if not np.isfinite(corr) or corr < min_cor:
            log.debug(
                "%s trait-col %d: correlation %.3f < %.3f – skipping",
                block_key, j, corr if np.isfinite(corr) else float("nan"), min_cor,
            )
            n_skipped_cor += 1
            continue

        fills.append((missing, j, z_adj[missing].astype(np.float32)))
        n_imputed += int(missing.sum())

    return fills, n_imputed, n_skipped_size, n_skipped_cor


def impute_z_block(
    z_block: np.ndarray,
    variants: np.ndarray,
    eaf_arr: np.ndarray,
    block_index: dict[str, dict],
    thresh: float = 0.9,
    min_cor: float = 0.7,
    out_mask: np.ndarray | None = None,
    workers: int = 1,
    vcf_paths: "list[tuple[str, str | None]] | None" = None,
    vcf_threads: int = 8,
) -> None:
    """Impute missing z-scores in *z_block* using LD-based elastic-net regression.

    Modifies *z_block* (shape V×B, float32) in-place: NaN cells that can be
    imputed with block-level correlation ≥ *min_cor* are filled with rescaled
    predictions. Sets the corresponding cells of *out_mask* to True when
    supplied.

    Parameters
    ----------
    z_block : float32 array (V, B)
        z-score matrix for a trait batch; NaN = missing.
    variants : structured array with fields 'id', 'chrom', 'pos', 'a1', 'a2'
    eaf_arr : float32 array (V,)
        Effect allele frequencies for all variants; used for SE prediction.
    block_index : dict from build_block_index()
    thresh : float
        Fraction of LD variance retained in the PCA step.
    min_cor : float
        Minimum block-level Pearson correlation between imputed and observed
        z-scores to accept imputed values for that block.
    out_mask : optional bool array (V, B)
        Positions that were imputed are set True.
    workers : int
        Number of parallel worker processes for LD block processing.
    vcf_paths : list of (vcf_path, vcf_build) per trait column, or None.
        When provided, each block worker reads VCF regions for traits with
        missing data and trains the elastic net on the dense z-score set
        (all LD panel variants in the block).  None = sparse mode only.
    vcf_threads : int
        Inner thread count for parallel VCF region reads within each block
        worker process.
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed

    total_imputed = 0
    blocks_skipped_size = 0
    blocks_skipped_cor = 0
    # each entry: (vi_global, missing_local_mask, j, fill_z)
    all_fills: list = []

    items = list(block_index.items())

    def _submit_one(block_key: str, info: dict):
        vi = np.array(info["variant_indices"])
        z_mini = z_block[vi, :].copy()
        block_eaf = eaf_arr[vi]
        tsv_path_str = str(info["tsv_path"]) if info.get("tsv_path") else ""
        ldeig = info.get("ldeig_path")
        ldeig_path_str = str(ldeig) if ldeig is not None else ""

        # Traits with at least one NaN in this block and a valid VCF path
        block_vcf: list[tuple[int, str, "str | None"]] = []
        if vcf_paths is not None and tsv_path_str:
            for j, (vp, vb) in enumerate(vcf_paths):
                if vp and not np.isfinite(z_mini[:, j]).all():
                    block_vcf.append((j, vp, vb))

        args = (
            block_key,
            str(info["ld_path"]),
            info["ld_row_indices"],
            z_mini,
            block_eaf,
            thresh,
            min_cor,
            tsv_path_str,
            ldeig_path_str,
            block_vcf,
            vcf_threads,
        )
        return vi, args

    if workers <= 1:
        for block_key, info in items:
            vi, args = _submit_one(block_key, info)
            fills, n_imp, n_ss, n_sc = _impute_block_process(args)
            all_fills.extend((vi, m, j, fz) for m, j, fz in fills)
            total_imputed += n_imp
            blocks_skipped_size += n_ss
            blocks_skipped_cor += n_sc
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {}
            for block_key, info in items:
                vi, args = _submit_one(block_key, info)
                fut = pool.submit(_impute_block_process, args)
                futures[fut] = vi
            for fut in as_completed(futures):
                vi = futures[fut]
                fills, n_imp, n_ss, n_sc = fut.result()
                all_fills.extend((vi, m, j, fz) for m, j, fz in fills)
                total_imputed += n_imp
                blocks_skipped_size += n_ss
                blocks_skipped_cor += n_sc

    # Apply all fills after parallel phase (safe: non-overlapping rows per block)
    for vi, missing_mask, j, fill_z in all_fills:
        z_block[vi[missing_mask], j] = fill_z
        if out_mask is not None:
            out_mask[vi[missing_mask], j] = True

    log.info(
        "Imputation complete: %d cells filled; %d blocks skipped (size/fit), "
        "%d block×trait combos skipped (low correlation)",
        total_imputed, blocks_skipped_size, blocks_skipped_cor,
    )
