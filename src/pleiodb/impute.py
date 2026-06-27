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
                index[key] = {
                    "ld_path": ld_path,
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


def _load_ld_submatrix(ld_path: Path, ld_row_indices: list[int]) -> np.ndarray:
    """Load the full LD matrix from *ld_path* and return the submatrix.

    *ld_path* points directly to the ``.unphased.vcor1.gz`` file.
    The submatrix rows and columns correspond to *ld_row_indices*.
    """
    import pandas as pd

    full = pd.read_csv(ld_path, header=None, delimiter="\t").values.astype(np.float64)
    idx = np.array(ld_row_indices)
    return full[np.ix_(idx, idx)]


# ---------------------------------------------------------------------------
# Issue 035 – public impute_z_block API
# ---------------------------------------------------------------------------

def impute_z_block(
    z_block: np.ndarray,
    variants: np.ndarray,
    eaf_arr: np.ndarray,
    block_index: dict[Path, dict],
    thresh: float = 0.9,
    min_cor: float = 0.7,
    out_mask: np.ndarray | None = None,
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
    """
    V, B = z_block.shape
    total_imputed = 0
    blocks_skipped_size = 0
    blocks_skipped_cor = 0

    for block_key, info in block_index.items():
        vi = np.array(info["variant_indices"])   # pleiodb row indices
        li = info["ld_row_indices"]              # rows in LD matrix
        ld_path = info["ld_path"]

        # EAF for this block's variants
        block_eaf = eaf_arr[vi]
        nan_eaf = ~np.isfinite(block_eaf)
        if nan_eaf.any():
            log.warning(
                "%s: %d variant(s) have NaN EAF; SE-outlier step skipped for them",
                block_key, nan_eaf.sum(),
            )

        try:
            ld_sub = _load_ld_submatrix(ld_path, li)
        except Exception as exc:  # noqa: BLE001
            log.warning("Could not load LD matrix for %s: %s – skipping", block_key, exc)
            blocks_skipped_size += 1
            continue

        try:
            eigenvalues, eigenvectors = _ld_pca(ld_sub, thresh)
        except Exception as exc:  # noqa: BLE001
            log.debug("PCA failed for %s: %s – skipping", block_key, exc)
            blocks_skipped_size += 1
            continue

        n_comp = len(eigenvalues)

        for j in range(B):
            z_col = z_block[vi, j].astype(np.float64)

            # SE-based outlier detection: treat unstable betas as missing
            with np.errstate(divide="ignore", invalid="ignore"):
                se_hat = np.where(
                    np.isfinite(block_eaf) & (block_eaf > 0) & (block_eaf < 1),
                    1.0 / np.sqrt(2.0 * block_eaf * (1.0 - block_eaf)),
                    np.nan,
                )
            se_obs = np.full(len(vi), np.nan)
            # We only have z-scores, not SE, in z_block; SE outlier step
            # is best-effort based on z-score magnitude vs theoretical SE.
            # Skip full SE outlier detection here; it requires SE values
            # which are not stored after ingestion.  The _se_outliers helper
            # is still available for callers that supply SE arrays.

            missing = ~np.isfinite(z_col)
            if missing.sum() == 0:
                continue  # nothing to impute for this block × trait

            if (~missing).sum() < 2:
                blocks_skipped_size += 1
                continue

            z_pred = _elastic_net_impute(z_col, eigenvectors, n_comp)
            if z_pred is None:
                blocks_skipped_size += 1
                continue

            z_adj, corr = _poly_rescale(z_col, z_pred, npoly=3)

            if not np.isfinite(corr) or corr < min_cor:
                log.debug(
                    "%s trait-col %d: correlation %.3f < %.3f – skipping",
                    block_key, j, corr if np.isfinite(corr) else float("nan"), min_cor,
                )
                blocks_skipped_cor += 1
                continue

            # Fill missing positions only
            fill_idx = vi[missing]
            z_block[fill_idx, j] = z_adj[missing].astype(np.float32)
            if out_mask is not None:
                out_mask[fill_idx, j] = True
            total_imputed += int(missing.sum())

    log.info(
        "Imputation complete: %d cells filled; %d blocks skipped (size/fit), "
        "%d block×trait combos skipped (low correlation)",
        total_imputed, blocks_skipped_size, blocks_skipped_cor,
    )
