## Parent PRD

`issues/prd.md`

## What to build

Implement the public entry point `impute_z_block()` in `src/pleiodb/impute.py`, which
composes the kernel helpers from `issues/033-imputation-kernel.md` and the block index from
`issues/034-ld-block-index.md` into a single call that the build pipeline can invoke.

**Signature:**
```python
def impute_z_block(
    z_block: np.ndarray,          # float32, shape (V, B), modified in-place
    variants: np.ndarray,         # structured array with 'id', 'chrom', 'pos', 'a1', 'a2'
    eaf_arr: np.ndarray,          # float32, shape (V,)
    block_index: dict,            # output of build_block_index()
    thresh: float = 0.9,
    min_cor: float = 0.7,
    out_mask: np.ndarray | None = None,  # bool, shape (V, B), set True for imputed cells
) -> None:
```

**Per-block, per-trait logic:**

1. Extract z-scores for the matched variants (`variant_indices`) from the current trait column.
2. Compute theoretical SE from EAF; call `_se_outliers` to find unstable beta positions and
   add them to the missing set.
3. Load the LD submatrix via `_load_ld_submatrix`; call `_ld_pca` to obtain eigenvectors.
4. Call `_elastic_net_impute` with the z-scores and eigenvectors.
5. Call `_poly_rescale`; if the returned correlation is below `min_cor`, log at DEBUG level
   and skip this block × trait (write nothing).
6. Fill only the missing positions in `z_block` with rescaled predictions; set the
   corresponding cells in `out_mask` to True.

**Logging:** emit an `INFO`-level summary per trait batch: number of blocks processed, number
skipped (too few observed variants, low correlation), and total cells imputed. Emit a
`WARNING` if `eaf_arr` has NaN values for any variant in the block (EAF is required for SE
estimation).

## Acceptance criteria

- [ ] `impute_z_block` is exported from `impute.py`.
- [ ] Calling it with a `z_block` that has deliberate NaN gaps and a synthetic block index
      (constructed in-memory, using `_load_ld_submatrix` patched to return a fixture matrix)
      fills those NaN positions with non-NaN floats and sets `out_mask` to `True` at the same
      positions.
- [ ] Positions that were already non-NaN before the call are unchanged.
- [ ] When the block-level Pearson correlation is below `min_cor` the block is skipped: the
      NaN values remain NaN and `out_mask` is not updated.
- [ ] When a block has fewer than 2 observed (non-NaN) z-scores for a given trait, it is
      skipped gracefully (no exception).
- [ ] `out_mask=None` is valid (no mask tracking performed).
- [ ] Log output (captured via `caplog`) contains an INFO-level summary after a batch
      completes, and a DEBUG message for each skipped block.
- [ ] `pytest tests/test_impute.py` (including the new tests for this function) passes
      without requiring the real LD panel on disk.

## Blocked by

- `issues/033-imputation-kernel.md` (needs `_ld_pca`, `_elastic_net_impute`, `_poly_rescale`,
  `_se_outliers`)

## User stories addressed

- User story 1 (imputation runs during build without separate post-processing)
- User story 3 (min_cor threshold skips poorly-imputed blocks)
- User story 4 (thresh parameter controls PCA variance fraction)
- User story 5 (log messages report block processing summary)
- User story 7 (hg38 warning emitted when block coordinate mismatch is suspected)
- User story 12 (imputation logic isolated in a dedicated module)
