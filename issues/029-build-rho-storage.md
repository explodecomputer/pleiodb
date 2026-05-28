## Parent PRD

`issues/prd-rho-matrix.md`

## What to build

A `build_rho(db_path, z_null_thresh=1.0, min_nulls=500, workers=8, chunk_size=512)`
function that computes the full T×T rho matrix from the z-scores already stored in the
database and writes it to `rho.bin` / `rho.cidx` alongside the existing database files.

The function:
- Reads the stored z-score matrix chunk by chunk (no VCF re-reads).
- For each upper-triangle pair (j, k) filters variants to `|z_j| < z_null_thresh AND
  |z_k| < z_null_thresh` and calls `estimate_rho_cml` (from issue
  `028-estimate-rho-cml-pure-function.md`).
- Sets the diagonal to 1.0 explicitly; mirrors the upper triangle to the lower triangle.
- Stores the result as a symmetric float16 chunked matrix at `rho.bin` / `rho.cidx`,
  consistent with the z-score matrix storage format.
- Parallelises pair computation across `workers` threads/processes.
- Records `rho_chunk_shape` and `rho_null_z_thresh` in `meta.json` (replacing the old
  `lambda_chunk_shape` / `lambda_null_z_thresh` keys if present).

The `min_nulls` parameter defaults to 500 but is passable as an argument so that
integration tests using the ~500-variant test dataset can set it to a lower value
(e.g. 50) and still receive finite rho estimates rather than NaN for every pair.

This issue also retires the old `build_lambda` function: rename / replace it entirely.
The `lambda.bin` / `lambda.cidx` storage format is not migrated; existing databases
will simply lack `rho.bin` (which `pleiodb info` will warn about — see
`030-gwas-database-rename-info-warning.md`).

See PRD §"Variant source", §"Symmetry and diagonal", and §"Modules to build or modify"
for the full specification.

## Acceptance criteria

- [ ] `build_rho(db_path)` runs to completion on the test database (built from
      `tests/test_data/`).  Pass `min_nulls=50` (or similar) to avoid NaN-everywhere
      with the small test dataset.
- [ ] After `build_rho`, `rho.bin` and `rho.cidx` exist inside the database directory.
- [ ] The diagonal of the resulting matrix is 1.0 for all traits.
- [ ] All off-diagonal finite values are in the open interval (−1, 1).
- [ ] The matrix is symmetric: `rho[j, k] == rho[k, j]` for all j, k.
- [ ] Pairs with fewer than `min_nulls` null variants store NaN (verifiable by setting
      `min_nulls` to a value larger than the available variant count).
- [ ] `meta.json` contains `rho_chunk_shape` and `rho_null_z_thresh`; the old
      `lambda_*` keys are absent.
- [ ] `build_lambda` is removed (or re-exports as a deprecated alias that raises a
      deprecation warning).

## Blocked by

- `issues/028-estimate-rho-cml-pure-function.md`

## User stories addressed

- User story 1 (rename lambda → rho)
- User story 3 (separate post-build step)
- User story 10 (NaN for < min_nulls null variants, configurable threshold)
- User story 11 (parallelised across trait pairs)
- User story 14 (float16 symmetric chunked storage)
