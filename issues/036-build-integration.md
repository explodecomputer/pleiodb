## Parent PRD

`issues/prd.md`

## What to build

Wire imputation into the full `pleiodb build` pipeline end-to-end:

1. **`build_database()` changes** — add `ld_dir`, `ld_ancestry`, `ld_thresh`, `ld_min_cor`
   parameters. Before the trait loop, call `build_block_index()` once if `ld_dir` is set.
   Inside each trait batch, after VCF parsing fills `z_block`, call `impute_z_block()` and
   accumulate the resulting `imputed_mask` into running COO lists. After all batches, write
   `imputed.coo.zst`.

2. **CLI** — add four new options to the `build` sub-command:
   `--ld-dir`, `--ld-ancestry` (default `EUR`), `--ld-thresh` (default `0.9`),
   `--ld-min-cor` (default `0.7`). When `--ld-dir` is supplied without
   `--variants-build hg38`, emit a `WARNING`.

3. **`imputed.coo.zst`** — written to the root of the output directory alongside
   `zscore.bin`. Format: zstd-compressed uint32 pairs `(v_idx, t_idx)`, sorted
   lexicographically by `(v_idx, t_idx)`. An empty file (0 pairs) is written when imputation
   is disabled so consumers can always probe for its existence.

4. **`traits.tsv` extension** — add an `n_variants_imputed` column (integer count of imputed
   z-scores per trait) after the existing `n_variants` column. The `n_variants` column
   continues to count *observed* (non-imputed) z-scores only.

The `n_variants_arr` counter (observed) is unaffected; a new `n_imputed_arr` accumulates
imputed counts in the same t_block loop.

## Acceptance criteria

- [ ] `build_database()` accepts `ld_dir=None` and behaves identically to the pre-imputation
      build when `ld_dir` is not set (no regressions).
- [ ] `pleiodb build --help` lists `--ld-dir`, `--ld-ancestry`, `--ld-thresh`, `--ld-min-cor`.
- [ ] Running `pleiodb build` with `--ld-dir` pointing at a synthetic 2-block LD fixture and
      a traits TSV that has deliberate missing variants produces an `imputed.coo.zst` file
      whose decoded pairs match the expected imputed positions.
- [ ] `imputed.coo.zst` is also written (as an empty file) when `--ld-dir` is not supplied.
- [ ] The `(v_idx, t_idx)` pairs in `imputed.coo.zst` are sorted lexicographically.
- [ ] `traits.tsv` contains an `n_variants_imputed` column; its values equal the number of
      imputed cells in `imputed.coo.zst` for each trait.
- [ ] `n_variants` in `traits.tsv` counts only observed (pre-imputation) z-scores; its value
      does not increase after imputation fills cells.
- [ ] Decoded z-scores at imputed positions round-trip through `encode_z` / `decode_z`
      without corruption (spot-check in the integration test).
- [ ] A WARNING is logged when `--ld-dir` is passed without `--variants-build hg38`.
- [ ] `pytest tests/test_integration.py` passes with no regressions.
- [ ] A new test in `tests/test_integration.py` (or a companion file) builds a database with
      the synthetic LD fixture and asserts on `imputed.coo.zst` contents, `n_variants`, and
      `n_variants_imputed` in `traits.tsv`.

## Blocked by

- `issues/034-ld-block-index.md`
- `issues/035-impute-z-block-api.md`

## User stories addressed

- User story 1 (opt-in via `--ld-dir` flag)
- User story 3 (min_cor passed through CLI)
- User story 4 (thresh passed through CLI)
- User story 5 (build log messages)
- User story 6 (partial LD panel doesn't abort build)
- User story 7 (warning for non-hg38 build)
- User story 8 (imputed z-scores stored at same precision)
- User story 9 (`imputed.coo.zst` distinguishes observed from imputed)
- User story 10 (`n_variants` counts observed only)
- User story 11 (`n_variants_imputed` column)
