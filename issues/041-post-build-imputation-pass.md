## Parent PRD

`issues/prd-imputation-performance.md`

## What to build

Remove `impute_z_block` from the trait batch loop in `build_database()`. After all trait
batches are written and both the z-score and Neff files are closed, perform a single
imputation pass:

1. Read the full V×T z-score matrix from the just-written chunks.
2. Read the full V×T Neff matrix from the just-written chunks.
3. Run `impute_z_block` once on the full z-score matrix — each LD block is loaded exactly
   once regardless of how many trait batches the build used.
4. For every imputed position `(v, t)`, set `neff[v, t] = neff_base[t]` (the per-trait
   median Neff already computed during ingestion and stored in `neff_base`).
5. Rewrite both the z-score and Neff chunked matrix files atomically (write to temp paths,
   then replace originals).
6. Populate `n_imputed_arr` and adjust `n_variants_arr` from the imputation mask.
7. Build the COO imputed-positions mask from the full V×T `imputed_mask`.

The batch loop continues to handle VCF ingestion, var_y/Neff estimation, and chunk writing
unchanged. Only the imputation call and the per-batch `imputed_mask_block` construction are
removed from the loop.

## Acceptance criteria

- [ ] `impute_z_block` is no longer called inside the trait batch loop
- [ ] After a build with `--ld-dir`, imputed z-score values are present in the database
- [ ] Imputed positions have finite Neff equal to the per-trait median Neff (`neff_base[t]`),
      not NaN
- [ ] Observed z-scores (non-imputed) are identical before and after the architectural change
- [ ] `n_imputed` and `n_variants` per trait in `traits.tsv` are correct
- [ ] The COO imputed-positions mask reflects the full V×T imputation pass
- [ ] A build crash after `close_write()` but before the rewrite completes leaves the
      pre-imputation files intact (temp → replace atomicity)
- [ ] Integration test in `test_integration.py` passes: checks imputed z-scores non-NaN,
      imputed Neff non-NaN, observed z-scores unchanged

## Blocked by

None — can start immediately.

## User stories addressed

- User story 1
- User story 2
- User story 4
- User story 6
- User story 7
- User story 9
- User story 10
- User story 11
- User story 12
