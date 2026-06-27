## Parent PRD

`issues/prd-imputation-performance.md`

## What to build

Create a 5-trait trial build and a comparison script that validates the new post-build
imputation pass produces results consistent with the expected behaviour, and measures
the wall-clock speedup.

1. Create `scratch/build2/traits_trial5.tsv` — 5 traits drawn from the existing
   `traits_trial100.tsv` that are known to have good LD panel coverage (i.e. variants
   that match LD blocks after liftover).
2. Run a build with the new code using the 5-trait list and the hg38 LD panel.
3. Write a comparison script `scratch/build2/validate_trial5.py` (or add a section to
   the existing `validate_trial.py`) that checks:
   - Imputed cell count is non-zero
   - Imputed positions have finite z-scores
   - Imputed positions have finite Neff (not NaN)
   - Observed (non-imputed) z-scores match an independent reference (e.g. re-reading
     from a reference build run without imputation)
   - `n_imputed` per trait in `traits.tsv` is consistent with the COO mask
4. Record per-block LD load timing from the build log for use in the Quarto doc.

## Acceptance criteria

- [ ] `scratch/build2/traits_trial5.tsv` exists with exactly 5 traits
- [ ] 5-trait build completes without error using `--ld-dir` and `--variants-build hg19`
- [ ] Validation script reports: imputed cells > 0, all imputed z-scores finite, all
      imputed Neff finite
- [ ] No regression in observed z-scores compared to a build without imputation
- [ ] Build log contains per-block timing lines usable for the performance section

## Blocked by

- Blocked by `issues/041-post-build-imputation-pass.md`

## User stories addressed

- User story 3
- User story 5
