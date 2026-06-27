## Parent PRD

`issues/prd-imputation-performance.md`

## What to build

Add a performance section to `scratch/imputation_benchmark.qmd` documenting the before/after
improvement from moving imputation to a single post-build pass.

The section should include:

1. **Timing data from the 5-trait trial** — wall-clock time for the imputation step with
   the new single-pass architecture, compared to the projected time under the old per-batch
   approach (extrapolated from the per-block load times recorded in the build log).
2. **Projected speedup table** — estimated total imputation time for representative build
   sizes (e.g. 100 traits / 1 batch, 512 traits / 1 batch, 4159 traits / 9 batches) under
   old and new architectures, using the measured per-block timing.
3. **Neff correctness note** — brief explanation that imputed positions now carry a finite
   Neff (per-trait median), enabling beta recovery.

## Acceptance criteria

- [ ] New section added to `scratch/imputation_benchmark.qmd` after the existing build-level
      imputation section
- [ ] Section includes at least one timing figure or table with old vs new comparison
- [ ] Projected speedup for a 9-batch (4159-trait) build is shown
- [ ] Document renders without error (`quarto render` or equivalent check)

## Blocked by

- Blocked by `issues/042-five-trait-consistency-trial.md`

## User stories addressed

- User story 8
