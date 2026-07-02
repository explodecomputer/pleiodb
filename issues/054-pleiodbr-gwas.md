## Parent PRD

`issues/prd-pleiodbr.md`

## What to build

`gwas(db, trait_id)` — pull all V variants for a single trait. Look up the
trait column index, read the full z-score column, reconstruct beta/SE, compute
p-values, join variant metadata (ALID, EAF) and per-variant Neff, and mark
imputed status. Return the standard nine-column tibble with one row per
variant where z is not NA.

## Acceptance criteria

- [ ] `gwas(db, "ukb-b-19953")` returns a tibble with ~V rows and all nine
      columns
- [ ] Rows with z = NA are dropped
- [ ] `imputed` column is correct for the requested trait
- [ ] Unrecognised trait ID raises an informative error

## Blocked by

- `issues/052-pleiodbr-chunk-reader.md`

## User stories addressed

- User story 4
