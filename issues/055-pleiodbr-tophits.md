## Parent PRD

`issues/prd-pleiodbr.md`

## What to build

`tophits(db, traits, pval = 5e-8)` — retrieve significant hits for one or
more specified traits. `traits` is a required character vector of OpenGWAS
trait IDs.

Internally, convert `pval` to a z-score threshold and check whether a
pre-built COO mask exists at `masks/{thr}.coo.zst` (the database stores masks
for `5e-8` and `1e-5`). If a matching mask exists, load it and filter to the
requested traits. If no mask exists for the requested threshold, fall back to
scanning the relevant trait columns and applying the threshold directly.

Decompress the COO pairs (uint32 v_idx, t_idx), look up z/beta/SE/EAF/Neff
for each pair, join metadata, mark imputed status, return the standard tibble.

## Acceptance criteria

- [ ] `tophits(db, "ukb-b-19953", pval = 5e-8)` returns only rows with
      `pval <= 5e-8` for that trait
- [ ] `traits` argument is required; omitting it raises an error
- [ ] Works for a vector of multiple trait IDs
- [ ] Falls back gracefully when no pre-built mask matches `pval`
- [ ] `imputed` column is correct

## Blocked by

- `issues/052-pleiodbr-chunk-reader.md`

## User stories addressed

- User story 5
