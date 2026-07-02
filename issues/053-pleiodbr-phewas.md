## Parent PRD

`issues/prd-pleiodbr.md`

## What to build

`phewas(db, variant)` — first end-to-end query returning the full tibble
schema. `variant` is either an ALID string (`"1:103574777_C_G"`) or a region
string (`"1:103e6-104e6"`). Region strings are parsed by splitting on `:` then
`-`; positions accept scientific notation.

For a single ALID: look up the variant row index via the variants table, read
the full z-score row, reconstruct beta/SE (see formula in PRD), compute
two-sided p-values, attach EAF and Neff, join trait metadata, and mark
imputed status from `imputed.coo.zst`.

For a region string: find all variant indices whose `chrom` matches and `pos`
falls within `[start, end]`, then union their rows into a single tibble (one
row per variant × trait pair where z is not NA).

## Acceptance criteria

- [ ] `phewas(db, "16:53800954_C_T")` returns a tibble with T rows (one per
      trait) and all nine columns populated
- [ ] `phewas(db, "16:53e6-54e6")` returns a tibble covering all variants in
      that window across all traits
- [ ] Rows with z = NA are dropped
- [ ] `imputed` column is TRUE for pairs present in `imputed.coo.zst`
- [ ] Unrecognised ALID raises an informative error
- [ ] Malformed region string raises an informative error

## Blocked by

- `issues/052-pleiodbr-chunk-reader.md`

## User stories addressed

- User story 2
- User story 3
