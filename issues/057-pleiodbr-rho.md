## Parent PRD

`issues/prd-pleiodbr.md`

## What to build

`rho(db, traits_1, traits_2)` — retrieve pairwise phenotypic correlations from
the `rho.bin/.cidx` matrix (T×T, float16). `traits_1` and `traits_2` are
character vectors of OpenGWAS trait IDs; all pairs from the cross-product are
returned.

The rho matrix uses the same chunk layout as zscore (512×512, ti-major). Read
only the chunks needed for the requested trait pairs. Decode float16 values
(read as uint16, interpret as IEEE 754 half-precision via bit manipulation or
the `float16` helper). Return a long-format tibble with columns `trait_id_1`,
`trait_id_2`, `rho`.

## Acceptance criteria

- [ ] `rho(db, "ukb-b-19953", "ebi-a-GCST006867")` returns a 1-row tibble
      with the correct rho value
- [ ] `rho(db, c("ukb-b-19953", "ebi-a-GCST90018961"), c("ebi-a-GCST006867"))` returns a 2-row tibble
- [ ] NA rho values (missing pairs) are dropped or returned as NA (consistent
      with other functions)
- [ ] Unrecognised trait IDs raise informative errors

## Blocked by

- `issues/052-pleiodbr-chunk-reader.md`

## User stories addressed

- User story 7
