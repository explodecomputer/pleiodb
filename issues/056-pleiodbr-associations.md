## Parent PRD

`issues/prd-pleiodbr.md`

## What to build

`associations(db, variants, traits)` — arbitrary V×T block query. `variants`
is a character vector of ALIDs; `traits` is a character vector of OpenGWAS
trait IDs. Look up row and column indices, read the minimal set of chunks
covering the requested block, decode z/beta/SE/EAF/Neff, mark imputed status,
and return a long-format tibble (one row per variant × trait pair where z is
not NA).

## Acceptance criteria

- [ ] `associations(db, c("16:53800954_C_T", "19:45412079_C_T"), c("ukb-b-19953", "ebi-a-GCST006867"))` returns a tibble with up to 4 rows
- [ ] Result is long format (not a matrix)
- [ ] Rows with z = NA are dropped
- [ ] Unrecognised variant or trait IDs raise informative errors
- [ ] `imputed` column is correct

## Blocked by

- `issues/052-pleiodbr-chunk-reader.md`

## User stories addressed

- User story 6
