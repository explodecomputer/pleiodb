## Parent PRD

`issues/prd-hg38-liftover-variants.md`

## What to build

Extend the integration test suite with two new scenarios: (a) an end-to-end build using a
synthetic hg19 variant list against the existing hg38 LD fixture, asserting that the
liftover fires, `alid_hg38` is written, and more cells are imputed than a no-liftover
baseline; (b) a backwards-compatibility check that opens a database built without liftover
and confirms `db.variants["id_hg38"]` returns empty strings without error. Tests follow
the existing pattern in `test_integration.py` (small fixture data, temp directory, call
`build_database()` directly).

## Acceptance criteria

- [ ] New test builds a database with `variants_build="hg19"` and the existing hg38 LD
      fixture; asserts `alid_hg38` column is present in `variants.tsv`
- [ ] Same test asserts `meta.json` contains `"variants_hg38_stored": true`
- [ ] Same test asserts imputed cell count is greater than a baseline build without
      `--ld-dir` (i.e. liftover + LD matching produced real imputation)
- [ ] New backwards-compat test opens an hg38-built database (no liftover); asserts
      `db.variants["id_hg38"]` returns an array of empty strings with no exception
- [ ] All existing tests in `test_integration.py` and `test_impute.py` continue to pass

## Blocked by

- Blocked by `issues/039-store-hg38-alids-in-variants-tsv.md`

## User stories addressed

- User story 10 (meaningful reduction in missingness verified)
- User story 14 (backwards compatibility verified)
- User story 15 (integration tests cover end-to-end imputation path)
