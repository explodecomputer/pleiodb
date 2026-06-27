## Parent PRD

`issues/prd-dense-vcf-imputation.md`

## What to build

Store `vcf_path` and `vcf_build` columns in the `traits.tsv` written inside every `.pleiodb` database directory. When a database is opened, these columns are read back and available on the trait metadata. Databases built without these columns continue to work (empty strings for both fields).

See "VCF path persistence" in the parent PRD.

## Acceptance criteria

- [ ] `_write_traits_tsv` writes `vcf_path` and `vcf_build` as the last two columns of `traits.tsv`
- [ ] `_load_traits_tsv` reads `vcf_path` and `vcf_build` when present; returns empty string for both when absent (backward-compatible)
- [ ] A database built with the new code contains the correct paths when inspected with `pleiodb info` or by reading `traits.tsv` directly
- [ ] A database built with old code (no `vcf_path` column) can still be opened without error
- [ ] Unit test: build a small fixture database and assert `traits.tsv` contains `vcf_path` matching the input traits TSV

## Blocked by

None — can start immediately.

## User stories addressed

- User story 4
- User story 9
