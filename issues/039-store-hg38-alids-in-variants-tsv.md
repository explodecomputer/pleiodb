## Parent PRD

`issues/prd-hg38-liftover-variants.md`

## What to build

Extend the `variants.tsv` schema to store an optional `alid_hg38` column alongside the
existing `alid` and `eaf` columns, written when a liftover was performed during build.
Extend the database reader so that `db.variants["id_hg38"]` returns the hg38 ALID for
each variant, or an empty string for databases built without liftover. Record a
`variants_hg38_stored` flag in `meta.json`. All changes are backwards-compatible: existing
databases with a two-column `variants.tsv` are read without error.

## Acceptance criteria

- [ ] When liftover was performed, `variants.tsv` contains a third column `alid_hg38`
- [ ] When no liftover was performed, `variants.tsv` has the existing two-column format
- [ ] `db.variants["id_hg38"]` returns the hg38 ALID string for each variant when stored
- [ ] `db.variants["id_hg38"]` returns `""` for all variants when no liftover was stored
- [ ] `meta.json` contains `"variants_hg38_stored": true` when liftover ALIDs are present
- [ ] `meta.json` does not gain new keys when no liftover was performed
- [ ] Opening a pre-existing database (no `alid_hg38` column) raises no exception and
      returns empty strings for `id_hg38`

## Blocked by

- Blocked by `issues/038-build-pipeline-liftover-wiring.md`

## User stories addressed

- User story 5 (both ALIDs stored in database)
- User story 6 (column absent when no liftover)
- User story 7 (`db.variants["id_hg38"]` accessible)
- User story 8 (empty string for databases without liftover)
- User story 9 (`meta.json` flag)
- User story 14 (backwards compatibility)
