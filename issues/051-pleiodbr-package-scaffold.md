## Parent PRD

`issues/prd-pleiodbr.md`

## What to build

Scaffold the `pleiodbr` R package: `DESCRIPTION`, `NAMESPACE`, directory
structure, and `open_pleiodb()`. The function reads `meta.json`,
`variants.tsv`, and `traits.tsv` from a `.pleiodb` directory and returns an
S3 object of class `"pleiodb"`. A `print.pleiodb` method shows V, T, format
version, and path. Raise an informative error if `format_version` in
`meta.json` is newer than the package supports (current: 3).

No query functions yet — just open, inspect, and close.

## Acceptance criteria

- [ ] `install_github("explodecomputer/pleiodbr")` succeeds with no Python
      dependency
- [ ] `open_pleiodb("/path/to/main.pleiodb")` returns a `"pleiodb"` S3 object
- [ ] `print(db)` displays V, T, format version, and path
- [ ] Opening a dataset with `format_version > 3` raises a clear error
- [ ] `DESCRIPTION` declares `zstd`, `tibble`, `dplyr` as Imports

## Blocked by

None — can start immediately.

## User stories addressed

- User story 1
- User story 8
