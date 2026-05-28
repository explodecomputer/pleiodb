## Parent PRD

`issues/prd-rho-matrix.md`

## What to build

Minor updates to `GWASDatabase` to reflect the lambda → rho rename and to surface the
presence (or absence) of the rho matrix in `pleiodb info` output.

Changes:
- Rename the `lambda_matrix` property to `rho_matrix`; rename `get_lambda_block` to
  `get_rho_block`.
- Update `info()` to include a `rho_present` boolean key in its output dict.
- When the rho matrix is absent, `info()` must include a `warnings` list containing an
  entry such as `"rho matrix absent — run 'pleiodb rho <db>' to compute"`.
- Remove (or deprecate) any remaining `lambda`-named attributes, methods, or meta.json
  key references from `GWASDatabase`.

No new computation logic goes here — this issue is purely renaming and the info warning.

See PRD §"GWASDatabase — minor update" and user story 4.

## Acceptance criteria

- [ ] `db.rho_matrix` returns the rho matrix (or raises a clear error when absent).
- [ ] `db.lambda_matrix` no longer exists (or raises `AttributeError` / deprecation
      warning).
- [ ] `db.info()` returns a dict containing `"rho_present": true` when `rho.bin` exists
      in the database directory.
- [ ] `db.info()` returns `"rho_present": false` and `"warnings"` containing a rho-
      absent message when `rho.bin` is missing.
- [ ] `pleiodb info <db>` (CLI) emits the warning text when rho is absent, visible in
      the JSON output printed to stdout.
- [ ] Integration test: build a test database without running `build_rho`, call
      `pleiodb info`, assert `rho_present` is false and the warnings list is non-empty.
- [ ] Integration test: after running `build_rho`, call `pleiodb info`, assert
      `rho_present` is true and no rho warning appears.

## Blocked by

- `issues/029-build-rho-storage.md`

## User stories addressed

- User story 1 (rename)
- User story 4 (`pleiodb info` warns when rho absent)
