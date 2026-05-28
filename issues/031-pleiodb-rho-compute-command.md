## Parent PRD

`issues/prd-rho-matrix.md`

## What to build

The `pleiodb rho <db>` CLI command in **compute mode**: a thin wrapper over `build_rho`
that is invoked when neither `--traits` nor `--traits-file` is supplied.

The command replaces `pleiodb lambda`.  The old `pleiodb lambda` command is removed (or
kept as a hidden alias that prints a deprecation notice and delegates to `pleiodb rho`).

Options for compute mode:
- `--workers / -j` (default 8) — parallel workers for pair computation.
- `--null-thresh` (default 1.0) — z-score threshold for null variant selection.
- `--min-nulls` (default 500) — minimum null variants required per pair; pairs below
  this produce NaN.  Exposed so users can lower the threshold for small databases or
  during testing.
- `--chunk-size` (default 512) — storage chunk size for `rho.bin`.

The dual-mode detection logic (compute vs query) is implemented here but the query
branch remains a stub (or raises `NotImplementedError`) until
`032-pleiodb-rho-query-mode.md` is complete.

See PRD §"Dual-mode `pleiodb rho` command" (compute mode rows) and §"Modules to build
or modify — `pleiodb rho` CLI command".

## Acceptance criteria

- [ ] `pleiodb rho <db>` (no `--traits` / `--traits-file`) calls `build_rho` and
      exits 0.
- [ ] `--workers`, `--null-thresh`, `--min-nulls`, and `--chunk-size` are forwarded
      to `build_rho`.
- [ ] After the command completes, `rho.bin` and `rho.cidx` are present in the
      database directory.
- [ ] `pleiodb lambda` is removed or emits a deprecation warning and delegates to
      `pleiodb rho`.
- [ ] The command is exercised in the integration test suite (use `CliRunner` as in
      `TestHg19Build`): build test db, run `pleiodb rho <db> --min-nulls 50`, assert
      exit code 0 and `rho.bin` present.

## Blocked by

- `issues/029-build-rho-storage.md`
- `issues/030-gwas-database-rename-info-warning.md`

## User stories addressed

- User story 3 (separate post-build command)
- User story 11 (`--workers` parallelism)
