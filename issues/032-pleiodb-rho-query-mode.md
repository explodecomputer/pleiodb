## Parent PRD

`issues/prd-rho-matrix.md`

## What to build

The **query mode** of `pleiodb rho <db>`: activated when `--traits` or `--traits-file`
is supplied.  Reads the stored rho matrix (via `GWASDatabase.rho_matrix`), filters to
the requested traits, and emits the results.

Options for query mode:
- `--traits t1,t2,...` — comma-separated list of trait IDs.
- `--traits-file path` — file with one trait ID per line.
- `--format tsv|json` (default `tsv`) — output format, consistent with `pleiodb query`.
- `--matrix` — pivot to a square matrix (trait IDs as row and column headers) instead
  of the default pairwise list.
- `--output / -o` — write to file instead of stdout.

**Default output (pairwise list)**:
```
trait_id_1   trait_id_2   rho
t1           t2           0.0412
t1           t3           0.1823
t2           t3           0.0091
```
Only the upper triangle; no duplicate pairs; diagonal excluded.

**`--matrix` output**:
```
          t1      t2      t3
t1        1.0     0.0412  0.1823
t2        0.0412  1.0     0.0091
t3        0.1823  0.0091  1.0
```

Trait IDs must not contain commas (commas are the `--traits` delimiter; this constraint
is already documented in `CONTEXT.md`).

See PRD §"Query output format" and §"Dual-mode `pleiodb rho` command" (query mode rows).

## Acceptance criteria

- [ ] `pleiodb rho <db> --traits t1,t2` outputs exactly one data row (the single pair),
      with columns `trait_id_1`, `trait_id_2`, `rho`.
- [ ] `pleiodb rho <db> --traits t1,t2,t3` outputs exactly three data rows (all pairs
      of 3 traits).
- [ ] Diagonal (`rho[t, t]`) is not emitted in pairwise list output.
- [ ] `--matrix` output has shape T×T with all diagonal values = 1.0 and trait IDs
      as headers.
- [ ] `--traits-file <file>` produces the same output as the equivalent `--traits`
      invocation.
- [ ] `--format json` output parses as valid JSON; each object has keys `trait_id_1`,
      `trait_id_2`, `rho`.
- [ ] `--output <file>` writes the result to the specified file rather than stdout.
- [ ] CLI tests use `CliRunner` (as in `TestQueryOutput`) and parse TSV/JSON output
      rather than inspecting internal objects.

## Blocked by

- `issues/030-gwas-database-rename-info-warning.md`
- `issues/031-pleiodb-rho-compute-command.md`

## User stories addressed

- User story 5 (query by `--traits` comma-separated list)
- User story 6 (query by `--traits-file`)
- User story 7 (default pairwise list output)
- User story 8 (`--matrix` flag for square output)
- User story 9 (diagonal always 1.0)
- User story 12 (`--format tsv/json`)
