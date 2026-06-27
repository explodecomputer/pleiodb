## Parent PRD

`issues/prd-dense-vcf-imputation.md`

## What to build

Add `--ld-vcf-threads INTEGER` to the `pleiodb build` CLI (default 8). This controls the number of threads used within each block worker process for parallel VCF region reads. Wire it through `build_database` → `impute_z_block` → each block worker's inner `ThreadPoolExecutor`.

See "CLI" in the parent PRD.

## Acceptance criteria

- [ ] `pleiodb build --ld-vcf-threads 16` is accepted without error
- [ ] The value is passed through to the inner `ThreadPoolExecutor` inside `_impute_block_process`
- [ ] Default of 8 is used when the flag is omitted
- [ ] `pleiodb build --help` shows the new option with a description
- [ ] `build_database()` Python API gains a `ld_vcf_threads: int = 8` parameter

## Blocked by

- Blocked by `issues/048-wire-dense-imputation-pass.md`

## User stories addressed

- User story 5
