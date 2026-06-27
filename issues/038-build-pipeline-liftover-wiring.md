## Parent PRD

`issues/prd-hg38-liftover-variants.md`

## What to build

Wire `lift_variants()` into `build_database()` so that when `--ld-dir` is supplied and the
variant list is not already in hg38, the variant coordinates are lifted to hg38 before
`build_block_index()` is called. The hg38 variant array is passed to `build_block_index()`
in place of the original array, so LD panel ALID matching uses the correct coordinate
system. The existing warning about `variants_build != hg38` when an LD dir is supplied is
removed, as the liftover step now handles this case correctly.

No schema changes to `variants.tsv` or `meta.json` in this slice — that is covered in 039.

## Acceptance criteria

- [ ] When `--ld-dir` is set and `variants_build` is `hg19`, `lift_variants()` is called
      before `build_block_index()`
- [ ] `build_block_index()` receives the hg38-lifted variant array, not the original
- [ ] When `variants_build` is already `hg38`, no liftover is performed (existing behaviour)
- [ ] When `--ld-dir` is not set, no liftover is performed (existing behaviour)
- [ ] The warning `"--variants-build is hg19 (not hg38); variant matching may be incorrect"`
      is removed
- [ ] Build log shows substantially more LD block matches than the previous 839/95,378
      when run against the production LD panel with hg19 variants

## Blocked by

- Blocked by `issues/037-lift-variants-utility.md`

## User stories addressed

- User story 1 (imputation works correctly with hg19 variant list)
- User story 2 (build auto-lifts without user pre-lifting)
- User story 3 (log message on lift success count)
- User story 4 (failed liftovers handled gracefully)
- User story 13 (warning removed, superseded by liftover step)
