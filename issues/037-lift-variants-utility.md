## Parent PRD

`issues/prd-hg38-liftover-variants.md`

## What to build

A new standalone `lift_variants(variants, from_build, to_build)` function in the existing
`liftover.py` module. It accepts a pleiodb variants structured array and two build strings,
lifts each variant's coordinates using `pyliftover.LiftOver`, and returns a new array of
the same length with updated `chrom`, `pos`, and `id` (ALID) fields. Variants that fail
liftover are kept with their original values (they will be silently skipped at LD matching
time). The function logs how many variants lifted successfully and how many failed.

This function is independent of the build pipeline and can be tested in isolation with a
mocked liftover backend.

## Acceptance criteria

- [ ] `lift_variants(variants, from_build, to_build)` exists in `liftover.py`
- [ ] Returns a structured array of identical length and dtype to the input
- [ ] Successfully lifted variants have updated `chrom`, `pos`, and `id` (ALID) fields
- [ ] Variants where liftover returns no result retain their original field values
- [ ] Logs count of successfully lifted variants and count of failures at INFO level
- [ ] Build strings are normalised via the existing `normalise_build()` helper
- [ ] Unit tests cover: successful lift, liftover failure fallback, array-length invariant
- [ ] Unit tests mock `pyliftover.LiftOver` — no network calls or chain files required

## Blocked by

None — can start immediately.

## User stories addressed

- User story 3 (log message reporting lift success count)
- User story 4 (variants that fail liftover handled gracefully)
- User story 11 (liftover as a standalone, testable function)
- User story 12 (existing `make_lifted_lookup` unchanged)
