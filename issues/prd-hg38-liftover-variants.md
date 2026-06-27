# PRD: Lift Variant Coordinates to hg38 at Build Time

## Problem Statement

When building a pleiodb database from hg19 GWAS-VCF files, the user can optionally supply
an LD reference panel (which is in hg38 coordinates) to impute missing z-scores. However,
the LD block index is built by matching variant ALIDs (e.g. `1:12345_A_C`) against the LD
panel's variant list. Because the variant list is in hg19 and the LD panel is in hg38, the
positional identifiers do not match: only ~839 of 95,378 variants (0.9%) are found in the
panel, so imputation fires on almost no variants and the missing-data rate is barely reduced.

The build log already emits a warning about this — `--variants-build is hg19 (not hg38);
variant matching may be incorrect` — but does not resolve it.

## Solution

When `--variants-build hg19` (or any non-hg38 build) is supplied alongside `--ld-dir`,
automatically lift the variant list to hg38 before building the LD block index. Store both
the original (hg19) ALID and the lifted hg38 ALID in the database's `variants.tsv` file so
that the database is self-contained and consumers do not need to re-lift coordinates
externally.

If variants are already in hg38, no liftover is performed and the schema is unchanged.

## User Stories

1. As a pleiodb user, I want LD-based imputation to work correctly when my variant list is
   in hg19, so that missing z-scores are filled at all LD-matched positions rather than
   almost none.

2. As a pleiodb user, I want the build command to automatically lift my hg19 variant
   coordinates to hg38 when an LD panel is supplied, so that I do not need to pre-lift the
   variant list myself.

3. As a pleiodb user, I want the lift to happen silently and correctly, with a log message
   reporting how many variants were successfully lifted, so I can verify the lift without
   inspecting intermediate files.

4. As a pleiodb user, I want variants that fail liftover to be handled gracefully (kept with
   their original coordinates, excluded from LD matching), so that the build does not crash
   or silently corrupt the output.

5. As a pleiodb user, I want both the original hg19 ALID and the hg38 ALID stored in the
   database `variants.tsv`, so that downstream tools can query by either coordinate system
   without re-lifting.

6. As a pleiodb user, I want the hg38 ALID column to be absent from `variants.tsv` when the
   variant list was already in hg38 (no liftover performed), so that the file format stays
   minimal for the common case.

7. As a pleiodb database reader, I want `db.variants["id_hg38"]` to return the hg38 ALID
   for each variant (or an empty string when no liftover was stored), so I can use it
   without parsing meta.json.

8. As a pleiodb database reader, I want `db.variants["id_hg38"]` to return an empty string
   for databases built without liftover, so that existing code reading `id` still works
   unchanged.

9. As a pleiodb user, I want `meta.json` to record whether hg38 ALIDs are stored, so that
   tooling can detect the presence of dual coordinates without reading the variants file.

10. As a pleiodb user, I want the imputation fidelity to be substantially better after this
    fix, so that the validation report shows a meaningful reduction in per-variant
    missingness compared to no imputation.

11. As a pleiodb developer, I want the liftover of the variant array to be a standalone,
    testable function separate from the build pipeline, so it can be tested in isolation and
    reused in other contexts.

12. As a pleiodb developer, I want the existing `make_lifted_lookup` function (used for VCF
    coordinate matching) to remain unchanged, so that VCF matching behaviour is unaffected
    by this change.

13. As a pleiodb developer, I want the build to remove the existing warning about
    `variants_build != hg38` when an LD dir is supplied, replacing it with the liftover
    step, so that the warning no longer appears for a case that is now handled correctly.

14. As a pleiodb developer, I want existing databases built without liftover (hg38 variant
    lists, or no LD dir) to be read without modification, so that backwards compatibility
    is maintained.

15. As a pleiodb developer, I want the tests to exercise liftover with a small synthetic
    variant set against a mock LD panel, so that the end-to-end imputation path is covered
    without relying on the production LD panel.

## Implementation Decisions

### Liftover utility function (`liftover.py`)

A new function `lift_variants(variants, from_build, to_build) -> np.ndarray` is added to
the existing `liftover.py` module. It:
- Accepts a variants structured array (fields: id, chrom, pos, a1, a2) and two build
  strings normalised through the existing `normalise_build()` helper.
- Uses `pyliftover.LiftOver(from_build, to_build)` (auto-downloads chain file if not
  cached locally by pyliftover).
- Returns a new array of the same length and dtype with updated `chrom`, `pos`, and `id`
  (ALID) fields for variants that lift successfully.
- For variants that fail liftover, copies the original row unchanged (they will not match
  any LD block and are silently skipped).
- Logs the number of variants lifted successfully and the number that failed.

This function is deliberately separate from `make_lifted_lookup` (which serves a different
purpose: producing a positional lookup dict for VCF coordinate matching per-trait).

### Build pipeline changes (`build.py`)

In `build_database()`, after loading the variant array and before building the LD block
index:
- If `ld_dir` is provided and `canon_build != "hg38"`, call `lift_variants(variants,
  canon_build, "hg38")` to produce a parallel hg38 array.
- Pass the hg38 array (not the original) to `build_block_index()`.
- Store the hg38 array for writing to `variants.tsv`.
- If `ld_dir` is not provided, or `canon_build` is already hg38, skip the lift (no change
  to current behaviour).
- Remove the existing warning about non-hg38 `variants_build` when LD dir is supplied; it
  is superseded by the liftover step.

### Variants TSV schema change (`build.py`, `db.py`)

`_write_variants_tsv()` gains an optional `variants_hg38` parameter. When provided, a
third column `alid_hg38` is written after `eaf`. When absent, the two-column format is
unchanged.

`_load_variants_tsv()` in `db.py` detects the `alid_hg38` column from the header and
populates a new `id_hg38` field in the returned struct array. When the column is absent,
`id_hg38` is filled with empty strings. This is fully backwards-compatible.

`_VARIANTS_DT` in `db.py` gains a new field `("id_hg38", "U64")` defaulting to `""`.

### Database metadata (`meta.json`)

Two new keys are written when a liftover is performed:
- `"variants_hg38_stored": true` — signals that `alid_hg38` is present in `variants.tsv`.
- No new keys are written when no liftover is performed (existing databases unaffected).

### LD block index (`impute.py`)

`build_block_index()` requires no internal changes. It already matches on ALID strings;
passing it the hg38 variant array is sufficient.

## Testing Decisions

A good test exercises only the externally observable behaviour of a module — inputs and
outputs — not its internal implementation. Tests should not depend on the production LD
panel or live liftover chain files.

### `lift_variants` (unit tests in `test_liftover.py` or `test_impute.py`)

- Provide a small synthetic variants array with known hg19 positions.
- Mock `pyliftover.LiftOver` to return predictable hg38 positions.
- Assert that the returned array has updated `id`, `chrom`, `pos` for lifted variants.
- Assert that variants where liftover returns no result keep their original values.
- Assert array length is unchanged.

### Build integration (extend `test_integration.py`)

The existing integration tests in `test_integration.py` build small databases with fixture
data. Extend them to:
- Build a database with a synthetic hg19 variant list + synthetic LD panel (hg38) and
  `variants_build="hg19"`.
- Assert that `variants.tsv` in the output contains an `alid_hg38` column.
- Assert that `db.variants["id_hg38"]` is non-empty for lifted variants.
- Assert that `meta.json` contains `"variants_hg38_stored": true`.
- Assert that imputation produces more filled cells than a baseline without liftover
  (i.e., the block index matched more variants).

### Backwards compatibility (extend `test_integration.py`)

- Open an existing database built without liftover (hg38 variant list).
- Assert that `db.variants["id_hg38"]` returns empty strings for all rows.
- Assert no exception is raised.

### Prior art

The existing `test_integration.py` uses a small fixture LD panel (nested layout) and
builds a real database in a temp directory. New tests should follow the same pattern:
fixture data, temp dir, call `build_database()` directly.

## Out of Scope

- Liftover of the variant list when `--ld-dir` is not supplied. The hg38 ALIDs are only
  needed for LD panel matching; there is no reason to perform liftover otherwise.
- Liftover for builds other than hg19 → hg38. Only this direction is currently needed.
- Modifying the query API, rho computation, or significance mask logic.
- Storing hg38 positions in the binary matrices (`zscore.bin`, `neff.bin`); the V-axis
  ordering is unchanged.
- Providing a CLI command to re-lift an existing database post-build.

## Further Notes

- `pyliftover` auto-downloads chain files from UCSC on first use and caches them in
  `~/.pyliftover`. If the build environment has no internet access, the chain file
  `hg19ToHg38.over.chain.gz` will need to be pre-cached or provided at a known path.
  This is an operational concern, not a code change.
- The fraction of variants failing liftover from hg19 to hg38 is typically <1% for
  well-curated variant lists; the log message is sufficient for monitoring.
- The hg38 ALID stored in the database uses the same ALID convention as the rest of
  pleiodb (`CHROM:POS_A1_A2` with alphabetically ordered alleles). No new format is
  introduced.
