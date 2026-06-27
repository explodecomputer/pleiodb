## Problem Statement

LD-based z-score imputation is called once per trait batch during `build_database()`. With the
default batch size of 512 traits, a 4159-trait full build requires 9 separate imputation passes —
each of which loads all matched LD blocks from disk. A single large LD block (e.g. chr10
100331627–104382250, 433 MB compressed) takes ~22 seconds to load with the current pandas CSV
reader. At 1345 matched blocks per pass and 9 passes, the total imputation time is approximately
76 hours — making a full-scale imputed build infeasible.

Even the 100-trait trial build takes ~8.5 hours for imputation despite having only one batch,
because the per-block loading cost (~22 s × 1345 blocks) is dominated by reading and parsing
matrices far larger than the matched submatrix actually used.

There is also a correctness bug in the current code: imputed z-score positions have NaN stored
for Neff. Neff is derived from SE values (available only from the VCF during ingestion), so
variants that were not in the VCF — the imputed positions — get no Neff. Without a finite Neff,
downstream tools cannot recover betas from imputed z-scores via
`SE = sqrt(var_y / (Neff × 2·EAF·(1−EAF)))`.

The fix is to assign the per-trait median Neff (`neff_base[t]`) to all imputed positions.
This is appropriate because the stored per-variant Neff formula —
`neff[v] = var_y / (SE[v]² × 2·EAF[v]·(1−EAF[v]))` — is approximately EAF-independent:
the EAF terms in numerator and denominator cancel under the standard GWAS model, leaving
`neff[v] ≈ N_study` for all variants. The median across observed variants is therefore a valid
approximation for imputed positions regardless of their allele frequency.

## Solution

Move imputation out of the per-batch ingestion loop and into a single post-build pass:

1. During batch ingestion, write z-scores to disk as usual (with the missing-value sentinel
   for unobserved positions). Do not run imputation inside the batch loop.
2. After all trait batches are written and closed, perform a single imputation pass:
   - Read the complete V×T z-score matrix back from the written chunks.
   - Run `impute_z_block` once against that full matrix — each LD block is loaded exactly once
     regardless of how many trait batches the build used.
   - For each imputed position `(v, t)`, also set `neff[v, t] = neff_base[t]` (the per-trait
     median Neff computed during ingestion), so that betas can be recovered downstream.
   - Rewrite both the z-score and Neff files with imputed values filled in.
   - Update per-trait imputed-count metadata.
3. Validate correctness and performance with a small 5-trait trial that compares the old and
   new code paths for identical imputed-cell counts, imputed values, and observed z-scores.
4. Document the before/after timing in the existing Quarto imputation benchmark file.

## User Stories

1. As a database builder, I want the full 4159-trait imputed build to complete in a practical
   timeframe (hours, not days), so that I can run builds routinely.
2. As a database builder, I want imputation to load each LD block exactly once per build,
   regardless of how many trait batches the build uses, so that build time scales with the
   number of LD blocks rather than (n_batches × n_blocks).
3. As a database builder, I want the post-build imputation pass to produce results numerically
   identical to the current per-batch approach, so that I can trust the migration does not
   silently change any imputed values.
4. As a database builder, I want the batch ingestion loop (VCF reading, Neff estimation, chunk
   writing) to behave exactly as before, so that the non-imputation steps are unaffected.
5. As a database builder, I want a 5-trait trial build to serve as a fast regression test that
   can complete in minutes, so that I can validate correctness and measure speed improvement
   without waiting for a full build.
6. As a database builder, I want per-trait imputed-count and observed-count metadata in
   `traits.tsv` to remain correct after the architectural change, so that downstream analyses
   are not affected.
7. As a database builder, I want the imputed-positions COO mask (which records which cells were
   imputed) to be built from the single post-build pass, so that it reflects the full V×T matrix.
8. As a researcher, I want to see a before/after performance comparison in the Quarto benchmark
   document, so that the improvement is reproducible and auditable.
9. As a developer, I want the post-build imputation pass to reuse existing functions
   (`impute_z_block`, `get_block`, `encode_z`, `ChunkedMatrix`) with no new storage formats or
   schemas, so that the change is minimal and reviewable.
10. As a developer, I want the z-score rewrite step to be atomic (write temp → replace) so that
    a crashed build leaves either the pre- or post-imputation file, never a corrupt hybrid.
11. As a researcher querying an imputed database, I want imputed positions to have a finite Neff
    value stored, so that I can recover betas via `SE = sqrt(var_y / (Neff × 2·EAF·(1−EAF)))`.
12. As a developer, I want the Neff assigned to imputed positions to be the per-trait median
    Neff (already computed during ingestion), so that no additional data is needed and the
    approximation is EAF-independent.

## Implementation Decisions

- **Imputation removed from the batch loop**: The call to `impute_z_block` and the per-batch
  construction of `imputed_mask_block` are removed from `build_database()`'s main loop.
  The loop continues to handle VCF ingestion, Neff estimation, and chunk writing unchanged.

- **n_variants_arr computed pre-imputation**: During the batch loop, `n_variants_arr[t]` is
  set to the count of finite z-scores (observed only, before imputation). After the post-build
  pass, `n_imputed_arr[t]` is populated from the imputation mask, and `n_variants_arr` is
  adjusted accordingly (subtract imputed count). This matches the current semantics.

- **Post-build pass reads V×T into RAM**: After `zscore_mat.close_write()`, the full V×T
  float32 z-score matrix is loaded via `get_block(0, V, 0, T)` — the same call already used by
  `_build_rho()`. For the 100-trait trial this is ~38 MB; for the 4159-trait full build it is
  ~1.6 GB. The batch loop arrays (z_block, se_block, neff_block) are freed before this step, so
  peak RAM during the post-build pass is ~1.6 GB.

- **Neff assigned at imputed positions**: The per-trait median Neff (`neff_base[t]`, already
  computed during the batch loop) is written to `neff[v, t]` for all imputed positions. This is
  valid because the stored neff formula `neff[v] = var_y / (SE[v]² × 2·EAF[v]·(1−EAF[v]))`
  is approximately EAF-independent (EAF cancels under the standard GWAS model), so the median
  across observed variants approximates N_study and is appropriate for any imputed position.

- **Both z-score and Neff files rewritten atomically**: After imputation, two new
  `ChunkedMatrix` writers are opened at temp paths (z-score and Neff), all chunks are
  re-encoded and written, then both temp files replace their originals. This is the same
  pattern used for significance masks.

- **`impute_z_block` signature unchanged**: The function receives a (V, T) float32 array (same
  shape as before, just now the full V×T rather than V×B). All internal logic is unchanged.

- **`_load_ld_submatrix` unchanged**: Streaming the submatrix (rather than loading the full
  gzip file) is a secondary optimisation deferred to a later issue. The architectural change
  (one load per block vs. N_batches loads) already reduces total LD I/O by a factor of N_batches.

- **5-trait validation trial**: A small build using 5 traits is run twice — once with the
  current code (on a branch without the change) and once with the new code — and the resulting
  databases are compared: imputed cell positions, imputed z-score values, observed z-score
  values, and n_imputed per trait. Both builds also produce per-block timing logs.

- **Quarto performance documentation**: A new section is added to
  `scratch/imputation_benchmark.qmd` reporting the per-block LD load time, per-batch imputation
  time (old), single-pass time (new), and projected speedup for 1-batch and 9-batch builds.

## Testing Decisions

- **Good tests test observable behaviour, not internals**: Tests should assert on the database
  contents after a build (z-score values, n_imputed, COO mask) rather than on call counts or
  internal intermediate arrays.
- **`test_integration.py`** (prior art: `TestLiftoverImputation`): Add a test that runs
  `build_database()` with a small fixture (5 traits, a small LD panel), checks that imputed
  cell positions match the expected set, and checks that observed z-scores are unchanged.
- **Consistency test**: Run build with and without imputation on the same fixture, verify the
  set of imputed positions is non-empty and z-scores differ only at imputed positions.
- **No unit tests needed for the restructured batch loop**: The loop changes are mechanical
  (removal of imputation call, deferral of n_imputed update). Behaviour is verified by the
  integration test above.

## Out of Scope

- Streaming `_load_ld_submatrix` (replacing `pd.read_csv` with a gzip line-by-line parser).
  This is a complementary optimisation that can be addressed in a follow-up issue.
- Binary cache of LD panels (converting `.vcor1.gz` to seekable float32 files). Also deferred.
- Parallelising imputation across LD blocks. Not needed once the architectural fix is in place.
- Changes to the query API, rho computation, or any downstream consumers.

## Further Notes

The `get_block(0, V, 0, T)` call for reading back V×T z-scores is already proven in production
via `_build_rho()` (build.py). The pattern of writing a temp ChunkedMatrix and replacing the
original is also used in the significance mask and rho computation steps.

The 5-trait trial build is the fastest path to validating correctness before the full 100-trait
or 4159-trait rebuild. A trait list file of 5 traits should be created in `scratch/` pointing
at existing VCF files from the trial100 set.
