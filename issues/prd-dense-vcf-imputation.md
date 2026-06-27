## Problem Statement

LD-based z-score imputation in pleiodb trains its elastic net on the ~95k curated pleiodb variants that fall within each LD block. The LD reference panel blocks typically contain 5,000–22,000 variants, but the pleiodb variant list covers only **0.3–1.5% of them** (median 58 pleiodb variants out of 7,057 LD panel variants per block). The GWAS VCF files used to build the database contain z-scores for the vast majority of those LD panel variants — the information exists but is never extracted. As a result, the elastic net fits on a tiny fraction of the available signal, severely limiting imputation quality regardless of how much compute is applied.

## Solution

Split imputation into two clean phases:

1. **Build phase (unchanged):** Read VCFs for the 95k pleiodb positions as now, write the V×T z-score matrix. No imputation happens here.

2. **Imputation phase (new dense design):** For each LD block, identify which traits have at least one missing pleiodb variant in that block. For each such trait, read the VCF for just that block's genomic region to obtain dense z-scores for all LD panel variants in the block. Fit the elastic net on this dense training set, predict for the missing pleiodb positions, and write back only those predictions. No scratch files — VCF regions are read on-demand per block.

The final V×T matrix is unchanged in format. Only the elastic net training data improves.

## User Stories

1. As a database builder, I want imputed z-scores trained on hundreds of dense LD panel variants per block rather than dozens of sparse pleiodb variants, so that elastic net fits are well-conditioned and imputation quality is high.
2. As a database builder, I want the build phase to remain exactly as-is, so that the dense imputation is an additive improvement with no regression risk to existing behaviour.
3. As a database builder, I want VCF reads during the imputation phase to happen only for (block, trait) pairs where data is actually missing, so that fully-observed traits incur no unnecessary I/O.
4. As a database builder, I want VCF paths and build metadata stored inside the database's traits.tsv, so that the imputation phase can locate the original VCFs without re-supplying the traits TSV.
5. As a database builder, I want VCF reads within a block to run in parallel (via threads), and blocks to run in parallel (via processes), so that imputation completes in minutes rather than hours even for thousands of traits.
6. As a database builder, I want liftover to be handled automatically when a VCF's coordinate build differs from the LD panel (hg38), so that dense region queries return the correct variants regardless of VCF build.
7. As a researcher, I want imputed z-scores to be closer to the true GWAS value, so that downstream colocalization and fine-mapping analyses are more reliable.
8. As a developer, I want the dense z-score assembly (VCF reads → LD-row-indexed matrix) and the elastic net fitting to be independently testable, so that correctness can be verified without a full build.
9. As a database builder, I want the imputation phase to log clearly when a VCF is inaccessible for a given block, and fall back to sparse-only imputation for that block, rather than crashing the entire run.
10. As a database builder, I want the elastic net training set to be the dense VCF z-scores and the prediction targets to be the missing pleiodb positions — never trained on what it predicts — so there is no data leakage.

## Implementation Decisions

### VCF path persistence

- Add `vcf_path` and `vcf_build` columns to the `traits.tsv` written inside the `.pleiodb` directory during `_write_traits_tsv`.
- `_load_traits_tsv` (used by `GWASDatabase`) reads these columns if present; absent columns yield empty strings (backward-compatible).
- VCF paths are stored as written (absolute paths). No path rewriting.

### Imputation phase — block-level design

For each LD block during the imputation pass:

1. Extract the pleiodb z-scores for variants in this block (`z_pleiodb`, shape N_pleiodb_in_block × T) from the stored V×T matrix.
2. Identify `missing_traits`: the set of trait indices j where any pleiodb variant in this block has a NaN z-score.
3. If `missing_traits` is empty, skip this block entirely.
4. For each trait j in `missing_traits`, read the VCF region `{chrom}:{block_start}-{block_end}` from `vcf_paths[j]`, obtaining z-scores keyed by allele ID for all variants in the region.
5. Match VCF allele IDs to LD panel row indices using the block's TSV (which lists LD panel variants with their positions and alleles). Assemble `z_dense` (shape N_ld_block × |missing_traits|).
6. For each trait j in `missing_traits`:
   - Training set: finite rows of `z_dense[:, j]` (LD panel variants observed in the VCF).
   - Fit elastic net on LD eigenvectors × training z-scores.
   - Predict for all N_ld_block positions.
   - Select predictions at the LD row indices corresponding to missing pleiodb positions.
   - Apply `min_cor` filter; if correlation too low, skip this trait for this block.
7. Return fills (global pleiodb variant indices, trait indices, fill z-scores).

### Parallelism — nested two-level

- **Outer level (ProcessPoolExecutor):** blocks are independent, dispatched to worker processes. Each worker process owns one block at a time. `workers` parameter controls the pool size.
- **Inner level (ThreadPoolExecutor within each worker):** VCF region reads for the `missing_traits` of a block are dispatched to threads within the worker process. `vcf_threads` parameter (default 8) controls thread count per worker. bcftools subprocess calls release the GIL, so threading is effective here.
- The LD matrix is loaded and decomposed once per block per worker process (unchanged from current design).

### VCF region reading — new function

- New function `read_vcf_region(vcf_path, chrom, start, end, vcf_build, target_build)` in `vcf.py`.
  - Uses `bcftools view -r {chrom}:{start}-{end}` to extract the region.
  - Returns a dict `{allele_id → z_score}` for all variants found (no filtering to pleiodb list).
  - Allele ID format matches the LD panel TSV (`CHR:POS_OA_EA` or similar).
  - When `vcf_build != target_build`, lifts the block coordinates from hg38 to the VCF build before querying, using the existing liftover machinery.
- The existing `read_vcf` function (index-based, returns aligned arrays) is unchanged.

### LD panel variant matching within block

- The LD block TSV (`{block}.tsv`) contains columns `CHR`, `SNP`, `OA`, `EA`, `BP`. This is already read by `build_block_index` to get `ld_row_indices`. 
- During the imputation phase, the TSV is re-read to build a `{chrom:pos_oa_ea → ld_row_index}` lookup (both allele orderings, for flip handling). This lookup is used to match VCF variants returned by `read_vcf_region` to their row in the LD matrix.
- Variants returned by the VCF but not in the LD panel TSV are ignored (they have no LD row index).

### Liftover for region queries

- Block boundaries (`{start}-{end}`) are in hg38 (the LD panel coordinate system).
- When a VCF has `vcf_build != hg38`, liftover the block start/end coordinates from hg38 to the VCF's build before issuing the `bcftools view -r` call.
- A simple single-interval liftover is sufficient (not per-variant); the region just needs to encompass all LD panel variants in the block.

### No scratch files

- The imputation phase reads VCF regions on-demand; nothing is written to disk beyond the final updated V×T matrix (same atomic replace as the current post-build pass).
- No changes to the database directory structure.

### CLI

- Add `--ld-vcf-threads INTEGER` to `pleiodb build` (default 8): threads per block worker for parallel VCF reads.
- The existing `--ld-workers` (ProcessPoolExecutor block workers) is unchanged in meaning.

## Testing Decisions

Good tests assert external outcomes, not internal mechanics.

- **Dense training improves fit**: build a small database with a toy LD panel and fixtures where the pleiodb variant list covers <5% of a block's LD variants; the remaining 95% have z-scores in the VCF fixtures. Assert that dense imputation produces Pearson r > 0.8 between imputed and held-out true values, while sparse-only (current method) produces r < 0.5 for the same positions.
- **Observed z-scores unchanged**: assert that pleiodb positions that had a finite z-score before imputation still have the same value after dense imputation (no overwrite of observed data).
- **Only missing traits queried**: assert (via mock on `read_vcf_region`) that VCF reads are only issued for traits where `missing_mask` is non-empty for that block. No reads for fully-observed traits.
- **VCF path round-trip**: build a database and assert that `traits.tsv` inside the output directory contains the correct `vcf_path` and `vcf_build` columns.
- **Graceful VCF failure**: when `read_vcf_region` raises an exception, assert the block falls back to the existing sparse z-scores (no crash, logged warning).
- **Liftover applied**: when `vcf_build=hg19` and block is in hg38, assert the region coordinates passed to bcftools are in hg19 (verifiable via mock on `read_vcf_region`).

Prior art: `tests/test_integration.py` — `TestPostBuildImputation` provides the fixture-based imputation test pattern.

## Out of Scope

- Storing dense LD z-scores permanently in the database.
- Re-imputing an already-built database without access to the original VCF files.
- Using an external reference z-score panel (e.g. UK Biobank) instead of the study VCFs.
- Imputing variants outside the 95k pleiodb set (final matrix stays V×T).
- Changing the LD panel, block structure, or eigenvector decomposition.

## Further Notes

- **Scale**: with 256 cores as outer ProcessPoolExecutor workers and 8 VCF threads per worker, the effective parallelism is 256 blocks × 8 VCF reads = 2,048 concurrent VCF region queries. For T=4,159 traits and 1,345 blocks, estimated wall time for the imputation phase is under 20 minutes.
- **Coverage impact**: median training set grows from 58 to ~5,000 variants per block (86× more signal). Blocks currently skipped due to sparse `min_cor` failure should largely pass.
- **Memory**: peak per block worker = N_ld_block × |missing_traits| × 4 bytes. For N_ld=7,057 and T=4,159: 117 MB per worker process. With 256 workers: ~30 GB total — acceptable on a high-memory server.
- **bcftools requirement**: region queries require a CSI or TBI index on every VCF. Build already requires indexed VCFs; imputation inherits the same requirement.
