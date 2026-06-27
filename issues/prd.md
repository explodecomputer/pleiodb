# PRD: LD-Based Z-Score Imputation During Ingestion

## Problem Statement

When building a pleiodb database from thousands of GWAS-VCF files, a significant fraction of
(variant, trait) cells end up with NaN z-scores — not because the association is absent, but
because the GWAS VCF simply did not test that variant. In the test build this reached 27%
missingness. These gaps reduce the utility of pleiotropy queries and downstream analyses that
rely on complete coverage across traits.

The root cause is that different GWAS studies test different sets of variants (depending on
the genotyping array, imputation panel, and quality-control choices used by each study). A
variant present in the pleiodb variant list may be absent from many VCFs even though correlated
neighbouring variants in those VCFs carry information about the missing association.

## Solution

Integrate optional LD-based z-score imputation into the `pleiodb build` pipeline. After parsing
each batch of GWAS VCFs and obtaining the observed z-score matrix, use a pre-built LD reference
panel (organised by LD block) to predict z-scores for missing variants using an elastic-net
regression on the eigenvectors of the local LD matrix — the same statistical approach used in
the genotype-phenotype-map pipeline (`imputation_method.R`).

Imputed cells are stored identically to observed cells in the main z-score binary matrix, but
their positions are additionally recorded in a separate sparse flag file so downstream consumers
can distinguish observed from imputed values if they choose to.

Imputation is opt-in: users who do not supply `--ld-dir` get the current behaviour unchanged.

## User Stories

1. As a database builder, I want to pass a `--ld-dir` flag to `pleiodb build` pointing at a
   pre-built LD reference panel so that missing z-scores are automatically imputed during
   ingestion without a separate post-processing step.

2. As a database builder, I want the imputation to use the same LD reference data that the
   genotype-phenotype-map pipeline uses so that I do not need to maintain a second reference
   dataset.

3. As a database builder, I want to control the minimum block-level imputation quality (via
   `--ld-min-cor`) so that only well-imputed blocks contribute, and blocks with poor LD coverage
   or few observed variants are silently skipped.

4. As a database builder, I want to tune the fraction of LD variance retained in the PCA step
   (via `--ld-thresh`) to balance imputation completeness against overfitting.

5. As a database builder, I want a progress log that tells me how many blocks were processed,
   how many were skipped (too few variants or low correlation), and how many z-scores were
   imputed per trait batch so I can monitor the run and diagnose problems.

6. As a database builder, I want the build to succeed even if the LD reference panel does not
   cover some chromosomes or blocks, with those blocks simply skipped, so that a partial LD
   panel does not abort the entire build.

7. As a database builder, I want a clear warning (not an error) if `--ld-dir` is supplied but
   `--variants-build hg38` is not set, because the LD reference is in hg38 coordinates.

8. As a database analyst, I want the imputed z-scores stored at the same precision and encoding
   as observed z-scores so that all downstream queries work without modification.

9. As a database analyst, I want a separate sparse file (`imputed.coo.zst`) that records exactly
   which (variant, trait) cells were imputed so I can filter them out if my analysis requires
   only observed associations.

10. As a database analyst, I want the `n_variants` column in `traits.tsv` to count only
    *observed* (non-imputed) variants so I can compare coverage before and after imputation.

11. As a database analyst, I want an additional `n_variants_imputed` column in `traits.tsv`
    that reports how many z-scores were imputed per trait so I can assess per-study imputation
    yield.

12. As a developer, I want the imputation logic isolated in a dedicated module (`impute.py`)
    with a clean public API so that I can test imputation independently of the ingestion
    pipeline and later swap or upgrade the statistical method.

13. As a developer, I want the core elastic-net imputation function to accept plain numpy arrays
    (z-scores, LD submatrix, EAF) and return imputed z-scores plus a quality metric, so that it
    can be unit-tested with small synthetic inputs without touching the filesystem.

14. As a developer, I want the LD block index to be built once at the start of the build and
    reused across all trait batches so that the filesystem is not scanned repeatedly during the
    main ingestion loop.

15. As a developer, I want `scikit-learn` declared as an explicit package dependency so that
    users who install pleiodb get the imputation functionality without manual dependency
    installation.

## Implementation Decisions

### New module: impute.py

A standalone module in `src/pleiodb/` containing all imputation logic. Its public interface:

- **`build_block_index(variants, ld_dir, ancestry)`** — scans the LD panel directory tree once,
  matches pleiodb variant ALIDs to LD-panel SNP IDs (which are already in canonical
  alphabetical allele order matching the ALID format), and returns a mapping from block path to
  the paired lists of pleiodb variant indices and LD-matrix row indices. Only blocks with at
  least 2 matched variants are included.

- **`impute_z_block(z_block, variants, eaf_arr, block_index, thresh, min_cor, out_mask)`** —
  iterates over all blocks in the index, runs the elastic-net routine for each block × trait
  column, writes imputed values back into `z_block` in-place, and sets `True` in `out_mask`
  for every cell that was filled.

Internal helpers (not part of the public API, but extractable for testing):

- `_load_ld_submatrix(block_dir, ld_row_indices)` — reads `.unphased.vcor1.gz`, subsets
  rows/columns to the matched variant indices.
- `_ld_pca(ld_sub, thresh)` — calls `scipy.linalg.eigh`, flips to descending eigenvalue order,
  selects the minimum number of components whose cumulative variance fraction reaches `thresh`.
- `_elastic_net_impute(z, eigenvectors, n_comp)` — fits `ElasticNetCV(l1_ratio=0.5,
  fit_intercept=False, cv=5)` on the observed (non-NaN) z-scores using the first `n_comp`
  eigenvectors as features; predicts for all positions including missing ones; returns the raw
  predictions.
- `_poly_rescale(truth, predicted, npoly=3)` — removes Cook's-distance outliers from the
  observed (truth, predicted) pairs, fits a degree-`npoly` polynomial with `np.polyfit`,
  evaluates it on all predicted values; returns the rescaled array and the Pearson correlation
  between rescaled and observed values.
- `_se_outliers(se_obs, se_hat, outthresh=3)` — computes Cook's distance from a linear
  regression of `se_obs ~ se_hat` and flags positions beyond `outthresh` standard deviations;
  these positions are added to the missing set before imputation.

### SE-based outlier smoothing

Before fitting the elastic net, theoretical SE values are computed from EAF
(`sehat = 1 / sqrt(2 * af * (1 - af))`). Positions whose observed SE deviates from `sehat`
beyond the Cook's distance threshold are treated as missing and also imputed. This mirrors the
R implementation and prevents unstable BETAs from corrupting the imputation features.

### Quality filter

After polynomial rescaling, the Pearson correlation between the rescaled predictions and the
held-out observed z-scores is computed. If this correlation falls below `min_cor` (default 0.7),
no values are written for that block × trait combination and the block is skipped silently with
a debug-level log message.

### Integration into build.py

`build_database()` gains four new optional parameters: `ld_dir`, `ld_ancestry`, `ld_thresh`,
`ld_min_cor`. The block index is constructed once before the trait loop. Inside each trait
batch, imputation runs on the completed `z_block` (after VCF parsing, before encoding). A
running list of COO pairs (v_idx, t_idx) is accumulated across batches and written as
`imputed.coo.zst` after the main loop finishes. The `n_variants_arr` counter (used for
`n_variants` in `traits.tsv`) is updated only from observed z-scores; a separate
`n_imputed_arr` accumulates imputed counts.

### Neff for imputed variants

The `neff_block` is derived from observed SE values. For imputed positions (where SE is NaN),
Neff remains NaN. This is consistent with current behaviour and avoids fabricating uncertainty
estimates for imputed associations.

### Imputed mask storage

`imputed.coo.zst` is stored at the root of the `.pleiodb` output directory alongside
`zscore.bin` and `neff.bin`. Format: a zstd-compressed sequence of `uint32` pairs
`(v_idx, t_idx)`, sorted by `(v_idx, t_idx)`, identical in structure to the significance
masks in `masks/`. An empty file (0 pairs) is written when imputation is not used, so code
that reads the file can always expect it to exist after a build.

### CLI changes

Four new options added to the `build` sub-command:
- `--ld-dir PATH`: path to the LD reference panel root (e.g. `.../ld_reference_panel_hg38/EUR`)
- `--ld-ancestry TEXT` (default `EUR`)
- `--ld-thresh FLOAT` (default `0.9`): PCA variance threshold
- `--ld-min-cor FLOAT` (default `0.7`): minimum imputation correlation to accept a block

### Dependency

`scikit-learn` is added to the `dependencies` list in `pyproject.toml`. `numpy` and `scipy`
are already present.

### Coordinate build assumption

The LD reference panel uses hg38 coordinates. Pleiodb variant ALIDs encode the coordinate
of the variant list. If `--ld-dir` is supplied without `--variants-build hg38`, the build
emits a `WARNING` but continues; matching is attempted using the coordinates as-is.

## Testing Decisions

### What makes a good test

Tests should verify external behaviour through the public API — inputs in, outputs out — not
implementation details like internal loop structure or intermediate variable names. Tests for
the imputation module should work with small synthetic LD matrices and z-score arrays so they
run without filesystem access.

### Modules to test

**`impute.py` (unit tests, no filesystem)**

- `_poly_rescale`: given a known truth and predicted array, verify the output correlation and
  that rescaled values track truth closely; verify outlier removal doesn't crash on degenerate
  input.
- `_ld_pca`: given a small random positive-definite matrix, verify eigenvalues are in
  descending order and the selected number of components matches the `thresh` parameter.
- `_elastic_net_impute`: given a small synthetic dataset with some z-scores masked, verify
  that the function returns a full-length array with no NaNs (predictions for all positions).
- `impute_z_block` (integration-style, still no real LD files): construct a tiny block index
  pointing at an in-memory LD submatrix fixture; verify that NaN cells in `z_block` are filled
  and that the corresponding positions in `out_mask` are True.

**`build.py` integration tests (filesystem, existing pattern)**

The existing `tests/test_integration.py` tests build a small pleiodb from VCF fixtures and
inspect the resulting files. Add a test that:

1. Runs `build_database()` with a minimal synthetic LD reference directory (two small blocks,
   a handful of variants each) and a traits fixture that has deliberate gaps.
2. Verifies that `imputed.coo.zst` exists and contains the expected imputed (v, t) pairs.
3. Verifies that `traits.tsv` contains both `n_variants` (observed only) and
   `n_variants_imputed` columns with correct counts.
4. Verifies that the imputed z-scores round-trip through the quantize/decode path without
   corruption.

Prior art: `tests/test_integration.py` already builds a full database and asserts on
`traits.tsv` column values, `zscore.bin` decoded values, and significance masks — follow the
same fixture-building and assertion style.

## Out of Scope

- **hg19 LD reference panel**: only hg38 LD matrices are currently available; imputation for
  studies whose variant list is in hg19 is not supported in this version.
- **Per-variant imputation R² filtering**: the elastic-net method provides a block-level
  quality metric (correlation), not per-variant R². Per-variant filtering (as in the RAISS
  method) is not implemented.
- **Imputing Neff / SE for imputed variants**: only z-scores are imputed; effective sample
  size for imputed positions remains NaN.
- **Ancestries other than EUR**: the LD panel directory structure supports multiple ancestry
  subdirectories, but only EUR is validated and tested in this PRD.
- **Query-time imputation flag exposure**: the `pleiodb query` command does not gain any new
  filter for imputed vs. observed; that is left to future work once the storage format is
  proven.
- **Re-imputing an existing build**: there is no `pleiodb impute` sub-command to retrofit
  imputation onto a previously built database; imputation is build-time only.

## Further Notes

- The LD reference panel is already in use by the genotype-phenotype-map pipeline on the same
  machine at `/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/`.
  No data transfer or reformatting is required.
- LD panel SNP IDs (`chr:bp_A1_A2` with alphabetical allele order) already match the pleiodb
  ALID format, enabling direct string comparison for variant matching.
- The R implementation uses `glmnet` with `alpha=0.5`; the Python port uses `sklearn`'s
  `ElasticNetCV` with `l1_ratio=0.5`, which is the exact equivalent parameter mapping.
- The block-level PCA is computed on the fly from the `.unphased.vcor1.gz` matrix using
  `scipy.linalg.eigh`; the pre-computed `.ldeig.rds` R objects are not used, avoiding an
  `rpy2` dependency.
- Imputation is CPU-only; the GPU-accelerated RAISS implementation in the Python pipeline
  (`cupy`) is not used here.
