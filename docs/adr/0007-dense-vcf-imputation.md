# ADR-0007: Dense VCF imputation via precomputed LD eigenvectors

## Status

Accepted

## Context

The post-build imputation pass trains an elastic-net model on the z-scores of stored
variants that fall within each LD block.  The stored variant set covers only ~20–50 of
the ~5–10k reference-panel SNPs per block (0.4–1%).  Training on this sparse set gives
the model too little signal to reliably recover missing values.

The source GWAS-VCF files contain z-scores for every genotyped or imputed variant in the
study — often the full LD-panel variant set.  Reading those z-scores at build time (before
storing only the ~95k canonical positions) and using them as training data was the natural
fix, but the data were discarded before the imputation step ran.

## Decision

After the VCF ingestion loop completes, for each LD block the imputation worker:

1. Reads z-scores for **all** variants in the block region from the source GWAS-VCFs via
   `bcftools view -r` (with liftover when the VCF uses a different genome build).
2. Assembles a dense z-score vector of length N (all reference-panel SNPs in the block)
   using allele-orientation-aware matching.
3. Loads precomputed block eigenvectors from `.ldeig.rds` files that ship with the LD
   reference panel; caches them as `.ldeig.npz` after the first read.
4. Trains the elastic net on the dense z-vector and extracts predictions only at the
   stored (missing) positions.

The traits TSV gains two new optional columns — `vcf_path` and `build` — which are stored
in the database and supplied back to the imputation pass.  When these columns are absent
(old databases), the pass falls back to sparse mode using only stored variants.

## Eigenvector caching

Each `.ldeig.rds` file is ~662 MB (full N×N eigenvector matrix, R binary format).
Reading it via Rscript subprocess takes 3–5 s.  Only K ≈ 169–250 columns are needed
(90% variance threshold).  On first access the worker extracts and saves those columns as
a compressed `.ldeig.npz` cache (~7 MB); subsequent reads take ~62 ms.  The cache write
uses an atomic temp-file-then-rename pattern to avoid corruption under concurrent workers.

## Consequences

- Imputation quality improves because the model sees the full local LD signal.
- First-pass per-block cost increases slightly (VCF region reads via inner
  ThreadPoolExecutor), but the eigenvector cache amortises the RDS parse cost so overall
  wall-clock time is comparable: **0.71 s/block dense vs 0.94 s/block sparse**
  (chr1, 5 traits, 32 ProcessPoolExecutor workers, 8 VCF threads each).
- The `.ldeig.npz` cache files accumulate in the LD reference panel directory
  (~9 GB for 1345 blocks).  This is a one-time cost shared across all builds.
- `vcf_path` / `build` are optional — existing databases and builds without `--ld-dir`
  continue to work unchanged.
