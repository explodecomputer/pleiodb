# ADR 0006 — rho: CML estimator replaces lambda/Pearson correlation

## Status
Accepted

## Context

pleiodb stores a T×T matrix of pairwise trait correlations that captures the
combined effect of sample overlap and phenotypic correlation between studies.
This quantity is used downstream for z-score correction, multivariate models
(SEM, factor analysis), pleiotropy tests (GSMR, MR-RAPS), and conditional
F-statistics.

The original implementation (`build_lambda`) computed this matrix as the
**Pearson correlation of z-scores** across variants with |z| < 3.0 in both
traits, and stored the result in `lambda.bin/.cidx`.

Two problems were identified:

1. **Name collision**: "lambda" (λ) standardly denotes the per-trait genomic
   inflation factor (median(χ²)/0.456) in the GWAS literature. Using the same
   symbol for a T×T cross-trait matrix creates ambiguity.

2. **Estimator bias**: selecting variants on |z| < threshold and then computing
   Pearson correlation on that selected subset produces a biased estimate. The
   selection induces a truncated bivariate distribution; the raw Pearson r of
   truncated data is not an unbiased estimator of the underlying correlation.

## Decision

### Rename: lambda → rho

The matrix is renamed to `rho`. The canonical term is
**"sample-overlap-weighted phenotypic correlation"**, abbreviated `rho`.
File names change from `lambda.bin/.cidx` to `rho.bin/.cidx`.
The CLI command changes from `pleiodb lambda` to `pleiodb rho`.

### Estimator: conditional maximum likelihood (CML)

Replace the Pearson estimator with the conditional maximum likelihood
estimator from Forde et al. (`mr.simss::est_lambda`). Given null variants
(|z_j| < z_thresh AND |z_k| < z_thresh), the CML estimator maximises:

```
log L(ρ | data) ∝  −(A − 2ρB + C) / (2(1−ρ²))
                   − n · log P(|X| < z, |Y| < z ; ρ)
```

where A = Σz_j², B = Σz_j·z_k, C = Σz_k², n = number of null variants, and
the denominator is the probability mass of a bivariate normal with correlation
ρ falling in the square [−z, z]². This denominator is a function of ρ and
z_thresh only (not the data) and is precomputed on a grid to avoid repeated
2D integration.

### z-threshold: 1.0

The default null-variant threshold is set to z_thresh = 1.0. The original
code used 3.0 (retaining ~99.7% of variants); the reference implementation
uses 0.5 (too aggressive). z = 1.0 retains ~68% of variants per trait and
~46% per pair, giving the CML estimator sufficient data while excluding
most true associations. Accuracy requirements for this matrix are moderate
(it is an approximation used in downstream corrections, not a primary result).

### Variant basis: existing z-score matrix (Approach A)

The rho computation reads from the already-stored z-score matrix rather than
re-extracting variants from VCF files. The database variant list is typically
enriched for associations (variants known to associate with at least one trait
in the matrix), but for any specific pair (j, k) the majority of variants will
have low z in both traits, providing an adequate null set.

An alternative (Approach B) — extracting ~1M HapMap3 SNPs independently from
VCF files — was considered but rejected as unnecessarily complex given the
moderate accuracy requirement.

### Separate command, not folded into build

`pleiodb rho` remains a separate post-build step (not run automatically by
`pleiodb build`). This allows recomputation with different thresholds without
rebuilding the database, and avoids surprising users with a slow step at the
end of a long build. `pleiodb info` warns when the rho matrix is absent.

### Query interface

The `pleiodb rho` command doubles as a query interface. Presence of
`--traits` (comma-separated trait IDs) or `--traits-file` (one trait ID per
line) puts the command into query mode; absence triggers matrix computation.

Default query output is a **pairwise list**:

```
trait_id_1      trait_id_2      rho
ukb-b-19953     ieu-a-7         0.0412
ukb-b-19953     ukb-b-10787     0.1823
ieu-a-7         ukb-b-10787     0.0091
```

The `--matrix` flag pivots to a square matrix with trait IDs as row/column
headers — convenience for small T and direct input to R/Python matrix methods.

Trait IDs must not contain commas (commas delimit the `--traits` list).
This is already satisfied by all known trait ID conventions (e.g. `ieu-a-7`).

## Consequences

- `lambda.bin/.cidx` is replaced by `rho.bin/.cidx`; existing databases built
  with the old code will not have a `rho` matrix (they will trigger the
  `pleiodb info` warning).
- The CML estimator requires `scipy.optimize.minimize_scalar` and the bivariate
  normal CDF (`scipy.stats.multivariate_normal.cdf`); both are already
  dependencies.
- rho[j, j] is set explicitly to 1.0 rather than estimated.
- A minimum of 500 null variants is required per pair; pairs falling below this
  threshold are stored as NaN.
- The matrix is enforced to be symmetric: only the upper triangle is computed
  and mirrored.
