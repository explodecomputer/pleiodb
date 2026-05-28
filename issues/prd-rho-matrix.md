# PRD: rho — sample-overlap-weighted phenotypic correlation matrix

## Problem Statement

pleiodb currently computes a T×T pairwise matrix via `pleiodb lambda` that is
intended to capture sample overlap and phenotypic correlation between studies.
This matrix has two problems:

1. **Name collision**: "lambda" (λ) standardly denotes the per-trait genomic
   inflation factor (median(χ²)/0.456) in the GWAS literature. Naming the
   cross-trait matrix "lambda" creates ambiguity for anyone familiar with LDSC
   or standard GWAS QC terminology.

2. **Biased estimator**: the current implementation uses Pearson correlation of
   z-scores filtered to |z| < 3.0. Selecting variants on their z-score and
   then computing a correlation from that selected subset produces a biased
   estimate — the selection induces a truncated bivariate distribution, and the
   raw Pearson r of truncated data is not an unbiased estimator of the
   underlying correlation.

Additionally, there is no query interface: users cannot extract subsets of the
matrix from the command line.

## Solution

Replace the `lambda` matrix and `pleiodb lambda` command with a `rho` matrix
and `pleiodb rho` command that:

1. Renames the concept and storage files from `lambda` to `rho`.
2. Replaces the Pearson estimator with a conditional maximum-likelihood (CML)
   estimator that corrects for the truncation bias introduced by the z-score
   threshold.
3. Adds a query interface to extract subsets of the T×T rho matrix as either
   a pairwise list or a square matrix.

The quantity stored is `rho[j,k]` — an estimate of
`ρ_pheno × N_overlap[j,k] / √(N_j × N_k)`, where `ρ_pheno` is the
phenotypic correlation between traits j and k and `N_overlap[j,k]` is the
number of shared samples. This is a single combined quantity; the two
components cannot be separated without external information.

## User Stories

1. As a pleiodb user, I want the cross-trait matrix to be called `rho` rather
   than `lambda`, so that it does not conflict with the standard genomic
   inflation factor.
2. As a pleiodb user, I want the rho matrix to use an unbiased estimator, so
   that downstream corrections (z-score adjustment, conditional F-statistics,
   MR methods) are not systematically wrong.
3. As a pleiodb user, I want to compute the rho matrix as a separate post-build
   step, so that I can rebuild the database without recomputing rho, and
   recompute rho with different settings without rebuilding the database.
4. As a pleiodb user, I want `pleiodb info` to warn me when the rho matrix is
   absent from a database, so that I do not silently use a database that is
   missing a required component.
5. As a pleiodb user, I want to query rho for a specific set of traits by
   passing a comma-separated list on the command line, so that I can quickly
   extract values for a small number of traits without writing a file.
6. As a pleiodb user, I want to query rho for a specific set of traits by
   passing a file of trait IDs (one per line), so that I can extract values for
   a large number of traits programmatically.
7. As a pleiodb user, I want the default rho query output to be a pairwise list
   with columns `trait_id_1`, `trait_id_2`, `rho`, so that I can pipe it into
   tabular tools and join it with other outputs.
8. As a pleiodb user, I want a `--matrix` flag on rho queries, so that I can
   get a square matrix directly suitable for input to R or Python matrix
   operations.
9. As a pleiodb user, I want rho[t,t] to always be 1.0, so that the diagonal
   is correct without needing to be estimated.
10. As a pleiodb user, I want trait pairs with fewer than 500 shared null
    variants to return NaN rather than a noisy estimate, so that downstream
    methods do not receive unreliable values.
11. As a pleiodb user, I want the rho computation to be parallelised across
    trait pairs, so that it completes in reasonable time for large T.
12. As a pleiodb user, I want the rho output to support both TSV and JSON
    formats via a `--format` flag, consistent with `pleiodb query`.
13. As a developer building on pleiodb, I want the CML estimator to be exposed
    as a standalone Python function, so that I can call it directly with my own
    z-score arrays.
14. As a developer, I want the rho matrix to be stored as a symmetric float16
    chunked matrix at `rho.bin/.cidx`, so that it is compact and consistent
    with the z-score matrix storage format.
15. As a developer, I want the normalizing constant of the CML objective (the
    bivariate normal CDF over the truncation square) to be precomputed on a
    grid of rho values once per computation run, so that repeated evaluations
    during optimization are fast.

## Implementation Decisions

### Renamed concept and files

- The concept "lambda" is retired; the canonical term is **rho**.
- Storage files change from `lambda.bin/.cidx` to `rho.bin/.cidx`.
- The meta.json keys `lambda_chunk_shape` and `lambda_null_z_thresh` become
  `rho_chunk_shape` and `rho_null_z_thresh`.
- The `pleiodb lambda` CLI command becomes `pleiodb rho`.

### CML estimator

The estimator follows Forde et al. (`mr.simss::est_lambda`). Given variants
where |z_j| < z_thresh AND |z_k| < z_thresh, the objective is:

```
log L(ρ) ∝  −(A − 2ρB + C) / 2(1−ρ²)
            − n · log P(|X| < z, |Y| < z ; ρ)
```

where A = Σz_j², B = Σz_j·z_k, C = Σz_k², n = count of null variants, and
the denominator is the bivariate normal probability mass in the truncation
square. Optimization is over ρ ∈ (−1, 1) using scalar minimisation.

The normalizing constant `P(|X| < z, |Y| < z ; ρ)` is computed via the
bivariate normal CDF (scipy), which is exact and does not require numerical
2D integration. It is precomputed once on a fine grid of ρ values (e.g. 2000
points in (−1, 1)) and evaluated by interpolation during optimization.

The z-score threshold defaults to **1.0**. This is lower than the current
code's 3.0 but the CML estimator handles the truncation bias regardless of
threshold.

### Variant source

Rho is computed from the z-scores already stored in the database. No VCF files
are re-read. The database variant list is typically enriched for associations
(variants known to associate with at least one trait), but for any specific
pair (j, k) the majority of variants will have |z_j| < 1 AND |z_k| < 1 since
each variant typically associates with only a few traits.

### Symmetry and diagonal

Only the upper triangle of the T×T matrix is computed; values are mirrored to
fill the lower triangle before writing. The diagonal rho[t,t] is set to 1.0
explicitly without estimation.

### Dual-mode `pleiodb rho` command

The `pleiodb rho <db>` command has two modes determined by flags:

| Flags | Mode |
|---|---|
| No `--traits` or `--traits-file` | **Compute** full T×T matrix |
| `--traits t1,t2,...` or `--traits-file path` | **Query** matrix |

Compute mode accepts `--workers`, `--null-thresh` (z_thresh, default 1.0),
and `--chunk-size`.

Query mode accepts `--traits` (comma-separated), `--traits-file`, `--format`
(tsv/json), `--matrix` (pivot to square output), and `--output`.

### Query output format

**Default (pairwise list)**:
```
trait_id_1   trait_id_2   rho
t1           t2           0.0412
t1           t3           0.1823
t2           t3           0.0091
```
Only the upper triangle is emitted (unordered pairs, no duplicates).

**`--matrix` flag**:
```
          t1      t2      t3
t1        1.0     0.0412  0.1823
t2        0.0412  1.0     0.0091
t3        0.1823  0.0091  1.0
```

### Trait ID constraint

Trait IDs must not contain commas. This is satisfied by all known OpenGWAS
trait ID conventions and is documented in CONTEXT.md.

### Modules to build or modify

**`estimate_rho_cml(z_j, z_k, z_thresh)` — new deep module**
Pure function. Takes two float32 arrays of z-scores (already filtered to
|z_j| < z_thresh AND |z_k| < z_thresh), returns a scalar rho estimate.
Internally uses a precomputed normalizing constant grid. Stateless and
independently testable with synthetic bivariate normal samples.

**`build_rho(db_path, ...)` — replaces `build_lambda`**
Reads the stored z-score matrix chunk by chunk, applies the z-threshold
filter per pair, calls `estimate_rho_cml`, and writes the result to
`rho.bin/.cidx`. Parallelised over trait pairs.

**`GWASDatabase` — minor update**
Rename `lambda_matrix` → `rho_matrix` and `get_lambda_block` → `get_rho_block`.
Update `info()` to report rho matrix presence and emit a warning key when absent.

**`pleiodb rho` CLI command — replaces `pleiodb lambda`**
Dual-mode command. Compute mode: thin wrapper over `build_rho`. Query mode:
reads `rho_matrix`, filters to requested traits, formats output.

## Testing Decisions

Good tests verify observable behaviour through public interfaces only. They
should survive internal refactors; if a test breaks when you rename a variable
or restructure a function body without changing behaviour, it was testing
implementation.

**What makes a good test here**:
- Give `estimate_rho_cml` synthetic z-score pairs drawn from a known bivariate
  normal with a known ρ and assert the estimate is within a reasonable
  tolerance. Do not inspect intermediate values (A, B, C, grid).
- Test the CLI query output by parsing the TSV/JSON, not by inspecting internal
  matrix objects.
- Test NaN behaviour by constructing a pair with fewer than 500 null variants
  and asserting the output is NaN.

**Modules with tests**:

1. `estimate_rho_cml` — unit tests:
   - Known ρ = 0: z-scores from independent normals → estimate ≈ 0
   - Known ρ = 0.5: z-scores from correlated normals → estimate within ±0.05
   - Known ρ = −0.3: negative correlation recovered
   - n < 500: returns NaN
   - z_thresh parameter: different thresholds give consistent estimates

2. `build_rho` + `GWASDatabase.rho_matrix` — integration test via the existing
   `_build()` / `_open()` test helpers (see `TestHg19Build` as prior art):
   - Build a test database, run `build_rho`, check `rho_matrix` is present
   - Diagonal values are 1.0
   - All off-diagonal finite values are in (−1, 1)
   - Symmetry: rho[j,k] == rho[k,j]
   - `pleiodb info` output contains a rho-present key

3. `pleiodb rho` query mode — CLI integration tests:
   - `--traits t1,t2` returns exactly one row (the single pair)
   - `--matrix` output has correct shape and diagonal = 1.0
   - `--traits-file` gives same result as equivalent `--traits`
   - `--format json` parses as valid JSON with correct keys

Prior art: `TestQueryOutput` and `TestQueryIntersect` in `tests/test_integration.py`
demonstrate the pattern for CLI query tests using `click.testing.CliRunner` and
`io.StringIO` capture.

## Out of Scope

- **Decomposing rho into separate ρ_pheno and N_overlap components**: requires
  external data (known sample rosters or LD score regression).
- **HapMap3-based estimation (Approach B)**: re-extracting ~1M SNPs from VCF
  files as an independent reference set. Deferred pending evidence that
  Approach A is insufficient.
- **On-demand / lazy rho computation**: all pairs are computed eagerly at build
  time. Lazy per-pair computation deferred pending performance evidence.
- **Per-variant rho**: the matrix is trait×trait only; per-variant sample
  overlap estimates are not in scope.
- **Liftover**: rho estimation assumes z-scores are already in the database;
  coordinate issues were resolved at build time.

## Further Notes

- The reference implementation is `mr.simss::est_lambda` (Forde et al.):
  https://github.com/amandaforde/mr.simss/blob/main/R/est_lambda.R
- Design decisions (naming rationale, estimator choice, Approach A vs B) are
  documented in `docs/adr/0006-rho-cml-estimator-replaces-lambda-pearson.md`.
- Domain terminology is in `CONTEXT.md` under "rho" and "rho query interface".
- Accuracy requirements are moderate — rho is an approximation used in
  downstream corrections, not a primary reported result.
