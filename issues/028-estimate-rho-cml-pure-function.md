## Parent PRD

`issues/prd-rho-matrix.md`

## What to build

A pure, stateless Python function `estimate_rho_cml(z_j, z_k, z_thresh, min_nulls=500)` that
implements the conditional maximum-likelihood estimator for the sample-overlap-weighted
phenotypic correlation between two traits.

The function receives two float32 arrays of z-scores that have already been filtered to
`|z_j| < z_thresh AND |z_k| < z_thresh` (the "null variants" for this pair), plus the
threshold value and a minimum-count guard.  It returns a scalar rho estimate in (−1, 1),
or `NaN` when fewer than `min_nulls` null variants are available.

Internally the objective follows Forde et al. (`mr.simss::est_lambda`):

```
log L(ρ) ∝  −(A − 2ρB + C) / 2(1−ρ²)
            − n · log P(|X| < z, |Y| < z ; ρ)
```

where `A = Σz_j²`, `B = Σz_j·z_k`, `C = Σz_k²`, `n = len(z_j)`, and the
normalising constant `P(|X| < z, |Y| < z ; ρ)` is the bivariate normal probability
mass in the truncation square.  The normalising constant is precomputed on a fine
grid of ρ values (≥ 2000 points in (−1, 1)) once per call and used for interpolation
during scalar minimisation.

The `min_nulls` parameter defaults to 500 but must be settable to smaller values
(e.g. 50) so that integration tests that use the ~500-variant test dataset can
exercise the estimator without always hitting the NaN path.

See PRD §"CML estimator" and §"Modules to build or modify" for the full specification.

## Acceptance criteria

- [ ] `estimate_rho_cml(z_j, z_k, z_thresh)` is importable from `pleiodb` (or a
      submodule) as a public function.
- [ ] Given z-scores drawn from independent (ρ = 0) bivariate normals and filtered to
      `|z| < z_thresh`, the estimate is within ±0.05 of 0.
- [ ] Given z-scores drawn from a bivariate normal with ρ = 0.5 and filtered to
      `|z| < z_thresh`, the estimate is within ±0.05 of 0.5.
- [ ] Given z-scores drawn from a bivariate normal with ρ = −0.3 and filtered to
      `|z| < z_thresh`, the estimate is within ±0.05 of −0.3.
- [ ] When `len(z_j) < min_nulls`, the function returns `NaN` without raising.
- [ ] The `min_nulls` parameter defaults to 500 and can be overridden (e.g.
      `min_nulls=50`) to allow tests with small datasets to exercise the estimator.
- [ ] Different `z_thresh` values (0.5, 1.0, 2.0) give broadly consistent estimates
      on the same underlying data.
- [ ] The function is pure (no side-effects, no global mutable state) and can be
      called concurrently from multiple threads.

## Blocked by

None — can start immediately.

## User stories addressed

- User story 2 (unbiased estimator)
- User story 13 (standalone Python function)
- User story 15 (normalising constant precomputed on a grid)
