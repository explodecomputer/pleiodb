## Parent PRD

`issues/prd.md`

## What to build

Add `scikit-learn` to the package dependencies and create `src/pleiodb/impute.py` containing
the stateless mathematical core of the elastic-net imputation method. This slice has no
filesystem I/O — all functions accept and return plain numpy arrays — so it can be developed,
reviewed, and unit-tested in complete isolation from the ingestion pipeline.

The functions to implement mirror the internal helpers described in the PRD:

- **`_ld_pca(ld_matrix, thresh)`** — wraps `scipy.linalg.eigh`, flips eigenvalues/vectors to
  descending order, and returns `(eigenvalues, eigenvectors)` truncated to the minimum number
  of components whose cumulative variance fraction reaches `thresh`.

- **`_elastic_net_impute(z, eigenvectors, n_comp)`** — fits
  `sklearn.linear_model.ElasticNetCV(l1_ratio=0.5, fit_intercept=False, cv=5)` on the
  observed (non-NaN) z-scores using the first `n_comp` eigenvector columns as features,
  then predicts z-scores for *all* positions (observed and missing). Returns the full
  prediction array.

- **`_poly_rescale(truth, predicted, npoly=3)`** — removes Cook's-distance outliers from
  `(truth, predicted)` pairs (two-pass SD threshold), fits a degree-`npoly` polynomial
  through the origin via `np.polyfit`, evaluates it on all predicted values, and returns
  `(rescaled_array, pearson_correlation)` between rescaled predictions and observed truth.

- **`_se_outliers(se_obs, se_hat, outthresh=3)`** — fits `se_obs ~ se_hat` via OLS, computes
  Cook's distances, and returns a boolean mask of outlier positions using the same two-pass
  SD threshold as `_poly_rescale`.

These four helpers are the only deliverable of this issue. The public-facing `impute_z_block`
function is built in a later issue once the block index is also available.

## Acceptance criteria

- [ ] `scikit-learn` appears in the `dependencies` list of `pyproject.toml`.
- [ ] `src/pleiodb/impute.py` exists and imports without error.
- [ ] `_ld_pca`: given a small random positive-definite matrix and `thresh=0.9`, the returned
      eigenvalues are in strictly descending order and the number of components satisfies
      `cumsum(values) / sum(values) >= 0.9`.
- [ ] `_elastic_net_impute`: given a length-20 z-score array with 5 randomly masked NaN values
      and a matching synthetic eigenvector matrix, the returned array has length 20, no NaNs,
      and the values at non-masked positions are close to the values predicted by the fitted
      model (not just pass-through of the original observations).
- [ ] `_poly_rescale`: given truth and predicted arrays whose relationship is a known polynomial
      plus small noise, the returned rescaled values correlate with truth at r > 0.95; the
      function handles the degenerate case where all outlier removal would leave fewer than 2
      points (should return gracefully, not raise).
- [ ] `_se_outliers`: given arrays where one element is a clear outlier (10× the expected
      ratio), that element is flagged `True` and non-outlier elements are `False`.
- [ ] All four helpers are covered by unit tests in `tests/test_impute.py`; the tests use only
      synthetic numpy arrays and require no filesystem access.
- [ ] `pytest tests/test_impute.py` passes.

## Blocked by

None — can start immediately.

## User stories addressed

- User story 12 (imputation logic in a dedicated, testable module)
- User story 13 (core function accepts plain numpy arrays)
- User story 15 (scikit-learn as an explicit dependency)
