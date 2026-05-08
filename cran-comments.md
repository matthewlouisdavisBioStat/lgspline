# CRAN Comments

## Resubmission

This is a resubmission of lgspline, now version 1.1.0. The revision addresses issues found during early JCGS review and package testing:

1. Inequality-constrained fitting now uses the active-set method across the main fitting paths, including correlated models.
2. Tuning avoids unnecessary parallel work on small vector operations and parallelizes the initial grid by default; BFGS trial-step parallelism is available as an opt-in.
3. Built-in inequality constraints can now be enforced on different observation subsets through the named-list form of `qp_observations`.
4. Smoothness constraints now support `add_first_and_second_derivative_constraints = NULL`, which chooses the default constraint layout from the spline-expanded predictors.

Documentation and comments were updated. The changes are fully backward compatible.

---

## Main Changes

### Inequality Constraints

- active-set fitting is now used for uncorrelated fits, correlated Gaussian GEE, correlated non-Gaussian GEE, and `blockfit_solve()`
- structured-correlation fits use the Woodbury form when it is efficient; otherwise they fall back to dense QP/SQP as before
- derivative constraints, active-set behavior, and fallback rules are now documented in the help files
- `qp_observations` accepts either the legacy shared vector or a named list such as `Time:qp_negative_derivative`
- equality and inequality QR reduction can now be controlled separately with `qr_pivot_smoothing_constraints` and `qr_pivot_inequality_constraints`
- `parallel_qr_qp` adds worker-level QR reduction for derivative and curvature inequality blocks; range and monotonicity constraints remain exact

### Tuning

- added outer-level tuning controls `parallel_grideval` and `parallel_bfgs`; `parallel_grideval` defaults to `TRUE`, while `parallel_bfgs` is opt-in
- tuning now avoids nested use of the same cluster by disabling inner parallel kernels inside outer tuning workers
- inner tuning-gradient algebra and main-process tuning evaluations are kept serial, so `parallel_eigen` does not trigger nested worker dispatch while grid or BFGS tuning is already using the cluster
- `stats::optim()` tuning now uses the same safe criterion-evaluation guard as grid search and custom BFGS
- tuning gradients now compute the derivative with respect to each partition penalty matrix once and reuse it for wiggle, flat ridge, predictor-specific, and partition-specific penalties
- default GLM score and working-weight functions now use the submitted family object's `mu.eta` and `variance`, so common non-canonical links work without extra user-supplied functions
- cached fixed correlated-tuning quantities to avoid rebuilding whitening objects during repeated LOO/GCV evaluations
- penalty tuning now warm-starts constrained coefficient solves from the last accepted active set
- relaxed the structured-correlation Woodbury rank gate from `P/3` to `2P/3` to avoid unnecessary dense fallback

### Other Updates

- changed the small-sample tuning adjustment from `[(N + 2)/(N - 2)]^2` to `(N + 1)/(N - 1)`
- small C++ helper updates improve the symmetric-positive-definite inversion path, return exactly symmetric Gram matrices, and skip structurally zero constraint columns in the parallel `A'GA` chunk helper
- custom covariance helpers can now use `correlation_id`, `spacetime`, and matching extra arguments while fitted objects store the one-argument wrappers used by posterior methods
- if `VhalfInv_fxn` is supplied and `VhalfInv_par_init` is omitted, optimization now starts at `1e-2`
- `logLik.lgspline()` now accepts fixed coefficients through `B_predict` and fixed dispersions through `sigmasq_predict`
- `prior_loglik()` now uses the same `B_predict` and `sigmasq_predict` convention
- updated examples for nested re-fitting and Toeplitz correlation intervals
- removed development-only benchmark artifacts from the submitted source tree
- added focused tests for keyed `qp_observations`, active-set/correlation behavior, QR controls, tuning, non-canonical GLM defaults, and custom likelihood inputs

## Verification

- package reinstall from source completed successfully
- all tests passing

## R CMD check results

0 errors | 0 warnings | 0 notes

There are no downstream dependencies for this package.
