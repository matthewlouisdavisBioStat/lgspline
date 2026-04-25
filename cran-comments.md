# CRAN Comments

## Resubmission

This is a resubmission, lgspline version 1.0.3. Following some early feedback upon JCGS submission, 3 major issues were resolved:
1) The active set method for imposing inequality constraints, which the paper highlights, was avoided when modelling correlation structures
2) Some vector-vector and matrix-vector operations were being parallelized unecessarily in the tuning script, which actually slowed down parallel tuning by 2-5x
3) Inequality constraints are only enforced at certain observations, and
   users can now optionally choose different observation subsets for
   different built-in constraint types / variables via a named-list form
   of `qp_observations`

In addition, documentation and code comments were improved throughout. 

There are no issues with backwards compatibility - code will not break for users who wrote code using a previous version of the package for any feature. 

---

### Changes

#### 1) Inequality Constraints
- inequality-constrained fitting now uses a shared active-set wrapper around equality-only re-solves across the main fitting paths, including uncorrelated fits, correlated Gaussian GEE, correlated non-Gaussian GEE, and `blockfit_solve()`
- for structured correlation, the package now uses the Woodbury decomposition when the low-rank gate is met, so the constrained correlated solve can work through the low-rank `G_on^(1/2) F^(1/2)` representation rather than immediately forming dense `P x P` matrices (this was a previous bug)
- when the Woodbury gate is not met, when the Woodbury square-root is not usable, or when the Woodbury active-set route does not converge, the code falls back to the existing dense correlated solver; for Gaussian correlated response this is dense QP, and for iterative GLM paths this is dense SQP
- derivative-constraint handling and correlated active-set behavior were reviewed and documented more explicitly in the package help files so the solver flow and fallback rules now match the implementation
- active-set refinement can be slower in some simple bounded examples, especially when users disable tuning with `opt = FALSE`; for example, a range-bounded fit with bounds `(-150, 150)` may take longer because the active-set loop is still attempted before the dense fallback
- users can now choose which observations built-in QP constraints are
  enforced at, either with the legacy shared vector form or with a new
  named-list form keyed by entries such as
  `Time:qp_negative_derivative` and `Dose:qp_positive_derivative`
- package help for `lgspline()`, `process_qp()`, and the details page was
  updated to document the new `qp_observations` list form, the known key
  types, and the warning behavior for unknown keys
- focused tests were added for legacy shared-subset compatibility,
  per-constraint keyed subsets, per-variable keyed derivative subsets,
  and unknown-key warnings
- users can now control QR-pivot reduction separately for equality and
  inequality constraints with
  `qr_pivot_smoothing_constraints` and
  `qr_pivot_inequality_constraints`
- a new `parallel_qr_qp` option allows partition-wise QR pivoting of
  derivative/curvature inequality blocks across workers; built-in range
  and monotonicity constraints are intentionally left exact rather than
  reduced
- package help for `lgspline()`, `process_qp()`, and the details page was
  updated accordingly, including the grouped argument-list interface

#### Parallelism
- a small issue in the auxiliary parallelism benchmark materials was corrected during revision; this concerned development-time timing/benchmarking content rather than the core CRAN-facing fitting routines
- tuning now has two additional outer-level parallel controls:
  `parallel_grideval` and `parallel_bfgs`, both defaulting to `TRUE`
- `parallel_grideval` evaluates the initial tuning grid across the
  submitted cluster and, when more than six workers are available,
  augments the starting grid with additional random raw-scale penalty
  pairs drawn from an expanded range
- `parallel_bfgs` evaluates batches of damped BFGS trial steps across the
  submitted cluster and keeps the best improving candidate from each
  batch
- these outer tuning-parallel stages deliberately disable the inner
  per-criterion parallel kernels within each worker, avoiding nested use
  of the same cluster and improving compatibility with existing solver
  code
- the base `stats::optim()` tuning fallback now uses the same safe
  criterion-evaluation guard as the grid search and custom BFGS path, so
  isolated failed finite-difference probes are treated as poor candidate
  points rather than aborting the whole tuning run
- a regression in `blockfit_solve()` introduced while threading the new
  QR controls was corrected so `parallel_qr` is now passed through the
  Gaussian and weighted blockfit helper paths instead of falling back
  with `object 'parallel_qr' not found`
- related non-package benchmark artifacts used during development were
  removed from the submitted source tree so the CRAN submission now
  includes only package-relevant files; this cleanup also removes an
  unneeded tracked PDF artifact from `tests/testthat`
- this cleanup does not change the numerical fitting results reported by the package examples or tests; it only clarifies and streamlines the materials associated with the parallelism discussion

#### Miscellaneous
- Documentation updated throughout for clarity and consistency
- Small sample-size adjustment to tuning penalties changed from [(N+2)/(N-2)]^2 to (N+1)/(N-1) which is less dramatic, less controversial to reviewers.
- focused tests were cleaned so the current QR / QP checks no longer rely
  on `diffobj` for boolean diagnostics
- custom covariance helper functions submitted through `VhalfInv_fxn`,
  `Vhalf_fxn`, `VhalfInv_logdet`, `REML_grad`, and
  `custom_VhalfInv_loss` are now wrapped so they can use
  `correlation_id`, `spacetime`, and matching extra arguments while the
  fitted object still stores the simple one-argument wrappers used by
  posterior methods
- if `VhalfInv_fxn` is supplied and `VhalfInv_par_init` is omitted, the
  optimizer now initializes at `1e-2` rather than requiring the user to
  pass an explicit starting value

## Focused verification for this revision

- package reinstall from source completed successfully
- focused tests passing:
  - `tests/testthat/test-fitting.R`
  - `tests/testthat/test-advanced.R`
  - `tests/testthat/test_correlation_structure.R`
  - `tests/testthat/test-tuning-criterion.R`
  - `tests/testthat/test-parallel-qr.R`
  - `tests/testthat/test-glm-sqp.R`
  - `tests/testthat/test_qp_observations_list.R`
- focused verification script passing:
  - `scripts/codex_verify.R`
- additional direct canary fit passing:
  - `lgspline(y ~ spl(t), data.frame(t, y), family = binomial())`
  - clustered Gaussian tuning canary with
    `parallel_grideval = TRUE` and `parallel_bfgs = TRUE`
  - blockfit Gaussian no-correlation canary with `parallel_qr = FALSE`
    and a supplied cluster
  - custom-covariance canary with omitted `VhalfInv_par_init`
  - keyed `qp_observations` canary with bare `qp_range_lower` and
    derivative-specific observation subsets

## R CMD check results

0 errors | 0 warnings | 0 notes

There are no downstream dependencies for this package.
