# CRAN Comments

## Resubmission
This is a resubmission addressing numerical stability, interface consistency, and correctness issues identified during JSS review preparation. The changes span penalty parameterization, matrix operations, constrained estimation, S3 methods, and correlation-aware inference. 
It also introduces novel functionality, including a logLik S3 method, 2nd derivative constraints, an anti-derivative/integration function, and correlation-aware posterior sampling that propagates uncertainty in correlation parameters.

---

### Changes

#### 1. Numerical Stability & Performance
- Replaced `compute_trace_UGXX_wrapper` with `compute_trace_H`, which is more computationally stable.
- Introduced a three-tier fallback strategy for matrix inversion: Cholesky (`chol2inv`), Armadillo (`armaInv`), and eigendecomposition with ridge penalty for near-singular cases.
- Replaced custom `%**%` operator with `crossprod()` and `tcrossprod()` in critical hot-paths when possible to leverage BLAS/LAPACK optimization.
- Stabilized variance-covariance matrix computation via outer-product form $(UG^{1/2})(UG^{1/2})^\top$; corrected `each = K + 1` to `times = K + 1` in unscaling logic.
- Constraint matrix $\mathbf{A}$ now selects only linearly independent columns via pivoted QR to prevent redundancy and rank deficiency.
- Fixed division-by-zero in `compute_G_eigen` (`GhalfInv`) by replacing division by zero-valued `sqrt_inv_eigen_values` with multiplication by `sqrt(pmax(eig$values, 0))`.
- Consolidated three near-identical code paths in `compute_G_eigen` into a single `compute_one()` closure, fixing inconsistent `GhalfInv` computation for remainder chunks under non-Gaussian families with identity link.
- Removed incorrect max-norm scaling in `compute_trace_correction` (C++).
- Replaced all `gramMatrix()` calls (C++) with `crossprod()`. This was found to improve performance.
- Ensured `symmetric = TRUE` is passed to `eigen()` in `matsqrt()` and `matinvsqrt()`.

#### 2. Penalty Parameterization (tune_Lambda)
- Replaced softplus parameterization with exponential (`raw_penalty = exp(theta)`) throughout `tune_Lambda`. The softplus chain rule factor $\lambda / (1 + \lambda)$ was approximately 20% incorrect at moderate values; the exponential factor $\lambda$ is exact, verified by finite-difference agreement to machine precision.
- Updated gradient chain rule at 4 sites ($\lambda_w$, $\lambda_r$, predictor-specific, and partition-specific penalties) from $\lambda / (1 + \lambda)$ to $\lambda$.
- Updated regularizer gradient at 2 sites from $c \times (\lambda - 1)$ to $c \times (\lambda - 1) \times \lambda$.
- Replaced `softplus(new_lambda)` with `exp(new_lambda)` in 2 BFGS validity checks.
- Renamed three arguments to accept penalties on the raw (non-negative) scale, removing the requirement for callers to apply inverse-softplus transforms: `invsoftplus_initial_wiggle` → `initial_wiggle`, `invsoftplus_initial_flat` → `initial_flat`, `invsoftplus_penalty_vec` → `penalty_vec`.
- Corrected upstream penalty vector construction: output variable renamed to `penalty_vec`; default initialization changed from `rnorm(1, 0, 0.00001)` to `exp(rnorm(1, 0, 0.00001))`; removed `log(exp(x) - 1)` inverse-softplus transform for user-supplied penalties.
- Fixed bug where `penalty_vec < 0` (a comparison) was used instead of `penalty_vec <- c()` (an assignment) in the empty-vector fallback.
- Fixed inverted observation weights condition in `gcvu_fxn`: changed `==` to `!=` to prevent double-counting weights already absorbed into pre-scaled Gaussian/identity data.

#### 3. Main Interface (lgspline.R)
- Replaced deprecated `expansions_only` argument with `dummy_fit`, which runs the full pipeline but sets coefficients to zero and skips expensive fitting steps.
- Added `auto_encode_factors` (default `TRUE`) for automatic one-hot encoding of factors and character vectors via the formula interface.
- Fixed formula parsing for 3-way interactions (`spl(a, b, c)`). Previously a bug was preventing 3-way interactions from appearing occasionally.
- Resolved parallel cluster scoping issues via explicit `shared_env` export to worker nodes.
- `predict()` now returns first derivatives as a named list per predictor variable rather than a concatenated vector.
- Added `offset()` term detection, extraction, and storage in the formula interface.
- Added conditional dispatch to `blockfit_solve` before falling back to `get_B`, with guard conditions: `blockfit = TRUE`, non-empty `flat_cols`, $K > 0$.
- When `VhalfInv` is supplied but `Vhalf` is not, `Vhalf` is now unconditionally computed as `invert(VhalfInv)` for all family/link combinations.
- Added `return_lagrange_multipliers` option.
- `generate_posterior` closure updated with `enforce_constraints` and `max_rejection_draws` arguments; the standalone `generate_posterior()` S3-like wrapper forwards both new arguments. See Section 4 for the constraint-handling design.

#### 4. Constrained Estimation (get_B & blockfit_solve)
- Removed blockfit path from `get_B`. Blockfitting is handled by the new standalone `blockfit_solve` function which performs backfitting.
- Replaced recursive IRLS scheme in `get_B` with an explicit `for` loop (max 100 iterations), eliminating deep call-stack growth. Signature parameters `prevB`, `prevUnconB`, `iter_count`, and `prev_diff` are retained for call-site compatibility but ignored internally.
- Fixed operator precedence bug in blockfit MSE fallback: `w*(y - mu^2)` corrected to `w*(y - mu)^2`.
- `blockfit_solve` implements spline/flat column separation, Gaussian identity backfitting, GLM IRLS outer loop, and optional SQP refinement for QP inequality constraints.
- Added `qp_positive_2ndderivative` and `qp_negative_2ndderivative` to enforce convexity or concavity using SQP constraints, similar to first derivative constraints.
- Both `get_B` and `blockfit_solve` now return a `qp_info` list containing the raw `solve.QP` output components at convergence: `solution`, `lagrangian` (Lagrange multipliers for all active constraints), `active_constraints` (indices into the combined equality/inequality constraint matrix), `iact` (raw 1-based active inequality indices from `solve.QP`), `info_matrix` (unscaled penalised information matrix $\mathbf{D}$), `Amat_combined`, `bvec_combined`, `meq_combined`, `converged`, and `final_deviance`. These are used downstream to construct the active-constraint projection matrix $\mathbf{A}_{\mathrm{active}}$ for $\mathbf{U}$ and for Lagrange multiplier retrieval, providing correct variance-covariance estimates under binding inequality constraints.
- `generate_posterior` updated to use `qp_info$Amat_active` when constructing $\mathbf{U}$ for posterior draws, so the projection accounts for both equality and binding inequality constraints at the MAP. An accept/reject option (`enforce_constraints = TRUE`) is provided with a `max_rejection_draws` safety valve; the default remains the unconstrained Gaussian approximation with a warning, consistent with standard practice for constrained MAP spline posteriors.

#### 5. Correlation-Aware Inference
- Under `VhalfInv`, inference quantities ($\mathbf{G}$, $\mathbf{G}^{1/2}$, $\mathbf{U}$, `trace_XUGX`, `varcovmat`, Lagrange multipliers, and posterior draws) are now computed from whitened Gram matrices $\mathbf{X}_k^\top \mathbf{V}^{-1} \mathbf{X}_k$ rather than the unwhitened $\mathbf{X}_k^\top \mathbf{X}_k$.
- `varcovmat`, `trace_XUGX`, and `generate_posterior` now correctly incorporate GLM working weights ($\tilde{\mathbf{W}}$) and user observation weights ($\mathbf{D}$) into the whitened Gram matrix under `VhalfInv`.
- `sigmasq_tilde` (Gaussian + identity): fixed incorrect squaring form `mean((w * V^{-1/2}(y-ŷ))^2)` to `mean(w * (V^{-1/2}(y-ŷ))^2)`.
- Fixed `leave_one_out`: (1) pre-existing bug where `matmult_U` was called with $\mathbf{G}$ instead of $\mathbf{G}^{1/2}$, computing $\mathrm{diag}(\mathbf{X} \mathbf{U} \mathbf{G}^2 \mathbf{X}^\top)$ rather than the correct $\mathrm{diag}(\mathbf{X} \mathbf{U} \mathbf{G} \mathbf{X}^\top)$; (2) under `VhalfInv`, hat diagonal now uses $\mathbf{G}_{\mathrm{correct}}^{1/2}$; (3) added `leverage_threshold` argument (default 0.9999) to set high-leverage LOO entries to `NA` with a warning.
- Fixed `REML_objective`: whitened residuals were missing $\sqrt{\mathbf{W}}$ division; corrected to apply $(\tilde{\mathbf{W}}\mathbf{D})^{1/2}$ before $\mathbf{V}^{-1/2}$.
- Fixed `REML_grad`: `t(t(M)*v)` pattern (column scaling) corrected to `glm_weights * VhalfInvX` (row scaling).
- Corrected pre-whitening of $\mathbf{X}$ and $\mathbf{y}$ in `lgspline.fit`: inputs are now preserved in unwhitened form; whitening is applied internally within `get_B` and `blockfit_solve` where the full matrix structure is available.
- Added `generate_posterior_correlation()`: a standalone function that propagates correlation parameter uncertainty into the coefficient posterior. For each Monte Carlo draw, the correlation parameter vector $\boldsymbol{\rho}$ is sampled from its approximate normal posterior $\mathcal{N}(\hat{\boldsymbol{\rho}}_{\mathrm{REML}}, \mathbf{H}^{-1}_{\mathrm{BFGS}})$ on the unbounded working scale, the coefficients are re-estimated with $\mathbf{V}^{-1/2}(\boldsymbol{\rho}^{(m)})$ held fixed using the already-expanded design matrices $\mathbf{X}_k$, constraint matrix $\mathbf{A}$, and tuned penalty matrices $\boldsymbol{\Lambda}$ from the original fit (avoiding re-partitioning, re-expansion, or penalty re-tuning), and then a single coefficient draw is obtained from the re-estimated quantities. Draws producing non-positive-definite correlation matrices are rejected and redrawn (up to 50 attempts, then point estimate with warning). Re-estimation is performed at the level of the whitened GLS Gram matrix $\mathbf{X}^\top \mathbf{V}^{-1} \mathbf{X} + \boldsymbol{\Lambda}$, constrained projection $\mathbf{U}$, and post-fit inference, operating entirely on stored model object fields without requiring access to the original raw predictor matrix.
- Updated standalone `generate_posterior()` wrapper with `draw_correlation` argument (default `FALSE`) plus supporting arguments (`correlation_param_mean`, `correlation_param_vcov`, `correlation_VhalfInv_fxn`, `correlation_Vhalf_fxn`) that are forwarded to `generate_posterior_correlation()` when `draw_correlation = TRUE`.

#### 6. S3 Methods & Inference (Methods.R)
- `print.summary.lgspline` now uses `stats::printCoefmat` for p-value formatting and significance stars.
- Refactored `wald_univariate` to return an S3 object of class `wald_lgspline` with dedicated `print`, `summary`, and `plot` methods.
- Added `logLik.lgspline`: computes exact log-likelihood for Gaussian/identity; GLS log-likelihood under `VhalfInv`; falls back to `family$aic` or deviance-based approximation otherwise. Includes `include_prior` argument.
- Added `confint.lgspline`.
- Predict function can now automatically accept both one-hot encoded and character/factor versions of variables. The `find_extremum` function was adapted for this change as well.

#### 7. Weibull AFT
- Fixed score formula in `weibull_qp_score_function`: `exp(z) - sigma * status * w` corrected to `sigma * w * (exp(z) - status)` so both terms are scaled identically.
- Fixed double-indexing of `status` in `weibull_dispersion_function` and `weibull_family()$custom_dev.resids`: `status[order_indices]` was applied twice; corrected to use `status` after initial subsetting.
- Widened Brent optimization bounds for scale parameter from $[\mathrm{init\_scale}/5,\, \mathrm{init\_scale} \times 5]$ to $[\mathrm{init\_scale}/10,\, \mathrm{init\_scale} \times 10]$.
- Applied `crossprod()` substitutions in `unconstrained_fit_weibull` and `unconstrained_fit_default`.
- Updated Roxygen documentation across all Weibull functions to unambiguously distinguish scale ($\sigma$) from dispersion ($\sigma^2$).

#### 8. Anti-Derivative Function
- Constructed a novel, comprehensive function for analytically and numerically evaluated indefinite and definite integrals for fitted models.

#### 9. Helper Functions & Miscellaneous
- `weibull_family()`: added `$aic`, `$loglik`, and `$dev.resids` members for correct `logLik.lgspline` dispatch.
- `create_onehot()`: added `drop_first = FALSE` parameter.
- `prior_loglik()`: added option to return the full multivariate normal prior log-likelihood.
- Introduced `.compute_dist_block`, `.rank_dists`, and `.reml_grad_from_dV` helpers to reduce redundancy in correlation structure code; switched to exponential parameterization for correlation structure fitting.
- `dispersion_function` signature updated to include `VhalfInv`.
- Fixed character-to-numeric resolution for `do_not_cluster_on_these` in both formula and direct-matrix paths.
- Fixed `find_extremum`: added column names to `predictors_vals` after `transf()`/`inv_transf()` cycles; subsetted L-BFGS-B bounds to the variables being optimized when `select_vars_fl = TRUE`.

#### 10. Documentation
- Corrected sign in `@param initial_wiggle` default: `exp(c(-25, 14, ...))` → `exp(c(-25, -14, ...))`.
- Updated `@param auto_encode_factors` to document default `TRUE`.
- Updated `@param dummy_fit` to replace deprecated `expansions_only` language.
- Updated `@param blockfit`, `@param return_lagrange_multipliers`, `@param do_not_cluster_on_these`, `@param VhalfInv`, `@param Vhalf`.
- Updated `@return` entries for `G`, `Ghalf`, `varcovmat`, `trace_XUGX`, `sigmasq_tilde`, `lagrange_multipliers`, `generate_posterior`, `find_extremum`, `predict`, and `A`.
- Added `@details` sections documenting correlation-aware inference, blockfit dispatch, and exponential penalty parameterization, and added more information elsewhere.
- Added `@seealso` entries for `logLik.lgspline`, `confint.lgspline`, `leave_one_out`, and `blockfit_solve`.
- Rewrote `get_B` Roxygen header with updated path descriptions and structured section comments.
- Updated `generate_posterior` documentation with full coefficient-draw and dispersion-draw equations, constraint-handling rationale, `enforce_constraints`/`max_rejection_draws` parameter descriptions, and correlation parameter posterior equations.
- Added full Roxygen documentation for `generate_posterior_correlation()` with: display equations for the correlation parameter draw $\boldsymbol{\rho}^{(m)} \sim \mathcal{N}(\hat{\boldsymbol{\rho}}_{\mathrm{REML}}, \mathbf{H}^{-1}_{\mathrm{BFGS}})$, the whitened Gram matrix $\mathbf{X}^\top \mathbf{V}^{-1} \mathbf{X} + \boldsymbol{\Lambda}$, and the constraint projection $\mathbf{U} = \mathbf{I} - \mathbf{G}_{\mathrm{correct}} \mathbf{A} (\mathbf{A}^\top \mathbf{G}_{\mathrm{correct}} \mathbf{A})^{-1} \mathbf{A}^\top$; enumerated algorithm description; computational cost analysis; normal approximation quality caveats; and a complete worked example with exchangeable correlation.

---

## R CMD check results

0 errors | 0 warnings | 0 notes

---

## Reverse dependencies

There are currently no downstream dependencies for this package.
