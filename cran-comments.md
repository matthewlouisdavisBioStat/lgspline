# CRAN Comments

## Resubmission
This is a resubmission addressing numerical stability, interface consistency, and correctness issues identified during JSS review preparation. 

In addition to computational and correctness issues found during ongoing package development, this revision also incorporates detailed editorial feedback received during a Journal of Statistical Software (JSS) pre-review on academic software standards, including interface design, S3 method completeness, documentation, examples, and replication-readiness.

It also introduces novel functionality, including a logLik S3 method, 2nd derivative constraints, a numerical integration interface for fitted models, correlation-aware posterior sampling that propagates uncertainty in correlation parameters, equation printing, and Cox Proportional Hazards and Negative-Binomial regression models and tests.

I understand large overhauls such as this aren't preferred; although the changes are expansive, there should not be any compatibility issues with previous versions of the package, and many of the changes were necessary for ensuring the mathematical rigor and correctness of the package.

I do not expect future changes to the package to be this large in the future, I expect only documentation changes and targeted bug fixes. Hence, this is version 1.0.1.

(A bug of too-long runtimes and mis-spellings accounted for v 1.0.0)

Thank you for your time.

---

### Changes

#### 1. Numerical Stability & Performance
- Replaced `compute_trace_UGXX_wrapper` with `compute_trace_H`, which is more computationally stable.
- Introduced a three-tier fallback strategy for matrix inversion: Cholesky (`chol2inv`), Armadillo (`armaInv`), and eigendecomposition with ridge penalty for near-singular cases.
- Replaced custom `%**%` operator with `crossprod()` and `tcrossprod()` when possible to leverage BLAS/LAPACK optimization.
- Stabilized variance-covariance matrix computation via outer-product form $(UG^{1/2})(UG^{1/2})^\top$; corrected `each = K + 1` to `times = K + 1` in unscaling logic.
- Constraint matrix $\mathbf{A}$ now selects only linearly independent columns via pivoted QR to prevent redundancy and rank deficiency.
- Fixed division-by-zero in `compute_G_eigen` (`GhalfInv`) by replacing division by zero-valued `sqrt_inv_eigen_values` with multiplication by `sqrt(pmax(eig$values, 0))`.
- Consolidated three near-identical code paths in `compute_G_eigen` into a single `compute_one()` closure, fixing inconsistent `GhalfInv` computation for remainder chunks under non-Gaussian families with identity link.
- Removed incorrect max-norm scaling in `compute_trace_correction` (C++).
- Replaced remaining internal `gramMatrix()` call sites with `crossprod()` in the package code.
- Ensured `symmetric = TRUE` is passed to `eigen()` in `matsqrt()` and `matinvsqrt()`, which were made more efficient and computationally stable as well.
- For `nr_iterate`, replaced `crossprod(t(invert(H)), g)` with `solve(H, g)` (direct linear solve, no explicit inverse); scaling changed from `sqrt(mean(abs(H)))` to `sqrt(max(abs(diag(H))))` for tighter condition-number reduction on diagonal-dominant information matrices; eigendecomposition fallback retained for near-singular cases.
- For `damped_newton_r, Newton direction is now computed once per outer iteration and reused across damping half-steps, eliminating redundant gradient/Hessian/factorization evaluations (up to `max_dmp_steps` per iteration). Damping loop no longer accepts a worsening step when all half-steps fail; instead returns the current iterate immediately.

#### 2. Penalty Parameterization (tune_Lambda)
- Replaced softplus parameterization with exponential (`raw_penalty = exp(theta)`) throughout `tune_Lambda`. The softplus chain rule factor $\lambda / (1 + \lambda)$ was approximately 20% incorrect at moderate values; the exponential factor $\lambda$ is exact, verified by finite-difference agreement to machine precision.
- Updated gradient chain rule at 4 sites ($\lambda_w$, $\lambda_r$, predictor-specific, and partition-specific penalties) from $\lambda / (1 + \lambda)$ to $\lambda$.
- Updated regularizer gradient at 2 sites from $c \times (\lambda - 1)$ to $c \times (\lambda - 1) \times \lambda$.
- Replaced `softplus(new_lambda)` with `exp(new_lambda)` in 2 BFGS validity checks.
- Renamed three arguments to accept penalties on the raw (non-negative) scale, removing the requirement for callers to apply inverse-softplus transforms: `invsoftplus_initial_wiggle` -> `initial_wiggle`, `invsoftplus_initial_flat` -> `initial_flat`, `invsoftplus_penalty_vec` -> `penalty_vec`.
- Corrected upstream penalty vector construction: output variable renamed to `penalty_vec`; default initialization changed from `rnorm(1, 0, 0.00001)` to `exp(rnorm(1, 0, 0.00001))`; removed `log(exp(x) - 1)` inverse-softplus transform for user-supplied penalties.
- Fixed bug where `penalty_vec < 0` (a comparison) was used instead of `penalty_vec <- c()` (an assignment) in the empty-vector fallback.
- Fixed inverted observation weights condition in `gcvu_fxn`: changed `==` to `!=` to prevent double-counting weights already absorbed into pre-scaled Gaussian/identity data.
- The construction of the cubic smoothing spline penalty matrix was rigorously verified and corrected as needed, see `?lgspline::get_2ndDerivPenalty`.

#### 3. Main Interface (lgspline.R)
- Several of the following changes respond directly to JSS editor feedback on R interface conventions and user-facing complexity.
- Replaced deprecated `expansions_only` argument with `dummy_fit`, which runs the full pipeline but sets coefficients to zero and skips expensive fitting steps.
- Added `auto_encode_factors` (default `TRUE`) for automatic one-hot encoding of factors and character vectors via the formula interface.
- Added optional grouped argument lists (`penalty_args`, `tuning_args`, `expansion_args`, `constraint_args`, `qp_args`, `parallel_args`, `covariance_args`, `return_args`, `glm_args`) so related controls can be supplied more compactly while preserving backwards-compatible direct arguments.
- Refactored formula/data preprocessing into a standalone exported helper `process_input()`. This now centralizes formula parsing, factor auto-encoding, offset extraction, factor-group bookkeeping, and backward-compatible argument translation. In particular, deprecated `expansions_only` inputs are still accepted there with warnings and mapped to `dummy_fit`, while misspelled `shur_correction_function` inputs trigger a warning.
- Fixed formula parsing for 3-way interactions (`spl(a, b, c)`). Previously a bug was preventing 3-way interactions from appearing occasionally.
- Resolved parallel cluster scoping issues via explicit `shared_env` export to worker nodes.
- `predict()` now returns first derivatives as a named list per predictor variable rather than a concatenated vector.
- Added `offset()` term detection, extraction, and storage in the formula interface.
- Added conditional dispatch to `blockfit_solve` before falling back to `get_B`, with guard conditions: `blockfit = TRUE`, non-empty `flat_cols`, $K > 0$.
- When `VhalfInv` is supplied but `Vhalf` is not, `Vhalf` is now unconditionally computed as `invert(VhalfInv)` for all family/link combinations, and the built-in correlation structures now also provide direct `Vhalf_fxn` constructors in the Gaussian case.
- Added `return_lagrange_multipliers` option.
- `estimate_dispersion` is no longer overridden to be TRUE when `estimate_varcovmat` is TRUE, which was adjusting the variance-covariance matrix by estimated exponential dispersion, even for families like ordinary binomial/Poisson which have fixed dispersion at 1.
- Quadratic programming arguments for `qp_negative/positive_derivative/2ndderivative` accept variable names and/or TRUE and FALSE, ensuring backwards compatibility while allowing unique QP constraints for each predictor.
- Allow for formulas like `y ~ s(t)` rather than just `y ~ spl(t)` so users will not have their workflows disrupted if they switch from the mgcv package, fitting similar models.
- Default `initial_wiggle` and `initial_flat` candidate grids were updated from earlier versions.

#### 4. Partitioning & Knot Placement (make_partitions.R / lgspline.fit)
- Unified 1-D and multi-D partitioning: the previous 1-D path used equally-spaced quantile knots on the standardized $[0, 1]$ scale; both paths now use kmeans clustering via `make_partitions()`. Knot locations in 1-D now concentrate in data-dense regions rather than being equally spaced, which is more appropriate for skewed or multimodal predictors.
- Internalized all standardization for clustering inside `make_partitions()`. Previously, `lgspline.fit` standardized the predictor matrix in-place, called `make_partitions`, then back-transformed; this exposed standardized-scale coordinates in the return list (centers, knots, `assign_partition` closure). Now `make_partitions` receives raw-scale predictors, standardizes a local copy, clusters on that copy, and back-transforms centers and knots before returning. The returned `assign_partition` closure accepts raw-scale input and standardizes internally.
- As a result, `make_partition_list` (returned and accepted by `lgspline.fit`) now contains raw-scale centers and knots, making it straightforward for users to inspect, modify, and reuse partition structures from a `dummy_fit`.
- `assign_partition` throughout the pipeline (training, prediction, dummy-fit early return) now consistently accepts raw-scale predictors. The special-cased 1-D closure `function(x) rowMeans(cbind(x))`, which silently required standardized input, has been removed.
- `partition_codes` and `partition_bounds` now use the integer-like scheme ($k - 0.5$ / $1{:}(K+1)$) for both 1-D and multi-D, eliminating the mixed continuous/integer encoding that previously differed by dimensionality.
- `transf` and `inv_transf` are still defined in `lgspline.fit` and exposed in the return list as `knot_scale_transf` / `knot_scale_inv_transf` for backward compatibility, but are no longer applied to predictors in the main fitting or prediction pipeline. All `inv_transf(knot_values)` and `inv_transf(knot_values_chunk)` calls in the constraint matrix construction have been removed; knot values are raw-scale throughout.
- `get_centers()` gains a `data_already_processed` parameter (default `FALSE`). When `TRUE`, the binary-column zeroing step inside `get_centers` is skipped, preventing double-zeroing when `make_partitions` has already prepared the clustering data.
- `make_partitions()` gains `standardize`, `standardize_mode`, `dummy_adder`, and `dummy_dividor` parameters. `standardize_mode = "auto"` selects minmax scaling for one effective clustering dimension and z-score for multiple dimensions, matching the previous per-path behavior. The function now additionally returns `standardize_transf`, `standardize_inv_transf`, and `centers_std` for diagnostic use.
- Custom 1-D knots (`custom_knots` with `q_predictors == 1`) are handled by a dedicated branch that builds a lightweight partition list using raw-scale breakpoints and `findInterval`-based assignment, bypassing kmeans while preserving full compatibility with `knot_expand_list`.

#### 5. Sum-to-Zero Constraints for Encoded Factors (process_input.R / lgspline.fit)
- When `auto_encode_factors = TRUE` and factor or character columns are one-hot encoded (without dropping a reference level), identifiability requires that the sum of indicator-level coefficients is zero within each partition. `process_input` now constructs a `factor_groups` named list mapping each original factor column name to the integer column positions of its one-hot indicators within the predictor matrix, and returns this in its output.
- `lgspline.fit` accepts `factor_groups` as a new parameter. For each group with at least two resolved indicator positions, one equality constraint column per partition is appended to `constraint_vectors` before A is assembled: each column places a $1$ at every expansion row corresponding to a group indicator within that partition's coefficient block, enforcing $\sum_{i} \beta_{ji,k} = 0$ for partition $k$. Groups with fewer than two resolved positions are silently ignored.
- These factor constraints are appended before the intercept and offset equality constraints, and before the pivoted QR linear-independence reduction is applied to A, so any redundancy with user-supplied constraints is handled automatically.
- Added `@param factor_groups` Roxygen documentation to `lgspline.fit`.

#### 6. Constrained Estimation (get_B & blockfit_solve)
- Removed blockfit path from `get_B`. Blockfitting is handled by the new standalone `blockfit_solve` function which performs backfitting.
- Replaced recursive call in `get_B` with an explicit `for` loop (max 100 iterations), eliminating deep call-stack growth. Signature parameters `prevB`, `prevUnconB`, `iter_count`, and `prev_diff` are retained for call-site compatibility but ignored internally.
- Fixed operator precedence bug in blockfit MSE fallback: `w*(y - mu^2)` corrected to `w*(y - mu)^2`.
- `blockfit_solve` implements spline/flat column separation, Gaussian identity backfitting, GLM newton-raphson outer loop, and optional SQP refinement for QP inequality constraints.
- Added `qp_positive_2ndderivative` and `qp_negative_2ndderivative` to enforce convexity or concavity using SQP constraints, similar to first derivative constraints.
- Added Woodbury-accelerated GEE solver paths inside `get_B` for correlated-data fits when the effective correction from `V^{-1} - I` is low-rank enough to exploit. Dense GEE solving remains as the fallback when that acceleration is not advantageous.
- Added automatic routing for inequality constraints in `get_B`: block-separable constraints are handled by a partition-wise active-set refinement, while genuinely cross-partition/global constraints fall back to dense `quadprog::solve.QP`.
- Refactored built-in QP-constraint assembly out of `lgspline.fit` into a standalone exported helper `process_qp()`, which now handles range, monotonicity, first-derivative, and second-derivative constraints in one place and can be tested independently.
- Both `get_B` and `blockfit_solve` now return a `qp_info` list containing the raw `solve.QP` output components at convergence: `solution`, `lagrangian` (Lagrange multipliers for all active constraints), `active_constraints` (indices into the combined equality/inequality constraint matrix), `iact` (raw 1-based active inequality indices from `solve.QP`), `info_matrix` (unscaled penalised information matrix $\mathbf{D}$), `Amat_combined`, `bvec_combined`, `meq_combined`, `converged`, and `final_deviance`. These are used downstream to construct the active-constraint projection matrix $\mathbf{A}_{\mathrm{active}}$ for $\mathbf{U}$ and for Lagrange multiplier retrieval, providing correct variance-covariance estimates under binding inequality constraints.
- The method of fitting SQPs and covariance structures was made much more efficient.

#### 7. Constrained Posterior Sampling (generate_posterior)
- Replaced the SQP-based projection step with elliptical slice sampling (Murray, Adams & MacKay 2010), implemented entirely in base R.
- The elliptical slice sampler proposes on ellipses $\beta(\theta) = \beta_{\mathrm{cur}} \cos\theta + \nu \sin\theta$ where $\nu \sim N(0, \sigma^2 L L^\top)$ and shrinks the angle bracket until a feasible point is found. Equality constraints (smoothness at knots) are enforced by the $\mathbf{U}$ projection; inequality constraints (monotonicity, range bounds, custom QP) are enforced by the slice sampler. This produces exact draws from the truncated multivariate normal posterior.
- For `num_draws > 1` with active inequality constraints, draws form a Markov chain (the chain state is maintained across draws); unconstrained draws remain i.i.d.
- $L_{\mathrm{post}} = (1 / \mathrm{sd}_y) \, \mathbf{U} \, \mathbf{G}^{1/2}$ is precomputed once outside the draw loop (it does not depend on $\sigma^2$). This eliminates a redundant recomputation that was present in both the single-draw and multi-draw paths.
- Dispersion draws, back-transformation, and output packaging are factored into shared helpers (`.draw_sigmasq`, `.backtransform`, `.package`) to remove three instances of duplicated logic.
- Removed dead code: an unreachable `Ghalf_block <- Reduce(...)` path that drew random normals into a vector that was never used, and a vestigial `mu_mode` assignment.
- Added `override_*` arguments (`override_B_raw`, `override_U`, `override_Ghalf_correct`, `override_Ghalf`, `override_VhalfInv`, `override_sigmasq_tilde`, `override_trace_XUGX`) to the closure so that `generate_posterior_correlation` can inject updated covariance-side quantities directly, bypassing the broken shallow-copy pattern (see Section 8).
- When `override_Ghalf_correct` is provided, the GLS path skips rebuilding the whitened Gram matrix and computes $L_{\mathrm{post}}$ directly from the supplied half-inverse, saving one `collapse_block_diagonal` + `gramMatrix` + `matinvsqrt` call per correlation draw.

#### 8. Correlation-Aware Inference
- Under `VhalfInv`, inference quantities ($\mathbf{G}$, $\mathbf{G}^{1/2}$, $\mathbf{U}$, `trace_XUGX`, `varcovmat`, Lagrange multipliers, and posterior draws) are now computed from whitened Gram matrices $\mathbf{X}_k^\top \mathbf{V}^{-1} \mathbf{X}_k$ rather than the unwhitened $\mathbf{X}_k^\top \mathbf{X}_k$.
- `varcovmat`, `trace_XUGX`, and `generate_posterior` now correctly incorporate GLM working weights ($\tilde{\mathbf{W}}$) and user observation weights ($\mathbf{D}$) into the whitened Gram matrix under `VhalfInv`.
- `sigmasq_tilde` (Gaussian + identity): fixed incorrect squaring form `mean((w * V^{-1/2}(y-\\tilde{y}))^2)` to `mean(w * (V^{-1/2}(y-\\tilde{y}))^2)`.
- Fixed `leave_one_out`: the hat diagonal is now computed from the correct factorizations of $\mathbf{X}\mathbf{U}\mathbf{G}\mathbf{X}^{\top}$ in both the uncorrelated and `VhalfInv` paths, rather than pairing $\mathbf{X}\mathbf{U}\mathbf{G}^{1/2}$ with raw $\mathbf{X}$. This also corrected the correlated-path typo in the whitened Gram calculation and aligned the weighting there with the fitted GLS system via $\sqrt{\mathbf{W}}$. Added `leverage_threshold` argument (default 100) to flag high-leverage observations with a warning.
- Fixed `REML_objective`: whitened residuals were missing $\sqrt{\mathbf{W}}$ division; corrected to apply $(\tilde{\mathbf{W}}\mathbf{D})^{1/2}$ before $\mathbf{V}^{-1/2}$.
- Fixed `REML_grad`: `t(t(M)*v)` pattern (column scaling) corrected to `glm_weights * VhalfInvX` (row scaling).
- Corrected pre-whitening of $\mathbf{X}$ and $\mathbf{y}$ in `lgspline.fit`: inputs are now preserved in unwhitened form; whitening is applied internally within `get_B` and `blockfit_solve` where the full matrix structure is available.
- Added `generate_posterior_correlation()`: a standalone function that propagates correlation parameter uncertainty into the coefficient posterior using a user-custom variance-covariance matrix, or by default, the inverse BFGS approximate Hessian used to estimate the correlation parameters.

#### 9. S3 Methods & Inference (Methods.R)
- These changes directly address JSS feedback requesting more conventional classes and methods for returned objects, including coefficient-style summaries, `confint()`, and `logLik()`.
- `print.summary.lgspline` now uses `stats::printCoefmat` for p-value formatting and significance stars.
- Refactored `wald_univariate` to return an S3 object of class `wald_lgspline` with dedicated `print`, `summary`, and `plot` methods.
- Added `coef.wald_lgspline` and `confint.wald_lgspline` methods for extracting estimates and intervals directly from `wald_lgspline` objects.
- Added `logLik.lgspline`: computes exact log-likelihood for Gaussian/identity; GLS log-likelihood under `VhalfInv`; falls back to `family$aic` or deviance-based approximation otherwise. Includes `include_prior` argument.
- Added `confint.lgspline`.
- `confint.lgspline` now appends working-scale Wald intervals for estimated correlation parameters when `VhalfInv_params_estimates` and `VhalfInv_params_vcov` are available.
- Predict function can now automatically accept both one-hot encoded and character/factor versions of variables. The `find_extremum` function was adapted for this change as well.
- Updated `generate_posterior` documentation so the public wrapper description matches the current implementation, including the present `enforce_qp_constraints` argument behavior and the underlying elliptical-slice-sampling closure used for constrained draws.
- Updated plotting to include `se.fit` argument for plotting standard errors in 1-dimensional plots, and automatically set non-selected predictors to 0.
- Included method for printing out the equation of the fitted function `equation`.
- Added exported `leave_one_out()` for fast leave-one-out cross-validated predictions in the Gaussian/identity case, using the constrained hat-matrix shortcut and accounting for observation weights and optional correlation structures.

#### 10. Weibull AFT
- Fixed score formula in `weibull_qp_score_function`: `exp(z) - sigma * status * w` corrected to `sigma * w * (exp(z) - status)` so both terms are scaled identically.
- Fixed double-indexing of `status` in `weibull_dispersion_function` and `weibull_family()$custom_dev.resids`: `status[order_indices]` was applied twice; corrected to use `status` after initial subsetting.
- Widened Brent optimization bounds for scale parameter from $[\mathrm{init\_scale}/5,\, \mathrm{init\_scale} \times 5]$ to $[\mathrm{init\_scale}/10,\, \mathrm{init\_scale} \times 10]$.
- Applied `crossprod()` substitutions in `unconstrained_fit_weibull` and `unconstrained_fit_default`.
- Updated Roxygen documentation across all Weibull functions to unambiguously distinguish scale ($\sigma$) from dispersion ($\sigma^2$).

#### 11. Numerical Integration
- Added an `integrate.lgspline` S3 method that applies Gauss-Legendre quadrature to fitted `lgspline` objects through the model's `predict()` method, on either the response scale or the link scale.
- Added an `integrate` S3 generic with `integrate.default` delegating to `stats::integrate`, so users can call `integrate(fit, ...)` for `lgspline` objects without changing the default behavior for ordinary functions.

#### 12. Helper Functions & Miscellaneous
- `weibull_family()`: added `$aic`, `$loglik`, and `$dev.resids` members for correct `logLik.lgspline` dispatch.
- Added correctly spelled `weibull_schur_correction()` and retained `weibull_shur_correction()` as an exported alias for backwards compatibility.
- `create_onehot()`: added `drop_first = FALSE` parameter.
- `prior_loglik()`: added option to return the full multivariate normal prior log-likelihood.
- Introduced `.compute_dist_block`, `.rank_dists`, and `.reml_grad_from_dV` helpers to reduce redundancy in correlation structure code; switched to exponential parameterization for correlation structure fitting.
- `dispersion_function` signature updated to include `VhalfInv`, and modified to by default divide by `family$variance(mu)`, yielding estimates equivalent to Method-of-Moment estimators.
- Fixed character-to-numeric resolution for `do_not_cluster_on_these` in both formula and direct-matrix paths.
- Fixed `find_extremum`: added column names to `predictors_vals` after `transf()`/`inv_transf()` cycles; subsetted L-BFGS-B bounds to the variables being optimized when `select_vars_fl = TRUE`.
- Added workaround for inverse link functions that yield invalid values of response for certain families, i.e. `1/(1+exp(-0)) = 0.5` is not valid for binomial data, so will return a warning using `family=binomial()`. Fix is to swap family with quasi version for these cases, for sake of `glm.fit()` and hot-starting the coefficients for full penalized Newton-Raphson optimization. See `unconstrained_fit_default`.
- Function `safe_replace_var` added to prevent a re-naming bug in plotting.

#### 13. Exact Variance-Covariance & Other
- Added optional (default FALSE) `exact_varcovmat` argument to `lgspline()` and `lgspline.fit()`. When TRUE, the returned `varcovmat` is replaced with the exact frequentist variance-covariance of the constrained penalized estimator rather than the default asymptotic (Bayesian posterior) version.

#### 14. Documentation
- A substantial documentation pass was motivated by JSS editorial feedback on academic presentation, implementation transparency, and example quality.
- Audited and tightened the main Roxygen headers, in-line comments, and `lgspline-details.R` so the argument descriptions and implementation notes match the current code paths, with particular attention to penalty parameterization, correlation-aware inference, posterior sampling, QP helper behavior, and leave-one-out diagnostics.
- Updated package references and `DESCRIPTION` citations to align with the long-form details page, including explicit discussion of the default k-means-based partitioning scheme and its citations.
- Updated `@param dummy_fit` to replace deprecated `expansions_only` language.
- Updated `@param blockfit`, `@param return_lagrange_multipliers`, `@param do_not_cluster_on_these`, `@param VhalfInv`, `@param Vhalf`, and `@param exact_varcovmat`.
- Added `@param factor_groups` to `lgspline.fit` documenting the named list structure, the sum-to-zero constraint mechanism, and the relationship to `process_input` and `auto_encode_factors`.
- Updated `@return` entries for `G`, `Ghalf`, `varcovmat`, `trace_XUGX`, `sigmasq_tilde`, `lagrange_multipliers`, `generate_posterior`, `find_extremum`, `predict`, and `A`.
- Added `@details` sections documenting correlation-aware inference, blockfit dispatch, exponential penalty parameterization, the unified kmeans partitioning scheme, and the current numerical integration interface.
- Added `@seealso` entries for `logLik.lgspline`, `confint.lgspline`, `leave_one_out`, and `blockfit_solve`.
- Rewrote `get_B` Roxygen header with updated path descriptions and structured section comments.
- Updated `generate_posterior` documentation with the current coefficient-draw and dispersion-draw equations, the underlying elliptical-slice-sampling mechanism for constrained draws, the current `enforce_qp_constraints` wrapper behavior, and correlation-parameter posterior equations.
- Added full Roxygen documentation for `generate_posterior_correlation()` with: display equations for the correlation parameter draw $\boldsymbol{\rho}^{(m)} \sim \mathcal{N}(\hat{\boldsymbol{\rho}}_{\mathrm{REML}}, \mathbf{H}^{-1}_{\mathrm{BFGS}})$, the whitened Gram matrix $\mathbf{X}^\top \mathbf{V}^{-1} \mathbf{X} + \boldsymbol{\Lambda}$, and the constraint projection $\mathbf{U} = \mathbf{I} - \mathbf{G}_{\mathrm{correct}} \mathbf{A} (\mathbf{A}^\top \mathbf{G}_{\mathrm{correct}} \mathbf{A})^{-1} \mathbf{A}^\top$; enumerated algorithm description; computational cost analysis; normal approximation quality caveats; and a complete worked example with exchangeable correlation.
- Updated `make_partitions` Roxygen documentation to reflect internalized standardization, raw-scale return values, new parameters (`standardize`, `standardize_mode`, `dummy_adder`, `dummy_dividor`), and new return list fields (`standardize_transf`, `standardize_inv_transf`, `centers_std`).
- Updated `get_centers` Roxygen documentation to reflect the new `data_already_processed` parameter.
- Decomposed `get_B` and `tune_Lambda` into sub-functions, with individual documentation, and placed in scripts separate from HelperFunctions.R.
- Replaced "nr" for number of rows with "N_obs", "qcols" and "q" predictors with "q_predictors", "nc" for number of columns of partition design matrix with "p_expansions", and "nca" for number of columns of A matrix with "R_constraints"; these match notation provided in `lgspline-details` and are more informative for a reader.
- Updated examples to use conditional `survival` usage for suggested-package compliance.

#### 15. Tests
- Expanded `testthat` coverage to exercise formula preprocessing (`process_input`), factor auto-encoding / `factor_groups`, offset handling, `integrate.lgspline`, `prior_loglik`, and `leave_one_out` against explicit refits in a Gaussian linear case.
- Added narrow regression tests around the posterior/correlation path, including forwarding of `enforce_qp_constraints` and availability of efficient Gaussian correlation square-root matrix constructors.
- The new tests were chosen to cover exported user-facing functionality efficiently, in addition to the existing fit, prediction, plotting, GLM/QP, parallel, and correlation-structure tests.

#### 16. Survival & Count Regression Helpers
- Moved Weibull AFT helper functions (`loglik_weibull`, `weibull_scale`, `weibull_family`, `weibull_dispersion_function`, `weibull_glm_weight_function`, `weibull_qp_score_function`, `weibull_schur_correction`, `unconstrained_fit_weibull`, `lgspline_weibull`) from `HelperFunctions.R` into a standalone `weibull_helpers.R` script.
- Added `cox_helpers.R`: Cox proportional hazards regression via Breslow partial log-likelihood, including `loglik_cox`, `score_cox`, `info_cox`, `cox_family`, `cox_glm_weight_function`, `cox_dispersion_function`, `cox_qp_score_function`, `cox_schur_correction`, `unconstrained_fit_cox`, and the convenience wrapper `lgspline_cox`. No external survival package dependencies are required for fitting.
- Added `negbin_helpers.R`: NB2 negative binomial regression with log link, including `loglik_negbin`, `score_negbin`, `info_negbin`, `negbin_theta`, `negbin_family`, `negbin_glm_weight_function`, `negbin_dispersion_function`, `negbin_qp_score_function`, `negbin_schur_correction`, `unconstrained_fit_negbin`, and the convenience wrapper `lgspline_negbin`. Shape parameter $\theta$ is profiled via Brent's method, with Schur complement correction for $\theta$ uncertainty analogous to the Weibull scale correction. The dispersion function uses `VhalfInv` when available to whiten Pearson residuals for moment-based initialization of $\theta$, supporting GEE paths.
- All three helper scripts follow identical interface conventions (`unconstrained_fit_fxn`, `glm_weight_function`, `dispersion_function`, `qp_score_function`, `schur_correction_function`) and slot directly into `get_B` Path 3 (no correlation) or Path 1b (GEE) without modifications to the core fitting machinery.
- Added `testthat` tests for Cox PH and negative binomial helpers covering log-likelihood correctness, score vanishing at the MLE, comparison against `survival::coxph` and `MASS::glm.nb`, information matrix positive definiteness, Newton ascent direction, unconstrained fit recovery, QP score interface, weighted likelihood, dispersion/weight functions, Schur correction semi-definiteness, and numerical gradient verification.

## R CMD check results

0 errors | 0 warnings | 0 notes

There are no downstream dependencies for this package.

