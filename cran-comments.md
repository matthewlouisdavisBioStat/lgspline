# CRAN Comments

## Resubmission
This is a resubmission addressing minor improvements and documentation to the package.

---

### Changes

#### 1 Modified GCV (`gcv_gamma`)
- Added a user-facing `gcv_gamma` argument to `lgspline()` and `lgspline.fit()`, defaulting to `1.4`, and threaded it through `tune_Lambda()` so automatic tuning now uses the modified GCV denominator
  \[
  N\{1 - \gamma\,\mathrm{tr}(H)/N\}^{2}
  \]
  rather than the ordinary `gamma = 1` form. Setting `gcv_gamma = 1` recovers standard GCV exactly.
- Updated `.compute_gcvu_gradient()` so the denominator derivative includes the same `gcv_gamma` factor, preserving consistency between the objective and the closed-form gradient used by the custom BFGS path.
- Documented the rationale in `lgspline-details.R` with the primary citation: Kim, Y.-J. and Gu, C. (2004), *Smoothing Spline Gaussian Regression: More Scalable Computation via Efficient Approximation*, *Journal of the Royal Statistical Society: Series B*, 66(2), 337-356, Section 4, equation 4.1. The default `gcv_gamma = 1.4` follows their recommended range for reducing occasional severe undersmoothing.

#### 2 Tuning criterion switch (`tuning_criterion`)
- Added a user-facing `tuning_criterion` argument to `lgspline()`, `lgspline.fit()`, and `tune_Lambda()`, with allowed values `"loo"` and `"gcv"`, defaulting to `"loo"`.
- The new default tuning path uses exact leave-one-out cross-validation on the same transformed or linearized tuning problem already used internally for GCV. For Gaussian identity-link tuning this is exact PRESS/LOO; for correlation structures and/or non-identity links it is computed on the same working scale already used by the package's tuning machinery.
- The exact LOO implementation computes the hat-matrix diagonal and its derivative from blockwise constrained-\(G\) quantities, without explicitly forming the full projection matrix \(U\) or the full hat matrix \(H\). This keeps the computation aligned with the existing efficient constrained-fitting path.
- `gcv_gamma` remains available for backward compatibility but is now used only when `tuning_criterion = "gcv"`; it is ignored when `tuning_criterion = "loo"`.
- Updated `lgspline-details.R` and `lgspline.R` to describe the new criterion choice, the LOO formula, the working-scale interpretation for non-Gaussian/correlated tuning, and the post-optimization sample-size adjustment that now divides tuned penalties by \(((N+2)/(N-2))^2\) for both GCV- and LOO-based tuning.

#### 3 Tuning-gradient documentation clarifications
- Clarified the package documentation to distinguish between the criterion itself and its optimizer-facing gradients.
- Documented that the shared wiggle and flat tuning directions are differentiated directly, whereas predictor- and partition-specific penalties continue to use the existing lower-cost ratio approximation.
- Documented an empirical caveat for exact LOO tuning: the observation-wise analytic leverage derivative can be numerically delicate in some datasets even when the LOO criterion itself is stable. Users who prefer a more conservative path can set `use_custom_bfgs = FALSE` to rely on finite-difference optimization instead of the native analytic-gradient BFGS path.

#### 4 Large-sample tuning guidance
- Added a short practical note in the tuning documentation recommending `tuning_criterion = "gcv"` as the more practical choice for very large samples, with a rough guideline of `N > 250,000`.
- This note is descriptive only: it does not change defaults or behavior, and is included to set expectations about the relative computational cost of exact LOO versus GCV at large `N`.

#### 5 Small-sample post-optimization penalty adjustment
- Updated `tune_Lambda()` so the sample-size correction is now applied after both `gcv` and `loo` tuning paths, rather than only after `gcv`.
- The correction now decreases tuned penalties by dividing by `((N+2)/(N-2))^2` (equivalently multiplying by `((N-2)/(N+2))^2`), matching the intended response when small samples show occasional over-smoothing rather than under-smoothing.
- Added a regression test covering both criteria to verify that the returned penalties equal the optimizer solution times the reciprocal factor and are therefore smaller than the unadjusted optimized values.

## R CMD check results

0 errors | 0 warnings | 0 notes

There are no downstream dependencies for this package.
