# CRAN Comments

## Resubmission

This is a resubmission v1.0.2 addressing minor improvements and documentation to the package in preparation for Journal of Computational and Graphical Statistics submission.

It was found empirically that generalized cross-validation, even when correctly implemented and modified, yielded subpar results vs. simple leave-one-out for tuning this specific method.

The most important change is switching from generalized cross-validation to leave-one-out cross-validation for tuning by default, and corrected some math involving the computation of gradients for optimizing penalties using gradient-based optimization.

No functionality is changed otherwise. There should be no disruption to previous workflows or usage.

---

### Changes

#### 1 Modified GCV (`gcv_gamma`)
- Added a user-facing `gcv_gamma` argument to `lgspline()` and `lgspline.fit()`, defaulting to `1.2`, and threaded it through `tune_Lambda()` so automatic tuning now uses the modified GCV denominator
  \[
  N\{1 - \gamma\,\mathrm{tr}(H)/N\}^{2}
  \]
  rather than the ordinary `gamma = 1` form. Setting `gcv_gamma = 1` recovers standard GCV exactly.
- Updated `.compute_gcvu_gradient()` so the denominator derivative includes the same `gcv_gamma` factor, preserving consistency between the objective and the closed-form gradient used by the custom BFGS path.
- Documented the rationale in `lgspline-details.R` with the primary citation: Kim, Y.-J. and Gu, C. (2004), *Smoothing Spline Gaussian Regression: More Scalable Computation via Efficient Approximation*, *Journal of the Royal Statistical Society: Series B*, 66(2), 337-35. The default `gcv_gamma = 1.2` follows their recommended range for reducing occasional severe undersmoothing.

#### 2 Tuning criterion switch (`tuning_criterion`)
- Added a user-facing `tuning_criterion` argument to `lgspline()`, `lgspline.fit()`, and `tune_Lambda()`, with allowed values `"loo"` and `"gcv"`, defaulting to `"loo"`.
- The new default tuning path uses exact leave-one-out cross-validation on the same transformed or linearized tuning problem already used internally for GCV. For Gaussian identity-link tuning this is exact PRESS/LOO; for correlation structures and/or non-identity links it is computed on the same working scale already used by the package's tuning machinery.
- The exact LOO implementation computes the hat-matrix diagonal and its derivative from blockwise constrained-\(G\) quantities, without explicitly forming the full projection matrix \(U\) or the full hat matrix \(H\). This keeps the computation aligned with the existing efficient constrained-fitting path.
- `gcv_gamma` remains available for backward compatibility but is now used only when `tuning_criterion = "gcv"`; it is ignored when `tuning_criterion = "loo"`.
- Updated `lgspline-details.R` and `lgspline.R` to describe the new criterion choice, the LOO formula, the working-scale interpretation for non-Gaussian/correlated tuning, and the post-optimization sample-size adjustment that now divides tuned penalties by \(((N+2)/(N-2))^2\) for both GCV- and LOO-based tuning.
- Fixed several computational issues with previous calculation of gradients for tuning. 
- Clarified the package documentation for GCV vs. LOO.

#### 3 Small-sample post-optimization penalty adjustment
- Updated `tune_Lambda()` so the sample-size correction is now applied after both `gcv` and `loo` tuning paths, rather than only after `gcv`.
- The correction now decreases tuned penalties by dividing by `((N+2)/(N-2))^2` (equivalently multiplying by `((N-2)/(N+2))^2`), matching the intended response when small samples show occasional over-smoothing rather than under-smoothing.

#### 4 Plot formula display option for marginal plots
- Added a new plotting argument `include_all_terms_in_formulas` to `plot.lgspline()` and the internal plotting wrapper in `R/lgspline.R`.
- The default remains unchanged: when plotting a subset of predictors via `vars`, formula legends/hover text continue to show only the terms involving the plotted predictor(s), which is the more practical display for high-dimensional models.
- When `include_all_terms_in_formulas = TRUE` and `show_formulas = TRUE`, the plotted marginal curve is still drawn exactly as before, but the displayed formula now retains the full fitted partition equation, including non-plotted standalone linear terms and other fixed-variable terms.
- This change was made to improve interpretability for low-dimensional marginal plots without changing existing defaults or fitted values.

## R CMD check results

0 errors | 0 warnings | 0 notes

There are no downstream dependencies for this package.
