# CRAN Comments

## Resubmission

This release v1.2 adds additive spline-group fitting. Separate formula terms such as
`y ~ spl(t1) + spl(t2)` now fit separate univariate smooths by default, each
with its own partitioning, penalties, smoothness constraints, and built-in
quadratic-programming constraints. Joined terms such as `spl(t1, t2)` keep the
previous shared-partition behavior. The low-level interface exposes the same
control through `spline_groups`. GAM-style formulas can now use `s()` as a
direct alias for `spl()` through the default `use_s_alias = TRUE` option.
Additive fits now include additive-aware post-fit wrappers for prediction,
plotting, equations, integration, likelihoods, Wald inference, posterior
simulation, extremum search, and leave-one-out calculations. Covariance,
posterior simulation, fitted standard errors, correlation-aware posterior
draws, and leave-one-out now use the full appended additive design with
block-diagonal penalties and active constraints. 

Additional documentation and user examples added.

## R CMD check results

0 errors | 0 warnings | 0 notess

There are no downstream dependencies for this package.
