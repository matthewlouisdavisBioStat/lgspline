# CRAN Comments

## Resubmission

This is a resubmission, lgspline version 1.0.3. Following some early feedback upon JCGS submission, 3 major issues were resolved:
1) The active set method for imposing inequality constraints, which the paper highlights, was avoided when modelling correlation structures
2) Some vector-vector and matrix-vector operations were being parallelized unecessarily in the tuning script, which actually slowed down parallel tuning by 2-5x
3) Inequality constraints are only enforced at certain observations 

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
- Users can now choose which observations to enforce specific QP constraints at

#### Parallelism
- a small issue in the auxiliary parallelism benchmark materials was corrected during revision; this concerned development-time timing/benchmarking content rather than the core CRAN-facing fitting routines
- related non-package benchmark artifacts used during development were removed from the submitted source tree so the CRAN submission now includes only package-relevant files
- this cleanup does not change the numerical fitting results reported by the package examples or tests; it only clarifies and streamlines the materials associated with the parallelism discussion

#### Miscellaneous
- Documentation updated throughout for clarity and consistency
- Small sample-size adjustment to tuning penalties changed from [(N+2)/(N-2)]^2 to (N+1)/(N-1) which is less dramatic, less controversial to reviewers.

## R CMD check results

0 errors | 0 warnings | 0 notes

There are no downstream dependencies for this package.
