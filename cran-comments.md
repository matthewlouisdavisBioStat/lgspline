## Resubmission

This is a resubmission. In this version I have addressed all the points raised by CRAN reviewers (1-5) plus additional documentation issues (6). Two minor changes to functionality (7) were also made.

### 1. Acronym explanation
* Expanded "AFT" to "accelerated failure time (AFT)" throughout the package
* Updated all Weibull AFT-related helper function documentation

### 2. References in DESCRIPTION
* Added properly formatted references with DOIs and ISBNs:
  - Ezhov et al. (2018) <doi:10.1515/jag-2017-0029> for Lagrangian multipliers approach
  - Searle et al. (2009) <ISBN:978-0470009598> for correlation structure estimation
  - Nocedal & Wright (2006) <doi:10.1007/978-0-387-40065-5> for quadratic programming
  - Wahba (1990) <doi:10.1137/1.9781611970128> for smoothing splines
  - Wood (2006) <ISBN:978-1584884743> for generalized additive models
* Added WORDLIST to /inst/ folder and incorporated spell checks to avoid notes on author names and et al

### 3. Added missing \value tags
* Added detailed return value descriptions to:
  - lgspline.fit.Rd: Comprehensive list of returned components
  - print.lgspline.Rd: Documented invisible return of original object
  - print.summary.lgspline.Rd: Documented invisible return of summary object
  - summary.lgspline.Rd: Detailed list of summary components returned

### 4. Fixed example code issues
* Replaced "dontrun" with "donttest" for examples that take longer than 5 seconds
* Unwrapped examples executable in < 5 seconds 
* Applied these changes to:
  - lgspline.Rd main example documentation
  - Methods: find_extremum, generate_posterior, plot.lgspline, predict.lgspline, coef.lgspline, wald.univariate
  - Additional functions: leave_one_out, prior_loglik

### 5. Fixed par() reset in examples
* Added par(oldpar) at the end of the predict.lgspline.Rd example
* Improved layout with ncol=2 instead of nrow=2 for better visualization in predict.lgspline

### 6. Additional Documentation Improvements
* Standardized notation across all examples (using 't' for predictors instead of 'x') in methods.R, prior_loglik.R, and leave_one_out.R.
* Fixed function calls in examples:
  - Added log=TRUE parameter to dnorm() in prior_loglik example
  - Corrected PRESS calculation in leave_one_out example to use "(y-loo)^2" correctly in place of "loo^2"
* Removed internal notes from documentation of print.summary.lgspline
* Fixed parameter descriptions in print.summary.lgspline
* Moved first 3D plotting example into "donttest" for lgspline examples
* Modified description of lgspline to change "effecs" to "effects" and inserting "interpretable" as a qualifier for interaction/non-spline terms.
* Corrected description of shape/rate terms of "generate_posterior" to be consistent with actual implementation
* Clarified in wald_univariate that for Gaussian response with identity link, t-intervals/tests/statistics are default over Wald analogous. 
* Changed parallel example of lgspline to use only 1 core, and run with K=1 instead of K=15 for efficiency and compatibility
* Added inst/WORDLIST folder to prevent notes on misspelled words that are mostly author names and "et al". 
* Removed donttest for Weibull helper functions 
* provided detailed examples for Weibull helper functions


### 7. Function Improvements
* Caught error: changed "std(max_min_C)" to be "std_X(max_min_C)" (line 6249 in new code, line 6207 previous). This ensures that the scaling applied to polynomial expansions matches the scaling applied to corresponding penalties.
* Inserted condition such that if both newdata and new_predictors are left NULL for "predict.lgspline", that input used for fitting will be used instead by default instead of an error being returned. I clarified this in documentation. This is for user-friendliness. 

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are currently no downstream dependencies for this package.