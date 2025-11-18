# CRAN Comments

## Resubmission
This is a resubmission that addresses some previous documentation and example shortcomings in preparation for journal of statistical software submission.

Most importantly, detailed documentation on fitting default correlation structures and interpreting model output was provided. 

As well, the wiggle_penalty and flat_ridge_penalty were decreased in lgspline.R example when running on the volcano dataset.

No functional changes were made, just modifying documentation and improving the user-friendliness of an example. 


### Changes
1. **Documentation Issues**:
   * Included detailed documentation on default correlation structures and how to analyze + interpret model output, affects lgspline-details.R
   * Changed “list” to “vector” of penalties in documentation describing the return of optimized penalties Ie. “Optional list of custom penalties specified.” changed to “Optional vector of custom penalties specified”, for documentation of both "predictor_penalties" and "partition_penalties" arguments to "lgspline", affects lgspline.R
   * Dropped the "wiggle_penalty" and "flat_ridge_penalty" values for the example fitting to volcano dataset to 2e-7 and 1e-2 respectively for user friendliness



## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are currently no downstream dependencies for this package.
