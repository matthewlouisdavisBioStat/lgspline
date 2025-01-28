#' Find Extremum of Fitted Lagrangian Multiplier Smoothing Spline
#'
#' S3 method for finding maximum or minimum of a fitted lgspline model using various optimization strategies.
#'
#' @param object A fitted lgspline model object
#' @param quick_heuristic Logical; whether to search only the top-performing partition. Default TRUE
#' @param initial Optional initial values for optimization, useful for fixing binary predictors. Default NULL
#' @param parallel Logical; whether to run optimization in parallel across partitions. Default FALSE
#' @param cl Optional cluster object for parallel processing. Default NULL
#' @param B_predict Optional custom coefficient matrix for prediction. Default NULL
#' @param minimize Logical; whether to find minimum instead of maximum. Default FALSE
#' @param stochastic Logical; whether to add noise for stochastic optimization. Default FALSE
#' @param stochastic_draw function; function for generating random draws from a distribution, analogous to posterior_predictive_draw. Takes two arguments, first is predicted values (mu), second is the square-root of sigmasq_tilde (sigma).
#' @param sigmasq_tilde Numeric; Variance for stochastic optimization. Default (object$sigmasq_tilde)
#' @param custom_objective_function function; Optional custom objective function for maximization/minimization. Default NULL.
#'   Function should take arguments:
#'   - mu: Prediction of response
#'   - sigma: Standard deviation
#'   - y_best: Best observed response
#' @param custom_objective_gradient Optional gradient function for custom optimization objective. Default NULL.
#'   Function should take arguments:
#'   - mu: Prediction of response
#'   - sigma: Standard deviation
#'   - y_best: Best observed response thus far
#'   - d_mu: Gradient of fitted function (use for chain-rule computations)
#'
#' @return A list containing:
#' \itemize{
#'   \item t: Input values at extremum
#'   \item y: Function value at extremum
#' }
#'
#' @details
#' This method finds the extremum (maximum or minimum) of the fitted Lagrangian multiplier smoothing splinepline.
#' It supports several optimization strategies:
#' \itemize{
#'   \item Quick heuristic search within top-performing partition
#'   \item Stochastic optimization for exploration
#'   \item Custom objective functions for advanced optimization
#' }
#'
#' Parallelism is experimental, and should not be used unless the user knows what they are doing.
#'
#' @examples
#' \dontrun{
#' # Find maximum of fitted model
#' max_point <- find_extremum(model_fit)
#'
#' # Find minimum
#' min_point <- find_extremum(model_fit, minimize = TRUE)
#'
#' }
#'
#' @export
find_extremum <- function(object,
                          quick_heuristic = TRUE,
                          initial = NULL,
                          parallel = FALSE,
                          cl = NULL,
                          B_predict = NULL,
                          minimize = FALSE,
                          stochastic = FALSE,
                          stochastic_draw = function(mu,
                                                     sigma){
                           N <- length(mu)
                           rnorm(
                             N, mu, sigma
                           )},
                          sigmasq_tilde = object$sigmasq_tilde,
                          custom_objective_function = NULL,
                          custom_objective_gradient = NULL) {
  ## Delegate to internal method
  object$find_extremum(
    quick_heuristic = quick_heuristic,
    initial = initial,
    parallel = parallel,
    cl = cl,
    B_predict = B_predict,
    minimize = minimize,
    stochastic = stochastic,
    stochastic_draw = stochastic_draw,
    sigmasq_tilde = sigmasq_tilde,
    custom_objective_function = custom_objective_function,
    custom_objective_gradient = custom_objective_gradient
  )
}

#' Generate Posterior Samples from Fitted Lagrangian Multiplier Smoothing Spline
#'
#' S3 method for generating posterior samples from a fitted lgspline model.
#'
#' @param object A fitted lgspline model object
#' @param new_sigmasq_tilde Dispersion parameter for sampling. Default object$sigmasq_tilde
#' @param new_predictors New data for posterior predictive sampling. Default object$X[[1]]
#' @param theta_1 First parameter for prior gamma distribution of inverse-dispersion. Default 0
#' @param theta_2 Second parameter for prior gamma distribution of inverse-dispersion. Default 0
#' @param posterior_predictive_draw Random number generation function. Default standard normal generator
#' @param draw_dispersion Logical; whether to sample dispersion parameter. Default TRUE
#' @param include_posterior_predictive Logical; whether to generate predictive samples. Default FALSE.
#' @param num_draws Number of posterior draws. Default 1
#'
#' @return A list containing:
#' \itemize{
#'   \item post_draw_coefficients: Matrix of posterior coefficient draws
#'   \item post_draw_sigmasq: Vector of posterior dispersion parameter draws
#'   \item post_pred_draw: Matrix of posterior predictive draws (if include_posterior_predictive = TRUE)
#' }
#'
#' @details
#' This method implements several posterior sampling strategies:
#' \itemize{
#'   \item Coefficient sampling using model variance-covariance matrix
#'   \item Optional dispersion parameter sampling
#' }
#'
#' @examples
#' \dontrun{
#' # Generate standard posterior draws
#' post_draws <- generate_posterior(model_fit)
#'
#' # Generate posterior with custom random number generator
#' custom_draws <- generate_posterior(
#'   model_fit,
#'   posterior_predictive_draw = function(N, mean, sqrt_dispersion)
#'     rt(N, df = 3, ncp = mean)
#' )
#' }
#'
#' @export
generate_posterior <-  function(object,
                                new_sigmasq_tilde = object$sigmasq_tilde,
                                new_predictors = object$X[[1]],
                                theta_1 = 0,
                                theta_2 = 0,
                                posterior_predictive_draw = function(N,
                                                 mean,
                                                 sqrt_dispersion,
                                                 ...){
                                  rnorm(N, mean, sqrt_dispersion)
                                },
                                draw_dispersion = TRUE,
                                include_posterior_predictive = FALSE,
                                num_draws = 1,
                                        ...) {
  # Delegate to the model's internal posterior generation method
  object$generate_posterior(
    new_sigmasq_tilde = new_sigmasq_tilde,
    new_predictors = new_predictors,
    theta_1 = theta_1,
    theta_2 = theta_2,
    posterior_predictive_draw = posterior_predictive_draw,
    draw_dispersion = draw_dispersion,
    include_posterior_predictive = include_posterior_predictive,
    num_draws = num_draws,
    ...
  )
}


#' Plot Method for Lagrangian Multiplier Smoothing Spline Models
#'
#' Visualize the fitted Lagrangian multiplier smoothing splinepline model
#'
#' @param object lgspline model fit object to plot (default uses current model)
#' @param show_formulas Logical, whether to display model formulas (default FALSE)
#' @param digits Number of decimal places for formula display (default 4)
#' @param legend_pos Position of legend for 1D plots (default "topright")
#' @param custom_response_lab Label for response variable (default "y")
#' @param custom_predictor_lab Label for predictor variable, x-axis, for one-dimensional plots only (default NULL, which gets converted to column name of predictor if so)
#' @param custom_predictor_lab1 Label for predictor variable, x1-axis, for two-dimensional plots only (default NULL, which gets converted to corresponding column name of predictors if so)
#' @param custom_predictor_lab2 Label for predictor variable, x2-axis, for two-dimensional plots only (default NULL, which gets converted to corresponding column name of predictors if so)
#' @param custom_formula_lab Label for predicted/fitted response on link-fxn scale. When NULL (by default), is set equal to "link(E[custom_response_lab])".
#' @param custom_title Plot title (default "Fitted Function")
#' @param text_size_formula Passes into cex argument of legend() for 1-D plots, and into hover font size of 2-D plots. When NULL (by default), it gets converted to 0.8 for 1-D, 8 for 2-D.
#' @param new_predictors Passes a matrix of new predictors predictors to plot (default uses in-sample fit)
#' @param xlim Passes into xlim argument of plot() if not NULL, for 1-D only
#' @param ylim Passes into the ylim argument of plot() if not NULL, for 1-D only
#' @param ... Additional arguments passed to plot/plotly for 1-D/2-D respectively
#'
#' @return
#' Generates a plot of the fitted model:
#' \itemize{
#'   \item 1D models: Line plot with partitions
#'   \item 2D models: 3D scatter plot
#'   \item Optional formula annotations
#' }
#'
#' @details
#' Plotting strategies:
#' \itemize{
#'   \item Supports 1D and 2D model visualizations
#'   \item Optional smooth interpolation
#'   \item Customizable plot aesthetics
#' }
#'
#' @examples
#' \dontrun{
#' # Standard plot with formulas
#' plot(model_fit, show_formulas = TRUE)
#'
#' }
#'
#' @export
plot.lgspline <- function(object,
                          show_formulas = FALSE,
                          digits = 4,
                          legend_pos = "topright",
                          custom_response_lab = "y",
                          custom_predictor_lab = NULL,
                          custom_predictor_lab1 = NULL,
                          custom_predictor_lab2 = NULL,
                          custom_formula_lab = NULL,
                          custom_title = "Fitted Function",
                          text_size_formula = NULL,
                          new_predictors = NULL,
                          xlim = NULL,
                          ylim = NULL,
                          ...) {
     object$plot(model_fit_in = object,
                 show_formulas = show_formulas,
                 digits = digits,
                 legend_pos = legend_pos,
                 custom_response_lab = custom_response_lab,
                 custom_predictor_lab = custom_predictor_lab,
                 custom_predictor_lab1 = custom_predictor_lab1,
                 custom_predictor_lab2 = custom_predictor_lab2,
                 custom_formula_lab = custom_formula_lab,
                 custom_title = custom_title,
                 text_size_formula = text_size_formula,
                 new_predictors = new_predictors,
                 xlim = xlim,
                 ylim = ylim,
                 ...)
}

#' Predict Method for Fitted Lagrangian Multiplier Smoothing Spline
#'
#' S3 method for generating predictions from a fitted lgspline model.
#'
#' @param object A fitted lgspline model object
#' @param newdata New data for prediction. Default NULL uses training data
#' @param parallel Logical; whether to use parallel processing. Default FALSE
#' @param cl Optional cluster object for parallel processing. Default NULL
#' @param chunk_size Size of computational chunks for parallel processing. Default NULL
#' @param num_chunks Number of chunks for parallel processing. Default NULL
#' @param rem_chunks Number of remainder chunks for parallel processing. Default NULL
#' @param B_predict Optional custom coefficient matrix for prediction. Default object$B
#' @param take_first_derivatives Logical; whether to compute first derivatives. Default FALSE
#' @param take_second_derivatives Logical; whether to compute second derivatives. Default FALSE
#' @param just_expansions Logical; whether to return only basis expansions. Default FALSE
#'
#' @return One of:
#' \itemize{
#'   \item Vector of predicted values (default)
#'   \item List containing predictions and derivatives (if derivatives requested)
#'   \item List of basis expansions (if just_expansions = TRUE)
#' }
#'
#' @details
#' This method implements several prediction strategies:
#' \itemize{
#'   \item Standard value prediction
#'   \item First and second derivative computation
#' }
#'
#' Parallelism for making predictions is experimental, and should not be used unless the user knows what they are doing.
#'
#' @examples
#' \dontrun{
#' ## Standard predictions
#' preds <- predict(model_fit, newdata)
#'
#' ## Compute first derivatives
#' deriv_preds1 <- predict(model_fit, newdata,
#'                       take_first_derivatives = TRUE)
#'
#' ## Compute second derivatives
#' deriv_preds2 <- predict(model_fit, newdata,
#'                       take_second_derivatives = TRUE)
#'
#' }
#'
#' @export
predict.lgspline <- function(object,
                             newdata = NULL,
                             parallel = FALSE,
                             cl = NULL,
                             chunk_size = NULL,
                             num_chunks = NULL,
                             rem_chunks = NULL,
                             B_predict = object$B,
                             take_first_derivatives = FALSE,
                             take_second_derivatives = FALSE,
                             just_expansions = FALSE) {
  # Delegate to the model's internal prediction method
  object$predict(
    new_predictors = newdata,
    parallel = parallel,
    cl = cl,
    chunk_size = chunk_size,
    num_chunks = num_chunks,
    rem_chunks = rem_chunks,
    B_predict = B_predict,
    take_first_derivatives = take_first_derivatives,
    take_second_derivatives = take_second_derivatives,
    just_expansions = just_expansions
  )
}


#' Extract Coefficients from Fitted Lagrangian Multiplier Smoothing Spline
#'
#' S3 method for extracting coefficients from a fitted lgspline model.
#'
#' @param object A fitted lgspline model object
#'
#' @return A list of coefficient vectors, one per partition, containing:
#' \itemize{
#'   \item Intercept terms
#'   \item Linear terms
#'   \item Quadratic terms
#'   \item Cubic terms
#'   \item Interaction terms (if present)
#' }
#'
#' @examples
#' \dontrun{
#' # Extract coefficients
#' coefficients <- coef(model_fit)
#'
#' # Print coefficients for first partition
#' print(coefficients[[1]])
#' }
#'
#' @export
coef.lgspline <- function(object) {
  object$B  # Return coefficients for each partition
}


#' Univariate Wald Tests for Lagrangian Multiplier Smoothing Spline
#'
#' S3 method for performing univariate Wald tests on individual coefficients.
#'
#' @param object A fitted lgspline model object
#' @param scale_vcovmat_by Scaling factor for variance-covariance matrix. Default 1
#' @param cv Critical value of confidence interval bounds. Default to the same as passed through the lgspline() function when fitting.
#'
#' @return A list containing:
#' \itemize{
#'   \item est: Coefficient estimates
#'   \item se: Standard errors
#'   \item stat: Test statistics
#'   \item interval_lb: Lower confidence bounds
#'   \item interval_ub: Upper confidence bounds
#'   \item pval: Two-sided p-values
#' }
#'
#' @details
#' Performs individual Wald tests for each coefficient.
#' Requires variance-covariance matrix (return_varcovmat = TRUE).
#'
#' @examples
#' \dontrun{
#'   wald_univariate(model_fit)
#' }
#'
#' @export
wald_univariate <- function(object, scale_vcovmat_by = 1, cv) {
  if(is.null(object$varcovmat)) {
    stop("Wald tests require return_varcovmat = TRUE during model fitting")
  }
  object$wald_univariate(scale_vcovmat_by, cv)
}


