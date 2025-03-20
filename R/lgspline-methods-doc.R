#' Find Extremum of Fitted Lagrangian Multiplier Smoothing Spline
#'
#' Finds global extrema of a fitted lgspline model using deterministic or stochastic
#' optimization strategies. Supports custom objective functions for advanced applications like
#' Bayesian optimization acquisition functions.
#'
#' @param object A fitted lgspline model object containing partition information and fitted values
#' @param vars Vector; A vector of numeric indices (or character variable names) of predictors to optimize for. If NULL (by default), all predictors will be optimized.
#' @param quick_heuristic Logical; whether to search only the top-performing partition. When TRUE (default),
#'        optimizes within the best partition. When FALSE, initiates searches from all partition local maxima.
#' @param initial Numeric vector; Optional initial values for optimization. Useful for fixing binary
#'        predictors or providing starting points. Default NULL
#' @param B_predict Matrix; Optional custom coefficient list for prediction. Useful for posterior
#'        draws in Bayesian optimization. Default NULL
#' @param minimize Logical; whether to find minimum instead of maximum. Default FALSE
#' @param stochastic Logical; whether to add noise for stochastic optimization. Enables better
#'        exploration of the function space. Default FALSE
#' @param stochastic_draw Function; Generates random noise/modifies predictions for stochastic optimization, analogous to
#'        posterior_predictive_draw. Takes three arguments:
#'        \itemize{
#'          \item mu: Vector of predicted values
#'          \item sigma: Vector of standard deviations (square-root of sigmasq_tilde)
#'          \item ...: Additional arguments to pass through
#'        }
#'        Default \code{rnorm(length(mu), mu, sigma)}
#' @param sigmasq_predict Numeric; Variance parameter for stochastic optimization. Controls
#'        the magnitude of random perturbations. Defaults to object$sigmasq_tilde
#' @param custom_objective_function Function; Optional custom objective function for optimization.
#'        Takes arguments:
#'        \itemize{
#'          \item mu: Vector of predicted response values
#'          \item sigma: Vector of standard deviations
#'          \item y_best: Numeric; Best observed response value
#'          \item ...: Additional arguments passed through
#'        }
#'        Default NULL
#' @param custom_objective_derivative Function; Optional gradient function for custom optimization
#'        objective. Takes arguments:
#'        \itemize{
#'          \item mu: Vector of predicted response values
#'          \item sigma: Vector of standard deviations
#'          \item y_best: Numeric; Best observed response value
#'          \item d_mu: Gradient of fitted function (for chain-rule computations)
#'          \item ...: Additional arguments passed through
#'        }
#'        Default NULL
#'
#' @details
#' This method finds extrema (maxima or minima) of the fitted function or composite functions
#' of the fit. The optimization process can be customized through several approaches:
#' \itemize{
#'   \item Partition-based search: Either focuses on the top-performing partition (quick_heuristic = TRUE)
#'         or searches across all partition local maxima
#'   \item Stochastic optimization: Adds random noise during optimization for better exploration
#'   \item Custom objectives: Supports user-defined objective functions and gradients for
#'         specialized optimization tasks like Bayesian optimization
#' }
#'
#' @return A list containing:
#' \itemize{
#'   \item{t}{Numeric vector of input values at the extremum}
#'   \item{y}{Numeric value of the objective function at the extremum}
#' }
#'
#' @examples
#' \donttest{
#' ## Basic usage with simulated data
#' set.seed(1234)
#' x <- runif(1000, -10, 10)
#' y <- 2*sin(x) + -0.06*x^2 + rnorm(length(x))
#' model_fit <- lgspline(x, y)
#' plot(model_fit)
#'
#' ## Find global maximum and minimum
#' max_point <- find_extremum(model_fit)
#' min_point <- find_extremum(model_fit, minimize = TRUE)
#' abline(v = max_point$t, col = 'blue')  # Add maximum point
#' abline(v = min_point$t, col = 'red')   # Add minimum point
#'
#' ## Advanced usage: custom objective functions
#' # expected improvement acquisition function
#' ei_custom_objective_function = function(mu, sigma, y_best, ...) {
#'   d <- y_best - mu
#'   d * pnorm(d/sigma) + sigma * dnorm(d/sigma)
#' }
#' # derivative of ei
#' ei_custom_objective_derivative = function(mu, sigma, y_best, d_mu, ...) {
#'   d <- y_best - mu
#'   z <- d/sigma
#'   d_z <- -d_mu/sigma
#'   pnorm(z)*d_mu - d*dnorm(z)*d_z + sigma*z*dnorm(z)*d_z
#' }
#'
#' ## Single iteration of Bayesian optimization
#' post_draw <- generate_posterior(model_fit)
#' acq <- find_extremum(model_fit,
#'                      stochastic = TRUE,  # Enable stochastic exploration
#'                      B_predict = post_draw$post_draw_coefficients,
#'                      sigmasq_predict = post_draw$post_draw_sigmasq,
#'                      custom_objective_function = ei_custom_objective_function,
#'                      custom_objective_derivative = ei_custom_objective_derivative)
#' abline(v = acq$t, col = 'green')  # Add acquisition point
#' }
#'
#' @seealso
#' \code{\link{lgspline}} for fitting the model,
#' \code{\link{generate_posterior}} for generating posterior draws
#'
#' @export
find_extremum <- function(object,
                          vars = NULL,
                          quick_heuristic = TRUE,
                          initial = NULL,
                          B_predict = NULL,
                          minimize = FALSE,
                          stochastic = FALSE,
                          stochastic_draw = function(mu,
                                                     sigma){
                            N <- length(mu)
                            rnorm(
                              N, mu, sigma
                            )},
                          sigmasq_predict = object$sigmasq_tilde,
                          custom_objective_function = NULL,
                          custom_objective_derivative = NULL,
                          ...) {
  ## Delegate to internal method
  object$find_extremum(
    vars = vars,
    quick_heuristic = quick_heuristic,
    initial = initial,
    B_predict = B_predict,
    minimize = minimize,
    stochastic = stochastic,
    stochastic_draw = stochastic_draw,
    sigmasq_predict = sigmasq_predict,
    custom_objective_function = custom_objective_function,
    custom_objective_derivative = custom_objective_derivative
  )
}

#' Generate Posterior Samples from Fitted Lagrangian Multiplier Smoothing Spline
#'
#' Draws samples from the posterior distribution of model parameters and optionally generates
#' posterior predictive samples. Uses Laplace approximation for non-Gaussian responses.
#'
#' @param object A fitted lgspline model object containing model parameters and fit statistics
#' @param new_sigmasq_tilde Numeric; Dispersion parameter for sampling. Controls variance of
#'        posterior draws. Default object$sigmasq_tilde
#' @param new_predictors Matrix; New data matrix for posterior predictive sampling. Should match
#'        structure of original predictors. Default = predictors as input to \code{lgspline}.
#' @param theta_1 Numeric; Shape parameter for prior gamma distribution of inverse-dispersion.
#'        Default 0 implies uniform prior
#' @param theta_2 Numeric; Rate parameter for prior gamma distribution of inverse-dispersion.
#'        Default 0 implies uniform prior
#' @param posterior_predictive_draw Function; Random number generator for posterior predictive
#'        samples. Takes arguments:
#'        \itemize{
#'          \item N: Integer; Number of samples to draw
#'          \item mean: Numeric vector; Predicted mean values
#'          \item sqrt_dispersion: Numeric vector; Square root of dispersion parameter
#'          \item ...: Additional arguments passed through
#'        }
#' @param draw_dispersion Logical; whether to sample the dispersion parameter from its
#'        posterior distribution. When FALSE, uses point estimate. Default TRUE
#' @param include_posterior_predictive Logical; whether to generate posterior predictive
#'        samples for new observations. Default FALSE
#' @param num_draws Integer; Number of posterior draws to generate. Default 1
#'
#' @details
#' Implements posterior sampling using the following approach:
#' \itemize{
#'   \item Coefficient posterior: Assumes sqrt(N)B ~ N(Btilde, sigma^2UG)
#'   \item Dispersion parameter: Sampled from inverse-gamma distribution with user-specified
#'         prior parameters (theta_1, theta_2) and model-based sufficient statistics
#'   \item Posterior predictive: Generated using custom sampling function, defaulting to
#'         Gaussian for standard normal responses
#' }
#'
#' For the dispersion parameter, the sampling process follows:
#' \preformatted{
#' shape <- theta_1 + 0.5 * model_fit$N
#' rate <- theta_2 + 0.5 * sum((model_fit$y - model_fit$ytilde)^2)
#' post_draw_sigmasq <- 1/rgamma(1, shape, rate)
#' }
#'
#' Users can modify sufficient statistics by adjusting theta_1 and theta_2 relative to
#' the default model-based values.
#'
#' @return A list containing:
#' \itemize{
#'   \item{post_draw_coefficients}{List of length num_draws containing posterior coefficient samples}
#'   \item{post_draw_sigmasq}{List of length num_draws containing posterior dispersion parameter
#'         samples (or repeated point estimate if draw_dispersion = FALSE)}
#'   \item{post_pred_draw}{List of length num_draws containing posterior predictive samples
#'         (only if include_posterior_predictive = TRUE)}
#' }
#'
#' @examples
#' \donttest{
#' ## Generate example data
#' x <- runif(1000, -10, 10)
#' true_y <- 2*sin(x) + -0.06*x^2
#' y <- rnorm(length(true_y), true_y, 1)
#'
#' ## Fit model (using unstandardized expansions for consistent inference)
#' model_fit <- lgspline(x, y,
#'                      K = 7,
#'                      standardize_expansions_for_fitting = FALSE)
#'
#' ## Compare Wald (= t-intervals here) to Monte Carlo credible intervals
#' # Get Wald intervals
#' wald <- wald_univariate(model_fit,
#'                         cv = qt(0.975, df = model_fit$trace_XUGX))
#' wald_bounds <- cbind(wald[,"lower_ci"], wald[,"upper_ci"])
#'
#' ## Generate posterior samples (uniform prior)
#' post_draws <- generate_posterior(model_fit,
#'                                theta_1 = -1,
#'                                theta_2 = 0,
#'                                num_draws = 2000)
#'
#' ## Convert to matrix and compute credible intervals
#' post_mat <- Reduce('cbind',
#'                   lapply(post_draws$post_draw_coefficients,
#'                         function(x) Reduce("rbind", x)))
#' post_bounds <- t(apply(post_mat, 1, quantile, c(0.025, 0.975)))
#'
#' ## Compare intervals
#' print(round(cbind(wald_bounds, post_bounds), 4))
#' }
#'
#' @seealso
#' \code{\link{lgspline}} for model fitting,
#' \code{\link{wald_univariate}} for Wald-type inference
#'
#' @export
generate_posterior <- function(object,
                               new_sigmasq_tilde = object$sigmasq_tilde,
                               new_predictors = object$X[[1]],
                               theta_1 = 0,
                               theta_2 = 0,
                               posterior_predictive_draw = function(N, mean,
                                                                    sqrt_dispersion, ...) {
                                 rnorm(N, mean, sqrt_dispersion)
                               },
                               draw_dispersion = TRUE,
                               include_posterior_predictive = FALSE,
                               num_draws = 1,
                               ...) {

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
#' Creates visualizations of fitted spline models, supporting both 1D line plots and 2D surface
#' plots with optional formula annotations and customizable aesthetics.
#'
#' @param object A fitted lgspline model object containing the model fit to be plotted
#' @param show_formulas Logical; whether to display analytical formulas for each partition.
#'        Default FALSE
#' @param digits Integer; Number of decimal places for coefficient display in formulas.
#'        Default 4
#' @param legend_pos Character; Position of legend for 1D plots ("top", "bottom", "left",
#'        "right", "topleft", etc.). Default "topright"
#' @param custom_response_lab Character; Label for response variable axis. Default "y"
#' @param custom_predictor_lab Character; Label for predictor axis in 1D plots. If NULL
#'        (default), uses predictor column name
#' @param custom_predictor_lab1 Character; Label for first predictor axis (x1) in 2D plots.
#'        If NULL (default), uses first predictor column name
#' @param custom_predictor_lab2 Character; Label for second predictor axis (x2) in 2D plots.
#'        If NULL (default), uses second predictor column name
#' @param custom_formula_lab Character; Label for fitted response on link function scale.
#'        If NULL (default), uses "link(E[custom_response_lab])" for non-Gaussian models
#'        with non-identity link, otherwise uses custom_response_lab
#' @param custom_title Character; Main plot title. Default "Fitted Function"
#' @param text_size_formula Numeric; Text size for formula display. Passed to cex in legend()
#'        for 1D plots and hover font size for 2D plots. If NULL (default), uses 0.8 for 1D
#'        and 8 for 2D
#' @param legend_args List; Additional arguments passed to legend() for 1D plots
#' @param new_predictors Matrix; Optional new predictor values for prediction. If NULL
#'        (default), uses original fitting data
#' @param xlim Numeric vector; Optional x-axis limits for 1D plots. Default NULL
#' @param ylim Numeric vector; Optional y-axis limits for 1D plots. Default NULL
#' @param color_function Function; Returns colors for plotting by partition, must return K+1 vector of valid colors. Defaults to NULL, in which case \code{rainbow(K+1)} is used for 1D and \code{color_function <- colorRampPalette(RColorBrewer::brewer.pal(8, "Spectral"))(K+1)} used in multiple.
#' @param add Logical; If TRUE, adds to existing plot (1D only). Similar to add in
#'        \code{\link[graphics]{hist}}. Default FALSE
#' @param vars Numeric or character vector; Optional indices for selecting variables to plot. Can either be numeric (the column indices of "predictors" or "data") or character (the column names, if available from "predictors" or "data")
#' @param ... Additional arguments passed to underlying plot functions:
#'        \itemize{
#'          \item 1D: Passed to \code{\link[graphics]{plot}}
#'          \item 2D: Passed to \code{\link[plotly]{plot_ly}}
#'        }
#'
#' @details
#' Produces different visualizations based on model dimensionality:
#' \itemize{
#'   \item 1D models: Line plot showing fitted function across partitions, with optional
#'         data points and formula annotations
#'   \item 2D models: Interactive 3D surface plot using plotly, with hover text showing
#'         predicted values and optional formula display
#' }
#'
#' Partition boundaries are indicated by color changes in both 1D and 2D plots.
#'
#' When plotting using "select_vars" option, it is recommended to use the
#' "new_predictors" argument to set all terms not involved with plotting to 0
#' to avoid non-sensical results. But for some cases, it may be useful to set
#' other predictors fixed at certain values. By default, observed values in the
#' data set are used.
#'
#' The function relies on linear expansions being present - if (for example) a
#' user includes the argument "_1_" or "_2_" in "exclude_these_expansions", then
#' this function will not be able to extract the predictors needed for plotting.
#'
#' For this case, try constraining the effects of these terms to 0 instead using
#' "constraint_vectors" and "constraint_values" argument, so they are kept in
#' the expansions but their corresponding coefficients will be 0.
#'
#' @return
#' Returns
#' \itemize{
#'   \item 1D: Base R plot object showing fitted curve and optional formulas
#'   \item 2D: Plotly object showing interactive surface plot
#' }
#'
#' @examples
#' \donttest{
#' ## Generate example data
#' set.seed(1234)
#' x <- runif(1000, -10, 10)
#' y <- 2*sin(x) + -0.06*x^2 + rnorm(length(x))
#'
#' ## Fit model with 10 partitions
#' model_fit <- lgspline(x, y, K = 9)
#'
#' ## Basic plot
#' plot(model_fit)
#'
#' ## Customized plot with formulas
#' plot(model_fit,
#'      show_formulas = TRUE,           # Show partition formulas
#'      custom_response_lab = 'Price',   # Custom axis labels
#'      custom_predictor_lab = 'Size',
#'      custom_title = 'Price vs Size',  # Custom title
#'      digits = 2,                      # Round coefficients
#'      legend_pos = 'bottom',           # Move legend
#'      text_size_formula = 0.375,       # Adjust formula text size
#'      pch = 16,                        # Point style
#'      cex.main = 1.25)                 # Title size
#' }
#'
#' @seealso
#' \code{\link{lgspline}} for model fitting,
#' \code{\link[graphics]{plot}} for additional 1D plot parameters,
#' \code{\link[plotly]{plot_ly}} for additional 2D plot parameters
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
                          legend_args = list(),
                          new_predictors = NULL,
                          xlim = NULL,
                          ylim = NULL,
                          color_function = NULL,
                          add = FALSE,
                          vars = c(),
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
              legend_args = legend_args,
              new_predictors = new_predictors,
              xlim = xlim,
              ylim = ylim,
              color_function = color_function,
              add = add,
              vars = vars,
              ...)
}

#' Predict Method for Fitted Lagrangian Multiplier Smoothing Spline
#'
#' Generates predictions, derivatives, and basis expansions from a fitted lgspline model.
#' Supports both in-sample and out-of-sample prediction with optional parallel processing.
#'
#' @param object A fitted lgspline model object containing model parameters and fit
#' @param newdata Matrix or data.frame; New predictor values for out-of-sample prediction.
#'        If NULL (default), uses training data
#' @param parallel Logical; whether to use parallel processing for prediction computations.
#'        Experimental feature - use with caution. Default FALSE
#' @param cl Optional cluster object for parallel processing. Required if parallel=TRUE.
#'        Default NULL
#' @param chunk_size Integer; Size of computational chunks for parallel processing.
#'        Default NULL
#' @param num_chunks Integer; Number of chunks for parallel processing. Default NULL
#' @param rem_chunks Integer; Number of remainder chunks for parallel processing.
#'        Default NULL
#' @param B_predict Matrix; Optional custom coefficient matrix for prediction.
#'        Default object$B
#' @param take_first_derivatives Logical; whether to compute first derivatives of the
#'        fitted function. Default FALSE
#' @param take_second_derivatives Logical; whether to compute second derivatives of the
#'        fitted function. Default FALSE
#' @param expansions_only Logical; whether to return only basis expansions without
#'        computing predictions. Default FALSE
#'
#' @details
#' Implements multiple prediction capabilities:
#' \itemize{
#'   \item Standard prediction: Returns fitted values for new data points
#'   \item Derivative computation: Calculates first and/or second derivatives
#'   \item Basis expansion: Returns design matrix of basis functions
#'   \item Correlation structures: Supports non-Gaussian GLM correlation via
#'         variance-covariance matrices
#' }
#'
#' The parallel processing feature is experimental and should be used with caution.
#' When enabled, computations are split across chunks and processed in parallel,
#' which may improve performance for large datasets.
#'
#' @return Depending on the options selected, returns one of:
#' \itemize{
#'   \item Numeric vector of predicted values (default)
#'   \item List containing:
#'     \itemize{
#'       \item predictions: Numeric vector of fitted values
#'       \item first_deriv: Numeric vector of first derivatives (if requested)
#'       \item second_deriv: Numeric vector of second derivatives (if requested)
#'     }
#'   \item List of basis expansions (if expansions_only = TRUE)
#' }
#'
#' @examples
#' \donttest{
#' ## Generate example data
#' set.seed(1234)
#' x <- runif(1000, -10, 10)
#' y <- 2*sin(x) + -0.06*x^2 + rnorm(length(x))
#'
#' ## Fit model
#' model_fit <- lgspline(x, y)
#'
#' ## Generate predictions for new data
#' newdata <- sort(rnorm(10000))
#' preds <- predict(model_fit, newdata)
#'
#' ## Compute derivatives
#' deriv1 <- predict(model_fit, newdata,
#'                  take_first_derivatives = TRUE)
#' deriv2 <- predict(model_fit, newdata,
#'                  take_second_derivatives = TRUE)
#'
#' ## Visualize results
#' layout(matrix(c(1,1,2,2,3,3,1,1,2,2,3,3),
#'              byrow = TRUE, nrow = 2))
#'
#' ## Plot function
#' plot(newdata, preds,
#'      main = 'Fitted Function',
#'      xlab = 't',
#'      ylab = "f(t)")
#'
#' ## Plot first derivative
#' plot(newdata,
#'      deriv1$first_deriv,
#'      main = 'First Derivative',
#'      xlab = 't',
#'      ylab = "f'(t)")
#'
#' ## Plot second derivative
#' plot(newdata,
#'      deriv2$second_deriv,
#'      main = 'Second Derivative',
#'      xlab = 't',
#'      ylab = "f''(t)")
#' }
#'
#' @seealso
#' \code{\link{lgspline}} for model fitting,
#' \code{\link{plot.lgspline}} for visualizing predictions
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
                             expansions_only = FALSE) {
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
    expansions_only = expansions_only
  )
}


#' Extract Coefficients from Fitted Lagrangian Multiplier Smoothing Spline
#'
#' Extracts polynomial coefficients for each partition from a fitted lgspline model. Returns
#' either a list of partition-specific coefficients.
#'
#' @param object A fitted lgspline model object containing coefficient vectors.
#'
#' @details
#' For each partition, coefficients represent a polynomial expansion of the predictor(s) by column index, for example:
#' \itemize{
#'   \item intercept: Constant term
#'   \item _1_: Linear term
#'   \item _1_^2: Quadratic term
#'   \item _1_^3: Cubic term
#'   \item Additional interaction terms for multivariate models
#' }
#'
#' If column/variable names are present, indices will be replaced with column/variable names.
#'
#' Coefficients can be accessed either as separate vectors per partition or combined into
#' a single matrix using \code{Reduce('cbind', coef(model_fit))}.
#'
#' @return
#' A list of matrices, one per partition, where:
#' \itemize{
#'   \item Each vector contains coefficient values as a single-column matrix
#'   \item Row names indicate the term type (intercept, linear, quadratic, cubic)
#'   \item List elements are named "partition1", "partition2", etc.
#'   \item Each matrix has an attribute "names" identifying the terms in addition to "rownames"
#' }
#'
#' @examples
#' \donttest{
#'
#' ## Simulate some data and fit using default settings
#' set.seed(1234)
#' x <- runif(1000, -10, 10)
#' y <- 2*sin(x) + -0.06*x^2 + rnorm(length(x))
#' model_fit <- lgspline(x, y)
#'
#' ## Extract coefficients
#' coefficients <- coef(model_fit)
#'
#' ## Print coefficients for first partition
#' print(coefficients[[1]])
#'
#' ## Compare coefficients across all partitions
#' print(Reduce('cbind', coefficients))
#' }
#'
#' @export
coef.lgspline <- function(object) {
  object$B  # Return coefficients for each partition
}


#' Univariate Wald Tests and Confidence Intervals for Lagrangian Multiplier Smoothing Splines
#'
#' Performs coefficient-specific Wald tests and constructs confidence intervals for fitted
#' lgspline models, providing inference for each partition and term within partitions.
#'
#' @param object A fitted lgspline model object containing coefficient estimates and
#'        variance-covariance matrix (requires return_varcovmat = TRUE in fitting)
#' @param scale_vcovmat_by Numeric; Scaling factor for variance-covariance matrix.
#'        Adjusts standard errors and test statistics. Default 1
#' @param cv Numeric; Critical value for confidence interval construction. Defaults to
#'        value specified in lgspline() fit. Common choices:
#'        \itemize{
#'          \item qnorm(0.975) for normal-based 95% intervals
#'          \item qt(0.975, df) for t-based 95% intervals, where df = N - trace(XUGX)
#'        }
#'
#' @details
#' For each coefficient, provides:
#' \itemize{
#'   \item Point estimates
#'   \item Standard errors from the model's variance-covariance matrix
#'   \item Two-sided Wald test statistics and p-values
#'   \item Confidence intervals using specified critical values
#' }
#'
#' @return A data frame with rows for each coefficient (across all partitions) and columns:
#' \itemize{
#'   \item{estimate}{Numeric; Coefficient estimate}
#'   \item{std_error}{Numeric; Standard error}
#'   \item{statistic}{Numeric; Wald test statistic (estimate/std_error)}
#'   \item{p_value}{Numeric; Two-sided p-value based on normal distribution}
#'   \item{lower_ci}{Numeric; Lower confidence bound (estimate - cv*std_error)}
#'   \item{upper_ci}{Numeric; Upper confidence bound (estimate + cv*std_error)}
#' }
#'
#' @examples
#' \donttest{
#'
#' ## Simulate some data and fit using default settings
#' set.seed(1234)
#' x <- runif(1000, -10, 10)
#' y <- 2*sin(x) + -0.06*x^2 + rnorm(length(x))
#' model_fit <- lgspline(x, y)
#'
#' ## The critical value defaults to qnorm(0.975),
#' # but with the cv (critical value) argument, we can set this to the
#' # critical value of a t-distribution, and obtain a 95% t-interval.
#' # Effective degrees of freedom = N - trace(XUGX^T)
#' wald_univariate(
#'   model_fit,
#'   cv = qt(0.975, model_fit$N - model_fit$trace_XUGX)
#' )
#' }
#'
#' @export
wald_univariate <- function(object, scale_vcovmat_by = 1, cv) {
  if(is.null(object$varcovmat)) {
    stop("Wald tests require return_varcovmat = TRUE during model fitting")
  }
  object$wald_univariate(scale_vcovmat_by, cv)
}


