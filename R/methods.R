#' lgspline: Lagrangian Multiplier Smoothing Splines
#'
#' @description
#' Allows for common S3 methods including print, summary, coef, plot, predict,
#' confint, and logLik, with additional inference methods provided.
#'
#' @details
#' This package implements various methods for working with lgspline models,
#' including printing, summarizing, plotting, and conducting statistical inference.
#'
#' @docType methods
#' @keywords internal
#' @name lgspline-methods
#' @rdname lgspline-methods
#' @aliases lgspline-methods
NULL

#' Print Method for lgspline Objects
#'
#' @description
#' Provides a standard print method for lgspline model objects to display
#' key model characteristics.
#'
#' @return Invisibly returns the original \code{lgspline} object \code{x}. This
#' function is called for printing a concise summary of the fitted model's key
#' characteristics (family, link, N, predictors, partitions, basis functions) to
#' the console.
#'
#' @param x An lgspline model object
#' @param ... Additional arguments (not used)
#'
#' @export
#' @method print lgspline
print.lgspline <- function(x, ...) {
  cat("Lagrangian Multiplier Smoothing Spline Model\n")
  cat("============================================\n")
  cat("Model Family (Link Function):",
      paste0(paste0(x$family)[1],
             " (",
             paste0(x$family)[2],
             ")", collapse = ""),
      "\n")
  cat("Number of Observations:", x$N, "\n")
  cat("Number of Predictors:", x$q, "\n")
  cat("Number of Partitions:", x$K + 1, "\n")
  cat("Basis Functions per Partition:", x$p, "\n")
  invisible(x)
}

#' Summary method for lgspline Objects
#' @param object An lgspline model object
#' @param ... Not used.
#' @return An object of class \code{summary.lgspline}. This object is a list
#' containing detailed information from \code{lgspline} fit, prepared for
#' display. Its main components are:
#' \describe{
#'  \item{model_family}{The \code{\link[stats]{family}} object or custom list specifying the distribution and link.}
#'  \item{observations}{The number of observations (N) used in the fit.}
#'  \item{predictors}{The number of original predictor variables (q) supplied.}
#'  \item{knots}{The number of partitions (K+1) minus 1.}
#'  \item{basis_functions}{The number of basis functions (coefficients) estimated per partition (p).}
#'  \item{estimate_dispersion}{A character string ("Yes" or "No") indicating if the dispersion parameter was estimated.}
#'  \item{cv}{The critical value (\code{critical_value} from the fit) used by the \code{print.summary.lgspline} method for confidence intervals.}
#'  \item{coefficients}{A matrix summarizing univariate inference results. Columns typically include 'Estimate', 'Std. Error', test statistic ('t value' or 'z value'), 'CI LB', 'CI UB', and 'Pr(>|t|)' or 'Pr(>|z|)'. This table is fully populated only if \code{return_varcovmat=TRUE} was set in the original \code{lgspline} call. Otherwise, it defaults to a single column of estimates.}
#'  \item{sigmasq_tilde}{The estimated (or fixed) dispersion parameter, \eqn{\tilde{\sigma}^2}.}
#'  \item{trace_XUGX}{The calculated trace term \eqn{\mathrm{trace}(\mathbf{X}\mathbf{U}\mathbf{G}\mathbf{X}^{\top})}, related to effective degrees of freedom.}
#'  \item{N}{Number of observations (N), re-included for convenience and printing.}
#' }
#' @export
summary.lgspline <- function(object, ...) {
  ## Create a brief summary
  summary_list <- list(
    model_family = object$family,
    observations = object$N,
    predictors = object$q,
    knots = object$K,
    basis_functions = object$p,
    estimate_dispersion = ifelse(!is.null(object$estimate_dispersion) && object$estimate_dispersion &&
                                   !is.null(object$sigmasq_tilde) && object$sigmasq_tilde != 1,
                                 'Yes',
                                 'No'),
    cv = object$critical_value,
    coefficients = NULL,
    sigmasq_tilde = object$sigmasq_tilde,
    trace_XUGX = object$trace_XUGX,
    N = object$N
  )

  ## Typical summaries for Wald inference, like lm() or glm()
  # [Change 2026-02-11] ensure coefficient matrix column names are
  # compatible with stats::printCoefmat() for proper p-value formatting.
  # Use exported wald_univariate() wrapper instead of
  # internal method to get properly normalized column ordering.
  if(!is.null(object$return_varcovmat) &&
     object$return_varcovmat &&
     !is.null(object$wald_univariate)){
    tr <- try({
      wald_obj <- wald_univariate(object)
      wald_obj$coefficients
    }, silent = TRUE)

    if(!inherits(tr, 'try-error') && !is.null(tr)){
      summary_list$coefficients <- tr
    } else {
      summary_list$coefficients <- cbind(Estimate = unlist(object$B))
    }
  } else {
    summary_list$coefficients <- cbind(Estimate = unlist(object$B))
  }

  class(summary_list) <- "summary.lgspline"
  return(summary_list)
}


#' Print Method for lgspline Object Summaries
#'
#' @description
#' Displays a formatted summary of the fitted \code{lgspline} model to the
#' console. Uses \code{\link[stats]{printCoefmat}} for coefficient tables to
#' obtain standard formatting of p-values and significance codes.
#'
#' @param x A summary.lgspline object, the result of calling \code{summary()} on an \code{lgspline} object.
#' @param ... Not used.
#' @return Invisibly returns the original \code{summary.lgspline} object \code{x}.
#' Like other print methods, this function is called to display a formatted
#' summary of the fitted \code{lgspline} model to the console.
#' This includes model dimensions, family information, dispersion estimate,
#' effective degrees of freedom, and a coefficient table for univariate inference
#' (if available) analogous to output from \code{\link[stats]{summary.glm}}.
#'
#' @seealso \code{\link[stats]{printCoefmat}}
#'
#' @export
print.summary.lgspline <- function(x, ...) {
  ## [Change 2026-02-11] use stats::printCoefmat() for coefficient table display
  # per reviewer recommendation, for standard p-value formatting.
  # printCoefmat's has.Pvalue = TRUE convention.
  cat("Lagrangian Multiplier Smoothing Spline Model Summary\n")
  cat("====================================================\n")
  cat("Model Family:", x$model_family[[1]], "\n")
  cat("Model Family (Link Function):",
      paste0(paste0(x$model_family)[1],
             " (",
             paste0(x$model_family)[2],
             ")", collapse = ""),
      "\n")
  cat("Observations:", x$observations, "\n")
  cat("Predictors:", x$predictors, "\n")
  cat("Partitions:", x$knots + 1, "\n")
  cat("Basis Functions per Partition:", x$basis_functions, "\n")
  if(length(unlist(x$coefficients)) > 1 && ncol(x$coefficients) > 1){
    cat("----------------------------------------------------\n")
    cat("Univariate Inference: \n")
    ## Identify which columns hold p-values for printCoefmat
    pval_col <- grep("^Pr\\(", colnames(x$coefficients))
    if(length(pval_col) > 0){
      ## printCoefmat expects the p-value column to be last;
      #  reorder so all non-p-value columns come first
      other_cols <- setdiff(seq_len(ncol(x$coefficients)), pval_col)
      print_mat <- x$coefficients[, c(other_cols, pval_col), drop = FALSE]
      pval_col_new <- ncol(print_mat)
      ## Re-detect test statistic column in reordered matrix
      tst_col <- grep("value$", colnames(print_mat))
      tst_col <- tst_col[!tst_col %in% pval_col_new]
      stats::printCoefmat(print_mat,
                          cs.ind = 1:2,
                          tst.ind = if(length(tst_col) > 0) tst_col[1] else NULL,
                          P.values = TRUE,
                          has.Pvalue = TRUE,
                          signif.stars = TRUE)
    } else {
      ## Fallback: no p-value column found, print as-is
      print(x$coefficients)
    }
    cat('\n')
    cat("Dispersion:", x$sigmasq_tilde, "\n")
    if(!is.null(x$N) && !is.null(x$trace_XUGX)){
      cat("Effective degrees of freedom:", x$N - x$trace_XUGX, "\n")
    } else {
      cat("Effective degrees of freedom: Not Available\n")
    }
    cat("Critical value for confidence intervals: ", x$cv, "\n")
    cat("----------------------------------------------------\n")
  }
  invisible(x)
}


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
#' @param ... Additional arguments passed to internal optimization routines.
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
#' @return A list containing the following components:
#' \describe{
#'   \item{t}{Numeric vector of input values at the extremum.}
#'   \item{y}{Numeric value of the objective function at the extremum.}
#' }
#'
#' @examples
#'
#' ## Basic usage with simulated data
#' set.seed(1234)
#' t <- runif(1000, -10, 10)
#' y <- 2*sin(t) + -0.06*t^2 + rnorm(length(t))
#' model_fit <- lgspline(t, y)
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
#'
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
                                                     sigma, ...){ # Added ...
                            N <- length(mu)
                            rnorm(
                              N, mu, sigma
                            )},
                          sigmasq_predict = object$sigmasq_tilde,
                          custom_objective_function = NULL,
                          custom_objective_derivative = NULL,
                          ...) {
  ## Delegate to internal method
  internal_find_extremum_func <- object$find_extremum
  if (!is.null(internal_find_extremum_func) && is.function(internal_find_extremum_func)) {
    return(internal_find_extremum_func(
      vars = vars,
      quick_heuristic = quick_heuristic,
      initial = initial,
      B_predict = B_predict,
      minimize = minimize,
      stochastic = stochastic,
      stochastic_draw = stochastic_draw,
      sigmasq_predict = sigmasq_predict,
      custom_objective_function = custom_objective_function,
      custom_objective_derivative = custom_objective_derivative,
      ...
    ))
  } else {
    stop("Internal find_extremum method not found or not a function.")
  }
}

#' Generate Posterior Samples from a Fitted Lagrangian Multiplier Smoothing
#' Spline with Optional Correlation Parameter Uncertainty
#'
#' Draws samples from the posterior distribution of model coefficients and
#' optionally generates posterior predictive samples. Uses a Laplace
#' approximation centred at the MAP estimate for non-Gaussian responses.
#'
#' When \code{draw_correlation = TRUE} and the fitted model contains an
#' estimated correlation structure, each draw first samples the correlation
#' parameters from their approximate normal posterior (centred at the BFGS
#' point estimate with covariance given by the inverse Hessian), re-estimates
#' the coefficients with the drawn correlation structure held fixed
#' (reusing all pre-computed design matrices, constraints, and penalties),
#' and then draws coefficients and (optionally) posterior predictive
#' realisations from the re-estimated quantities. This propagates
#' uncertainty in the correlation parameters through to all downstream
#' quantities without requiring access to the original raw predictor
#' matrix.
#'
#' @param object A fitted \code{lgspline} model object.
#' @param new_sigmasq_tilde Numeric scalar; dispersion parameter
#'        \eqn{\tilde{\sigma}^{2}} used as the point estimate when
#'        \code{draw_dispersion = FALSE}, and as the sufficient statistic
#'        when sampling \eqn{\sigma^{2}} from its inverse-gamma posterior.
#'        Default: \code{object$sigmasq_tilde}.
#' @param new_predictors Matrix; predictor matrix for posterior predictive
#'        sampling. Must be coercible to a numeric matrix with the same
#'        column structure as the original predictor input to
#'        \code{\link{lgspline}}. Default: in-sample predictors.
#' @param theta_1 Numeric scalar; shape increment for the inverse-gamma
#'        prior on \eqn{\sigma^{2}}. Setting \code{theta_1 = 0} (default)
#'        with \code{theta_2 = 0} implies a (improper) uniform prior.
#'        See Details.
#' @param theta_2 Numeric scalar; rate increment for the inverse-gamma
#'        prior on \eqn{\sigma^{2}}. See Details.
#' @param posterior_predictive_draw Function; sampler for posterior
#'        predictive realisations. Must accept arguments
#'        \code{N} (integer), \code{mean} (numeric vector),
#'        \code{sqrt_dispersion} (numeric scalar), and \code{...}.
#'        Defaults to \code{rnorm}.
#' @param draw_dispersion Logical; if \code{TRUE} (default), \eqn{\sigma^{2}}
#'        is drawn from its inverse-gamma posterior before sampling
#'        coefficients. If \code{FALSE}, \code{new_sigmasq_tilde} is used
#'        as a fixed point estimate throughout.
#' @param include_posterior_predictive Logical; if \code{TRUE}, posterior
#'        predictive realisations are generated at \code{new_predictors}
#'        for each draw. Default \code{FALSE}.
#' @param num_draws Positive integer; number of posterior draws. Default 1.
#' @param enforce_constraints Logical; if \code{TRUE} and inequality
#'        constraints were active during MAP estimation, an accept/reject
#'        loop is used to ensure each draw satisfies those constraints.
#'        \strong{Warning:} acceptance probability can be extremely low
#'        in high-dimensional or tightly constrained settings, causing
#'        the sampler to fall back to the MAP estimate after
#'        \code{max_rejection_draws} attempts. Default \code{FALSE}.
#'        See Details.
#' @param max_rejection_draws Positive integer; maximum number of
#'        accept/reject attempts per draw when
#'        \code{enforce_constraints = TRUE}. If exhausted, the MAP
#'        estimate is returned for that draw with a warning.
#'        Default \code{50L}.
#' @param draw_correlation Logical; if \code{TRUE} and the fitted model
#'        contains an estimated correlation structure (i.e.,
#'        \code{VhalfInv_fxn} and \code{VhalfInv_params_estimates} are
#'        present), each posterior draw first samples the correlation
#'        parameters from their approximate normal posterior, re-estimates
#'        the coefficients with the drawn correlation held fixed, and
#'        then draws coefficients from the re-estimated model. This
#'        propagates uncertainty in the correlation parameters through
#'        to the coefficient posterior. Default \code{FALSE}.
#' @param correlation_param_mean Numeric vector; mean of the approximate
#'        normal posterior for the correlation parameters on the
#'        unbounded (working) scale. Default: point estimates from the
#'        fitted model (\code{object$VhalfInv_params_estimates}). When
#'        supplied together with \code{correlation_param_vcov}, allows
#'        the user to override the model's estimates or to supply
#'        estimates when the model was fit with a fixed (non-optimised)
#'        correlation structure.
#' @param correlation_param_vcov Matrix; variance-covariance matrix of
#'        the approximate normal posterior for the correlation parameters
#'        on the unbounded (working) scale. Default: inverse Hessian from
#'        BFGS (\code{object$VhalfInv_params_vcov}). Must be symmetric
#'        positive semi-definite.
#' @param correlation_VhalfInv_fxn Function; maps the correlation
#'        parameter vector to \eqn{\mathbf{V}^{-1/2}}. Default:
#'        \code{object$VhalfInv_fxn}. Required when
#'        \code{draw_correlation = TRUE}.
#' @param correlation_Vhalf_fxn Function or \code{NULL}; maps the
#'        correlation parameter vector to \eqn{\mathbf{V}^{1/2}}.
#'        Default: \code{object$Vhalf_fxn}. When \code{NULL}, the
#'        inverse of \code{VhalfInv_fxn(par)} is computed internally
#'        (expensive for large \eqn{N}).
#' @param correlation_param_vcov_scale \code{NULL}; when NULL, will default to
#'        scaling the \code{correlation_param_vcov} by 1/(N-trace(H)).
#' @param include_warnings Logical; whether to emit warnings for
#'        constraint violations, degenerate draws, etc. Default
#'        \code{TRUE}.
#' @param ... Additional arguments forwarded to the GLM weight function,
#'        dispersion function, and \code{posterior_predictive_draw}.
#'
#' @details
#' \strong{Dispersion posterior.}
#' When \code{draw_dispersion = TRUE}, \eqn{\sigma^{2}} is drawn from
#'
#' \deqn{
#'   \sigma^{2} \mid \mathbf{y}
#'   \sim \mathrm{InvGamma}(\alpha_1, \alpha_2),
#' }
#'
#' where
#'
#' \deqn{
#'   \alpha_1 = \theta_1
#'   + \tfrac{1}{2}
#'     \bigl(N - s \cdot \mathrm{tr}(\mathbf{H})\bigr),
#'   \qquad
#'   \alpha_2 = \theta_2
#'   + \tfrac{1}{2}
#'     \bigl(N - s \cdot \mathrm{tr}(\mathbf{H})\bigr)
#'     \tilde{\sigma}^{2},
#' }
#'
#' \eqn{\mathbf{H} = \mathbf{X}\mathbf{U}\mathbf{G}\mathbf{X}^{\top}}
#' is the hat matrix, \eqn{s = 1} when \code{unbias_dispersion = TRUE}
#' and \eqn{s = 0} otherwise, and
#' \eqn{\tilde{\sigma}^{2}} = \code{new_sigmasq_tilde}.
#' Setting \eqn{\theta_1 = \theta_2 = 0} recovers the improper uniform
#' prior on \eqn{\sigma^{2}}.
#'
#' \strong{Correlation parameter posterior.}
#' When \code{draw_correlation = TRUE}, the correlation parameter vector
#' \eqn{\boldsymbol{\rho}} (on the unbounded working scale) is drawn as
#'
#' \deqn{
#'   \boldsymbol{\rho}^{(m)}
#'   \sim \mathcal{N}\bigl(
#'     \hat{\boldsymbol{\rho}}_{\mathrm{REML}},\,
#'     \mathbf{H}^{-1}_{\mathrm{BFGS}}
#'   \bigr),
#' }
#'
#' where \eqn{\hat{\boldsymbol{\rho}}_{\mathrm{REML}}} is the REML
#' point estimate and \eqn{\mathbf{H}^{-1}_{\mathrm{BFGS}}} is the
#' inverse Hessian from the BFGS optimiser (stored in
#' \code{VhalfInv_params_vcov}). For each draw of
#' \eqn{\boldsymbol{\rho}^{(m)}}, the coefficients are re-estimated
#' with \eqn{\mathbf{V}^{-1/2}(\boldsymbol{\rho}^{(m)})} held fixed.
#'
#' A single coefficient draw is then obtained from the re-estimated quantities. This yields
#' a Monte Carlo sample from the marginal posterior
#' \eqn{p(\boldsymbol{\beta} \mid \mathbf{y})} that integrates over
#' correlation parameter uncertainty.
#'
#' Re-estimation reuses all pre-computed structures (design matrices
#' \eqn{\mathbf{X}_k}, constraint matrix \eqn{\mathbf{A}}, penalty
#' matrices \eqn{\boldsymbol{\Lambda}}, partition assignments) from the
#' original fit. Only the GLS Gram matrices, coefficient solve, and
#' post-fit inference quantities (\eqn{\mathbf{G}}, \eqn{\mathbf{U}},
#' trace, dispersion, variance-covariance) are recomputed for each
#' draw. This avoids needing access to the original raw predictor
#' matrix.
#'
#' The normal approximation is exact asymptotically but may be poor for
#' small samples or when the likelihood surface for
#' \eqn{\boldsymbol{\rho}} is highly non-Gaussian. Draws that produce
#' non-positive-definite correlation matrices are rejected and redrawn
#' (up to 50 attempts per draw).
#'
#' \strong{Inequality constraints and approximate posteriors.}
#' When quadratic programming constraints are active, the MAP estimate
#' \eqn{\hat{\boldsymbol{\beta}}_{\mathrm{MAP}}} lies on or inside the
#' feasible region by construction. However, the unconstrained Gaussian
#' perturbation \eqn{\mathbf{U}\mathbf{G}^{1/2}\boldsymbol{z}} has
#' support over the entire coefficient space, so draws may violate the
#' constraints. Setting \code{enforce_constraints = TRUE} enables
#' accept/reject filtering at the cost of potentially low acceptance
#' rates.
#'
#' \strong{Correlation.}
#' When \code{draw_correlation = FALSE} (the default) and a correlation
#' structure is present, correlation parameters are fixed at their
#' estimated values, which ignores uncertainty in their estimation.
#' Setting \code{draw_correlation = TRUE} propagates that uncertainty at
#' the cost of re-estimating the coefficients for each draw.
#' Computational cost per draw is dominated by forming the whitened Gram
#' matrix and the constrained projection; knot placement, polynomial expansion, and GCV
#' penalty tuning are skipped entirely.
#'
#' @return
#' When \code{num_draws = 1}, a named list with elements:
#' \describe{
#'   \item{post_draw_coefficients}{List of length \eqn{K+1}; each element
#'         is a named \eqn{p \times 1} coefficient vector for one partition,
#'         on the original (unstandardised) scale.}
#'   \item{post_draw_sigmasq}{Numeric scalar; the drawn (or fixed)
#'         dispersion parameter \eqn{\sigma^{2(m)}}.}
#'   \item{post_pred_draw}{Numeric vector of length \eqn{N_{\mathrm{new}}};
#'         posterior predictive realisations (only when
#'         \code{include_posterior_predictive = TRUE}).}
#'   \item{post_draw_correlation_params}{Numeric vector; the drawn
#'         correlation parameters on the working scale (only when
#'         \code{draw_correlation = TRUE}).}
#' }
#' When \code{num_draws > 1}, each element above becomes a list of
#' length \code{num_draws}, and \code{post_pred_draw} (if requested)
#' is an \eqn{N_{\mathrm{new}} \times \code{num_draws}} matrix.
#'
#' @examples
#' ## Generate correlated data
#' set.seed(42)
#' n_blocks <- 100
#' block_size <- 5
#' N <- n_blocks * block_size
#' rho_true <- 0.3
#'
#' t <- seq(-5, 5, length.out = N)
#' true_mean <- sin(t)
#'
#' errors <- Reduce("rbind",
#'   lapply(1:n_blocks, function(i) {
#'     sigma <- diag(block_size) + rho_true *
#'       (matrix(1, block_size, block_size) - diag(block_size))
#'     matsqrt(sigma) %*% rnorm(block_size)
#'   })
#' )
#'
#' y <- true_mean + errors * 0.5
#'
#' ## Fit model with correlation structure
#' model_fit <- lgspline(t, y,
#'   K = 3,
#'   correlation_id = rep(1:n_blocks, each = block_size),
#'   correlation_structure = "exchangeable",
#'   include_warnings = FALSE
#' )
#'
#' ## Draw from posterior with correlation uncertainty propagated
#' post <- generate_posterior(model_fit,
#'   draw_correlation = TRUE,
#'   num_draws = 50,
#'   include_warnings = FALSE
#' )
#'
#' ## Compare to draws without correlation uncertainty
#' post_fixed <- generate_posterior(model_fit,
#'   draw_correlation = FALSE,
#'   num_draws = 50
#' )
#'
#' ## Posterior draws of correlation parameter (on working scale)
#' corr_draws <- unlist(post$post_draw_correlation_params)
#' rho_draws <- exp(-exp(corr_draws))
#' print(summary(rho_draws))
#'
#' @seealso
#' \code{\link{lgspline}} for model fitting,
#' \code{\link{generate_posterior_correlation}} for the standalone
#' correlation-aware sampler,
#' \code{\link{wald_univariate}} for Wald-type inference
#'
#' @export
generate_posterior <- function(object,
                               new_sigmasq_tilde = object$sigmasq_tilde,
                               new_predictors = NULL,
                               theta_1 = 0,
                               theta_2 = 0,
                               posterior_predictive_draw =
                                 function(N, mean, sqrt_dispersion, ...) {
                                   rnorm(N, mean, sqrt_dispersion)
                                 },
                               draw_dispersion = TRUE,
                               include_posterior_predictive = FALSE,
                               num_draws = 1,
                               enforce_constraints = FALSE,
                               max_rejection_draws = 50L,
                               draw_correlation = FALSE,
                               correlation_param_mean = NULL,
                               correlation_param_vcov = NULL,
                               correlation_VhalfInv_fxn = NULL,
                               correlation_Vhalf_fxn = NULL,
                               correlation_param_vcov_scale = NULL,
                               include_warnings = TRUE,
                               ...) {

  ## Warn when inequality constraints are active
  has_qp <- !is.null(object$quadprog_list) &&
    !identical(object$quadprog_list, list(NA)) &&
    !is.null(object$quadprog_list$qp_Amat)

  if(has_qp && include_warnings){
    warning(
      "\n\t Inequality constraints were active during MAP estimation. ",
      "Posterior draws are from the unconstrained Gaussian approximation ",
      "and may violate those constraints. ",
      "Set enforce_constraints = TRUE to use accept/reject sampling, ",
      "but note this can be extremely slow or degenerate when many ",
      "constraints are active or the feasible region is small.\n"
    )
  }

  ## Route to correlation-aware sampler when requested
  if(draw_correlation){

    ## Scale variance covariance matrix of correlation parameters by
    if(is.null(correlation_param_vcov_scale)){
      edf <- ifelse(is.null(model_fit$trace_XUGX) ||
                      is.na(model_fit$trace_XUGX),
                    0,
                    model_fit$trace_XUGX)
      scale_by <- model_fit$N - edf
    } else{
      scale_by <- correlation_param_vcov_scale
    }


    return(
      generate_posterior_correlation(
        object                     = object,
        new_sigmasq_tilde          = new_sigmasq_tilde,
        new_predictors             = new_predictors,
        theta_1                    = theta_1,
        theta_2                    = theta_2,
        posterior_predictive_draw  = posterior_predictive_draw,
        draw_dispersion            = draw_dispersion,
        include_posterior_predictive = include_posterior_predictive,
        num_draws                  = num_draws,
        enforce_constraints        = enforce_constraints,
        max_rejection_draws        = max_rejection_draws,
        correlation_param_mean     = correlation_param_mean,
        correlation_param_vcov_sc  = correlation_param_vcov,
        correlation_VhalfInv_fxn   = correlation_VhalfInv_fxn,
        correlation_Vhalf_fxn      = correlation_Vhalf_fxn,
        include_warnings           = include_warnings,
        ...
      )
    )
  }

  ## Standard path: delegate to the internal closure stored in the object
  internal_genpost_func <- object$generate_posterior
  if(!is.null(internal_genpost_func) && is.function(internal_genpost_func)){
    if(is.null(new_predictors)){
      return(internal_genpost_func(
        new_sigmasq_tilde          = new_sigmasq_tilde,
        theta_1                    = theta_1,
        theta_2                    = theta_2,
        posterior_predictive_draw  = posterior_predictive_draw,
        draw_dispersion            = draw_dispersion,
        include_posterior_predictive = include_posterior_predictive,
        num_draws                  = num_draws,
        enforce_constraints        = enforce_constraints,
        max_rejection_draws        = max_rejection_draws,
        ...
      ))
    } else {
      return(internal_genpost_func(
        new_sigmasq_tilde          = new_sigmasq_tilde,
        new_predictors             = new_predictors,
        theta_1                    = theta_1,
        theta_2                    = theta_2,
        posterior_predictive_draw  = posterior_predictive_draw,
        draw_dispersion            = draw_dispersion,
        include_posterior_predictive = include_posterior_predictive,
        num_draws                  = num_draws,
        enforce_constraints        = enforce_constraints,
        max_rejection_draws        = max_rejection_draws,
        ...
      ))
    }
  } else {
    stop("Internal generate_posterior method not found in lgspline object.")
  }
}

#' Generate Posterior Samples Propagating Correlation Parameter Uncertainty
#'
#' Performs a full Bayesian routine that accounts for uncertainty in
#' correlation structure parameters. For each Monte Carlo draw, the
#' correlation parameter vector is sampled from its approximate normal
#' posterior, the coefficient estimates are recomputed with the drawn
#' correlation held fixed (reusing all pre-computed design matrices,
#' constraints, and penalties from the original fit), and then a single
#' coefficient (and optionally posterior predictive) draw is obtained.
#'
#' This function is called internally by \code{\link{generate_posterior}}
#' when \code{draw_correlation = TRUE}, but can also be used directly for
#' finer control over the correlation sampling step.
#'
#' @param object A fitted \code{lgspline} model object that was fit with a
#'        correlation structure (i.e., \code{VhalfInv_fxn} and
#'        \code{VhalfInv_params_estimates} are present, or user supplies
#'        them via the override arguments below).
#' @param new_sigmasq_tilde Numeric scalar; dispersion parameter
#'        \eqn{\tilde{\sigma}^{2}} used for the coefficient posterior.
#'        Default: \code{object$sigmasq_tilde}. Note that for each
#'        re-estimation a new dispersion is obtained; this argument
#'        serves as the starting value / override when
#'        \code{draw_dispersion = FALSE}.
#' @param new_predictors Matrix or \code{NULL}; predictor matrix for
#'        posterior predictive sampling. Default: in-sample predictors.
#' @param theta_1 Numeric scalar; shape increment for the inverse-gamma
#'        prior on \eqn{\sigma^{2}}. Default 0.
#' @param theta_2 Numeric scalar; rate increment for the inverse-gamma
#'        prior on \eqn{\sigma^{2}}. Default 0.
#' @param posterior_predictive_draw Function; sampler for posterior
#'        predictive realisations. Default: \code{rnorm}.
#' @param draw_dispersion Logical; if \code{TRUE} (default),
#'        \eqn{\sigma^{2}} is drawn from its inverse-gamma posterior
#'        within each re-estimated model.
#' @param include_posterior_predictive Logical; if \code{TRUE}, posterior
#'        predictive realisations are generated for each draw.
#'        Default \code{FALSE}.
#' @param num_draws Positive integer; number of posterior draws. Each
#'        draw involves one correlation parameter sample and one
#'        coefficient re-estimation. Default 1.
#' @param enforce_constraints Logical; passed to the coefficient
#'        posterior sampler within each re-estimated model.
#'        Default \code{FALSE}.
#' @param max_rejection_draws Positive integer; passed to the coefficient
#'        posterior sampler. Default \code{50L}.
#' @param correlation_param_mean Numeric vector or \code{NULL}; mean of
#'        the approximate normal posterior for the correlation parameters
#'        on the working (unbounded) scale. Default: point estimates from
#'        the fitted model (\code{object$VhalfInv_params_estimates}).
#'        When supplied, overrides the model's stored estimates. This
#'        allows the user to provide external estimates or to enable
#'        correlation draws for a model that was fit with a fixed (not
#'        optimised) correlation structure.
#' @param correlation_param_vcov_sc Matrix or \code{NULL}; variance-covariance
#'        matrix of the approximate normal posterior on the working scale.
#'        Default: \code{object$VhalfInv_params_vcov}. Must be symmetric
#'        positive semi-definite. When supplied, overrides the model's
#'        stored inverse Hessian pre-scaled. No more scaling is peformed here.
#' @param correlation_VhalfInv_fxn Function or \code{NULL}; maps the
#'        correlation parameter vector to \eqn{\mathbf{V}^{-1/2}}.
#'        Default: \code{object$VhalfInv_fxn}. Required for
#'        constructing the correlation matrix from drawn parameters.
#' @param correlation_Vhalf_fxn Function or \code{NULL}; maps the
#'        correlation parameter vector to \eqn{\mathbf{V}^{1/2}}.
#'        Default: \code{object$Vhalf_fxn}. When \code{NULL}, the
#'        inverse of \code{VhalfInv_fxn(par)} is computed internally
#'        (expensive for large \eqn{N}).
#' @param correlation_param_vcov_sc \code{NULL}; when NULL, will default to
#'        scaling the \code{correlation_param_vcov} by 1/(N-trace(H)).
#' @param include_warnings Logical; whether to emit warnings.
#'        Default \code{TRUE}.
#' @param ... Additional arguments forwarded to downstream functions
#'        (GLM weight function, dispersion function, posterior
#'        predictive draw, etc.).
#'
#' @details
#' The algorithm proceeds as follows for each of the \code{num_draws}
#' iterations:
#'
#' \enumerate{
#'   \item \strong{Draw correlation parameters.}
#'     \eqn{\boldsymbol{\rho}^{(m)} \sim
#'     \mathcal{N}(\hat{\boldsymbol{\rho}},\,
#'     \mathbf{H}^{-1}_{\mathrm{BFGS}})} on the unbounded working
#'     scale. If the draw produces a non-positive-definite correlation
#'     matrix (i.e., \code{VhalfInv_fxn} fails), the draw is rejected
#'     and redrawn up to 50 times. If all attempts fail, the point
#'     estimate is used with a warning.
#'
#'   \item \strong{Re-estimate coefficients.}
#'     Using the already-expanded design matrices
#'     \eqn{\mathbf{X}_k}, constraint matrix \eqn{\mathbf{A}}, and
#'     tuned penalty matrices \eqn{\boldsymbol{\Lambda}} from the
#'     original fit, recompute the whitened Gram matrices with the
#'     drawn \eqn{\mathbf{V}^{-1/2}(\boldsymbol{\rho}^{(m)})}, solve
#'     for new MAP coefficients via constrained GLS, and update the
#'     post-fit inference quantities (\eqn{\mathbf{G}},
#'     \eqn{\mathbf{U}}, trace, dispersion, variance-covariance).
#'     This avoids repeating knot placement, partitioning, polynomial
#'     expansion, or penalty tuning. Specifically, the whitened
#'     penalised Gram matrix is formed as
#'
#'     \deqn{
#'       \mathbf{X}^{\top}\mathbf{V}^{-1}(\boldsymbol{\rho}^{(m)})
#'       \mathbf{X} + \boldsymbol{\Lambda},
#'     }
#'
#'     the corrected covariance is
#'     \eqn{\mathbf{G}_{\mathrm{correct}} =
#'     (\mathbf{X}^{\top}\mathbf{V}^{-1}\mathbf{X} +
#'     \boldsymbol{\Lambda})^{-1}},
#'     and the constraint projection is
#'     \eqn{\mathbf{U} = \mathbf{I} - \mathbf{G}_{\mathrm{correct}}
#'     \mathbf{A}
#'     (\mathbf{A}^{\top}\mathbf{G}_{\mathrm{correct}}\mathbf{A})^{-1}
#'     \mathbf{A}^{\top}}.
#'
#'   \item \strong{Draw coefficients.}
#'     From the re-estimated model, draw a single set of coefficients
#'     (and optionally dispersion and posterior predictive values)
#'     using the standard \code{generate_posterior} machinery with
#'     \code{draw_correlation = FALSE}.
#' }
#'
#' \strong{Normal approximation quality.}
#' The normal approximation for \eqn{\boldsymbol{\rho}} is based on
#' the BFGS inverse Hessian at the REML optimum. This is sometimes
#' asymptotically valid but may be poor for small samples, near
#' boundary estimates (e.g., correlation near 0 or 1), or when the
#' REML surface is multimodal. The BFGS hessian is not gauranteed to converge to
#' true Hessian and thus, the true observed information. Users should inspect
#' \code{VhalfInv_params_vcov} and consider whether the
#' approximation is reasonable for their application.
#'
#' @return
#' When \code{num_draws = 1}, a named list with elements:
#' \describe{
#'   \item{post_draw_coefficients}{List of length \eqn{K+1}; each element
#'         is a named coefficient vector for one partition on the original
#'         scale.}
#'   \item{post_draw_sigmasq}{Numeric scalar; the drawn dispersion.}
#'   \item{post_pred_draw}{Numeric vector; posterior predictive
#'         realisations (only when
#'         \code{include_posterior_predictive = TRUE}).}
#'   \item{post_draw_correlation_params}{Numeric vector; the drawn
#'         correlation parameters on the working scale.}
#' }
#' When \code{num_draws > 1}:
#' \describe{
#'   \item{post_draw_coefficients}{List of \code{num_draws} lists, each
#'         containing \eqn{K+1} coefficient vectors.}
#'   \item{post_draw_sigmasq}{List of \code{num_draws} numeric scalars.}
#'   \item{post_pred_draw}{\eqn{N_{\mathrm{new}} \times \code{num_draws}}
#'         matrix (only when \code{include_posterior_predictive = TRUE}).}
#'   \item{post_draw_correlation_params}{List of \code{num_draws} numeric
#'         vectors.}
#' }
#'
#' @examples
#' ## See ?generate_posterior for a complete example with
#' ## draw_correlation = TRUE, which calls this function internally.
#'
#' @seealso
#' \code{\link{generate_posterior}} for the unified interface,
#' \code{\link{lgspline}} for model fitting,
#' \code{\link{lgspline.fit}} for the low-level fitting interface
#'
#' @export
generate_posterior_correlation <- function(
    object,
    new_sigmasq_tilde = object$sigmasq_tilde,
    new_predictors = NULL,
    theta_1 = 0,
    theta_2 = 0,
    posterior_predictive_draw =
      function(N, mean, sqrt_dispersion, ...) {
        rnorm(N, mean, sqrt_dispersion)
      },
    draw_dispersion = TRUE,
    include_posterior_predictive = FALSE,
    num_draws = 1,
    enforce_constraints = FALSE,
    max_rejection_draws = 50L,
    correlation_param_mean = NULL,
    correlation_param_vcov_sc = NULL,
    correlation_VhalfInv_fxn = NULL,
    correlation_Vhalf_fxn = NULL,
    include_warnings = TRUE,
    ...
) {

  ## ## 1. Resolve correlation parameter mean, vcov, and functions.
  if(is.null(correlation_param_mean)){
    correlation_param_mean <- object$VhalfInv_params_estimates
    if(is.null(correlation_param_mean)){
      stop(
        "\n\t Cannot draw correlation parameters: no point estimates ",
        "found (VhalfInv_params_estimates is NULL) and no ",
        "'correlation_param_mean' was supplied.\n"
      )
    }
  }
  correlation_param_mean <- c(correlation_param_mean)
  n_corr_par <- length(correlation_param_mean)

  if(is.null(correlation_param_vcov_sc)){
    correlation_param_vcov_sc <- object$VhalfInv_params_vcov /
                                 (model_fit$N - ifelse(
                                   is.null(model_fit$trace_XUGX),
                                   0,
                                   model_fit$trace_XUGX))
    if(is.null(correlation_param_vcov_sc)){
      stop(
        "\n\t Cannot draw correlation parameters: no inverse Hessian ",
        "found (VhalfInv_params_vcov is NULL) and no ",
        "'correlation_param_vcov_sc' was supplied.\n"
      )
    }
  }
  correlation_param_vcov_sc <- as.matrix(correlation_param_vcov_sc)
  if(nrow(correlation_param_vcov_sc) != n_corr_par ||
     ncol(correlation_param_vcov_sc) != n_corr_par){
    stop(
      "\n\t 'correlation_param_vcov_sc' must be a ",
      n_corr_par, " x ", n_corr_par, " matrix.\n"
    )
  }

  if(is.null(correlation_VhalfInv_fxn)){
    correlation_VhalfInv_fxn <- object$VhalfInv_fxn
    if(is.null(correlation_VhalfInv_fxn)){
      stop(
        "\n\t No VhalfInv_fxn found and no ",
        "'correlation_VhalfInv_fxn' was supplied.\n"
      )
    }
  }

  if(is.null(correlation_Vhalf_fxn)){
    correlation_Vhalf_fxn <- object$Vhalf_fxn
  }

  ## ## 2. Cholesky factor of correlation parameter vcov for sampling.
  vcov_chol <- tryCatch({
    chol(correlation_param_vcov_sc)
  }, error = function(e){
    if(include_warnings){
      warning(
        "\n\t correlation_param_vcov_sc is not positive definite; ",
        "adding a small ridge to the diagonal.\n"
      )
    }
    chol(correlation_param_vcov_sc +
           diag(sqrt(.Machine$double.eps), n_corr_par))
  })

  ## ## 3. Extract components from fitted object for re-estimation.
  ##    No raw predictor matrix required.
  K  <- object$K
  nc <- object$p
  N  <- object$N
  A  <- object$A
  nca <- if(!is.null(A)) ncol(A) else 0

  ## Standardized design matrices per partition
  X_std <- lapply(object$X, object$std_X)

  y_og <- object$y
  order_list <- object$order_list
  og_order <- object$og_order
  family <- object$family
  penalties <- object$penalties
  expansion_scales <- object$expansion_scales
  mean_y <- object$mean_y
  sd_y <- object$sd_y
  backtransform_coefficients <- object$backtransform_coefficients
  raw_expansion_names <- object$raw_expansion_names

  ## Response in partition order, standardized scale
  y_std <- lapply(order_list, function(inds){
    cbind((y_og[inds] - mean_y) / sd_y)
  })

  ## Observation weights
  wts_og <- object$weights
  if(is.list(wts_og)){
    obs_wts_vec <- unlist(wts_og)[og_order]
  } else {
    obs_wts_vec <- wts_og
  }

  ## Unique penalty per partition flag
  unique_pp <- (length(penalties$L_partition_list) == (K + 1))

  ## Constraint values (for equality constraints)
  has_constraint_values <- length(object$constraint_values) > 0

  ## ## 4. Helper: re-estimate coefficients given a new VhalfInv.
  #     Mirrors the post-tuning fitting + post-fit inference block
  #     of lgspline.fit for the VhalfInv path, operating entirely
  #     on stored, already-expanded components.
  .reestimate_with_VhalfInv <- function(VhalfInv_new){

    ## Full N x P block-diagonal design matrix (standardized scale)
    X_full <- collapse_block_diagonal(X_std)

    ## Reorder VhalfInv to partition-based row order
    po <- unlist(order_list)
    VhalfInv_po <- VhalfInv_new[po, po]

    ## Whitened design and response
    VhalfInvX <- VhalfInv_po %**% X_full
    VhalfInvy <- VhalfInv_po %**% cbind(unlist(y_std))

    ## Full penalty matrix (block diagonal)
    Lambda_full <- collapse_block_diagonal(
      lapply(1:(K + 1), function(k){
        if(unique_pp){
          penalties$Lambda + penalties$L_partition_list[[k]]
        } else {
          penalties$Lambda
        }
      })
    )

    ## Penalized GLS Gram: X^T V^{-1} X + Lambda
    gram_gls <- crossprod(VhalfInvX) + Lambda_full

    ## G_correct = (X^T V^{-1} X + Lambda)^{-1}
    G_correct <- invert(gram_gls)
    Ghalf_correct <- matinvsqrt(gram_gls)

    ## Unconstrained penalized GLS estimate
    Xy_gls <- crossprod(VhalfInvX, VhalfInvy)
    B_unc <- G_correct %**% Xy_gls

    ## Constraint projection U = I - G A (A^T G A)^{-1} A^T
    P_total <- nc * (K + 1)
    if(!is.null(A) && nca > 0){
      GA <- G_correct %**% A
      AGA_inv <- invert(crossprod(A, GA))
      U <- diag(P_total) - GA %**% tcrossprod(AGA_inv, A)
    } else {
      U <- diag(P_total)
    }

    ## Constrained MAP coefficients
    if(has_constraint_values && !is.null(A) && nca > 0){
      c_vec <- unlist(object$constraint_values)
      B_vec <- c(U %**% B_unc +
                   G_correct %**% A %**% AGA_inv %**% cbind(c_vec))
    } else {
      B_vec <- c(U %**% B_unc)
    }

    ## Partition into per-partition lists (standardized scale)
    B_raw <- lapply(1:(K + 1), function(k){
      cbind(B_vec[1:nc + (k - 1) * nc])
    })

    ## Back-transform to original scale
    B <- lapply(B_raw, function(b){
      b_scaled <- b * sd_y
      b_scaled[1] <- b_scaled[1] + mean_y
      b_out <- backtransform_coefficients(b_scaled)
      names(b_out) <- raw_expansion_names
      b_out
    })
    names(B) <- paste0("partition", 1:(K + 1))

    ## Fitted values on original response scale.
    #  object$X stores the unstandardized design matrices per partition.
    ytilde_parts <- lapply(1:(K + 1), function(k){
      if(nrow(object$X[[k]]) == 0) return(numeric(0))
      c(object$X[[k]] %**% B[[k]])
    })
    ytilde <- family$linkinv(unlist(ytilde_parts)[og_order])

    ## Trace: ||V^{-1/2} X U Ghalf_correct||_F^2
    UGhalf <- U %**% Ghalf_correct
    trace_XUGX <- sum((VhalfInvX %**% UGhalf)^2)

    ## Dispersion estimate
    if(paste0(family)[1] == "gaussian" &&
       paste0(family)[2] == "identity"){
      resid_w <- c(VhalfInv_new %**% cbind(y_og - ytilde))
      scale_by <- if(object$unbias_dispersion){
        N / (N - trace_XUGX)
      } else {
        1
      }
      sigmasq_tilde <- mean(obs_wts_vec * resid_w^2) * scale_by
    } else {
      sigmasq_tilde <- object$sigmasq_tilde
    }

    ## Variance-covariance matrix (backtransformed scale)
    d <- rep(c(1, 1 / expansion_scales), times = K + 1)
    varcovmat <- tcrossprod(UGhalf) * sigmasq_tilde
    varcovmat <- t(t(varcovmat * d) * d)

    return(list(
      B = B,
      B_raw = B_raw,
      ytilde = ytilde,
      sigmasq_tilde = sigmasq_tilde,
      trace_XUGX = trace_XUGX,
      G_correct = G_correct,
      Ghalf_correct = Ghalf_correct,
      U = U,
      varcovmat = varcovmat,
      VhalfInv = VhalfInv_new
    ))
  }

  ## ## 5. Single-draw function.
  .one_corr_draw <- function(){

    ## 5a. Draw correlation parameters from N(mean, vcov)
    max_corr_reject <- 50L
    n_corr_reject <- 0L
    phi_draw <- NULL
    VhalfInv_draw <- NULL

    repeat {
      z <- rnorm(n_corr_par)
      phi_candidate <- correlation_param_mean + c(crossprod(vcov_chol, z))

      ## Validate: VhalfInv_fxn must succeed and return conformable matrix
      tr <- try({
        VhalfInv_candidate <- correlation_VhalfInv_fxn(phi_candidate)
        stopifnot(
          is.matrix(VhalfInv_candidate),
          all(is.finite(VhalfInv_candidate)),
          nrow(VhalfInv_candidate) == N,
          ncol(VhalfInv_candidate) == N
        )
      }, silent = TRUE)

      if(!inherits(tr, "try-error")){
        phi_draw <- phi_candidate
        VhalfInv_draw <- VhalfInv_candidate
        break
      }

      n_corr_reject <- n_corr_reject + 1L
      if(n_corr_reject >= max_corr_reject){
        if(include_warnings){
          warning(
            "\n\t Failed to draw a valid correlation parameter after ",
            max_corr_reject, " attempts. Using point estimate.\n"
          )
        }
        phi_draw <- correlation_param_mean
        VhalfInv_draw <- correlation_VhalfInv_fxn(phi_draw)
        break
      }
    }

    ## 5b. Re-estimate coefficients with drawn VhalfInv
    reest <- try(
      .reestimate_with_VhalfInv(VhalfInv_draw),
      silent = TRUE
    )

    if(inherits(reest, "try-error")){
      if(include_warnings){
        warning(
          "\n\t Re-estimation with drawn correlation parameters failed. ",
          "Falling back to original model for this draw.\n"
        )
      }
      one_draw <- object$generate_posterior(
        new_sigmasq_tilde          = new_sigmasq_tilde,
        theta_1                    = theta_1,
        theta_2                    = theta_2,
        posterior_predictive_draw  = posterior_predictive_draw,
        draw_dispersion            = draw_dispersion,
        include_posterior_predictive = include_posterior_predictive,
        num_draws                  = 1,
        enforce_constraints        = enforce_constraints,
        max_rejection_draws        = max_rejection_draws,
        ...
      )
      one_draw$post_draw_correlation_params <- phi_draw
      return(one_draw)
    }

    ## 5c. Shallow copy of the object with re-estimated fields.
    #      The internal generate_posterior closure reads B, B_raw,
    #      ytilde, G, Ghalf, U, sigmasq_tilde, trace_XUGX,
    #      varcovmat, VhalfInv from the model_fit list. Overwriting
    #      these on a copy makes the closure use re-estimated values.
    tmp_obj <- object
    tmp_obj$B <- reest$B
    tmp_obj$B_raw <- reest$B_raw
    tmp_obj$ytilde <- reest$ytilde
    tmp_obj$sigmasq_tilde <- reest$sigmasq_tilde
    tmp_obj$trace_XUGX <- reest$trace_XUGX
    tmp_obj$U <- reest$U
    tmp_obj$varcovmat <- reest$varcovmat
    tmp_obj$VhalfInv <- VhalfInv_draw

    ## Set per-partition G and Ghalf from the full G_correct
    tmp_obj$Ghalf <- lapply(1:(K + 1), function(k){
      idx <- 1:nc + (k - 1) * nc
      reest$Ghalf_correct[idx, idx, drop = FALSE]
    })
    tmp_obj$G <- lapply(1:(K + 1), function(k){
      idx <- 1:nc + (k - 1) * nc
      reest$G_correct[idx, idx, drop = FALSE]
    })

    refit_sigmasq <- reest$sigmasq_tilde
    if(is.na(refit_sigmasq) || !is.finite(refit_sigmasq)){
      refit_sigmasq <- new_sigmasq_tilde
    }

    ## 5d. Draw coefficients from the re-estimated model
    if(!is.null(new_predictors)){
      one_draw <- tmp_obj$generate_posterior(
        new_sigmasq_tilde          = refit_sigmasq,
        new_predictors             = new_predictors,
        theta_1                    = theta_1,
        theta_2                    = theta_2,
        posterior_predictive_draw  = posterior_predictive_draw,
        draw_dispersion            = draw_dispersion,
        include_posterior_predictive = include_posterior_predictive,
        num_draws                  = 1,
        enforce_constraints        = enforce_constraints,
        max_rejection_draws        = max_rejection_draws,
        ...
      )
    } else {
      one_draw <- tmp_obj$generate_posterior(
        new_sigmasq_tilde          = refit_sigmasq,
        theta_1                    = theta_1,
        theta_2                    = theta_2,
        posterior_predictive_draw  = posterior_predictive_draw,
        draw_dispersion            = draw_dispersion,
        include_posterior_predictive = include_posterior_predictive,
        num_draws                  = 1,
        enforce_constraints        = enforce_constraints,
        max_rejection_draws        = max_rejection_draws,
        ...
      )
    }

    one_draw$post_draw_correlation_params <- phi_draw
    return(one_draw)
  }

  ## Execute draws and collect results.
  results <- lapply(seq_len(num_draws), function(m) .one_corr_draw())

  if(num_draws == 1){
    return(results[[1]])
  }

  post_draw_coefficients <- lapply(results, `[[`, "post_draw_coefficients")
  post_draw_sigmasq <- lapply(results, `[[`, "post_draw_sigmasq")
  post_draw_correlation_params <- lapply(
    results, `[[`, "post_draw_correlation_params"
  )

  if(include_posterior_predictive){
    post_pred_draw <- Reduce(
      "cbind",
      lapply(results, `[[`, "post_pred_draw")
    )
    return(list(
      post_draw_coefficients       = post_draw_coefficients,
      post_draw_sigmasq            = post_draw_sigmasq,
      post_pred_draw               = post_pred_draw,
      post_draw_correlation_params = post_draw_correlation_params
    ))
  }

  return(list(
    post_draw_coefficients       = post_draw_coefficients,
    post_draw_sigmasq            = post_draw_sigmasq,
    post_draw_correlation_params = post_draw_correlation_params
  ))
}


#' Plot Method for Lagrangian Multiplier Smoothing Spline Models
#'
#' Creates visualizations of fitted spline models, supporting both 1D line plots and 2D surface
#' plots with optional formula annotations and customizable aesthetics.
#' (Wrapper for internal plot method)
#'
#' @param x A fitted lgspline model object containing the model fit to be plotted
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
#' @param color_function Function; Returns colors for plotting by partition, must return K+1 vector of valid colors. Defaults to NULL, in which case \code{grDevices::rainbow(K+1)} is used for 1D and \code{grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Spectral"))(K+1)} used in multiple.
#' @param add Logical; If TRUE, adds to existing plot (1D only). Similar to add in
#'        \code{\link[graphics]{hist}}. Default FALSE
#' @param vars Numeric or character vector; Optional indices for selecting variables to plot. Can either be numeric (the column indices of "predictors" or "data") or character (the column names, if available from "predictors" or "data")
#' @param legend_order Numeric specifying the re-arranged default order of partitions in the legend.
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
#' @return Returns
#' \describe{
#'   \item{1D}{Invisibly returns NULL (base R plot is drawn to device).}
#'   \item{2D}{Plotly object showing interactive surface plot.}
#' }
#'
#' @examples
#'
#' ## Generate example data
#' set.seed(1234)
#' t_data <- runif(1000, -10, 10)
#' y_data <- 2*sin(t_data) + -0.06*t_data^2 + rnorm(length(t_data))
#'
#' ## Fit model with 10 partitions
#' model_fit <- lgspline(t_data, y_data, K = 9)
#'
#' ## Basic plot
#' plot(model_fit)
#'
#' ## Customized plot with formulas
#' plot(model_fit,
#'      show_formulas = TRUE,         # Show partition formulas
#'      custom_response_lab = 'Price',  # Custom axis labels
#'      custom_predictor_lab = 'Size',
#'      custom_title = 'Price vs Size', # Custom title
#'      digits = 2,                    # Round coefficients
#'      legend_pos = 'bottom',         # Move legend
#'      text_size_formula = 0.375,     # Adjust formula text size
#'      pch = 16,                      # Point style
#'      cex.main = 1.25)               # Title size
#'
#' @seealso
#' \code{\link{lgspline}} for model fitting,
#' \code{\link[graphics]{plot}} for additional 1D plot parameters,
#' \code{\link[plotly]{plot_ly}} for additional 2D plot parameters
#'
#' @export
plot.lgspline <- function(x,
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
                          legend_order = NULL, # [Change 2026-02-14] Include
                          ...) {
  # Use the model's internal plotting function
  internal_plot_func <- x$plot
  if (!is.null(internal_plot_func) && is.function(internal_plot_func)) {
    plot_result <- internal_plot_func(model_fit_in = x,
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
    if(inherits(plot_result, "plotly")) {
      return(plot_result)
    } else {
      invisible(NULL)
    }
  } else {
    stop("Internal plot method not found or not a function.")
  }
}


#' Predict Method for Fitted Lagrangian Multiplier Smoothing Spline
#'
#' Generates predictions, derivatives, and basis expansions from a fitted lgspline model.
#' Supports both in-sample and out-of-sample prediction with optional parallel processing.
#' (Wrapper for internal predict method)
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
#' @param B_predict List; Optional custom per-partition coefficient list for prediction,
#'        e.g. from generate_posterior(). Default NULL (uses object$B).
#' @param take_first_derivatives Logical; whether to compute first derivatives of the
#'        fitted function. Default FALSE
#' @param take_second_derivatives Logical; whether to compute second derivatives of the
#'        fitted function. Default FALSE
#' @param expansions_only Logical; whether to return only basis expansions without
#'        computing predictions. Default FALSE
#' @param new_predictors Matrix or data frame; overrides 'newdata' if provided.
#' @param ... Additional arguments passed to internal prediction methods.
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
#' If newdata and new_predictor are left NULL, default input used for model fitting
#' will be used. Priority will be awarded to new_predictor over newdata when
#' both are not NULL.
#'
#' To obtain fitted values, users may also call model_fit$predict() or
#' model_fit$ytilde for an lgspline object "model_fit".
#'
#' The parallel processing feature is experimental and should be used with caution.
#' When enabled, computations are split across chunks and processed in parallel,
#' which may improve performance for large datasets.
#'
#' @return Depending on the options selected, returns the following:
#' \describe{
#'   \item{predictions}{Numeric vector of predicted values (default case, or if derivatives requested).}
#'   \item{first_deriv}{Numeric vector of first derivatives (if take_first_derivatives = TRUE).}
#'   \item{second_deriv}{Numeric vector of second derivatives (if take_second_derivatives = TRUE).}
#'   \item{expansions}{List of basis expansions (if expansions_only = TRUE).}
#' }
#'
#' With derivatives included, output is in the form of a list with elements
#' "preds", "first_deriv", and "second_deriv" for the vector of predictions,
#' first derivatives, and second derivatives respectively.
#'
#' Important, make sure the input new_predictors/newdata matches the input
#' structure of the data used to fit the model - do not include additional
#' predictors or columns that weren't originally included.
#'
#' @examples
#'
#' ## Generate example data
#' set.seed(1234)
#' t <- runif(1000, -10, 10)
#' y <- 2*sin(t) + -0.06*t^2 + rnorm(length(t))
#'
#' ## Fit model
#' model_fit <- lgspline(t, y)
#'
#' ## Generate predictions for new data
#' newdata <- matrix(sort(rnorm(10000)), ncol = 1) # Ensure matrix format
#' preds <- predict(model_fit, newdata)
#'
#' ## Compute derivative
#' deriv1_res <- predict(model_fit, newdata,
#'                       take_first_derivatives = TRUE)
#' deriv2_res <- predict(model_fit, newdata,
#'                       take_second_derivatives = TRUE)
#'
#' ## Visualize results
#' oldpar <- par(no.readonly = TRUE) # Save current par settings
#' layout(matrix(c(1,1,2,2,3,3), byrow = TRUE, ncol = 2))
#'
#' ## Plot function
#' plot(newdata[,1], preds,
#'        main = 'Fitted Function',
#'        xlab = 't',
#'        ylab = "f(t)", type = 'l')
#'
#' ## Plot first derivative
#' plot(newdata[,1],
#'        deriv1_res$first_deriv,
#'        main = 'First Derivative',
#'        xlab = 't',
#'        ylab = "f'(t)", type = 'l')
#'
#' ## Plot second derivative
#' plot(newdata[,1],
#'        deriv2_res$second_deriv,
#'        main = 'Second Derivative',
#'        xlab = 't',
#'        ylab = "f''(t)", type = 'l')
#'
#' par(oldpar) # Reset to original par settings
#'
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
                             B_predict = NULL,
                             take_first_derivatives = FALSE,
                             take_second_derivatives = FALSE,
                             expansions_only = FALSE,
                             new_predictors = NULL,
                             ...) {

  internal_predict_func <- object$predict
  if(is.null(internal_predict_func) || !is.function(internal_predict_func)){
    stop("Internal predict method not found or not a function.")
  }

  ## Resolve coefficient list. B_predict = NULL means use object$B.
  #  Explicit non-NULL B_predict (e.g. from generate_posterior) is passed
  #  through directly. Using is.null() rather than !missing() avoids the
  #  bug where missing() returns TRUE even when B_predict is supplied via
  #  a named argument in certain paths.
  B_predict_val <- if(!is.null(B_predict)) B_predict else object$B

  ## Resolve predictor source. new_predictors takes priority over newdata
  if(!is.null(new_predictors)){
    predictors_val <- new_predictors
  } else if(!is.null(newdata)){
    predictors_val <- newdata
  } else {
    predictors_val <- NULL
  }

  ## Unwrap nested data.frame: data.frame(new_predictors = data.frame(...))
  #  produces a single-column df whose sole column is itself a df.
  #  Flatten it before any further processing.
  if(inherits(predictors_val, "data.frame") && ncol(predictors_val) == 1 &&
     inherits(predictors_val[[1]], "data.frame")){
    predictors_val <- predictors_val[[1]]
  }

  ## Leave data.frames with non-numeric columns as-is so predict_function
  #  can run its factor-encoding logic before matrix coercion
  if(!is.null(predictors_val)){
    if(inherits(predictors_val, "data.frame") &&
       any(!sapply(predictors_val, function(x) is.numeric(x) ||
                   is.integer(x)))){
      ## has non-numeric columns, leave as data.frame
    } else {
      predictors_val <- try(
        methods::as(cbind(predictors_val), 'matrix'),
        silent = TRUE
      )
      if(inherits(predictors_val, 'try-error')){
        stop('\n \t newdata / new_predictors cannot be coerced to a matrix. \n')
      }
    }
  }

  ## In-sample path omits new_predictors so the internal function
  #  uses its default (the closed-over training predictors). Out-of-sample
  #  path passes new_predictors explicitly by name to avoid positional
  #  mismatch with the internal function signature.
  if(is.null(predictors_val)){
    internal_predict_func(
      parallel              = parallel,
      cl                    = cl,
      chunk_size            = chunk_size,
      num_chunks            = num_chunks,
      rem_chunks            = rem_chunks,
      B_predict             = B_predict_val,
      take_first_derivatives  = take_first_derivatives,
      take_second_derivatives = take_second_derivatives,
      expansions_only       = expansions_only,
      ...
    )
  } else {
    internal_predict_func(
      new_predictors        = predictors_val,
      parallel              = parallel,
      cl                    = cl,
      chunk_size            = chunk_size,
      num_chunks            = num_chunks,
      rem_chunks            = rem_chunks,
      B_predict             = B_predict_val,
      take_first_derivatives  = take_first_derivatives,
      take_second_derivatives = take_second_derivatives,
      expansions_only       = expansions_only,
      ...
    )
  }
}


#' Extract model coefficients
#'
#' Extracts polynomial coefficients for each partition from a fitted lgspline model.
#'
#' @param object A fitted lgspline model object containing coefficient vectors.
#' @param ... Not used.
#'
#' @details
#' For each partition, coefficients represent a polynomial expansion of the predictor(s) by column index, for example:
#' \itemize{
#'   \item intercept: Constant term
#'   \item v: Linear term
#'   \item v_^2: Quadratic term
#'   \item v^3: Cubic term
#'   \item _v_x_w_: Interaction between v and w
#' }
#'
#' If column/variable names are present, indices will be replaced with column/variable names.
#'
#' Coefficients can be accessed either as separate vectors per partition or combined into
#' a single matrix using \code{Reduce('cbind', coef(model_fit))}.
#'
#' @return
#' A list where each element corresponds to a partition and contains a single-column matrix
#' of coefficient values for that partition. Row names indicate the term type. Returns NULL if
#' coefficients are not found in the object.
#' \describe{
#'   \item{partition1, partition2, ...}{Matrices containing coefficients for each partition.}
#' }
#'
#' @examples
#'
#' ## Simulate some data and fit using default settings
#' set.seed(1234)
#' t <- runif(1000, -10, 10)
#' y <- 2*sin(t) + -0.06*t^2 + rnorm(length(t))
#' model_fit <- lgspline(t, y)
#'
#' ## Extract coefficients
#' coefficients <- coef(model_fit)
#'
#' ## Print coefficients for first partition
#' print(coefficients[[1]])
#'
#' ## Compare coefficients across all partitions
#' print(Reduce('cbind', coefficients))
#'
#' @seealso \code{\link{lgspline}}
#' @export
coef.lgspline <- function(object, ...) {
  if(is.null(object$B)) {
    warning("Coefficient component 'B' not found in object.", call. = FALSE)
    return(NULL)
  }
  object$B
}


## =========================================================================
## wald_univariate: Univariate Wald Tests and Confidence Intervals
## =========================================================================
## [Change 2026-02-11] Refactored per reviewer feedback.
##   - Return value is now a classed "wald_lgspline" object (a list) with
##     $coefficients containing a properly formatted matrix.
##   - Dedicated print, summary, and plot methods are provided.
##   - Backward-compatible list accessors ($estimate, $std_error, etc.)
##     are included for convenience.
##   - print.wald_lgspline uses stats::printCoefmat() for standard
##     p-value formatting and significance stars.

#' Univariate Wald Tests and Confidence Intervals for Lagrangian Multiplier Smoothing Splines
#'
#' Performs coefficient-specific Wald tests and constructs confidence intervals for fitted
#' lgspline models. (Wrapper for internal wald_univariate method). For Gaussian family
#' with identity-link, a t-distribution replaces a normal distribution (and t-intervals, t-tests etc.)
#' over Wald when mentioned.
#'
#' @param object A fitted lgspline model object containing coefficient estimates and
#'        variance-covariance matrix (requires return_varcovmat = TRUE in fitting).
#' @param scale_vcovmat_by Numeric; Scaling factor for variance-covariance matrix.
#'        Adjusts standard errors and test statistics. Default 1.
#' @param cv Numeric; Critical value for confidence interval construction. If missing,
#'        defaults to value specified in lgspline() fit (`object$critical_value`) or
#'        `qnorm(0.975)` as a fallback. Common choices:
#'        \itemize{
#'          \item qnorm(0.975) for normal-based 95\% intervals
#'          \item qt(0.975, df) for t-based 95\% intervals, where df = N - trace(XUGX)
#'        }
#' @param ... Additional arguments passed to the internal `wald_univariate` method.
#'
#' @details
#' For each coefficient, provides:
#' \itemize{
#'   \item Point estimates
#'   \item Standard errors from the model's variance-covariance matrix
#'   \item Two-sided test statistics and p-values
#'   \item Confidence intervals using specified critical values
#' }
#'
#' @return An object of class \code{"wald_lgspline"}, which is a list with components:
#' \describe{
#'   \item{coefficients}{Matrix with one row per coefficient and columns: \code{Estimate},
#'     \code{Std. Error}, test statistic (\code{t value} or \code{z value}),
#'     p-value (\code{Pr(>|t|)} or \code{Pr(>|z|)}), \code{CI LB}, \code{CI UB}.}
#'   \item{critical_value}{The critical value used for confidence intervals.}
#'   \item{family}{The GLM family from the fitted model.}
#'   \item{N}{Number of observations.}
#'   \item{trace_XUGX}{Effective degrees of freedom trace term.}
#'   \item{statistic_name}{Character: \code{"t value"} or \code{"z value"}.}
#'   \item{p_value_name}{Character: \code{"Pr(>|t|)"} or \code{"Pr(>|z|)"}.}
#' }
#'
#' Print, summary, and plot methods are provided for this class; see
#' \code{\link{print.wald_lgspline}}, \code{\link{summary.wald_lgspline}},
#' \code{\link{plot.wald_lgspline}}.
#'
#' @examples
#'
#' ## Simulate some data and fit using default settings
#' set.seed(1234)
#' t <- runif(1000, -10, 10)
#' y <- 2*sin(t) + -0.06*t^2 + rnorm(length(t))
#' # Ensure varcovmat is returned for Wald tests
#' model_fit <- lgspline(t, y, return_varcovmat = TRUE)
#'
#' ## Use default critical value (likely qnorm(0.975) if not set in fit)
#' wald_default <- wald_univariate(model_fit)
#' print(wald_default)
#'
#' ## Specify t-distribution critical value
#' eff_df <- NA
#' if(!is.null(model_fit$N) && !is.null(model_fit$trace_XUGX)) {
#'    eff_df <- model_fit$N - model_fit$trace_XUGX
#' }
#' if (!is.na(eff_df) && eff_df > 0) {
#'   wald_t <- wald_univariate(
#'     model_fit,
#'     cv = stats::qt(0.975, eff_df)
#'   )
#'   print(wald_t)
#' } else {
#'   warning("Effective degrees of freedom invalid.")
#' }
#'
#' ## Extract the coefficient matrix directly
#' coef_table <- wald_default$coefficients
#' print(coef_table)
#'
#' ## Plot coefficient estimates with confidence intervals
#' plot(wald_default)
#'
#'
#' @seealso
#' \code{\link{lgspline}}, \code{\link{confint.lgspline}},
#' \code{\link{print.wald_lgspline}}, \code{\link{summary.wald_lgspline}},
#' \code{\link{plot.wald_lgspline}}
#' @export
wald_univariate <- function(object, scale_vcovmat_by = 1, cv, ...) {
  if (is.null(object$varcovmat)) {
    stop("Wald tests require return_varcovmat = TRUE during model fitting")
  }

  if (missing(cv)) {
    if (!is.null(object$critical_value)) {
      cv <- object$critical_value
    } else {
      cv <- stats::qnorm(0.975)
      warning("Critical value 'cv' not provided, defaulting to qnorm(0.975).",
              call. = FALSE)
    }
  }

  ## Determine test statistic and p-value column labels
  if(object$family$family == 'gaussian' && object$family$link == 'identity'){
    stat_name <- 't value'
    p_name <- 'Pr(>|t|)'
  } else {
    stat_name <- 'z value'
    p_name <- 'Pr(>|z|)'
  }

  internal_wald_func <- object$wald_univariate
  if (is.null(internal_wald_func) || !is.function(internal_wald_func)) {
    stop("Internal wald_univariate method not found or not a function.")
  }

  res <- internal_wald_func(scale_vcovmat_by = scale_vcovmat_by, cv = cv, ...)

  ## Normalize result to a properly labelled matrix.
  # The internal method may return a list of vectors or a matrix/data.frame.
  coef_mat <- NULL
  if(is.list(res) && !is.data.frame(res) && !is.matrix(res)){
    ## Internal method returned a list -- cbind into matrix
    coef_mat <- tryCatch({
      mat <- Reduce('cbind', res)
      ## Expect 6 columns: estimate, se, stat, ci_lb, ci_ub, pval
      # (order from internal method based on summary.lgspline code)
      if(ncol(mat) >= 6){
        colnames(mat) <- c('Estimate', 'Std. Error', stat_name,
                           'CI LB', 'CI UB', p_name)
        mat <- mat[, c('Estimate', 'Std. Error', stat_name,
                       p_name, 'CI LB', 'CI UB'), drop = FALSE]
      }
      mat
    }, error = function(e) NULL)
    if(is.null(coef_mat)){
      ## Fallback: just use the raw cbind
      coef_mat <- tryCatch(Reduce('cbind', res), error = function(e) NULL)
    }
  } else if(is.matrix(res) || is.data.frame(res)){
    coef_mat <- as.matrix(res)
  }

  ## Final fallback
  if(is.null(coef_mat)){
    warning("Could not normalize wald_univariate output to matrix.",
            call. = FALSE)
    coef_mat <- cbind(Estimate = unlist(object$B))
  }

  ## Assign row names from coefficient names if not already present
  if(is.null(rownames(coef_mat)) && !is.null(object$B) && is.list(object$B) &&
     length(unlist(object$B)) == nrow(coef_mat)){
    rn <- tryCatch({
      unlist(lapply(seq_along(object$B), function(k) {
        part_names <- names(object$B[[k]])
        if(is.null(part_names))
          part_names <- paste0("Term", seq_len(length(object$B[[k]])))
        paste0("partition", k, "_", part_names)
      }))
    }, error = function(e) NULL)
    if(!is.null(rn) && length(rn) == nrow(coef_mat)){
      rownames(coef_mat) <- rn
    }
  }

  ## Build classed output object
  out <- list(
    coefficients = coef_mat,
    critical_value = cv,
    family = object$family,
    N = object$N,
    trace_XUGX = object$trace_XUGX,
    statistic_name = stat_name,
    p_value_name = p_name
  )
  class(out) <- "wald_lgspline"
  return(out)
}


#' Print Method for wald_lgspline Objects
#'
#' @description
#' Prints the coefficient table from \code{\link{wald_univariate}} using
#' \code{\link[stats]{printCoefmat}} for standard p-value formatting with
#' significance stars.
#'
#' @param x An object of class \code{"wald_lgspline"}, as returned by
#'   \code{\link{wald_univariate}}.
#' @param ... Additional arguments passed to \code{\link[stats]{printCoefmat}}.
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link{wald_univariate}}, \code{\link[stats]{printCoefmat}}
#' @method print wald_lgspline
#' @export
print.wald_lgspline <- function(x, ...) {
  ## [Change 2026-02-11] added for proper classed output
  # with printCoefmat-based display.
  # printCoefmat's has.Pvalue = TRUE convention.
  cat("Univariate Wald Tests\n")
  cat("---------------------\n")
  mat <- x$coefficients
  pval_col <- grep("^Pr\\(", colnames(mat))
  if(length(pval_col) > 0){
    ## printCoefmat expects the p-value column to be last;
    ## reorder so all non-p-value columns come first
    other_cols <- setdiff(seq_len(ncol(mat)), pval_col)
    print_mat <- mat[, c(other_cols, pval_col), drop = FALSE]
    pval_col_new <- ncol(print_mat)
    ## Re-detect test statistic column in reordered matrix
    tst_col <- grep("value$", colnames(print_mat))
    tst_col <- tst_col[!tst_col %in% pval_col_new]
    stats::printCoefmat(print_mat,
                        cs.ind = 1:2,
                        tst.ind = if(length(tst_col) > 0) tst_col[1] else NULL,
                        P.values = TRUE,
                        has.Pvalue = TRUE,
                        signif.stars = TRUE,
                        ...)
  } else {
    print(mat, ...)
  }
  cat("\nCritical value:", x$critical_value, "\n")
  invisible(x)
}


#' Summary Method for wald_lgspline Objects
#'
#' @description
#' Provides a more detailed summary of the Wald test results, including
#' effective degrees of freedom and dispersion information.
#'
#' @param object An object of class \code{"wald_lgspline"}, as returned by
#'   \code{\link{wald_univariate}}.
#' @param ... Not used.
#'
#' @return Invisibly returns \code{object}.
#'
#' @seealso \code{\link{wald_univariate}}, \code{\link{print.wald_lgspline}}
#' @method summary wald_lgspline
#' @export
summary.wald_lgspline <- function(object, ...) {
  ## [Change 2026-02-11] added per reviewer request.
  cat("Univariate Wald Tests Summary\n")
  cat("=============================\n")
  cat("Family:", object$family$family, "\n")
  cat("Link:", object$family$link, "\n")
  cat("N:", object$N, "\n")
  if(!is.null(object$N) && !is.null(object$trace_XUGX)){
    cat("Effective df:", object$N - object$trace_XUGX, "\n")
  }
  cat("Test statistic:", object$statistic_name, "\n")
  cat("Critical value:", object$critical_value, "\n")
  cat("-----------------------------\n")
  print(object)
  invisible(object)
}


#' Plot Method for wald_lgspline Objects
#'
#' @description
#' Produces a forest-style plot of coefficient estimates with confidence
#' intervals from \code{\link{wald_univariate}} results.
#'
#' @param x An object of class \code{"wald_lgspline"}, as returned by
#'   \code{\link{wald_univariate}}.
#' @param parm An optional vector of row indices or coefficient names
#'   selecting which parameters to plot. By default, all are shown.
#' @param main Character; plot title. Default \code{"Coefficient Estimates and CIs"}.
#' @param xlab Character; x-axis label. Default \code{"Estimate"}.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}.
#'
#' @return Invisibly returns \code{NULL}. A base R plot is drawn to the
#'   current graphics device.
#'
#' @seealso \code{\link{wald_univariate}}, \code{\link{confint.lgspline}}
#' @method plot wald_lgspline
#' @export
plot.wald_lgspline <- function(x,
                               parm = NULL,
                               main = "Coefficient Estimates and CIs",
                               xlab = "Estimate",
                               ...) {
  ## [Change 2026-02-11] added plot method.
  mat <- x$coefficients
  if(!is.null(parm)){
    if(is.character(parm)){
      mat <- mat[rownames(mat) %in% parm, , drop = FALSE]
    } else {
      mat <- mat[parm, , drop = FALSE]
    }
  }

  ## Extract estimates and CI bounds
  est_col <- which(colnames(mat) == "Estimate")
  lb_col <- which(colnames(mat) == "CI LB")
  ub_col <- which(colnames(mat) == "CI UB")

  if(length(est_col) == 0 || length(lb_col) == 0 || length(ub_col) == 0){
    warning("Cannot produce forest plot: missing Estimate, CI LB, or CI UB columns.",
            call. = FALSE)
    print(mat)
    return(invisible(NULL))
  }

  ests <- mat[, est_col]
  lbs <- mat[, lb_col]
  ubs <- mat[, ub_col]
  n_coef <- length(ests)
  idx <- seq_len(n_coef)

  ## Coefficient labels
  labs <- rownames(mat)
  if(is.null(labs)) labs <- paste0("Coef ", idx)

  ## Draw plot
  xlims <- range(c(lbs, ubs), na.rm = TRUE)
  xlims <- xlims + diff(xlims) * c(-0.05, 0.05)

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(mar = c(4.5, max(nchar(labs)) * 0.5 + 2, 3, 1))

  graphics::plot(ests, idx, xlim = xlims, yaxt = "n",
                 ylab = "", xlab = xlab, main = main,
                 pch = 19, ...)
  graphics::axis(2, at = idx, labels = labs, las = 1, cex.axis = 0.7)
  graphics::segments(lbs, idx, ubs, idx, lwd = 2)
  graphics::abline(v = 0, lty = 2, col = "grey50")
  invisible(NULL)
}


## =========================================================================
## confint.lgspline: Confidence Intervals for lgspline Model Coefficients
## =========================================================================
## [Change 2026-02-11] New method requesting a confint
#  method for the lgspline class.

#' Confidence Intervals for lgspline Model Coefficients
#'
#' @description
#' Computes confidence intervals for all or selected coefficients of a fitted
#' \code{lgspline} model, using the Wald-based approach from the estimated
#' variance-covariance matrix. If correlation parameters are present and their
#' variance-covariance matrix is available, confidence intervals for those
#' parameters (on the transformed scale) are also returned.
#'
#' @param object A fitted \code{lgspline} model object. Must have been fit
#'   with \code{return_varcovmat = TRUE}.
#' @param parm An optional specification of which regression parameters to give
#'   confidence intervals for, either a vector of numbers (indices) or a
#'   vector of names. If missing, all regression parameters are considered.
#' @param level The confidence level required. Default 0.95.
#' @param ... Additional arguments passed to \code{\link{wald_univariate}}.
#'
#' @details
#' For Gaussian family with identity link, t-distribution quantiles are used
#' with effective degrees of freedom \eqn{N - \mathrm{trace}(\mathbf{X}\mathbf{U}\mathbf{G}\mathbf{X}^{\top})}.
#' For all other families, normal quantiles are used.
#'
#' Regression parameters are handled via \code{\link{wald_univariate}}.
#' Correlation parameters (if present) are handled analogously using a Wald
#' interval on the working (un-transformed, unbounded) scale, if the inverse
#' approximate-Hessian of the BFGS algorithm was returned.
#'
#' @return A matrix with columns giving lower and upper confidence limits
#'   for each parameter. The column names encode the confidence level,
#'   e.g., \code{2.5 \%} and \code{97.5 \%} for 95\% intervals.
#'
#' @method confint lgspline
#' @export
confint.lgspline <- function(object, parm, level = 0.95, ...) {

  if(is.null(object$varcovmat)){
    stop("confint requires return_varcovmat = TRUE during model fitting")
  }

  alpha <- 1 - level

  ## Choose critical value based on family
  if(object$family$family == "gaussian" &&
     object$family$link == "identity"){

    eff_df <- object$N - object$trace_XUGX

    if(!is.null(eff_df) && !is.na(eff_df) &&
       is.finite(eff_df) && eff_df > 0){
      cv <- stats::qt(1 - alpha / 2, df = eff_df)
    } else {
      cv <- stats::qnorm(1 - alpha / 2)
      warning("Effective df non-positive; using normal quantiles.",
              call. = FALSE)
    }

  } else {
    cv <- stats::qnorm(1 - alpha / 2)
  }

  ## Regression parameter intervals
  wald_res <- wald_univariate(object, cv = cv, ...)

  coef_mat <- wald_res$coefficients
  lb_col <- which(colnames(coef_mat) == "CI LB")
  ub_col <- which(colnames(coef_mat) == "CI UB")

  if(length(lb_col) == 0 || length(ub_col) == 0){
    stop("Could not extract CI columns from wald_univariate output.")
  }

  ci <- coef_mat[, c(lb_col, ub_col), drop = FALSE]

  pct <- format(100 * c(alpha / 2, 1 - alpha / 2),
                digits = 3, trim = TRUE)
  colnames(ci) <- paste0(pct, " %")

  ## Subset regression parameters if requested
  if(!missing(parm)){
    if(is.character(parm)){
      ci <- ci[rownames(ci) %in% parm, , drop = FALSE]
    } else {
      ci <- ci[parm, , drop = FALSE]
    }
  }


  ## Correlation parameter intervals (on untransformed scale)
  if(!is.null(object$VhalfInv_params_estimates) &&
     !is.null(object$VhalfInv_params_vcov) &&
     all(!is.na(object$VhalfInv_params_estimates)) &&
     all(!is.na(object$VhalfInv_params_vcov))){

    par_est <- object$VhalfInv_params_estimates
    vcov_mat <- object$VhalfInv_params_vcov

    se <- sqrt(diag(vcov_mat)) /
          sqrt(object$N - object$trace_XUGX)

    for(i in seq_along(par_est)){

      work_ci <- par_est[i] + c(-cv, cv) * se[i]
      row_name <- paste0("Correlation parameter ", i)
      ci_corr <- matrix(work_ci, nrow = 1)
      colnames(ci_corr) <- colnames(ci)
      rownames(ci_corr) <- row_name
      ci <- rbind(ci, ci_corr)
    }
  }

  return(ci)
}


## =========================================================================
## logLik.lgspline: Extract Log-Likelihood from lgspline Objects
## =========================================================================
## [Change 2026-02-11] New method. Replaces the need
#  for a standalone loglik_weibull function by providing a generic
#  logLik method returning a proper "logLik" object. Can handle correlation
#   structures and prior penalizations

#' Extract Log-Likelihood from a Fitted lgspline Model
#'
#' @description
#' Returns the log-likelihood of a fitted \code{lgspline} model as a
#' \code{"logLik"} object, enabling use with \code{\link[stats]{AIC}},
#' \code{\link[stats]{BIC}}, and other model comparison tools.
#'
#' @param object A fitted \code{lgspline} model object.
#' @param include_prior Logical. Default \code{TRUE}. When \code{TRUE},
#'   the log-prior penalty
#'   \eqn{-\frac{1}{2\tilde{\sigma}^2}\sum_k \boldsymbol{\beta}_k^{\top}
#'   \boldsymbol{\Lambda}_k \boldsymbol{\beta}_k}
#'   is added to the marginal log-likelihood, giving the penalized
#'   (maximum a posteriori) log-likelihood that is coherent with the
#'   smoothing spline objective function. Set to \code{FALSE} to obtain
#'   the unpenalized marginal GLS log-likelihood, which is more
#'   appropriate when comparing models with different penalty structures,
#'   different numbers of knots, or when using external information
#'   criteria that account for smoothing degrees of freedom separately.
#' @param ... Not used.
#'
#' @details
#' For Gaussian family with identity link without a correlation structure,
#' the exact log-likelihood is computed as:
#' \deqn{-\frac{N}{2}\log(2\pi\tilde{\sigma}^2) -
#'   \frac{1}{2\tilde{\sigma}^2}\sum_{i=1}^{N}(y_i - \hat{y}_i)^2}
#'
#' When a correlation structure is present (\code{VhalfInv} is non-NULL),
#' the GLS log-likelihood is used instead:
#' \deqn{-\frac{N}{2}\log(2\pi\tilde{\sigma}^2)
#'   + \log|\mathbf{V}^{-1/2}|
#'   - \frac{1}{2\tilde{\sigma}^2}
#'     \|\mathbf{V}^{-1/2}(\mathbf{y} - \hat{\mathbf{y}})\|^2}
#'
#' where \eqn{\log|\mathbf{V}^{-1/2}|} is obtained from
#' \code{VhalfInv_logdet} if supplied (efficient path), or computed
#' directly from \code{VhalfInv} otherwise.
#'
#' When \code{include_prior = TRUE} (the default), the log-prior
#' contribution from \code{\link{prior_loglik}} is appended:
#' \deqn{-\frac{1}{2\tilde{\sigma}^2}
#'   \sum_{k=1}^{K+1}\boldsymbol{\beta}_k^{\top}\boldsymbol{\Lambda}_k
#'   \boldsymbol{\beta}_k}
#' This gives the joint penalized log-likelihood coherent with the
#' smoothing spline MAP objective. When \code{include_prior = FALSE}
#' the unpenalized marginal likelihood is returned.
#'
#' Note that this function returns the marginal (full) GLS
#' log-likelihood, not the REML log-likelihood. The REML objective
#' additionally subtracts
#' \eqn{\frac{1}{2}\log|\mathbf{X}^{\top}\mathbf{V}^{-1}\mathbf{X}
#' + \boldsymbol{\Lambda}|} to account for uncertainty in the fixed
#' effects. The marginal likelihood is the conventional choice for
#' AIC/BIC comparisons of fixed effects structure (different \code{K},
#' different formulas), and is consistent with \code{REML = FALSE} in
#' \code{lme} and \code{gls}.
#'
#' For other GLM families this method attempts to use
#' \code{family$aic()} if available (as in standard R families) to back
#' out the log-likelihood. When a correlation structure is present the
#' whitened residuals and fitted values are passed in place of the
#' originals. If \code{family$aic()} is not available, the
#' deviance-based approximation
#' \eqn{-0.5 \times \text{deviance} / \tilde{\sigma}^2} is returned,
#' which is valid for model comparison but omits a family-specific
#' constant.
#'
#' The effective degrees of freedom (\code{df} attribute) is set to
#' N-\eqn{\text{trace}(\mathbf{XUGX}^{\top})}, the smoothing spline analogue
#' of the number of parameters.
#'
#' @return An object of class \code{"logLik"} with attributes:
#' \describe{
#'   \item{df}{Effective degrees of freedom (trace of hat matrix).}
#'   \item{nobs}{Number of observations.}
#' }
#'
#' @examples
#'
#' ## Simulate data and fit model
#' set.seed(1234)
#' t <- runif(1000, -10, 10)
#' y <- 2*sin(t) + -0.06*t^2 + rnorm(length(t))
#' model_fit <- lgspline(t, y)
#'
#' ## Extract log-likelihood (penalized, default)
#' ll <- logLik(model_fit)
#' print(ll)
#'
#' ## Unpenalized marginal likelihood
#' ll_unpen <- logLik(model_fit, include_prior = FALSE)
#' print(ll_unpen)
#'
#' ## Use with AIC/BIC (penalized form recommended for smoothing splines)
#' AIC(model_fit)
#' BIC(model_fit)
#'
#' ## Compare models with different numbers of knots using unpenalized
#' ## likelihood (penalty structures differ so prior is not comparable)
#' fit_k3 <- lgspline(t, y, K = 3)
#' fit_k7 <- lgspline(t, y, K = 7)
#' AIC(fit_k3, fit_k7)
#'
#' @seealso \code{\link{lgspline}}, \code{\link{prior_loglik}},
#'   \code{\link[stats]{logLik}}, \code{\link[stats]{AIC}},
#'   \code{\link[stats]{BIC}}
#'
#' @method logLik lgspline
#' @export
logLik.lgspline <- function(object,
                            include_prior = TRUE,
                            new_weights = NULL, ...) {
  N  <- object$N
  mu <- object$ytilde
  y  <- object$y
  sigma2 <- object$sigmasq_tilde
  fam <- object$family
  if(!is.null(new_weights)){
    if(length(new_weights) %in% c(1, N)){
      wt <- new_weights
    } else {
      stop('\n\t weights argument must be a scalar or N-length vector\n')
    }
  } else {
    wt  <- object$weights
  }

  ## Effective degrees of freedom
  edf <- if(!is.null(object$trace_XUGX)) N - object$trace_XUGX else NA_real_

  ## Detect correlation structure
  has_corr <- !is.null(object$VhalfInv)

  ## log|V^{-1/2}|: prefer the cheap scalar path when available
  if(has_corr){
    logdet_Vhalfinv <- tryCatch({
      if(!is.null(object$VhalfInv_logdet) &&
         !is.null(object$VhalfInv_params_estimates)){
        object$VhalfInv_logdet(object$VhalfInv_params_estimates)
      } else {
        determinant(object$VhalfInv, logarithm = TRUE)$modulus[[1L]]
      }
    }, error = function(e) NA_real_)

    resid_w <- tryCatch(
      c(object$VhalfInv %*% (cbind(y - mu)*sqrt(wt))),
      error = function(e) NULL
    )
  }

  ll <- NA_real_

  if(fam$family == "gaussian" && fam$link == "identity"){

    if(has_corr && !is.null(resid_w) && is.finite(logdet_Vhalfinv)){
      ## Exact GLS log-likelihood
      ss_w <- sum(wt*(resid_w)^2)
      ll   <- -0.5 * N * log(2 * pi * sigma2) +
        logdet_Vhalfinv -
        0.5 * ss_w / sigma2
    } else {
      ## Standard uncorrelated Gaussian log-likelihood
      resid <- y - mu
      ss    <- sum(wt*(resid)^2)
      ll    <- -0.5 * N * log(2 * pi * sigma2) - 0.5 * ss / sigma2
    }

  } else if(!is.null(fam$aic)){

    ## When a correlation structure is present, pass whitened quantities
    #  to family$aic() so that the residual contribution reflects V^{-1}.
    #  For families where aic() depends on the raw scale (e.g. binomial
    #  with known n_i) this is an approximation; users requiring exact
    #  correlated-GLM likelihoods should supply custom_VhalfInv_loss.
    if(has_corr && !is.null(resid_w)){
      y_eval  <- tryCatch(
        c(object$VhalfInv %*% cbind(y)),
        error = function(e) y
      )
      mu_eval <- tryCatch(
        c(object$VhalfInv %*% cbind(mu)),
        error = function(e) mu
      )
    } else {
      y_eval  <- y
      mu_eval <- mu
    }

    ## Deviance resids
    dev_resids <- tryCatch(
      fam$dev.resids(y_eval, mu_eval, wt),
      error = function(e) NULL
    )

    if(!is.null(dev_resids)){
      dev <- sum(dev_resids)
      aic_val <- tryCatch(
        fam$aic(y_eval, rep(1, N), mu_eval, wt, dev),
        error = function(e) NULL
      )
      if(!is.null(aic_val) && is.finite(aic_val)){
        ll_base <- -0.5 * aic_val
        ll <- if(has_corr && is.finite(logdet_Vhalfinv)){
          ll_base + logdet_Vhalfinv
        } else {
          ll_base
        }
      }
    }
  }

  ## Fallback: deviance-based approximation
  if(!is.finite(ll)){

    if(has_corr && !is.null(resid_w)){
      y_fb  <- tryCatch(
        c(object$VhalfInv %*% cbind(y)),
        error = function(e) y
      )
      mu_fb <- tryCatch(
        c(object$VhalfInv %*% cbind(mu)),
        error = function(e) mu
      )
    } else {
      y_fb  <- y
      mu_fb <- mu
    }

    ## Deviance resids
    dev_resids <- tryCatch(
      fam$dev.resids(y_fb, mu_fb, wt),
      error = function(e) NULL
    )

    if(!is.null(dev_resids)){
      dev <- sum(dev_resids)
      ll_base <- -0.5 * dev / sigma2
      ll <- if(has_corr && is.finite(logdet_Vhalfinv)){
        ll_base + logdet_Vhalfinv
      } else {
        ll_base
      }
      warning(
        paste0(
          "Log-likelihood computed via deviance approximation",
          if(has_corr){
            " with correlation structure log-determinant correction"
          } else {
            ""
          },
          "; a family-specific constant may be omitted. ",
          "Valid for relative model comparison only."
        ),
        call. = FALSE
      )
    } else {
      warning(
        "Could not compute log-likelihood for this family.",
        call. = FALSE
      )
    }
  }

  ## Add log-prior contribution if requested and ll was successfully computed
  if(include_prior && is.finite(ll)){
    lp <- tryCatch(
      as.numeric(prior_loglik(object, sigmasq = sigma2)),
      error = function(e) {
        warning(
          "Could not evaluate prior_loglik; prior term omitted.",
          call. = FALSE
        )
        0
      }
    )
    ll <- ll + lp
  }

  attr(ll, "df")   <- edf
  attr(ll, "nobs") <- N
  class(ll)        <- "logLik"
  return(ll)
}
