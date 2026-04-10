#' lgspline: S3 Methods
#'
#' @description
#' S3 methods for lgspline objects: print, summary, coef, plot, predict,
#' confint, logLik, and inference helpers.
#'
#' @docType methods
#' @keywords internal
#' @name lgspline-methods
#' @rdname lgspline-methods
#' @aliases lgspline-methods
NULL

#' Print Method for lgspline Objects
#'
#' Prints a concise summary of the fitted model to the console.
#'
#' @param x An lgspline model object.
#' @param ... Not used.
#'
#' @return Invisibly returns \code{x}.
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


#' Summary Method for lgspline Objects
#'
#' @param object An lgspline model object.
#' @param ... Not used.
#'
#' @return An object of class \code{summary.lgspline}, a list containing:
#' \describe{
#'   \item{model_family}{The \code{\link[stats]{family}} object.}
#'   \item{observations}{Number of observations N.}
#'   \item{predictors}{Number of predictor variables q.}
#'   \item{knots}{Number of knots K.}
#'   \item{basis_functions}{Basis functions per partition p.}
#'   \item{estimate_dispersion}{"Yes" or "No".}
#'   \item{cv}{Critical value used for confidence intervals.}
#'   \item{coefficients}{Coefficient matrix from \code{\link{wald_univariate}},
#'         or a single-column estimate matrix if \code{return_varcovmat = FALSE}.}
#'   \item{sigmasq_tilde}{Estimated dispersion \eqn{\tilde{\sigma}^2}.}
#'   \item{trace_XUGX}{Trace of the hat matrix \eqn{\mathrm{trace}(\mathbf{XUGX}^\top)}.}
#'   \item{N}{Number of observations.}
#' }
#'
#' @export
summary.lgspline <- function(object, ...) {
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

  ## Wald inference table (requires return_varcovmat = TRUE)
  #  Column names are compatible with stats::printCoefmat()
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


#' Print Method for lgspline Summaries
#'
#' Displays a formatted model summary using \code{\link[stats]{printCoefmat}}
#' for the coefficient table.
#'
#' @param x A \code{summary.lgspline} object.
#' @param ... Not used.
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link[stats]{printCoefmat}}
#' @export
print.summary.lgspline <- function(x, ...) {
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
    ## Reorder so p-value column is last, as required by printCoefmat
    pval_col <- grep("^Pr\\(", colnames(x$coefficients))
    if(length(pval_col) > 0){
      other_cols <- setdiff(seq_len(ncol(x$coefficients)), pval_col)
      print_mat <- x$coefficients[, c(other_cols, pval_col), drop = FALSE]
      pval_col_new <- ncol(print_mat)
      tst_col <- grep("value$", colnames(print_mat))
      tst_col <- tst_col[!tst_col %in% pval_col_new]
      stats::printCoefmat(print_mat,
                          cs.ind = 1:2,
                          tst.ind = if(length(tst_col) > 0) tst_col[1] else NULL,
                          P.values = TRUE,
                          has.Pvalue = TRUE,
                          signif.stars = TRUE)
    } else {
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


#' Find the Extremum of a Fitted lgspline
#'
#' Finds the global maximum or minimum of a fitted lgspline using L-BFGS-B,
#' with options for partition-based heuristics, stochastic exploration, and
#' custom objective functions (e.g., acquisition functions for Bayesian
#' optimization).
#'
#' @param object A fitted lgspline model object.
#' @param vars Integer or character vector; indices or names of predictors to
#'        optimize over. Default NULL optimizes all predictors.
#' @param quick_heuristic Logical; if TRUE (default) searches only the
#'        best-performing partition. If FALSE, initiates searches from all
#'        partition local maxima.
#' @param initial Numeric vector; optional starting values. Useful for
#'        fixing binary predictors. Default NULL.
#' @param B_predict List; optional coefficient list for prediction, e.g.
#'        from \code{\link{generate_posterior}}. Default NULL uses
#'        \code{object$B}.
#' @param minimize Logical; find minimum instead of maximum. Default FALSE.
#' @param stochastic Logical; add noise during optimization for exploration.
#'        Default FALSE.
#' @param stochastic_draw Function; generates noise for stochastic
#'        optimization. Takes \code{mu}, \code{sigma}, and \code{...}.
#'        Default \code{rnorm(length(mu), mu, sigma)}.
#' @param sigmasq_predict Numeric; variance for stochastic draws.
#'        Default \code{object$sigmasq_tilde}.
#' @param custom_objective_function Function; optional custom objective.
#'        Takes \code{mu}, \code{sigma}, \code{y_best}, \code{...}.
#'        Default NULL.
#' @param custom_objective_derivative Function; optional gradient of
#'        \code{custom_objective_function}. Takes \code{mu}, \code{sigma},
#'        \code{y_best}, \code{d_mu}, \code{...}. Default NULL.
#' @param ... Additional arguments passed to internal optimization routines.
#'
#' @return A list with elements:
#' \describe{
#'   \item{t}{Numeric vector; predictor values at the extremum.}
#'   \item{y}{Numeric; objective value at the extremum.}
#' }
#'
#' @examples
#'
#' set.seed(1234)
#' t <- runif(1000, -10, 10)
#' y <- 2*sin(t) + -0.06*t^2 + rnorm(length(t))
#' model_fit <- lgspline(t, y)
#' plot(model_fit)
#'
#' max_point <- find_extremum(model_fit)
#' min_point <- find_extremum(model_fit, minimize = TRUE)
#' abline(v = max_point$t, col = 'blue')
#' abline(v = min_point$t, col = 'red')
#'
#' ## Expected improvement acquisition function
#' ei_obj <- function(mu, sigma, y_best, ...) {
#'   d <- y_best - mu
#'   d * pnorm(d/sigma) + sigma * dnorm(d/sigma)
#' }
#' ei_deriv <- function(mu, sigma, y_best, d_mu, ...) {
#'   d <- y_best - mu
#'   z <- d/sigma
#'   d_z <- -d_mu/sigma
#'   pnorm(z)*d_mu - d*dnorm(z)*d_z + sigma*z*dnorm(z)*d_z
#' }
#'
#' post_draw <- generate_posterior(model_fit)
#' acq <- find_extremum(model_fit,
#'                      stochastic = TRUE,
#'                      B_predict = post_draw$post_draw_coefficients,
#'                      sigmasq_predict = post_draw$post_draw_sigmasq,
#'                      custom_objective_function = ei_obj,
#'                      custom_objective_derivative = ei_deriv)
#' abline(v = acq$t, col = 'green')
#'
#' @seealso \code{\link{lgspline}}, \code{\link{generate_posterior}}
#' @export
find_extremum <- function(object,
                          vars = NULL,
                          quick_heuristic = TRUE,
                          initial = NULL,
                          B_predict = NULL,
                          minimize = FALSE,
                          stochastic = FALSE,
                          stochastic_draw = function(mu, sigma, ...){
                            N <- length(mu)
                            rnorm(N, mu, sigma)
                          },
                          sigmasq_predict = object$sigmasq_tilde,
                          custom_objective_function = NULL,
                          custom_objective_derivative = NULL,
                          ...) {
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


#' Generate Posterior Samples from a Fitted lgspline
#'
#' Draws from the posterior distribution of model coefficients, with optional
#' dispersion sampling, posterior predictive draws, and propagation of
#' uncertainty in estimated correlation parameters.
#'
#' Uses a Laplace approximation centred at the MAP estimate for non-Gaussian
#' responses.
#'
#' @details
#' \strong{Dispersion posterior.}
#' When \code{draw_dispersion = TRUE}, \eqn{\sigma^2} is drawn from
#' \deqn{
#'   \sigma^2 \mid \mathbf{y} \sim
#'   \mathrm{InvGamma}(\alpha_1, \alpha_2),
#' }
#' where
#' \deqn{
#'   \alpha_1 = \theta_1 + \tfrac{1}{2}(N - s \cdot \mathrm{tr}(\mathbf{H})),
#'   \quad
#'   \alpha_2 = \theta_2 + \tfrac{1}{2}(N - s \cdot \mathrm{tr}(\mathbf{H}))
#'              \tilde{\sigma}^2,
#' }
#' \eqn{\mathbf{H} = \mathbf{XUGX}^\top} is the hat matrix, \eqn{s = 1}
#' when \code{unbias_dispersion = TRUE} (else \eqn{s = 0}), and
#' \eqn{\theta_1 = \theta_2 = 0} recovers an improper uniform prior.
#'
#' \strong{Correlation parameter posterior.}
#' When \code{draw_correlation = TRUE} and the fitted model contains an
#' estimated correlation structure, each draw first samples
#' \eqn{\boldsymbol{\rho}} from
#' \deqn{
#'   \boldsymbol{\rho}^{(m)} \sim
#'   \mathcal{N}(\hat{\boldsymbol{\rho}}_{\mathrm{REML}},
#'               \mathbf{H}^{-1}_{\mathrm{BFGS}}),
#' }
#' rebuilds the posterior covariance under the drawn correlation
#' structure (reusing all pre-computed design matrices, constraints, and
#' penalty matrices) and then draws coefficients from the updated
#' posterior. Knot placement, partitioning, coefficient re-estimation,
#' and GCV tuning are skipped entirely.
#' Draws producing non-positive-definite correlation matrices are
#' rejected and redrawn (up to 50 attempts).
#'
#' When \code{draw_correlation = FALSE} (default), correlation parameters
#' are fixed at their estimated values.
#'
#' \strong{Inequality constraints.}
#' Active QP inequalities can be enforced during posterior sampling via
#' elliptical slice sampling, producing draws from the corresponding
#' truncated multivariate normal posterior on the coefficient scale.
#' The public \code{enforce_qp_constraints} argument is forwarded to the
#' stored sampler for both the standard and correlation-aware posterior
#' paths.
#'
#' @param object A fitted \code{lgspline} model object.
#' @param new_sigmasq_tilde Numeric; dispersion \eqn{\tilde{\sigma}^2} used
#'        as the point estimate when \code{draw_dispersion = FALSE}.
#'        Default \code{object$sigmasq_tilde}.
#' @param new_predictors Matrix; predictor matrix for posterior predictive
#'        sampling. Default uses in-sample predictors.
#' @param theta_1 Numeric; shape increment for the inverse-gamma prior on
#'        \eqn{\sigma^2}. Default 0.
#' @param theta_2 Numeric; rate increment for the inverse-gamma prior.
#'        Default 0.
#' @param posterior_predictive_draw Function; sampler for posterior predictive
#'        realisations. Must accept \code{N}, \code{mean},
#'        \code{sqrt_dispersion}, \code{...}. Default \code{rnorm}.
#' @param draw_dispersion Logical; sample \eqn{\sigma^2} from its posterior.
#'        Default TRUE.
#' @param include_posterior_predictive Logical; generate posterior predictive
#'        draws at \code{new_predictors}. Default FALSE.
#' @param num_draws Positive integer; number of draws. Default 1.
#' @param enforce_qp_constraints Logical; if TRUE, enforce active QP
#'        inequality constraints during posterior sampling via the stored
#'        elliptical-slice constrained sampler. Default TRUE.
#' @param draw_correlation Logical; propagate correlation parameter
#'        uncertainty. Requires \code{VhalfInv_fxn} and
#'        \code{VhalfInv_params_estimates} in the fitted object.
#'        Default FALSE.
#' @param correlation_param_mean Numeric vector; mean of the approximate
#'        normal posterior for correlation parameters on the unbounded
#'        (working) scale. Default: \code{object$VhalfInv_params_estimates}.
#' @param correlation_param_vcov Matrix; variance-covariance for correlation
#'        parameter draws. Default: inverse Hessian from BFGS
#'        (\code{object$VhalfInv_params_vcov}).
#' @param correlation_VhalfInv_fxn Function; maps correlation parameter
#'        vector to \eqn{\mathbf{V}^{-1/2}}. Default \code{object$VhalfInv_fxn}.
#' @param correlation_Vhalf_fxn Function or NULL; maps to
#'        \eqn{\mathbf{V}^{1/2}}. Passed through to
#'        \code{\link{generate_posterior_correlation}}; the current correlation-aware posterior path only requires \code{correlation_VhalfInv_fxn}.
#' @param correlation_param_vcov_scale NULL or numeric; if supplied,
#'        divides a user-supplied \code{correlation_param_vcov} by this
#'        value before passing it to
#'        \code{\link{generate_posterior_correlation}}. When NULL,
#'        no additional scaling is applied.
#' @param include_warnings Logical; emit warnings for degenerate draws,
#'        constraint violations, etc. Default TRUE.
#' @param ... Additional arguments forwarded to the GLM weight function,
#'        dispersion function, and \code{posterior_predictive_draw}.
#'
#' @return When \code{num_draws = 1}, a named list:
#' \describe{
#'   \item{post_draw_coefficients}{List of length K+1; per-partition
#'         coefficient vectors on the original scale.}
#'   \item{post_draw_sigmasq}{Drawn (or fixed) dispersion.}
#'   \item{post_pred_draw}{Posterior predictive vector (only when
#'         \code{include_posterior_predictive = TRUE}).}
#'   \item{post_draw_correlation_params}{Drawn correlation parameters
#'         on the working scale (only when \code{draw_correlation = TRUE}).}
#' }
#' When \code{num_draws > 1}, each element becomes a list of length
#' \code{num_draws}, and \code{post_pred_draw} (if requested) is an
#' \eqn{N_{\mathrm{new}} \times M} matrix, where \eqn{M = \mathrm{num\_draws}}.
#'
#' @examples
#'
#' \donttest{
#' set.seed(1234)
#' n_blocks <- 100; block_size <- 5; N <- n_blocks * block_size
#' rho_true <- 0.3
#' t <- seq(-5, 5, length.out = N)
#' true_mean <- sin(t)
#' errors <- Reduce("rbind",
#'   lapply(1:n_blocks, function(i) {
#'     sigma <- diag(block_size) + rho_true *
#'       (matrix(1, block_size, block_size) - diag(block_size))
#'     matsqrt(sigma) %*% rnorm(block_size)
#'   })
#' )
#' y <- true_mean + errors * 0.5
#'
#' model_fit <- lgspline(t, y,
#'   K = 3,
#'   correlation_id = rep(1:n_blocks, each = block_size),
#'   correlation_structure = "exchangeable",
#'   include_warnings = FALSE
#' )
#'
#' ## Propagate correlation uncertainty across 50 draws
#' post <- generate_posterior(model_fit,
#'   draw_correlation = TRUE, num_draws = 50,
#'   include_warnings = FALSE
#' )
#'
#' ## Fixed correlation parameters for comparison
#' post_fixed <- generate_posterior(model_fit, num_draws = 50)
#'
#' corr_draws <- unlist(post$post_draw_correlation_params)
#' rho_draws <- exp(-exp(corr_draws))
#' print(summary(rho_draws))
#' }
#'
#' @seealso
#' \code{\link{lgspline}},
#' \code{\link{generate_posterior_correlation}},
#' \code{\link{wald_univariate}}
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
                               enforce_qp_constraints = TRUE,
                               draw_correlation = FALSE,
                               correlation_param_mean = NULL,
                               correlation_param_vcov = NULL,
                               correlation_VhalfInv_fxn = NULL,
                               correlation_Vhalf_fxn = NULL,
                               correlation_param_vcov_scale = NULL,
                               include_warnings = TRUE,
                               ...) {

  has_qp <- !is.null(object$quadprog_list) &&
    !identical(object$quadprog_list, list(NA)) &&
    !is.null(object$quadprog_list$qp_Amat)

  ## Route to correlation-aware sampler when requested
  if(draw_correlation){

    ## Resolve the scaling divisor for the correlation parameter vcov.
    #  When the user supplies correlation_param_vcov_scale, use it directly.
    #  Otherwise default to 1.
    if(is.null(correlation_param_vcov_scale)){
      scale_by <- 1
    } else {
      scale_by <- correlation_param_vcov_scale
    }

    ## Apply scaling to the user-supplied vcov before handoff.
    #  If the user supplied correlation_param_vcov, scale it here.
    #  If NULL, generate_posterior_correlation builds it internally
    #  from object$VhalfInv_params_vcov with the same default scaling.
    if(!is.null(correlation_param_vcov)){
      scaled_vcov <- correlation_param_vcov / scale_by
    } else {
      scaled_vcov <- NULL
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
        enforce_qp_constraints     = enforce_qp_constraints,
        correlation_param_mean     = correlation_param_mean,
        correlation_param_vcov_sc  = scaled_vcov,
        correlation_VhalfInv_fxn   = correlation_VhalfInv_fxn,
        correlation_Vhalf_fxn      = correlation_Vhalf_fxn,
        include_warnings           = include_warnings,
        ...
      )
    )
  }

  ## Standard path: delegate to closure stored in the object
  internal_genpost_func <- object$generate_posterior
  if(!is.null(internal_genpost_func) && is.function(internal_genpost_func)){
    if(is.null(new_predictors)){
      return(internal_genpost_func(
        new_sigmasq_tilde            = new_sigmasq_tilde,
        theta_1                      = theta_1,
        theta_2                      = theta_2,
        posterior_predictive_draw    = posterior_predictive_draw,
        draw_dispersion              = draw_dispersion,
        include_posterior_predictive = include_posterior_predictive,
        num_draws                    = num_draws,
        enforce_qp_constraints       = enforce_qp_constraints,
        ...
      ))
    } else {
      return(internal_genpost_func(
        new_sigmasq_tilde            = new_sigmasq_tilde,
        new_predictors               = new_predictors,
        theta_1                      = theta_1,
        theta_2                      = theta_2,
        posterior_predictive_draw    = posterior_predictive_draw,
        draw_dispersion              = draw_dispersion,
        include_posterior_predictive = include_posterior_predictive,
        num_draws                    = num_draws,
        enforce_qp_constraints       = enforce_qp_constraints,
        ...
      ))
    }
  } else {
    stop("Internal generate_posterior method not found in lgspline object.")
  }
}


#' Generate Posterior Samples Propagating Correlation Parameter Uncertainty
#'
#' Called internally by \code{\link{generate_posterior}} when
#' \code{draw_correlation = TRUE}, but can be used directly for finer control.
#' For each draw, samples the correlation parameter vector from its approximate
#' normal posterior, rebuilds the posterior covariance under that drawn
#' correlation structure without re-solving for a new coefficient mode, then
#' draws coefficients from the updated posterior.
#'
#' @details
#' Each draw proceeds in three steps:
#' \enumerate{
#'   \item \strong{Draw correlation parameters.}
#'     \eqn{\boldsymbol{\rho}^{(m)} \sim
#'     \mathcal{N}(\hat{\boldsymbol{\rho}}_{\mathrm{REML}},
#'               \mathbf{H}^{-1}_{\mathrm{BFGS}})}
#'     on the unbounded working scale. Draws producing a non-PD correlation
#'     matrix are rejected and redrawn (up to 50 attempts); if all fail, the
#'     point estimate is used with a warning.
#'
#'   \item \strong{Rebuild posterior covariance.}
#'     Using the already-expanded \eqn{\mathbf{X}_k}, \eqn{\mathbf{A}}, and
#'     \eqn{\boldsymbol{\Lambda}} from the original fit, recompute only the
#'     covariance-side quantities implied by the drawn correlation structure:
#'     \deqn{
#'       \mathbf{G}_{\mathrm{correct}}^{(m)} =
#'       \left(\mathbf{X}^{\top} \mathbf{W}\mathbf{D}
#'       \mathbf{V}^{-1}(\boldsymbol{\rho}^{(m)}) \mathbf{X}
#'       + \boldsymbol{\Lambda}\right)^{-1},
#'     }
#'     and from this the updated constraint projection
#'     \eqn{\mathbf{U}^{(m)}} and effective degrees of freedom
#'     \eqn{\mathrm{trace}(\mathbf{H}^{(m)})}.
#'
#'     The coefficient mode \eqn{\hat{\boldsymbol{\beta}}_{\mathrm{raw}}}
#'     and fitted mean \eqn{\tilde{\mathbf{y}}} are held fixed at the
#'     original fit values. Knot placement, partitioning, polynomial
#'     expansion, penalty tuning, and coefficient re-estimation are all
#'     skipped entirely.
#'
#'   \item \strong{Draw coefficients.}
#'     Updated quantities (\code{U}, \code{Ghalf_correct}, \code{VhalfInv},
#'     \code{sigmasq_tilde}, \code{trace_XUGX}) are passed to the stored
#'     closure via \code{override_*} arguments so that the draw is centred
#'     at the original mode but uses the covariance implied by the drawn
#'     correlation structure. The stored mode \code{object$B_raw} is passed
#'     as \code{override_B_raw} so it is not recomputed.
#' }
#'
#' \strong{Why the mode is held fixed.}
#' Re-solving for a new MAP estimate under each drawn correlation structure
#' is expensive, requires iterative solvers, and risks convergence failures
#' on draws far from the REML estimate. The posterior draw is centred at the
#' original mode, which remains a reasonable approximation when the REML
#' surface is not sharply peaked. The covariance update captures the primary
#' effect of correlation uncertainty on posterior width and shape.
#'
#' \strong{BFGS inverse Hessian caveat.}
#' The BFGS inverse Hessian approximation for the correlation parameter
#' covariance is asymptotically valid but may be poor for small samples,
#' near-boundary estimates, or multimodal REML surfaces. It is not guaranteed
#' to converge to the observed information matrix. Users should inspect
#' \code{object$VhalfInv_params_vcov} before relying on these draws.
#'
#' @param object A fitted \code{lgspline} object with a correlation structure
#'        (i.e., \code{VhalfInv_fxn} and \code{VhalfInv_params_estimates}
#'        present, or supplied via override arguments).
#' @param new_sigmasq_tilde Numeric; dispersion starting value when
#'        \code{draw_dispersion = FALSE}. Default \code{object$sigmasq_tilde}.
#' @param new_predictors Matrix or NULL; predictor matrix for posterior
#'        predictive sampling. Default uses in-sample predictors.
#' @param theta_1 Numeric; shape increment for the inverse-gamma prior.
#'        Default 0.
#' @param theta_2 Numeric; rate increment for the inverse-gamma prior.
#'        Default 0.
#' @param posterior_predictive_draw Function; sampler for posterior predictive
#'        realisations. Default \code{rnorm}.
#' @param draw_dispersion Logical; sample \eqn{\sigma^2} within each
#'        draw. Default TRUE.
#' @param include_posterior_predictive Logical; generate posterior predictive
#'        draws. Default FALSE.
#' @param num_draws Positive integer; number of draws (each requires one
#'        correlation parameter sample and one covariance rebuild).
#'        Default 1.
#' @param enforce_qp_constraints Logical; if TRUE, enforce active QP
#'        inequality constraints during each coefficient draw via the
#'        stored elliptical-slice constrained sampler. Default TRUE.
#' @param correlation_param_mean Numeric vector or NULL; mean of the
#'        approximate normal posterior on the working scale. Default:
#'        \code{object$VhalfInv_params_estimates}. Supplying this allows
#'        correlation draws for models fit with a fixed (non-optimised)
#'        correlation structure.
#' @param correlation_param_vcov_sc Matrix or NULL; variance-covariance
#'        on the working scale. Default:
#'        \code{object$VhalfInv_params_vcov}. No further scaling is
#'        applied within this function.
#' @param correlation_VhalfInv_fxn Function or NULL; maps parameter vector
#'        to \eqn{\mathbf{V}^{-1/2}}. Default \code{object$VhalfInv_fxn}.
#' @param correlation_Vhalf_fxn Function or NULL; maps to
#'        \eqn{\mathbf{V}^{1/2}}. Not consumed in the current method body
#'        but resolved for interface consistency and potential future use.
#' @param include_warnings Logical; emit warnings. Default TRUE.
#' @param ... Additional arguments forwarded to the GLM weight function,
#'        dispersion function, and \code{posterior_predictive_draw}.
#'
#' @return When \code{num_draws = 1}, a named list:
#' \describe{
#'   \item{post_draw_coefficients}{List of length K+1; per-partition
#'         coefficient vectors on the original scale.}
#'   \item{post_draw_sigmasq}{Drawn dispersion.}
#'   \item{post_pred_draw}{Posterior predictive vector (only when
#'         \code{include_posterior_predictive = TRUE}).}
#'   \item{post_draw_correlation_params}{Drawn correlation parameters
#'         on the working scale.}
#' }
#' When \code{num_draws > 1}:
#' \describe{
#'   \item{post_draw_coefficients}{List of \code{num_draws} lists of K+1
#'         coefficient vectors.}
#'   \item{post_draw_sigmasq}{List of \code{num_draws} scalars.}
#'   \item{post_pred_draw}{\eqn{N_{\mathrm{new}} \times M} matrix, where
#'         \eqn{M = \mathrm{num\_draws}} (only when
#'         \code{include_posterior_predictive = TRUE}).}
#'   \item{post_draw_correlation_params}{List of \code{num_draws} vectors.}
#' }
#'
#' @examples
#' ## See ?generate_posterior for a complete worked example.
#'
#' @seealso
#' \code{\link{generate_posterior}},
#' \code{\link{lgspline}},
#' \code{\link{lgspline.fit}}
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
    enforce_qp_constraints = TRUE,
    correlation_param_mean = NULL,
    correlation_param_vcov_sc = NULL,
    correlation_VhalfInv_fxn = NULL,
    correlation_Vhalf_fxn = NULL,
    include_warnings = TRUE,
    ...
) {


  ## 1. Resolve correlation parameter mean, vcov, and functions.
  #     These are required for sampling phi. If the user did not supply
  #     overrides, pull from the fitted object.

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

  ## Default vcov: no scaling
  if(is.null(correlation_param_vcov_sc)){
    correlation_param_vcov_sc <- object$VhalfInv_params_vcov
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

  ## VhalfInv function: maps unbounded parameters to V^{-1/2}
  if(is.null(correlation_VhalfInv_fxn)){
    correlation_VhalfInv_fxn <- object$VhalfInv_fxn
    if(is.null(correlation_VhalfInv_fxn)){
      stop(
        "\n\t No VhalfInv_fxn found and no ",
        "'correlation_VhalfInv_fxn' was supplied.\n"
      )
    }
  }

  ## Vhalf function: resolved for interface consistency, not consumed.
  #  The covariance-only path only needs VhalfInv.
  if(is.null(correlation_Vhalf_fxn)){
    correlation_Vhalf_fxn <- object$Vhalf_fxn
  }


  ## 2. Cholesky factor of correlation parameter vcov for sampling.
  #     If the matrix is not positive definite (can happen with near-
  #     degenerate BFGS approximations), add a small ridge.

  vcov_chol <- tryCatch({
    chol(correlation_param_vcov_sc)
  }, error = function(e){
    ## The BFGS inverse Hessian can be indefinite when the optimizer
    #  converges poorly or the REML surface is nearly flat. A token
    #  ridge is not sufficient when eigenvalues are actually negative.
    #  Project to the nearest PD matrix by clamping eigenvalues, then
    #  add a small ridge for numerical headroom.
    eig <- eigen(correlation_param_vcov_sc, symmetric = TRUE)
    min_eval <- max(abs(eig$values)) * sqrt(.Machine$double.eps)
    eig$values <- pmax(eig$values, min_eval)
    vcov_pd <- crossprod(t(eig$vectors) * sqrt(eig$values))
    if(include_warnings){
      warning(
        "\n\t correlation_param_vcov_sc is not positive definite; ",
        "projecting to nearest PD matrix.\n"
      )
    }
    chol(vcov_pd)
  })

  ## 3. Extract stored components needed for the covariance rebuild.
  K  <- object$K
  nc <- object$p          # p_expansions: basis terms per partition
  N  <- object$N
  P_total <- nc * (K + 1) # total coefficient dimension

  A  <- object$A
  nca <- if(!is.null(A)) ncol(A) else 0

  ## Standardized design matrices (the scale on which the posterior
  #  sampler operates). These are the stored X, rescaled.
  X_std <- lapply(object$X, object$std_X)

  ## Full block-diagonal standardized design: N x P in partition order
  X_full_std <- collapse_block_diagonal(X_std)

  ## Permutation from original observation order to partition order.
  #  Under correlation, algebra must respect this mapping because
  #  VhalfInv is in original observation order while X is in
  #  partition order.
  order_list <- object$order_list
  og_order   <- object$og_order
  po <- unlist(order_list)  # partition order permutation

  y_og <- object$y         # response in original observation order
  family <- object$family
  penalties <- object$penalties
  expansion_scales <- object$expansion_scales
  mean_y <- object$mean_y
  sd_y   <- object$sd_y

  ## Observation weights D in original observation order.
  #  For the penalized information matrix we need sqrt(W * D) row-wise
  #  on the design, where W is the GLM working weight and D is the
  #  observation weight.
  wts_og <- object$weights
  if(is.list(wts_og)){
    obs_wts_vec <- unlist(wts_og)[og_order]
  } else {
    obs_wts_vec <- wts_og
  }

  ## GLM weight function.
  #  For the penalized information matrix we need the mean-variance
  #  relationship evaluated at the fitted values:
  #    gram_gls = X' (W D) V^{-1} X + Lambda
  #  where W = glm_weight_function(mu, y, ...).
  #  For Gaussian identity, W = 1 everywhere, so the weight is just D.
  #  For GLMs (logistic, Poisson, etc.), W carries the inverse
  #  variance function and is essential for correctness.
  #
  #  glm_weight_function is stored in .fit_call_args. If it is not
  #  available (e.g. objects created before .fit_call_args was introduced),
  #  fall back to the family$variance, which is correct for the default
  #  weight function.
  if(!is.null(object$.fit_call_args$glm_weight_function)){
    glm_weight_function <- object$.fit_call_args$glm_weight_function
  } else {
    ## Fallback: replicate the default glm_weight_function from lgspline()
    glm_weight_function <- function(mu, y, order_indices, family,
                                    dispersion, observation_weights, ...){
      if(any(!is.null(observation_weights))){
        family$variance(mu) * observation_weights
      } else {
        family$variance(mu)
      }
    }
  }

  ## Whether partition-specific penalties exist
  unique_pp <- (length(penalties$L_partition_list) == (K + 1))


  ## 4. Precompute the full block-diagonal penalty matrix Lambda_full.
  #     This is fixed across draws and can be computed once.
  #     Lambda_full = blockdiag(Lambda + L_partition_list[[k]])

  Lambda_full <- collapse_block_diagonal(
    lapply(1:(K + 1), function(k){
      if(unique_pp){
        penalties$Lambda + penalties$L_partition_list[[k]]
      } else {
        penalties$Lambda
      }
    })
  )


  ## 5. Determine whether we are Gaussian identity.
  #     This controls whether sigmasq_tilde is recomputed from
  #     whitened residuals under each drawn VhalfInv, or held fixed.

  is_gaussian_identity <- (paste0(family)[1] == "gaussian" &&
                             paste0(family)[2] == "identity")

  ## The fixed fitted mean and mode from the original fit.
  #  These are the quantities we explicitly do NOT recompute.
  ytilde_fixed <- object$ytilde
  B_raw_fixed  <- object$B_raw


  ## 5b. Precompute GLM-weighted design matrix.
  #      The penalized information matrix under correlation is:
  #        gram_gls = X' (W~ D) V^{-1} X + Lambda
  #      where W~ = glm_weight_function(mu, y, ...) is the GLM working
  #      weight (inverse mean-variance relationship) and D is the
  #      observation weight vector.
  #
  #      Since we hold ytilde fixed, the GLM weights W~ do not change
  #      across draws. We can precompute (W~ D)^{1/2} X once, so the
  #      per-draw helper only needs to apply V^{-1/2}.
  #
  #      For Gaussian identity, W~ = 1 everywhere and the combined
  #      weight reduces to sqrt(D). For GLMs (logistic, Poisson, etc.),
  #      W~ carries the inverse variance function and is essential for
  #      correctness of the information matrix.
  #
  #      We pass rep(1, N) for observation_weights inside
  #      glm_weight_function to avoid double-counting D, then fold D
  #      in separately. This mirrors the varcovmat block in lgspline.fit.
  W_glm <- c(glm_weight_function(
    ytilde_fixed[po],       # mu in partition order
    y_og[po],               # y in partition order
    1:N,                    # order_indices
    family,
    object$sigmasq_tilde,   # dispersion at the fixed estimate
    rep(1, N),              # unit obs weights to avoid double-counting D
    ...
  ))
  W_glm <- pmax(W_glm, .Machine$double.eps)

  ## Observation weights D in partition order
  D_po <- obs_wts_vec[po]

  ## Combined weight applied after whitening, matching the main GEE fit.
  combined_wt_po <- sqrt(W_glm * D_po)


  ## 6. Helper: rebuild posterior covariance given a drawn VhalfInv.
  #
  #     This replaces the previous .reestimate_with_VhalfInv() which
  #     re-solved for B_raw, B, and ytilde. Here we only rebuild:
  #       - Ghalf_correct  (dense correlated posterior square root)
  #       - U              (equality-constraint projection, G-dependent)
  #       - trace_XUGX     (effective df for dispersion draw)
  #       - sigmasq_tilde  (Gaussian identity: whitened residual-based;
  #                          non-Gaussian: held at original estimate)
  #
  #     The mode is held fixed at object$B_raw. The fitted mean is
  #     held fixed at object$ytilde.

  .rebuild_posterior_cov_with_VhalfInv <- function(VhalfInv_draw){

    ## Permute VhalfInv into partition order.
    #  VhalfInv_draw is N x N in original observation order.
    #  X_full_std is N x P in partition order.
    #  We need V^{-1/2}[po, po] to align with the design matrix.
    VhalfInv_po <- VhalfInv_draw[po, po]

    ## Whiten then apply the combined row weights, matching the solver path:
    ##   X' V^{-1/2} W~ D V^{-1/2} X.
    VhalfInvX <- VhalfInv_po %**% X_full_std
    VhalfInvX <- t(t(VhalfInvX) * combined_wt_po)

    ## Full penalized GLS Gram in the whitened system + Lambda  (P x P)
    #  This is the dense correlated information matrix. Under correlation,
    #  the exact posterior factor is based on this full-system Gram,
    #  NOT the stored diagonal blocks of Ghalf.
    gram_gls <- crossprod(VhalfInvX) + Lambda_full

    ## G_correct = gram_gls^{-1}     (needed for U)
    ## Ghalf_correct = gram_gls^{-1/2}  (posterior square root factor)
    G_correct     <- invert(gram_gls)
    Ghalf_correct <- matinvsqrt(gram_gls)

    ## Constraint projection U = I - G A (A' G A)^{-1} A'
    #  The feasible equality subspace itself is fixed by A, but the
    #  projection matrix U depends on the current G. If G changes
    #  with the drawn correlation structure, U changes too.
    #  If A is empty or NULL, U reduces to identity.
    if(!is.null(A) && nca > 0){
      GA      <- G_correct %**% A
      AGA_inv <- invert(crossprod(A, GA))
      U_draw  <- diag(P_total) - GA %**% tcrossprod(AGA_inv, A)
    } else {
      U_draw <- diag(P_total)
    }

    ## Effective df: trace(X' W~ D V^{-1} X U G_correct)
    #  = || V^{-1/2} (W~ D)^{1/2} X U G_correct^{1/2} ||_F^2
    #  This matters because the stored posterior closure uses it in
    #  the inverse-gamma dispersion draw. If the effective df changes
    #  under the drawn correlation structure, leaving it fixed is
    #  internally inconsistent.
    UGhalf_draw     <- U_draw %**% Ghalf_correct
    trace_XUGX_draw <- sum((VhalfInvX %**% UGhalf_draw)^2)

    ## Dispersion
    sigmasq_draw <- object$sigmasq_tilde

    ## Return only the covariance-side pieces. No B_raw, no B, no ytilde.
    list(
      Ghalf_correct = Ghalf_correct,
      U             = U_draw,
      trace_XUGX    = trace_XUGX_draw,
      sigmasq_tilde = sigmasq_draw,
      VhalfInv      = VhalfInv_draw
    )
  }


  ## 7. Single draw: sample phi, rebuild covariance, draw coefficients.

  .one_corr_draw <- function(){

    ## 7a. Draw correlation parameters from approximate normal posterior.
    #     Reject draws that produce non-finite or wrong-dimension V^{-1/2}.
    max_corr_reject <- 50L
    n_corr_reject   <- 0L
    phi_draw         <- NULL
    VhalfInv_draw    <- NULL

    repeat {
      z <- rnorm(n_corr_par)
      phi_candidate <- correlation_param_mean + c(crossprod(vcov_chol, z))

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
        phi_draw      <- phi_candidate
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
        phi_draw      <- correlation_param_mean
        VhalfInv_draw <- correlation_VhalfInv_fxn(phi_draw)
        break
      }
    }

    ## 7b. Rebuild posterior covariance under drawn VhalfInv.
    #     This is the covariance-only replacement for the old
    #     .reestimate_with_VhalfInv(). The coefficient mode stays
    #     fixed at object$B_raw.
    cov_rebuild <- try(
      .rebuild_posterior_cov_with_VhalfInv(VhalfInv_draw),
      silent = TRUE
    )

    if(inherits(cov_rebuild, "try-error")){
      if(include_warnings){
        warning(
          "\n\t Covariance rebuild with drawn correlation parameters ",
          "failed. Falling back to original model for this draw.\n"
        )
      }
      ## Fall back: draw from original posterior with no overrides
      one_draw <- object$generate_posterior(
        new_sigmasq_tilde            = new_sigmasq_tilde,
        theta_1                      = theta_1,
        theta_2                      = theta_2,
        posterior_predictive_draw    = posterior_predictive_draw,
        draw_dispersion              = draw_dispersion,
        include_posterior_predictive = include_posterior_predictive,
        num_draws                    = 1,
        enforce_qp_constraints       = enforce_qp_constraints,
        ...
      )
      one_draw$post_draw_correlation_params <- phi_draw
      return(one_draw)
    }

    ## 7c. Draw coefficients from updated posterior via overrides.
    #
    #     The closure captures model_fit by reference, so mutating a
    #     shallow copy of 'object' would have no effect. Instead we
    #     pass the rebuilt covariance-side quantities through
    #     override_* args:
    #
    #       override_B_raw          = object$B_raw  (FIXED mode)
    #       override_U              = U_draw
    #       override_Ghalf_correct  = Ghalf_correct_draw
    #       override_VhalfInv       = VhalfInv_draw
    #       override_sigmasq_tilde  = sigmasq_draw
    #       override_trace_XUGX     = trace_draw
    #
    #     The closure will then compute:
    #       L_post = (1/sd_y) * U_draw %*% Ghalf_correct_draw
    #     and draw:
    #       beta ~ N(B_raw_fixed, sigmasq * L_post L_post^T)
    #     (or use ESS when inequalities are enforced).

    ## Use the rebuilt dispersion for the InvGamma draw
    refit_sigmasq <- cov_rebuild$sigmasq_tilde
    if(is.na(refit_sigmasq) || !is.finite(refit_sigmasq)){
      refit_sigmasq <- new_sigmasq_tilde
    }

    override_args <- list(
      new_sigmasq_tilde            = refit_sigmasq,
      theta_1                      = theta_1,
      theta_2                      = theta_2,
      posterior_predictive_draw    = posterior_predictive_draw,
      draw_dispersion              = draw_dispersion,
      include_posterior_predictive = include_posterior_predictive,
      num_draws                    = 1,
      enforce_qp_constraints       = enforce_qp_constraints,
      override_B_raw               = B_raw_fixed,
      override_U                   = cov_rebuild$U,
      override_Ghalf_correct       = cov_rebuild$Ghalf_correct,
      override_VhalfInv            = cov_rebuild$VhalfInv,
      override_sigmasq_tilde       = refit_sigmasq,
      override_trace_XUGX          = cov_rebuild$trace_XUGX
    )
    if(!is.null(new_predictors)){
      override_args$new_predictors <- new_predictors
    }

    one_draw <- do.call(object$generate_posterior,
                        c(override_args, list(...)))
    one_draw$post_draw_correlation_params <- phi_draw
    return(one_draw)
  }


  ## 8. Execute draws and repack into the standard return shape.

  results <- lapply(seq_len(num_draws), function(m) .one_corr_draw())

  ## Single-draw return: flat list
  if(num_draws == 1){
    return(results[[1]])
  }

  ## Multi-draw return: repack per-draw results into the usual shape
  #  that matches the generate_posterior() multi-draw contract.
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


#' Plot Method for lgspline Objects
#'
#' Wrapper for the internal lgspline plot function. Produces a 1D line plot
#' (base R) or interactive 3D surface plot (plotly) depending on the number
#' of predictors, with optional formula overlays per partition. When plotting
#' a subset of variables via \code{vars}, non-plotted predictors are
#' automatically set to zero (or a user-specified value via
#' \code{fixed_values}).
#'
#' @details
#' Partition boundaries are indicated by color changes. For 1D models,
#' observation points can be overlaid. For 2D models, plotly is used.
#'
#' When using \code{vars} to plot a subset of predictors, the non-plotted
#' predictors are automatically set to zero. This can be overridden by
#' passing a named list to \code{fixed_values} (e.g.,
#' \code{fixed_values = list(Height = 75)}). The automatic zeroing replaces
#' the previous behavior where the user had to manually construct
#' \code{new_predictors} with non-plotted variables set to fixed values.
#'
#' When \code{se.fit = TRUE}, pointwise confidence bands are drawn around
#' the fitted function. These are Wald-type intervals constructed on the
#' link scale and back-transformed to the response scale, using
#' \code{cv} as the critical value to multiply se.fit by (default 1 for actual se).
#'
#' The function extracts predictor positions from linear expansion terms.
#' If linear terms are excluded (e.g., via \code{exclude_these_expansions}),
#' plotting will fail. As a workaround, constrain those terms to zero via
#' \code{constraint_vectors} / \code{constraint_values} so they remain in
#' the expansion but are zeroed out.
#'
#' @param x A fitted lgspline model object.
#' @param show_formulas Logical; display partition-level polynomial formulas.
#'        Default FALSE.
#' @param include_all_terms_in_formulas Logical; when \code{show_formulas = TRUE}
#'        and plotting only a subset of predictors via \code{vars}, include all
#'        fitted terms in the displayed formulas rather than only the terms
#'        involving the plotted predictor(s). Default FALSE retains the current
#'        marginal-only formula display.
#' @param digits Integer; decimal places for formula coefficients. Default 4.
#' @param legend_pos Character; legend position for 1D plots. Default
#'        \code{"topright"}.
#' @param custom_response_lab Character; response axis label. Default
#'        \code{"y"}.
#' @param custom_predictor_lab Character; predictor axis label (1D). Default
#'        NULL uses the column name.
#' @param custom_predictor_lab1 Character; first predictor axis label (2D).
#'        Default NULL.
#' @param custom_predictor_lab2 Character; second predictor axis label (2D).
#'        Default NULL.
#' @param custom_formula_lab Character; fitted response label on the link
#'        scale. Default NULL.
#' @param custom_title Character; plot title. Default \code{"Fitted Function"}.
#' @param text_size_formula Numeric; formula text size. Passed to \code{cex}
#'        (1D) or hover font size (2D). Default NULL (0.8 for 1D, 8 for 2D).
#' @param legend_args List; additional arguments passed to \code{legend()}
#'        (1D only).
#' @param new_predictors Matrix; optional predictor values for prediction.
#'        Default NULL. When \code{vars} is specified and
#'        \code{new_predictors} is NULL, a grid is automatically generated
#'        with non-plotted variables set to zero (or values from
#'        \code{fixed_values}).
#' @param xlim Numeric vector; x-axis limits (1D only). Default NULL.
#' @param ylim Numeric vector; y-axis limits (1D only). Default NULL.
#' @param color_function Function; returns K+1 colors, one per partition.
#'        Default NULL uses \code{grDevices::rainbow(K+1)} for 1D and a
#'        Spectral palette for 2D.
#' @param add Logical; add to an existing plot (1D only). Default FALSE.
#' @param vars Numeric or character vector; predictor indices or names to
#'        plot. Default \code{c()} plots all.
#' @param legend_order Numeric; re-ordered partition indices for the legend.
#' @param se.fit Logical; if TRUE, plot pointwise confidence bands. Default
#'        FALSE.
#' @param cv Numeric; critical value for confidence bands. Default 1.
#' @param band_col Character; color for confidence band fill. Default
#'        \code{"grey80"}.
#' @param band_border Character or NA; border color for confidence band
#'        polygon. Default NA (no border).
#' @param fixed_values Named list; fixed values for non-plotted predictors
#'        when \code{vars} is specified. Names should match predictor names.
#'        Default NULL sets non-plotted predictors to zero.
#' @param n_grid Integer; number of grid points for automatic grid
#'        generation when \code{vars} is specified and
#'        \code{new_predictors} is NULL. Default 200.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}
#'        (1D) or \code{\link[plotly]{plot_ly}} (2D).
#'
#' @return For 1D models: invisibly returns NULL (base R plot drawn to
#'   device). For 2D models: returns a plotly object.
#'
#' @examples
#'
#' set.seed(1234)
#' t_data <- runif(1000, -10, 10)
#' y_data <- 2*sin(t_data) + -0.06*t_data^2 + rnorm(length(t_data))
#' model_fit <- lgspline(t_data, y_data, K = 9)
#'
#' ## Basic plot
#' plot(model_fit)
#'
#' ## Plot with confidence bands
#' plot(model_fit,
#'      se.fit = TRUE,
#'      cv = 1.96,
#'      custom_title = 'Fitted Function with 95% CI')
#'
#' ## Multi-predictor: automatically zeros non-plotted variables
#' # plot(model_fit_2d, vars = 'x1', se.fit = TRUE)
#'
#' @seealso
#' \code{\link{lgspline}},
#' \code{\link[graphics]{plot}},
#' \code{\link[plotly]{plot_ly}}
#' @export
plot.lgspline <- function(x,
                          show_formulas = FALSE,
                          include_all_terms_in_formulas = FALSE,
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
                          legend_order = NULL,
                          se.fit = FALSE,
                          cv = 1,
                          band_col = "grey80",
                          band_border = NA,
                          fixed_values = NULL,
                          n_grid = 200,
                          ...) {

  ## When vars is supplied without new_predictors, build a plotting grid
  #  from the training predictors and hold the other variables fixed.

  og_cols <- x$.fit_call_args$og_cols
  q_pred <- x$q

  ## Only auto-generate new_predictors when vars is specified and
  ## new_predictors is NULL and we have >1 predictor.
  if(length(vars) > 0 && is.null(new_predictors) && q_pred > 1){

    ## Resolve vars to numeric indices
    if(is.character(vars)){
      if(is.null(og_cols)){
        stop('\n\t Character vars requires named predictor columns. ',
             'Use numeric column indices instead.\n')
      }
      vars_idx <- match(vars, og_cols)
      if(any(is.na(vars_idx))){
        stop('\n\t vars "', paste(vars[is.na(vars_idx)], collapse='", "'),
             '" not found in predictor names.\n')
      }
    } else {
      vars_idx <- as.integer(vars)
    }

    ## Recover training predictors from the design matrices.
    #  x$X is a list of unstandardized design matrices per partition.
    #  The linear columns (power1_cols and nonspline_cols) contain the
    #  original predictor values. We extract these in partition order
    #  and reorder to observation order.
    all_linear_cols <- sort(c(x$power1_cols, x$nonspline_cols))
    training_preds <- matrix(NA, x$N, q_pred)
    for(k in 1:(x$K + 1)){
      if(nrow(x$X[[k]]) == 0) next
      obs_idx <- x$order_list[[k]]
      ## Extract the q linear columns from the unstandardized expansion
      training_preds[obs_idx, ] <- x$X[[k]][, all_linear_cols, drop = FALSE]
    }

    ## Build the grid
    if(length(vars_idx) == 1){
      ## 1D plot: vary the single plotted variable
      var_range <- range(training_preds[, vars_idx], na.rm = TRUE)
      grid_vals <- seq(var_range[1], var_range[2], length.out = n_grid)

      ## Initialize all columns to their fixed values (default 0)
      new_pred_mat <- matrix(0, nrow = n_grid, ncol = q_pred)
      if(!is.null(og_cols)){
        colnames(new_pred_mat) <- og_cols
      }

      ## Set plotted variable
      new_pred_mat[, vars_idx] <- grid_vals

      ## Override non-plotted variables with fixed_values if provided
      if(!is.null(fixed_values)){
        for(nm in names(fixed_values)){
          if(!is.null(og_cols)){
            fix_idx <- match(nm, og_cols)
          } else {
            fix_idx <- as.integer(nm)
          }
          if(!is.na(fix_idx) && !(fix_idx %in% vars_idx)){
            new_pred_mat[, fix_idx] <- fixed_values[[nm]]
          }
        }
      }

    } else if(length(vars_idx) == 2){
      ## 2D plot: vary both plotted variables
      var1_range <- range(training_preds[, vars_idx[1]], na.rm = TRUE)
      var2_range <- range(training_preds[, vars_idx[2]], na.rm = TRUE)
      n_side <- ceiling(sqrt(n_grid))
      grid_expand <- expand.grid(
        seq(var1_range[1], var1_range[2], length.out = n_side),
        seq(var2_range[1], var2_range[2], length.out = n_side)
      )

      new_pred_mat <- matrix(0, nrow = nrow(grid_expand), ncol = q_pred)
      if(!is.null(og_cols)){
        colnames(new_pred_mat) <- og_cols
      }
      new_pred_mat[, vars_idx[1]] <- grid_expand[, 1]
      new_pred_mat[, vars_idx[2]] <- grid_expand[, 2]

      ## Override non-plotted variables with fixed_values if provided
      if(!is.null(fixed_values)){
        for(nm in names(fixed_values)){
          if(!is.null(og_cols)){
            fix_idx <- match(nm, og_cols)
          } else {
            fix_idx <- as.integer(nm)
          }
          if(!is.na(fix_idx) && !(fix_idx %in% vars_idx)){
            new_pred_mat[, fix_idx] <- fixed_values[[nm]]
          }
        }
      }
    }

    new_predictors <- new_pred_mat
  }

  ## For 1D plots, se.fit draws the confidence band first and then
  #  overlays the fitted curve.

  internal_plot_func <- x$plot
  if (!is.null(internal_plot_func) && is.function(internal_plot_func)) {

    ## Determine if this is a 1D plot scenario
    is_1d <- (x$q == 1) || (length(vars) == 1)

    if(se.fit && is_1d){
      ## Get predictions with standard errors
      #  First, determine which predictors to use for the SE computation
      if(!is.null(new_predictors)){
        se_predictors <- new_predictors
      } else {
        ## Use the training data; extract from X
        all_linear_cols <- sort(c(x$power1_cols, x$nonspline_cols))
        se_predictors <- matrix(NA, x$N, x$q)
        for(k in 1:(x$K + 1)){
          if(nrow(x$X[[k]]) == 0) next
          obs_idx <- x$order_list[[k]]
          se_predictors[obs_idx, ] <- x$X[[k]][, all_linear_cols, drop = FALSE]
        }
        if(!is.null(og_cols)){
          colnames(se_predictors) <- og_cols
        }
      }

      ## Get predictions with SE
      se_result <- x$predict(new_predictors = se_predictors,
                             se.fit = TRUE,
                             cv = cv)

      ## Order the fitted values by the plotted x-variable before drawing
      #  the confidence band polygon.
      if(length(vars) == 1){
        if(is.character(vars)){
          var_col <- match(vars, og_cols)
        } else {
          var_col <- as.integer(vars)
        }
      } else {
        var_col <- 1
      }
      x_vals <- se_predictors[, var_col]
      sort_order <- order(x_vals)
      x_sorted <- x_vals[sort_order]
      lower_sorted <- se_result$lower[sort_order]
      upper_sorted <- se_result$upper[sort_order]

      ## Set up the plot region first (if not adding to existing)
      if(!add){
        ## Set up a fresh plotting region, then draw the band underneath
        #  the fitted curves.
        ## Determine y-limits including the bands
        all_y <- c(se_result$lower, se_result$upper, se_result$fit, x$y)
        if(is.null(ylim)){
          ylim_use <- range(all_y, na.rm = TRUE)
        } else {
          ylim_use <- ylim
        }

        ## Determine axis label
        if(is.null(custom_predictor_lab)){
          if(!is.null(og_cols) && length(vars) == 1){
            if(is.character(vars)){
              custom_predictor_lab_use <- vars
            } else {
              custom_predictor_lab_use <- og_cols[var_col]
            }
          } else if(!is.null(og_cols)){
            custom_predictor_lab_use <- og_cols[1]
          } else {
            custom_predictor_lab_use <- "x"
          }
        } else {
          custom_predictor_lab_use <- custom_predictor_lab
        }

        ## Draw empty plot frame
        graphics::plot(range(x_sorted), ylim_use,
                       type = 'n',
                       xlab = custom_predictor_lab_use,
                       ylab = custom_response_lab,
                       main = custom_title,
                       xlim = xlim,
                       ...)

        ## Draw confidence band polygon
        graphics::polygon(
          x = c(x_sorted, rev(x_sorted)),
          y = c(lower_sorted, rev(upper_sorted)),
          col = band_col,
          border = band_border
        )

        ## Now overlay the fitted curves using add = TRUE
        plot_result <- internal_plot_func(
          model_fit_in = x,
          show_formulas = show_formulas,
          include_all_terms_in_formulas = include_all_terms_in_formulas,
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
          ylim = ylim_use,
          color_function = color_function,
          add = TRUE,  # overlay on the band
          vars = vars,
          legend_order = legend_order,
          ...)
      } else {
        ## When adding to an existing plot, just layer the band and curves.
        graphics::polygon(
          x = c(x_sorted, rev(x_sorted)),
          y = c(lower_sorted, rev(upper_sorted)),
          col = band_col,
          border = band_border
        )

        plot_result <- internal_plot_func(
          model_fit_in = x,
          show_formulas = show_formulas,
          include_all_terms_in_formulas = include_all_terms_in_formulas,
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
          add = TRUE,
          vars = vars,
          legend_order = legend_order,
          ...)
      }

      invisible(NULL)

    } else {
      ## Non-SE path: delegate directly to the internal plot function
      plot_result <- internal_plot_func(
        model_fit_in = x,
        show_formulas = show_formulas,
        include_all_terms_in_formulas = include_all_terms_in_formulas,
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
        legend_order = legend_order,
        ...)

      if(inherits(plot_result, "plotly")) {
        return(plot_result)
      } else {
        invisible(NULL)
      }
    }
  } else {
    stop("Internal plot method not found or not a function.")
  }
}


#' Predict Method for lgspline Objects
#'
#' Generates predictions, derivatives, and basis expansions from a fitted
#' lgspline model. Wrapper for the internal predict closure stored in the
#' object.
#'
#' @details
#' \code{new_predictors} takes priority over \code{newdata} when both are
#' supplied. When both are NULL, the training data is used.
#'
#' Fitted values are also accessible directly as \code{model_fit$ytilde} or
#' via \code{model_fit$predict()}.
#'
#' The parallel processing feature is experimental.
#'
#' Additional arguments passed through \code{...} include \code{se.fit}
#' and \code{cv} for pointwise interval summaries.
#'
#' Predictor input should use the original predictor columns. Named extra
#' columns are dropped when they can be identified as irrelevant to the
#' fitted expansions.
#'
#' @param object A fitted lgspline model object.
#' @param newdata Matrix or data.frame; new predictor values. Default NULL.
#' @param parallel Logical; use parallel processing (experimental).
#'        Default FALSE.
#' @param cl Cluster object for parallel processing. Default NULL.
#' @param chunk_size Integer; chunk size for parallel processing. Default NULL.
#' @param num_chunks Integer; number of chunks. Default NULL.
#' @param rem_chunks Integer; remainder chunks. Default NULL.
#' @param B_predict List; per-partition coefficient list for prediction, e.g.
#'        from \code{\link{generate_posterior}}. Default NULL uses
#'        \code{object$B}.
#' @param take_first_derivatives Logical; compute first derivatives.
#'        Default FALSE.
#' @param take_second_derivatives Logical; compute second derivatives.
#'        Default FALSE.
#' @param expansions_only Logical; return basis expansion matrix only.
#'        Default FALSE.
#' @param new_predictors Matrix or data.frame; overrides \code{newdata}.
#' @param ... Additional arguments passed to the internal predict method.
#'
#' @return A numeric vector of predictions, or a list when derivatives or
#'   interval summaries are requested:
#' \describe{
#'   \item{preds}{Numeric vector of predictions when derivatives are requested.}
#'   \item{fit}{Numeric vector of predictions when \code{se.fit = TRUE} and no
#'         derivatives are requested.}
#'   \item{first_deriv}{Numeric vector or named list of first derivatives (if
#'         requested).}
#'   \item{second_deriv}{Numeric vector or named list of second derivatives (if
#'         requested).}
#'   \item{se.fit}{Pointwise standard errors on the link scale (if requested).}
#'   \item{lower}{Pointwise lower interval bound (if requested).}
#'   \item{upper}{Pointwise upper interval bound (if requested).}
#'   \item{cv}{Critical value returned when \code{se.fit = TRUE} without
#'         derivative requests.}
#' }
#' If \code{expansions_only = TRUE}, returns a list of basis expansions.
#'
#' @examples
#'
#' set.seed(1234)
#' t <- runif(1000, -10, 10)
#' y <- 2*sin(t) + -0.06*t^2 + rnorm(length(t))
#' model_fit <- lgspline(t, y)
#'
#' newdata <- matrix(sort(rnorm(10000)), ncol = 1)
#' preds <- predict(model_fit, newdata)
#'
#' deriv1_res <- predict(model_fit, newdata, take_first_derivatives = TRUE)
#' deriv2_res <- predict(model_fit, newdata, take_second_derivatives = TRUE)
#'
#' oldpar <- par(no.readonly = TRUE)
#' layout(matrix(c(1,1,2,2,3,3), byrow = TRUE, ncol = 2))
#'
#' plot(newdata[,1], preds, main = 'Fitted Function',
#'      xlab = 't', ylab = "f(t)", type = 'l')
#' plot(newdata[,1], deriv1_res$first_deriv, main = 'First Derivative',
#'      xlab = 't', ylab = "f'(t)", type = 'l')
#' plot(newdata[,1], deriv2_res$second_deriv, main = 'Second Derivative',
#'      xlab = 't', ylab = "f''(t)", type = 'l')
#'
#' par(oldpar)
#'
#' @seealso \code{\link{lgspline}}, \code{\link{plot.lgspline}}
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

  ## Default to the fitted coefficients unless B_predict is supplied.
  B_predict_val <- if(!is.null(B_predict)) B_predict else object$B

  ## new_predictors takes priority over newdata
  if(!is.null(new_predictors)){
    predictors_val <- new_predictors
  } else if(!is.null(newdata)){
    predictors_val <- newdata
  } else {
    predictors_val <- NULL
  }

  ## Unwrap nested data.frame: data.frame(new_predictors = data.frame(...))
  #  produces a single-column df whose sole column is itself a df.
  if(inherits(predictors_val, "data.frame") && ncol(predictors_val) == 1 &&
     inherits(predictors_val[[1]], "data.frame")){
    predictors_val <- predictors_val[[1]]
  }

  ## Leave data.frames with non-numeric columns as-is so the internal
  #  predict function can handle factor encoding before matrix coercion.
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

  ## Omit new_predictors for in-sample predictions; otherwise pass the
  #  requested predictor matrix through explicitly.
  if(is.null(predictors_val)){
    internal_predict_func(
      parallel                = parallel,
      cl                      = cl,
      chunk_size              = chunk_size,
      num_chunks              = num_chunks,
      rem_chunks              = rem_chunks,
      B_predict               = B_predict_val,
      take_first_derivatives  = take_first_derivatives,
      take_second_derivatives = take_second_derivatives,
      expansions_only         = expansions_only,
      ...
    )
  } else {
    internal_predict_func(
      new_predictors          = predictors_val,
      parallel                = parallel,
      cl                      = cl,
      chunk_size              = chunk_size,
      num_chunks              = num_chunks,
      rem_chunks              = rem_chunks,
      B_predict               = B_predict_val,
      take_first_derivatives  = take_first_derivatives,
      take_second_derivatives = take_second_derivatives,
      expansions_only         = expansions_only,
      ...
    )
  }
}


#' Extract Coefficients from a Fitted lgspline
#'
#' Returns the per-partition polynomial coefficient lists from a fitted
#' lgspline model.
#'
#' @param object A fitted lgspline model object.
#' @param ... Not used.
#'
#' @details
#' Coefficient names reflect the polynomial expansion terms, e.g.:
#' \itemize{
#'   \item intercept
#'   \item v: linear term for predictor v
#'   \item v_^2: quadratic term
#'   \item v^3: cubic term
#'   \item _v_x_w_: two-way interaction
#' }
#' Column/variable names replace numeric indices when available.
#'
#' To get all coefficients as a single matrix:
#' \code{Reduce('cbind', coef(model_fit))}.
#'
#' @return A list of per-partition coefficient vectors. Returns NULL with a
#'   warning if \code{object$B} is not found.
#'
#' @examples
#'
#' set.seed(1234)
#' t <- runif(1000, -10, 10)
#' y <- 2*sin(t) + -0.06*t^2 + rnorm(length(t))
#' model_fit <- lgspline(t, y)
#'
#' coefficients <- coef(model_fit)
#' print(coefficients[[1]])
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


#' Univariate Wald Tests and Confidence Intervals for lgspline Coefficients
#'
#' Computes per-coefficient Wald tests and confidence intervals from a fitted
#' lgspline. For Gaussian identity-link models, t-statistics and t-intervals
#' are used; otherwise z-statistics.
#'
#' @param object A fitted lgspline object. Must have been fit with
#'        \code{return_varcovmat = TRUE}.
#' @param scale_vcovmat_by Numeric; scaling factor for the variance-covariance
#'        matrix. Default 1.
#' @param cv Numeric; critical value for confidence intervals. If missing,
#'        defaults to \code{object$critical_value} or \code{qnorm(0.975)}.
#' @param ... Additional arguments passed to the internal \code{wald_univariate}
#'        method.
#'
#' @return An object of class \code{"wald_lgspline"}, a list with:
#' \describe{
#'   \item{coefficients}{Matrix with columns: Estimate, Std. Error,
#'         t value or z value, Pr(>|t|) or Pr(>|z|), CI LB, CI UB.}
#'   \item{critical_value}{Critical value used.}
#'   \item{family}{GLM family from the fitted model.}
#'   \item{N}{Number of observations.}
#'   \item{trace_XUGX}{Effective df trace term.}
#'   \item{statistic_name}{"t value" or "z value".}
#'   \item{p_value_name}{"Pr(>|t|)" or "Pr(>|z|)".}
#'   \item{df.residual}{Residual degrees of freedom when supplied by the
#'         internal Wald method.}
#' }
#' Print, summary, and plot methods are available; see
#' \code{\link{print.wald_lgspline}}, \code{\link{summary.wald_lgspline}},
#' \code{\link{plot.wald_lgspline}}.
#'
#' @examples
#'
#' set.seed(1234)
#' t <- runif(1000, -10, 10)
#' y <- 2*sin(t) + -0.06*t^2 + rnorm(length(t))
#' model_fit <- lgspline(t, y, return_varcovmat = TRUE)
#'
#' wald_default <- wald_univariate(model_fit)
#' print(wald_default)
#'
#' ## t-distribution critical value
#' eff_df <- model_fit$N - model_fit$trace_XUGX
#' wald_t <- wald_univariate(model_fit, cv = qt(0.975, eff_df))
#' print(wald_t)
#'
#' coef_table <- wald_default$coefficients
#' plot(wald_default)
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

  ## t-tests for Gaussian identity; z-tests otherwise
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

  ## If internal method already returns a wald_lgspline object, rebuild the
  #  wrapper return for a consistent structure.
  if (inherits(res, "wald_lgspline")) {
    out <- list(
      coefficients   = res$coefficients,
      critical_value = cv,
      family         = object$family,
      N              = object$N,
      trace_XUGX     = object$trace_XUGX,
      statistic_name = stat_name,
      p_value_name   = p_name
    )
    class(out) <- "wald_lgspline"
    return(out)
  }

  ## Normalize to a labelled matrix.
  #  The internal method may return a list of vectors (old-style interface).
  coef_mat <- NULL
  if (is.list(res) && !is.data.frame(res) && !is.matrix(res)) {
    coef_mat <- tryCatch({
      mat <- cbind(
        Estimate      = res$est,
        `Std. Error`  = res$se,
        Statistic     = res$stat,
        `CI LB`       = res$interval_lb,
        `CI UB`       = res$interval_ub,
        p.value       = res$pval
      )
      colnames(mat)[3] <- stat_name
      colnames(mat)[6] <- p_name
      mat <- mat[, c("Estimate", "Std. Error", stat_name,
                     p_name, "CI LB", "CI UB"), drop = FALSE]
      mat
    }, error = function(e) NULL)
    if (is.null(coef_mat)) {
      coef_mat <- tryCatch(Reduce('cbind', res), error = function(e) NULL)
    }
  } else if (is.matrix(res) || is.data.frame(res)) {
    coef_mat <- as.matrix(res)
  }

  if (is.null(coef_mat)) {
    warning("Could not normalize wald_univariate output to matrix.",
            call. = FALSE)
    coef_mat <- cbind(Estimate = unlist(object$B))
  }

  ## Assign row names from coefficient names when missing
  if (is.null(rownames(coef_mat)) && !is.null(object$B) && is.list(object$B) &&
      length(unlist(object$B)) == nrow(coef_mat)) {
    rn <- tryCatch({
      unlist(lapply(seq_along(object$B), function(k) {
        part_names <- names(object$B[[k]])
        if (is.null(part_names))
          part_names <- paste0("Term", seq_len(length(object$B[[k]])))
        paste0("partition", k, "_", part_names)
      }))
    }, error = function(e) NULL)
    if (!is.null(rn) && length(rn) == nrow(coef_mat)) {
      rownames(coef_mat) <- rn
    }
  }

  out <- list(
    coefficients   = coef_mat,
    critical_value = cv,
    family         = object$family,
    N              = object$N,
    trace_XUGX     = object$trace_XUGX,
    statistic_name = stat_name,
    p_value_name   = p_name,
    df.residual    = res$df.residual
  )
  class(out) <- "wald_lgspline"
  return(out)
}


#' Print Method for wald_lgspline Objects
#'
#' Prints the coefficient table using \code{\link[stats]{printCoefmat}} with
#' significance stars.
#'
#' @param x A \code{"wald_lgspline"} object from \code{\link{wald_univariate}}.
#' @param digits Number of significant digits.
#' @param signif.stars Logical; show significance stars.
#' @param ... Additional arguments passed to \code{\link[stats]{printCoefmat}}.
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link{wald_univariate}}, \code{\link[stats]{printCoefmat}}
#' @method print wald_lgspline
#' @export
print.wald_lgspline <- function(x, digits = max(3, getOption("digits") - 3),
                                signif.stars = getOption("show.signif.stars"),
                                ...) {
  cat("\nWald Inference for lgspline Coefficients\n")
  cat("Family:", paste0(x$family)[1],
      " Link:", paste0(x$family)[2], "\n")
  if(!is.null(x$N)) cat("N =", x$N, "\n\n")

  mat <- x$coefficients
  pval_col <- grep("^Pr\\(", colnames(mat))
  if(length(pval_col) > 0){
    ## printCoefmat expects p-value column last
    other_cols <- setdiff(seq_len(ncol(mat)), pval_col)
    print_mat <- mat[, c(other_cols, pval_col), drop = FALSE]
    pval_col_new <- ncol(print_mat)
    tst_col <- grep("value$", colnames(print_mat))
    tst_col <- tst_col[!tst_col %in% pval_col_new]
    stats::printCoefmat(print_mat,
                        digits = digits,
                        signif.stars = signif.stars,
                        cs.ind = 1:2,
                        tst.ind = if(length(tst_col) > 0) tst_col[1] else NULL,
                        P.values = TRUE,
                        has.Pvalue = TRUE,
                        df = if(!is.null(x$df.residual) &&
                                length(x$df.residual) > 0 &&
                                is.finite(x$df.residual)) x$df.residual else NULL,
                        ...)
  } else {
    print(mat, ...)
  }
  cat("\nCritical value:", x$critical_value, "\n")
  if(!is.null(x$df.residual) && is.finite(x$df.residual)){
    cat("Residual degrees of freedom:", round(x$df.residual, 2), "\n")
  }
  invisible(x)
}


#' Summary Method for wald_lgspline Objects
#'
#' Prints a header with model info then delegates to
#' \code{\link{print.wald_lgspline}}.
#'
#' @param object A \code{"wald_lgspline"} object.
#' @param ... Passed to \code{\link{print.wald_lgspline}}.
#'
#' @return Invisibly returns \code{object}.
#'
#' @seealso \code{\link{wald_univariate}}, \code{\link{print.wald_lgspline}}
#' @method summary wald_lgspline
#' @export
summary.wald_lgspline <- function(object, ...) {
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
#' Forest-style plot of coefficient estimates with confidence intervals.
#'
#' @param x A \code{"wald_lgspline"} object.
#' @param parm Integer vector of coefficient indices or character vector of
#'        names to plot. Default NULL plots all.
#' @param which Integer vector of coefficient indices to plot (alternative to
#'        \code{parm}). Default NULL.
#' @param main Plot title. Default \code{"Coefficient Estimates and CIs"}.
#' @param xlab x-axis label. Default \code{"Estimate"}.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}.
#'
#' @return Invisibly returns NULL.
#'
#' @seealso \code{\link{wald_univariate}}, \code{\link{confint.lgspline}}
#' @method plot wald_lgspline
#' @export
plot.wald_lgspline <- function(x,
                               parm = NULL,
                               which = NULL,
                               main = "Coefficient Estimates and CIs",
                               xlab = "Estimate",
                               ...) {
  mat <- x$coefficients

  ## Support both parm (named/indexed) and which (index-only)
  if(!is.null(parm)){
    if(is.character(parm)){
      mat <- mat[rownames(mat) %in% parm, , drop = FALSE]
    } else {
      mat <- mat[parm, , drop = FALSE]
    }
  } else if(!is.null(which)){
    mat <- mat[which, , drop = FALSE]
  }

  est_col <- which(colnames(mat) == "Estimate")
  lb_col  <- grep("CI LB|^Lower", colnames(mat))
  ub_col  <- grep("CI UB|^Upper", colnames(mat))

  if(length(est_col) == 0 || length(lb_col) == 0 || length(ub_col) == 0){
    ## Try x$est path for old-style objects
    if(!is.null(x$est)){
      idx_sel <- if(!is.null(which)) which else seq_along(x$est)
      est <- x$est[idx_sel]
      lb  <- x$interval_lb[idx_sel]
      ub  <- x$interval_ub[idx_sel]
      nms <- names(est)
      if(is.null(nms)) nms <- paste0("beta[", idx_sel, "]")
      n <- length(est)
      xlim <- range(c(lb, ub), na.rm = TRUE)
      graphics::plot(est, n:1, xlim = xlim, yaxt = "n",
                     ylab = "", xlab = xlab, main = main,
                     pch = 16, ...)
      graphics::segments(lb, n:1, ub, n:1)
      graphics::axis(2, at = n:1, labels = nms, las = 1, cex.axis = 0.7)
      graphics::abline(v = 0, lty = 2, col = "grey50")
      return(invisible(x))
    }
    warning("Cannot produce forest plot: missing Estimate, CI LB, or CI UB columns.",
            call. = FALSE)
    print(mat)
    return(invisible(NULL))
  }

  ests <- mat[, est_col]
  lbs  <- mat[, lb_col]
  ubs  <- mat[, ub_col]
  n_coef <- length(ests)
  idx <- seq_len(n_coef)

  labs <- rownames(mat)
  if(is.null(labs)) labs <- paste0("Coef ", idx)

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


#' Extract Coefficients from a wald_lgspline Object
#'
#' @param object A \code{"wald_lgspline"} object.
#' @param ... Not used.
#'
#' @return Named numeric vector of coefficient estimates, or NULL if no
#'   estimate column is available.
#' @export
coef.wald_lgspline <- function(object, ...){
  if(!is.null(object$est)) return(object$est)
  if(!is.null(object$coefficients)){
    est_col <- which(colnames(object$coefficients) == "Estimate")
    if(length(est_col) > 0) return(object$coefficients[, est_col])
  }
  NULL
}


#' Extract Confidence Intervals from a wald_lgspline Object
#'
#' @param object A \code{"wald_lgspline"} object.
#' @param parm Parameter specification (ignored; all returned).
#' @param level Confidence level (ignored; uses the object's critical value).
#' @param ... Not used.
#'
#' @return Matrix with columns lower and upper, or NULL if confidence
#'   limits are not available.
#' @export
confint.wald_lgspline <- function(object, parm = NULL, level = NULL, ...){
  if(!is.null(object$interval_lb) && !is.null(object$interval_ub)){
    return(cbind(lower = object$interval_lb, upper = object$interval_ub))
  }
  if(!is.null(object$coefficients)){
    lb_col <- grep("CI LB|^Lower", colnames(object$coefficients))
    ub_col <- grep("CI UB|^Upper", colnames(object$coefficients))
    if(length(lb_col) > 0 && length(ub_col) > 0){
      return(cbind(lower = object$coefficients[, lb_col],
                   upper = object$coefficients[, ub_col]))
    }
  }
  NULL
}


#' Confidence Intervals for lgspline Coefficients
#'
#' Wald-based confidence intervals for regression coefficients and, when
#' available, correlation parameters (on the working scale).
#'
#' @details
#' For Gaussian identity-link models, t-distribution quantiles are used with
#' effective degrees of freedom
#' \eqn{N - \mathrm{trace}(\mathbf{XUGX}^\top)}.
#' All other families use normal quantiles.
#'
#' Correlation parameter intervals (if \code{VhalfInv_params_estimates} and
#' \code{VhalfInv_params_vcov} are present) are computed on the unbounded
#' working scale via a Wald interval.
#' @param object A fitted lgspline object with \code{return_varcovmat = TRUE}.
#' @param parm Optional vector of parameter indices or names. Default returns
#'        all regression parameters; working-scale correlation parameters are
#'        appended when available.
#' @param level Confidence level. Default 0.95.
#' @param ... Additional arguments passed to \code{\link{wald_univariate}}.
#'
#' @return A matrix with columns giving lower and upper confidence limits,
#'   named e.g. \code{2.5 \%} and \code{97.5 \%} for 95\% intervals.
#'   When available, rows for working-scale correlation parameters are
#'   appended after the regression coefficients.
#'
#' @method confint lgspline
#' @export
confint.lgspline <- function(object, parm, level = 0.95, ...) {
  if(is.null(object$varcovmat)){
    stop("confint requires return_varcovmat = TRUE during model fitting")
  }

  alpha <- 1 - level

  if(object$family$family == "gaussian" && object$family$link == "identity"){
    eff_df <- object$N - object$trace_XUGX
    if(!is.null(eff_df) && !is.na(eff_df) && is.finite(eff_df) && eff_df > 0){
      cv <- stats::qt(1 - alpha / 2, df = eff_df)
    } else {
      cv <- stats::qnorm(1 - alpha / 2)
      warning("Effective df non-positive; using normal quantiles.", call. = FALSE)
    }
  } else {
    cv <- stats::qnorm(1 - alpha / 2)
  }

  wald_res <- wald_univariate(object, cv = cv, ...)
  coef_mat <- wald_res$coefficients

  ## CI column detection
  lb_col <- grep("CI LB|Lower", colnames(coef_mat), ignore.case = TRUE)
  ub_col <- grep("CI UB|Upper", colnames(coef_mat), ignore.case = TRUE)
  if(length(lb_col) == 0 || length(ub_col) == 0){
    stop("Could not extract CI columns from wald_univariate output.")
  }
  lb_col <- lb_col[1]
  ub_col <- ub_col[1]

  ci <- coef_mat[, c(lb_col, ub_col), drop = FALSE]
  pct <- format(100 * c(alpha / 2, 1 - alpha / 2), digits = 3, trim = TRUE)
  colnames(ci) <- paste0(pct, " %")

  if(!missing(parm)){
    if(is.character(parm)){
      ci <- ci[rownames(ci) %in% parm, , drop = FALSE]
    } else {
      ci <- ci[parm, , drop = FALSE]
    }
  }

  ## Append correlation parameter intervals on the working scale
  if(!is.null(object$VhalfInv_params_estimates) &&
     !is.null(object$VhalfInv_params_vcov) &&
     all(!is.na(object$VhalfInv_params_estimates)) &&
     all(!is.na(object$VhalfInv_params_vcov))){

    par_est <- object$VhalfInv_params_estimates
    vcov_mat <- object$VhalfInv_params_vcov
    se <- sqrt(diag(vcov_mat))

    for(i in seq_along(par_est)){
      work_ci <- par_est[i] + c(-cv, cv) * se[i]
      ci_corr <- matrix(work_ci, nrow = 1)
      colnames(ci_corr) <- colnames(ci)
      rownames(ci_corr) <- paste0("Correlation parameter ", i)
      ci <- rbind(ci, ci_corr)
    }
  }

  return(ci)
}


#' Extract Log-Likelihood from a Fitted lgspline
#'
#' Returns the log-likelihood as a \code{"logLik"} object for use with
#' \code{\link[stats]{AIC}}, \code{\link[stats]{BIC}}, and other model
#' comparison tools.
#'
#' @details
#' \strong{Gaussian identity, no correlation.}
#' \deqn{
#'   \ell = -\frac{N}{2}\log(2\pi\tilde{\sigma}^2) -
#'   \frac{1}{2\tilde{\sigma}^2}\sum_{i}(y_i - \hat{y}_i)^2
#' }
#'
#' \strong{Gaussian identity, with correlation.}
#' GLS log-likelihood:
#' \deqn{
#'   \ell = -\frac{N}{2}\log(2\pi\tilde{\sigma}^2)
#'   + \log|\mathbf{V}^{-1/2}|
#'   - \frac{1}{2\tilde{\sigma}^2}
#'     \|\mathbf{V}^{-1/2}(\mathbf{y} - \hat{\mathbf{y}})\|^2
#' }
#' \eqn{\log|\mathbf{V}^{-1/2}|} is obtained from \code{VhalfInv_logdet}
#' when available, or computed directly from \code{VhalfInv}.
#'
#' \strong{Prior contribution.}
#' When \code{include_prior = TRUE} (default), the log-prior
#' \deqn{
#'   -\frac{1}{2\tilde{\sigma}^2}
#'   \sum_{k}\boldsymbol{\beta}_k^\top\boldsymbol{\Lambda}_k
#'   \boldsymbol{\beta}_k
#' }
#' is added, giving the penalised MAP log-likelihood coherent with the
#' smoothing spline objective. Set \code{include_prior = FALSE} for the
#' unpenalised marginal likelihood, which is more appropriate when comparing
#' models with different penalty structures or numbers of knots.
#'
#' \strong{Other GLM families.}
#' Uses \code{family$aic()} when available. For correlated models the
#' whitened residuals and fitted values are passed. When \code{family$aic()}
#' is unavailable, a deviance-based approximation is used (valid for
#' relative comparisons; a warning is emitted).
#'
#' This function returns the marginal (full) GLS log-likelihood, not the
#' REML log-likelihood. This is consistent with \code{REML = FALSE} in
#' \code{lme} and \code{gls}, and is the conventional choice for AIC/BIC
#' comparisons of fixed-effects structure.
#'
#' The \code{df} attribute is set to \eqn{N - \mathrm{trace}(\mathbf{XUGX}^\top)}.
#'
#' @param object A fitted lgspline model object.
#' @param include_prior Logical; add the log-prior penalty term. Default TRUE.
#' @param new_weights Numeric scalar or N-vector; optional observation weights
#'        overriding \code{object$weights}.
#' @param ... Not used.
#'
#' @return A \code{"logLik"} object with attributes \code{df} (effective
#'   degrees of freedom) and \code{nobs} (number of observations).
#'
#' @examples
#'
#' set.seed(1234)
#' t <- runif(1000, -10, 10)
#' y <- 2*sin(t) + -0.06*t^2 + rnorm(length(t))
#' model_fit <- lgspline(t, y)
#'
#' logLik(model_fit)
#' logLik(model_fit, include_prior = FALSE)
#'
#' AIC(model_fit)
#' BIC(model_fit)
#'
#' ## Compare models with different K using unpenalized likelihood
#' fit_k3 <- lgspline(t, y, K = 3)
#' fit_k7 <- lgspline(t, y, K = 7)
#' AIC(fit_k3, fit_k7)
#'
#' @seealso
#' \code{\link{lgspline}}, \code{\link{prior_loglik}},
#' \code{\link[stats]{logLik}}, \code{\link[stats]{AIC}},
#' \code{\link[stats]{BIC}}
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

  edf <- if(!is.null(object$trace_XUGX)) N - object$trace_XUGX else NA_real_
  has_corr <- !is.null(object$VhalfInv)

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
      c(object$VhalfInv %**% cbind(y - mu)),
      error = function(e) NULL
    )
  }

  ll <- NA_real_

  ## First try the exact Gaussian formulas when they apply.
  if(fam$family == "gaussian" && fam$link == "identity"){

    if(has_corr && !is.null(resid_w) && is.finite(logdet_Vhalfinv)){
      ss_w <- sum(wt*(resid_w)^2)
      ll   <- -0.5 * N * log(2 * pi * sigma2) +
        logdet_Vhalfinv -
        0.5 * ss_w / sigma2
    } else {
      resid <- y - mu
      ss    <- sum(wt*(resid)^2)
      ll    <- -0.5 * N * log(2 * pi * sigma2) - 0.5 * ss / sigma2
    }

  } else if(!is.null(fam$aic)){

    ## Otherwise ask the family for its AIC/log-likelihood contribution.
    ## For correlated models, pass whitened quantities to family$aic().
    #  For families where aic() depends on the raw scale (e.g. binomial)
    #  this is an approximation.
    if(has_corr && !is.null(resid_w)){
      y_eval  <- tryCatch(c(object$VhalfInv %*% cbind(y)), error = function(e) y)
      mu_eval <- tryCatch(c(object$VhalfInv %*% cbind(mu)), error = function(e) mu)
    } else {
      y_eval  <- y
      mu_eval <- mu
    }

    dev_resids <- tryCatch(fam$dev.resids(y_eval, mu_eval, wt), error = function(e) NULL)

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
      y_fb  <- tryCatch(c(object$VhalfInv %*% cbind(y)), error = function(e) y)
      mu_fb <- tryCatch(c(object$VhalfInv %*% cbind(mu)), error = function(e) mu)
    } else {
      y_fb  <- y
      mu_fb <- mu
    }

    dev_resids <- tryCatch(fam$dev.resids(y_fb, mu_fb, wt), error = function(e) NULL)

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
          if(has_corr) " with correlation structure log-determinant correction" else "",
          "; a family-specific constant may be omitted. ",
          "Valid for relative model comparison only."
        ),
        call. = FALSE
      )
    } else {
      warning("Could not compute log-likelihood for this family.", call. = FALSE)
    }
  }

  ## Add the penalty-induced prior term when requested.
  if(include_prior && is.finite(ll)){
    lp <- tryCatch(
      as.numeric(prior_loglik(object, sigmasq = sigma2)),
      error = function(e) {
        warning("Could not evaluate prior_loglik; prior term omitted.", call. = FALSE)
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


#' Print Closed-Form Fitted Equation from lgspline Model
#'
#' @description
#' Displays the closed-form polynomial equation for each partition of a fitted
#' lgspline model, along with partition boundary or cluster center information.
#' Optionally prints the first derivative, second derivative, or antiderivative
#' of the fitted equation with respect to a single specified variable.
#'
#' @param object A fitted lgspline model object.
#' @param x An object returned by \code{equation()} for printing.
#' @param digits Integer; decimal places for coefficient display. Default 4.
#' @param scientific Logical; use scientific notation for coefficients with
#'        absolute value < 1e-3 or > 1e4. Default FALSE.
#' @param show_bounds Logical; display partition bounds (1D) or knot midpoint
#'        boundaries (multi-D). Default TRUE.
#' @param predictor_names Character vector; custom names for predictor variables.
#'        If NULL (default), uses original column names or "_j_" labels.
#' @param response_name Character; label for response. If NULL (default), uses
#'        "y" for identity link Gaussian, or "link(E[y])" otherwise.
#' @param collapse_zero Logical; omit terms with coefficient exactly 0.
#'        Default TRUE.
#' @param first_derivative Default: NULL. Character name or integer index of
#'        the predictor variable with respect to which the first derivative
#'        is printed. Only one variable at a time is supported. When non-NULL,
#'        the printed equations show \eqn{df/dx_j} for each partition.
#' @param second_derivative Default: NULL. Character name or integer index of
#'        the predictor variable with respect to which the second derivative
#'        is printed. Only one variable at a time is supported. When non-NULL,
#'        the printed equations show \eqn{d^2f/dx_j^2} for each partition.
#'        Ignored if \code{first_derivative} is also non-NULL.
#' @param antiderivative Default: NULL. Character name or integer index of
#'        the predictor variable with respect to which the antiderivative
#'        (indefinite integral) is printed. Only one variable at a time is
#'        supported. When non-NULL, the printed equations show
#'        \eqn{\int f\, dx_j} for each partition, with an unspecified
#'        constant of integration \eqn{C}. Ignored if \code{first_derivative}
#'        or \code{second_derivative} is also non-NULL.
#' @param ... Not used.
#'
#' @details
#' For 1D models with K knots, partition boundaries are displayed as intervals
#' on the predictor scale. For multi-predictor models, partition boundaries are
#' computed as the midpoints between adjacent cluster centers along each
#' predictor dimension. When the model's \code{make_partition_list} contains
#' \code{knots} (midpoint boundaries between clusters), those are used directly.
#' Otherwise, cluster centers are displayed.
#'
#' Coefficients are displayed on the original (unstandardized) predictor scale.
#' For GLMs with non-identity link, the left-hand side shows the link function
#' applied to the expected response.
#'
#' \strong{Derivative and antiderivative modes.}
#' Only one of \code{first_derivative}, \code{second_derivative}, or
#' \code{antiderivative} may be non-NULL. If more than one is supplied, the
#' priority order is: first derivative, second derivative, antiderivative.
#'
#' Derivatives and antiderivatives are computed symbolically from the
#' polynomial coefficients. For a term \eqn{a x^n}, the first derivative is
#' \eqn{n a x^{n-1}}, the second derivative is \eqn{n(n-1) a x^{n-2}}, and
#' the antiderivative is \eqn{a x^{n+1}/(n+1)}. Cross-terms (interactions)
#' involving the target variable are differentiated or integrated with respect
#' to that variable only, treating all other variables as constants.
#'
#' A warning is emitted if the user attempts to differentiate or integrate with
#' respect to more than one variable simultaneously. Multi-variable calculus
#' operations should be performed one variable at a time by calling
#' \code{equation()} repeatedly.
#'
#' @return Invisibly returns a list with components:
#' \describe{
#'   \item{formulas}{Character vector of equation strings per partition.}
#'   \item{bounds}{Matrix or list of partition boundary information.}
#'   \item{link}{Character; link function name.}
#'   \item{mode}{Character; one of "equation", "first_derivative",
#'         "second_derivative", or "antiderivative".}
#'   \item{variable}{Character; the variable name for the calculus operation,
#'         or NULL if mode is "equation".}
#' }
#'
#' @examples
#'
#' ## 1D example
#' set.seed(1234)
#' t <- runif(500, -5, 5)
#' y <- 2*sin(t) + 0.1*t^2 + rnorm(length(t), 0, 0.5)
#' fit <- lgspline(t, y, K = 2)
#' equation(fit)
#' equation(fit, digits = 2, predictor_names = "time")
#'
#' ## First derivative with respect to predictor
#' equation(fit, first_derivative = 1)
#'
#' ## Second derivative
#' equation(fit, second_derivative = 1)
#'
#' ## Antiderivative
#' equation(fit, antiderivative = 1)
#'
#' ## 2D example with named predictors
#' x1 <- runif(300, 0, 10)
#' x2 <- runif(300, 0, 10)
#' y <- x1 + 0.5*x2 + 0.1*x1*x2 + rnorm(300)
#' fit2d <- lgspline(cbind(x1, x2), y, K = 3)
#' equation(fit2d, predictor_names = c("Length", "Width"))
#'
#' ## Derivative w.r.t. first variable only
#' equation(fit2d, first_derivative = "Length",
#'          predictor_names = c("Length", "Width"))
#'
#' ## GLM example
#' y_bin <- rbinom(500, 1, plogis(0.5*t))
#' fit_glm <- lgspline(t, y_bin, family = binomial(), K = 1)
#' equation(fit_glm)
#'
#' @seealso \code{\link{lgspline}}, \code{\link{plot.lgspline}},
#'   \code{\link{coef.lgspline}}
#'
#' @export
equation <- function(object, ...) {

  UseMethod("equation")
}

#' @rdname equation
#' @method equation lgspline
#' @export
equation.lgspline <- function(object,
                              digits = 4,
                              scientific = FALSE,
                              show_bounds = TRUE,
                              predictor_names = NULL,
                              response_name = NULL,
                              collapse_zero = TRUE,
                              first_derivative = NULL,
                              second_derivative = NULL,
                              antiderivative = NULL,
                              ...) {

  K  <- object$K
  nc <- object$p
  q  <- object$q
  b_names <- names(object$B[[1]])
  raw_nms <- object$raw_expansion_names

  ## Resolve predictor display names from B[[1]] at linear positions
  all_linear_cols <- sort(c(object$power1_cols, object$nonspline_cols))
  if (is.null(predictor_names)) {
    if (!is.null(b_names) && length(all_linear_cols) > 0) {
      predictor_names <- b_names[all_linear_cols]
    } else if (!is.null(raw_nms) && length(all_linear_cols) > 0) {
      predictor_names <- raw_nms[all_linear_cols]
    } else {
      predictor_names <- paste0("x", 1:q)
    }
  }
  if (length(predictor_names) < q) {
    predictor_names <- c(
      predictor_names,
      paste0("x", (length(predictor_names) + 1):q)
    )
  }

  ## Map from internal column index (the number inside "_j_") to
  #  1-based predictor display index into predictor_names.
  internal_idx_to_pred <- integer(0)
  if (!is.null(raw_nms) && length(all_linear_cols) > 0) {
    for (jj in seq_along(all_linear_cols)) {
      col_idx <- all_linear_cols[jj]
      raw_nm <- raw_nms[col_idx]
      idx_match <- regmatches(raw_nm, regexec("^_([0-9]+)_$", raw_nm))[[1]]
      if (length(idx_match) == 2) {
        internal_idx_to_pred[as.integer(idx_match[2])] <- jj
      }
    }
  }

  ## Calculus mode
  calc_mode <- "equation"
  calc_var_idx <- NULL
  calc_var_name <- NULL

  .resolve_calc_var <- function(spec, predictor_names, q, arg_name) {
    if (length(spec) > 1) {
      warning(
        "\n\t ", arg_name, " accepts only one variable at a time. ",
        "Only the first element will be used. Call equation() separately ",
        "for each variable.\n"
      )
      spec <- spec[1]
    }
    if (is.character(spec)) {
      idx <- match(spec, predictor_names)
      if (is.na(idx)) {
        stop(
          "\n\t ", arg_name, " variable '", spec,
          "' not found in predictor names: ",
          paste(predictor_names, collapse = ", "), ".\n"
        )
      }
      list(idx = idx, name = predictor_names[idx])
    } else {
      idx <- as.integer(spec)
      if (idx < 1 || idx > q) {
        stop(
          "\n\t ", arg_name, " index ", idx,
          " is out of range [1, ", q, "].\n"
        )
      }
      list(idx = idx, name = predictor_names[idx])
    }
  }

  if (!is.null(first_derivative)) {
    calc_mode <- "first_derivative"
    res <- .resolve_calc_var(first_derivative, predictor_names, q,
                             "first_derivative")
    calc_var_idx  <- res$idx
    calc_var_name <- res$name
  } else if (!is.null(second_derivative)) {
    calc_mode <- "second_derivative"
    res <- .resolve_calc_var(second_derivative, predictor_names, q,
                             "second_derivative")
    calc_var_idx  <- res$idx
    calc_var_name <- res$name
  } else if (!is.null(antiderivative)) {
    calc_mode <- "antiderivative"
    res <- .resolve_calc_var(antiderivative, predictor_names, q,
                             "antiderivative")
    calc_var_idx  <- res$idx
    calc_var_name <- res$name
  }

  ## Response name from link function
  fam <- object$family
  if (is.null(response_name)) {
    if (fam$family == "gaussian" && fam$link == "identity") {
      response_name <- "y"
    } else {
      response_name <- paste0(fam$link, "(E[y])")
    }
  }

  ## Parse a raw_expansion_name into structured factors.
  #  Raw names use "_j_" for variables and "x" (no underscores) as
  #  the interaction separator between complete "_j_" tokens.
  #  Examples: "_1_", "_1_^2", "_1_x_2_", "_2_x_1_^2", "_1_x_2_x_3_"
  .parse_term <- function(nm) {
    if (nm == "intercept" ||
        grepl("^intercept$", nm, ignore.case = TRUE)) {
      return(list(is_intercept = TRUE, factors = list()))
    }

    ## Extract all "_j_" and "_j_^n" tokens from the raw name
    #  Pattern: underscore, digits, underscore, optionally ^digits
    matches <- gregexpr("_([0-9]+)_(\\^[0-9]+)?", nm)
    tokens <- regmatches(nm, matches)[[1]]

    if (length(tokens) == 0) {
      return(list(is_intercept = FALSE,
                  factors = list(list(pred_idx = NA_integer_, power = 1L))))
    }

    factors <- lapply(tokens, function(tok) {
      ## Extract power if present
      pow_match <- regmatches(tok, regexec("\\^([0-9]+)$", tok))[[1]]
      if (length(pow_match) == 2) {
        power <- as.integer(pow_match[2])
        base_tok <- sub("\\^[0-9]+$", "", tok)
      } else {
        power <- 1L
        base_tok <- tok
      }
      ## Extract the column number
      idx_match <- regmatches(
        base_tok, regexec("^_([0-9]+)_$", base_tok)
      )[[1]]
      if (length(idx_match) == 2) {
        raw_col <- as.integer(idx_match[2])
        pred_idx <- internal_idx_to_pred[raw_col]
        if (is.na(pred_idx)) pred_idx <- NA_integer_
      } else {
        pred_idx <- NA_integer_
      }
      list(pred_idx = pred_idx, power = power)
    })
    list(is_intercept = FALSE, factors = factors)
  }

  ## Differentiate a single term w.r.t. variable at index var_idx
  .differentiate_term <- function(coef_val, parsed, var_idx) {
    if (parsed$is_intercept) return(NULL)
    hit <- NULL
    hit_pos <- 0L
    for (i in seq_along(parsed$factors)) {
      if (!is.na(parsed$factors[[i]]$pred_idx) &&
          parsed$factors[[i]]$pred_idx == var_idx) {
        hit <- parsed$factors[[i]]
        hit_pos <- i
        break
      }
    }
    if (is.null(hit)) return(NULL)
    new_coef <- coef_val * hit$power
    new_power <- hit$power - 1L
    new_factors <- parsed$factors
    if (new_power == 0L) {
      new_factors[[hit_pos]] <- NULL
    } else {
      new_factors[[hit_pos]]$power <- new_power
    }
    list(coef = new_coef,
         parsed = list(is_intercept = (length(new_factors) == 0),
                       factors = new_factors))
  }

  ## Second derivative: differentiate twice
  .differentiate2_term <- function(coef_val, parsed, var_idx) {
    first <- .differentiate_term(coef_val, parsed, var_idx)
    if (is.null(first)) return(NULL)
    if (first$parsed$is_intercept) return(NULL)
    .differentiate_term(first$coef, first$parsed, var_idx)
  }

  ## Antidifferentiate a single term w.r.t. variable at index var_idx
  .antidifferentiate_term <- function(coef_val, parsed, var_idx) {
    if (parsed$is_intercept) {
      new_factors <- list(list(pred_idx = var_idx, power = 1L))
      return(list(coef = coef_val,
                  parsed = list(is_intercept = FALSE,
                                factors = new_factors)))
    }
    hit_pos <- 0L
    for (i in seq_along(parsed$factors)) {
      if (!is.na(parsed$factors[[i]]$pred_idx) &&
          parsed$factors[[i]]$pred_idx == var_idx) {
        hit_pos <- i
        break
      }
    }
    new_factors <- parsed$factors
    if (hit_pos > 0) {
      old_power <- new_factors[[hit_pos]]$power
      new_power <- old_power + 1L
      new_coef <- coef_val / new_power
      new_factors[[hit_pos]]$power <- new_power
    } else {
      new_coef <- coef_val
      new_factors <- c(new_factors,
                       list(list(pred_idx = var_idx, power = 1L)))
    }
    list(coef = new_coef,
         parsed = list(is_intercept = FALSE, factors = new_factors))
  }

  ## Rebuild human-readable term from parsed structure
  .rebuild_term_name <- function(parsed, predictor_names) {
    if (parsed$is_intercept) return("intercept")
    parts <- vapply(parsed$factors, function(f) {
      if (is.na(f$pred_idx)) {
        nm <- "UNKNOWN"
      } else if (f$pred_idx <= length(predictor_names)) {
        nm <- predictor_names[f$pred_idx]
      } else {
        nm <- paste0("x", f$pred_idx)
      }
      if (f$power == 1) return(nm)
      paste0(nm, "^", f$power)
    }, character(1))
    paste(parts, collapse = "*")
  }

  ## Format one coefficient * term pair for display
  .fmt_coef <- function(val, term_nm, first_term = FALSE) {
    if (collapse_zero && abs(val) < .Machine$double.eps * 100) {
      return(NULL)
    }
    if (scientific && (abs(val) < 1e-3 || abs(val) > 1e4) && val != 0) {
      val_str <- format(val, digits = digits, scientific = TRUE)
    } else {
      val_str <- format(round(val, digits), nsmall = digits, trim = TRUE)
    }
    ## Intercept / constant
    if (term_nm == "intercept" || term_nm == "") {
      if (first_term) return(val_str)
      if (val >= 0) return(paste0(" + ", val_str))
      return(paste0(" - ", gsub("^-", "", val_str)))
    }
    ## Coefficient * term
    if (first_term) {
      if (abs(val - 1) < .Machine$double.eps * 100) return(term_nm)
      if (abs(val + 1) < .Machine$double.eps * 100) return(paste0("-", term_nm))
      return(paste0(val_str, "*", term_nm))
    }
    if (val >= 0) {
      if (abs(val - 1) < .Machine$double.eps * 100) return(paste0(" + ", term_nm))
      return(paste0(" + ", val_str, "*", term_nm))
    }
    if (abs(val + 1) < .Machine$double.eps * 100) return(paste0(" - ", term_nm))
    paste0(" - ", gsub("^-", "", val_str), "*", term_nm)
  }

  ## LHS label
  if (calc_mode == "first_derivative") {
    lhs_label <- paste0("d(", response_name, ")/d(", calc_var_name, ")")
  } else if (calc_mode == "second_derivative") {
    lhs_label <- paste0("d2(", response_name, ")/d(", calc_var_name, ")2")
  } else if (calc_mode == "antiderivative") {
    lhs_label <- paste0("integral ", response_name, " d(", calc_var_name, ")")
  } else {
    lhs_label <- response_name
  }

  ## Build formula string for each partition.
  #  Always parse from raw_expansion_names (reliable _j_ structure),
  #  display with predictor_names (human-readable).
  formulas <- character(K + 1)
  for (k in 1:(K + 1)) {
    coefs <- object$B[[k]]
    terms_out <- character(0)
    first <- TRUE

    for (i in seq_along(coefs)) {
      orig_coef <- coefs[i]
      parsed <- .parse_term(raw_nms[i])

      if (calc_mode == "equation") {
        display_nm <- .rebuild_term_name(parsed, predictor_names)
        term_str <- .fmt_coef(orig_coef, display_nm, first_term = first)
      } else if (calc_mode == "first_derivative") {
        result <- .differentiate_term(orig_coef, parsed, calc_var_idx)
        if (is.null(result)) next
        display_nm <- .rebuild_term_name(result$parsed, predictor_names)
        term_str <- .fmt_coef(result$coef, display_nm, first_term = first)
      } else if (calc_mode == "second_derivative") {
        result <- .differentiate2_term(orig_coef, parsed, calc_var_idx)
        if (is.null(result)) next
        display_nm <- .rebuild_term_name(result$parsed, predictor_names)
        term_str <- .fmt_coef(result$coef, display_nm, first_term = first)
      } else if (calc_mode == "antiderivative") {
        result <- .antidifferentiate_term(orig_coef, parsed, calc_var_idx)
        if (is.null(result)) next
        display_nm <- .rebuild_term_name(result$parsed, predictor_names)
        term_str <- .fmt_coef(result$coef, display_nm, first_term = first)
      }

      if (!is.null(term_str)) {
        terms_out <- c(terms_out, term_str)
        first <- FALSE
      }
    }

    if (length(terms_out) == 0) {
      formulas[k] <- paste0(lhs_label, " = 0")
    } else {
      rhs <- paste0(terms_out, collapse = "")
      if (calc_mode == "antiderivative") rhs <- paste0(rhs, " + C")
      formulas[k] <- paste0(lhs_label, " = ", rhs)
    }
  }

  ## Partition boundary info
  bounds_info <- NULL
  if (show_bounds && K > 0) {
    if (q == 1 && !is.null(object$knots)) {
      ## 1D interval bounds
      knot_vals <- sort(c(object$knots))
      bounds_info <- matrix(NA, nrow = K + 1, ncol = 2)
      colnames(bounds_info) <- c("lower", "upper")
      rownames(bounds_info) <- paste0("Partition ", 1:(K + 1))

      all_x <- unlist(lapply(object$X, function(Xk) {
        if (nrow(Xk) > 0 && length(object$power1_cols) > 0) {
          Xk[, object$power1_cols[1]]
        } else if (nrow(Xk) > 0 && length(object$nonspline_cols) > 0) {
          Xk[, object$nonspline_cols[1]]
        } else {
          NULL
        }
      }))
      if (length(all_x) > 0) {
        x_min <- min(all_x, na.rm = TRUE)
        x_max <- max(all_x, na.rm = TRUE)
      } else {
        x_min <- -Inf
        x_max <- Inf
      }

      bounds_info[1, 1] <- x_min
      bounds_info[1, 2] <- knot_vals[1]
      if (K >= 2) {
        for (kk in 2:K) {
          bounds_info[kk, 1] <- knot_vals[kk - 1]
          bounds_info[kk, 2] <- knot_vals[kk]
        }
      }
      bounds_info[K + 1, 1] <- knot_vals[K]
      bounds_info[K + 1, 2] <- x_max

    } else if (!is.null(object$make_partition_list)) {
      ## Multi-D bounding boxes from cluster centers
      centers <- object$make_partition_list$centers

      if (!is.null(centers)) {
        n_centers <- nrow(centers)
        n_dims <- ncol(centers)

        ## Only show spline predictor dimensions, not indicator columns
        n_spline <- length(object$power1_cols)
        if (n_spline > 0 && n_spline < n_dims) {
          show_dims <- seq_len(n_spline)
        } else {
          show_dims <- seq_len(n_dims)
        }

        center_colnames <- colnames(centers)
        if (is.null(center_colnames)) {
          if (!is.null(predictor_names) && n_dims <= length(predictor_names)) {
            center_colnames <- predictor_names[1:n_dims]
          } else {
            center_colnames <- paste0("x", 1:n_dims)
          }
        }

        ## Training predictor ranges per dimension
        training_ranges <- matrix(NA, 2, n_dims)
        rownames(training_ranges) <- c("min", "max")
        for (jj in seq_len(n_dims)) {
          col_vals <- unlist(lapply(object$X, function(Xk) {
            if (nrow(Xk) > 0 && jj <= length(all_linear_cols)) {
              Xk[, all_linear_cols[jj]]
            } else {
              NULL
            }
          }))
          if (length(col_vals) > 0) {
            training_ranges[1, jj] <- min(col_vals, na.rm = TRUE)
            training_ranges[2, jj] <- max(col_vals, na.rm = TRUE)
          }
        }

        ## Per-partition bounding boxes via midpoints between sorted centers
        bounds_lower <- matrix(NA, n_centers, n_dims)
        bounds_upper <- matrix(NA, n_centers, n_dims)
        colnames(bounds_lower) <- center_colnames
        colnames(bounds_upper) <- center_colnames

        for (jj in seq_len(n_dims)) {
          center_vals_j <- centers[, jj]
          ord_j <- order(center_vals_j)
          sorted_j <- center_vals_j[ord_j]
          if (length(sorted_j) > 1) {
            midpts_j <- (sorted_j[-length(sorted_j)] + sorted_j[-1]) / 2
          } else {
            midpts_j <- numeric(0)
          }
          for (rank in seq_along(ord_j)) {
            part_idx <- ord_j[rank]
            lo <- if (rank == 1 && !is.na(training_ranges[1, jj])) {
              training_ranges[1, jj]
            } else if (rank == 1) {
              -Inf
            } else {
              midpts_j[rank - 1]
            }
            hi <- if (rank == length(ord_j) && !is.na(training_ranges[2, jj])) {
              training_ranges[2, jj]
            } else if (rank == length(ord_j)) {
              Inf
            } else {
              midpts_j[rank]
            }
            bounds_lower[part_idx, jj] <- lo
            bounds_upper[part_idx, jj] <- hi
          }
        }

        bounds_info <- list(
          lower = bounds_lower[, show_dims, drop = FALSE],
          upper = bounds_upper[, show_dims, drop = FALSE],
          centers = centers,
          colnames = center_colnames[show_dims]
        )
      }
    }
  }

  ## Print
  cat("\n")
  div_major <- paste(rep("=", 72), collapse = "")
  div_minor <- paste(rep("-", 72), collapse = "")

  cat(div_major, "\n", sep = "")
  if (calc_mode == "equation") {
    cat("Fitted Equations: lgspline Model\n")
  } else if (calc_mode == "first_derivative") {
    cat("First Derivative w.r.t. ", calc_var_name,
        ": lgspline Model\n", sep = "")
  } else if (calc_mode == "second_derivative") {
    cat("Second Derivative w.r.t. ", calc_var_name,
        ": lgspline Model\n", sep = "")
  } else if (calc_mode == "antiderivative") {
    cat("Antiderivative w.r.t. ", calc_var_name,
        ": lgspline Model\n", sep = "")
  }
  cat(div_major, "\n", sep = "")
  cat("Family:", fam$family, "  Link:", fam$link, "\n")
  cat("Partitions:", K + 1, "  Basis functions per partition:", nc, "\n")
  cat(div_minor, "\n", sep = "")

  for (k in 1:(K + 1)) {
    cat("\n")
    cat("Partition", k, "\n")

    if (show_bounds && !is.null(bounds_info)) {
      if (q == 1 && is.matrix(bounds_info) && ncol(bounds_info) == 2) {
        lb <- bounds_info[k, 1]
        ub <- bounds_info[k, 2]
        lb_str <- if (is.finite(lb)) {
          format(round(lb, digits), nsmall = digits)
        } else { "-Inf" }
        ub_str <- if (is.finite(ub)) {
          format(round(ub, digits), nsmall = digits)
        } else { "+Inf" }
        lb_bracket <- if (k == 1) "[" else "("
        ub_bracket <- if (k == K + 1) "]" else ")"
        cat("  Bounds: ", predictor_names[1], " in ",
            lb_bracket, lb_str, ", ", ub_str, ub_bracket,
            "\n", sep = "")

      } else if (is.list(bounds_info) && !is.null(bounds_info$lower)) {
        n_dims <- ncol(bounds_info$lower)
        cnames <- bounds_info$colnames
        cat("  Bounds:\n")
        for (jj in seq_len(n_dims)) {
          lo <- bounds_info$lower[k, jj]
          hi <- bounds_info$upper[k, jj]
          lo_str <- if (is.finite(lo)) {
            format(round(lo, min(digits, 2)), nsmall = min(digits, 2))
          } else { "-Inf" }
          hi_str <- if (is.finite(hi)) {
            format(round(hi, min(digits, 2)), nsmall = min(digits, 2))
          } else { "+Inf" }
          lb_bracket <- if (lo == bounds_info$lower[
            which.min(bounds_info$lower[, jj]), jj]) "[" else "("
          ub_bracket <- if (hi == bounds_info$upper[
            which.max(bounds_info$upper[, jj]), jj]) "]" else ")"
          cat("    ", cnames[jj], " in ",
              lb_bracket, lo_str, ", ", hi_str, ub_bracket,
              "\n", sep = "")
        }
      }
    }

    cat("  ", formulas[k], "\n", sep = "")
  }

  cat("\n")
  cat(div_minor, "\n", sep = "")
  if (!is.null(object$sigmasq_tilde) && is.finite(object$sigmasq_tilde)) {
    cat("Dispersion (sigma^2):",
        format(object$sigmasq_tilde, digits = digits), "\n")
  }
  if (!is.null(object$trace_XUGX) && is.finite(object$trace_XUGX)) {
    cat("Effective df:",
        format(object$N - object$trace_XUGX, digits = 2), "\n")
  }
  cat(div_major, "\n", sep = "")
  cat("\n")

  invisible(list(
    formulas = formulas,
    bounds   = bounds_info,
    link     = fam$link,
    mode     = calc_mode,
    variable = calc_var_name
  ))
}


#' @rdname equation
#' @method print equation
#' @export
print.equation <- function(x, ...) {
  cat("\nFitted Equations")
  if (!is.null(x$mode) && x$mode != "equation") {
    mode_label <- switch(
      x$mode,
      "first_derivative"  = "First Derivative",
      "second_derivative" = "Second Derivative",
      "antiderivative"    = "Antiderivative"
    )
    cat(" (", mode_label, " w.r.t. ", x$variable, ")", sep = "")
  }
  cat(":\n")
  for (k in seq_along(x$formulas)) {
    cat("  Partition", k, ":", x$formulas[k], "\n")
  }
  if (!is.null(x$bounds)) {
    if (is.matrix(x$bounds)) {
      cat("\nPartition Bounds:\n")
      print(x$bounds)
    } else if (is.list(x$bounds) && !is.null(x$bounds$lower)) {
      cat("\nPartition Bounds (lower):\n")
      print(x$bounds$lower)
      cat("\nPartition Bounds (upper):\n")
      print(x$bounds$upper)
    }
  }
  invisible(x)
}


#' @rdname equation
#' @method equation lgspline
#' @export
equation.lgspline <- function(object,
                              digits = 4,
                              scientific = FALSE,
                              show_bounds = TRUE,
                              predictor_names = NULL,
                              response_name = NULL,
                              collapse_zero = TRUE,
                              first_derivative = NULL,
                              second_derivative = NULL,
                              antiderivative = NULL,
                              ...) {

  K  <- object$K
  nc <- object$p
  q  <- object$q
  b_names <- names(object$B[[1]])
  raw_nms <- object$raw_expansion_names

  ## Resolve predictor display names from B[[1]] at linear positions
  all_linear_cols <- sort(c(object$power1_cols, object$nonspline_cols))
  if (is.null(predictor_names)) {
    if (!is.null(b_names) && length(all_linear_cols) > 0) {
      predictor_names <- b_names[all_linear_cols]
    } else if (!is.null(raw_nms) && length(all_linear_cols) > 0) {
      predictor_names <- raw_nms[all_linear_cols]
    } else {
      predictor_names <- paste0("x", 1:q)
    }
  }
  if (length(predictor_names) < q) {
    predictor_names <- c(
      predictor_names,
      paste0("x", (length(predictor_names) + 1):q)
    )
  }

  ## Map from internal column index (the number inside "_j_") to
  #  1-based predictor display index into predictor_names.
  internal_idx_to_pred <- integer(0)
  if (!is.null(raw_nms) && length(all_linear_cols) > 0) {
    for (jj in seq_along(all_linear_cols)) {
      col_idx <- all_linear_cols[jj]
      raw_nm <- raw_nms[col_idx]
      idx_match <- regmatches(raw_nm, regexec("^_([0-9]+)_$", raw_nm))[[1]]
      if (length(idx_match) == 2) {
        internal_idx_to_pred[as.integer(idx_match[2])] <- jj
      }
    }
  }

  ## Calculus mode
  calc_mode <- "equation"
  calc_var_idx <- NULL
  calc_var_name <- NULL

  .resolve_calc_var <- function(spec, predictor_names, q, arg_name) {
    if (length(spec) > 1) {
      warning(
        "\n\t ", arg_name, " accepts only one variable at a time. ",
        "Only the first element will be used. Call equation() separately ",
        "for each variable.\n"
      )
      spec <- spec[1]
    }
    if (is.character(spec)) {
      idx <- match(spec, predictor_names)
      if (is.na(idx)) {
        stop(
          "\n\t ", arg_name, " variable '", spec,
          "' not found in predictor names: ",
          paste(predictor_names, collapse = ", "), ".\n"
        )
      }
      list(idx = idx, name = predictor_names[idx])
    } else {
      idx <- as.integer(spec)
      if (idx < 1 || idx > q) {
        stop(
          "\n\t ", arg_name, " index ", idx,
          " is out of range [1, ", q, "].\n"
        )
      }
      list(idx = idx, name = predictor_names[idx])
    }
  }

  if (!is.null(first_derivative)) {
    calc_mode <- "first_derivative"
    res <- .resolve_calc_var(first_derivative, predictor_names, q,
                             "first_derivative")
    calc_var_idx  <- res$idx
    calc_var_name <- res$name
  } else if (!is.null(second_derivative)) {
    calc_mode <- "second_derivative"
    res <- .resolve_calc_var(second_derivative, predictor_names, q,
                             "second_derivative")
    calc_var_idx  <- res$idx
    calc_var_name <- res$name
  } else if (!is.null(antiderivative)) {
    calc_mode <- "antiderivative"
    res <- .resolve_calc_var(antiderivative, predictor_names, q,
                             "antiderivative")
    calc_var_idx  <- res$idx
    calc_var_name <- res$name
  }

  ## Response name from link function
  fam <- object$family
  if (is.null(response_name)) {
    if (fam$family == "gaussian" && fam$link == "identity") {
      response_name <- "y"
    } else {
      response_name <- paste0(fam$link, "(E[y])")
    }
  }

  ## Parse a raw_expansion_name into structured factors.
  #  Raw names use "_j_" for variables and "x" (no underscores) as
  #  the interaction separator between complete "_j_" tokens.
  #  Examples: "_1_", "_1_^2", "_1_x_2_", "_2_x_1_^2", "_1_x_2_x_3_"
  #
  #  Key insight: strsplit on "_x_" eats the trailing/leading underscores

  #  of the _j_ tokens. Instead, use a regex that finds all _j_ or _j_^n
  #  tokens directly.
  .parse_term <- function(nm) {
    if (nm == "intercept" ||
        grepl("^intercept$", nm, ignore.case = TRUE)) {
      return(list(is_intercept = TRUE, factors = list()))
    }

    ## Extract all "_j_" and "_j_^n" tokens from the raw name
    #  Pattern: underscore, digits, underscore, optionally ^digits
    matches <- gregexpr("_([0-9]+)_(\\^[0-9]+)?", nm)
    tokens <- regmatches(nm, matches)[[1]]

    if (length(tokens) == 0) {
      return(list(is_intercept = FALSE,
                  factors = list(list(pred_idx = NA_integer_, power = 1L))))
    }

    factors <- lapply(tokens, function(tok) {
      ## Extract power if present
      pow_match <- regmatches(tok, regexec("\\^([0-9]+)$", tok))[[1]]
      if (length(pow_match) == 2) {
        power <- as.integer(pow_match[2])
        base_tok <- sub("\\^[0-9]+$", "", tok)
      } else {
        power <- 1L
        base_tok <- tok
      }
      ## Extract the column number
      idx_match <- regmatches(
        base_tok, regexec("^_([0-9]+)_$", base_tok)
      )[[1]]
      if (length(idx_match) == 2) {
        raw_col <- as.integer(idx_match[2])
        pred_idx <- internal_idx_to_pred[raw_col]
        if (is.na(pred_idx)) pred_idx <- NA_integer_
      } else {
        pred_idx <- NA_integer_
      }
      list(pred_idx = pred_idx, power = power)
    })
    list(is_intercept = FALSE, factors = factors)
  }

  ## Differentiate a single term w.r.t. variable at index var_idx
  .differentiate_term <- function(coef_val, parsed, var_idx) {
    if (parsed$is_intercept) return(NULL)
    hit <- NULL
    hit_pos <- 0L
    for (i in seq_along(parsed$factors)) {
      if (!is.na(parsed$factors[[i]]$pred_idx) &&
          parsed$factors[[i]]$pred_idx == var_idx) {
        hit <- parsed$factors[[i]]
        hit_pos <- i
        break
      }
    }
    if (is.null(hit)) return(NULL)
    new_coef <- coef_val * hit$power
    new_power <- hit$power - 1L
    new_factors <- parsed$factors
    if (new_power == 0L) {
      new_factors[[hit_pos]] <- NULL
    } else {
      new_factors[[hit_pos]]$power <- new_power
    }
    list(coef = new_coef,
         parsed = list(is_intercept = (length(new_factors) == 0),
                       factors = new_factors))
  }

  ## Second derivative: differentiate twice
  .differentiate2_term <- function(coef_val, parsed, var_idx) {
    first <- .differentiate_term(coef_val, parsed, var_idx)
    if (is.null(first)) return(NULL)
    if (first$parsed$is_intercept) return(NULL)
    .differentiate_term(first$coef, first$parsed, var_idx)
  }

  ## Antidifferentiate a single term w.r.t. variable at index var_idx
  .antidifferentiate_term <- function(coef_val, parsed, var_idx) {
    if (parsed$is_intercept) {
      new_factors <- list(list(pred_idx = var_idx, power = 1L))
      return(list(coef = coef_val,
                  parsed = list(is_intercept = FALSE,
                                factors = new_factors)))
    }
    hit_pos <- 0L
    for (i in seq_along(parsed$factors)) {
      if (!is.na(parsed$factors[[i]]$pred_idx) &&
          parsed$factors[[i]]$pred_idx == var_idx) {
        hit_pos <- i
        break
      }
    }
    new_factors <- parsed$factors
    if (hit_pos > 0) {
      old_power <- new_factors[[hit_pos]]$power
      new_power <- old_power + 1L
      new_coef <- coef_val / new_power
      new_factors[[hit_pos]]$power <- new_power
    } else {
      new_coef <- coef_val
      new_factors <- c(new_factors,
                       list(list(pred_idx = var_idx, power = 1L)))
    }
    list(coef = new_coef,
         parsed = list(is_intercept = FALSE, factors = new_factors))
  }

  ## Rebuild human-readable term from parsed structure
  .rebuild_term_name <- function(parsed, predictor_names) {
    if (parsed$is_intercept) return("intercept")
    parts <- vapply(parsed$factors, function(f) {
      if (is.na(f$pred_idx)) {
        nm <- "UNKNOWN"
      } else if (f$pred_idx <= length(predictor_names)) {
        nm <- predictor_names[f$pred_idx]
      } else {
        nm <- paste0("x", f$pred_idx)
      }
      if (f$power == 1) return(nm)
      paste0(nm, "^", f$power)
    }, character(1))
    paste(parts, collapse = "*")
  }

  ## Format one coefficient * term pair for display
  .fmt_coef <- function(val, term_nm, first_term = FALSE) {
    if (collapse_zero && abs(val) < .Machine$double.eps * 100) {
      return(NULL)
    }
    if (scientific && (abs(val) < 1e-3 || abs(val) > 1e4) && val != 0) {
      val_str <- format(val, digits = digits, scientific = TRUE)
    } else {
      val_str <- format(round(val, digits), nsmall = digits, trim = TRUE)
    }
    ## Intercept / constant
    if (term_nm == "intercept" || term_nm == "") {
      if (first_term) return(val_str)
      if (val >= 0) return(paste0(" + ", val_str))
      return(paste0(" - ", gsub("^-", "", val_str)))
    }
    ## Coefficient * term
    if (first_term) {
      if (abs(val - 1) < .Machine$double.eps * 100) return(term_nm)
      if (abs(val + 1) < .Machine$double.eps * 100) return(paste0("-", term_nm))
      return(paste0(val_str, "*", term_nm))
    }
    if (val >= 0) {
      if (abs(val - 1) < .Machine$double.eps * 100) return(paste0(" + ", term_nm))
      return(paste0(" + ", val_str, "*", term_nm))
    }
    if (abs(val + 1) < .Machine$double.eps * 100) return(paste0(" - ", term_nm))
    paste0(" - ", gsub("^-", "", val_str), "*", term_nm)
  }

  ## LHS label
  if (calc_mode == "first_derivative") {
    lhs_label <- paste0("d(", response_name, ")/d(", calc_var_name, ")")
  } else if (calc_mode == "second_derivative") {
    lhs_label <- paste0("d2(", response_name, ")/d(", calc_var_name, ")2")
  } else if (calc_mode == "antiderivative") {
    lhs_label <- paste0("integral ", response_name, " d(", calc_var_name, ")")
  } else {
    lhs_label <- response_name
  }

  ## Build formula string for each partition.
  #  Always parse from raw_expansion_names (reliable _j_ structure),
  #  display with predictor_names (human-readable).
  formulas <- character(K + 1)
  for (k in 1:(K + 1)) {
    coefs <- object$B[[k]]
    terms_out <- character(0)
    first <- TRUE

    for (i in seq_along(coefs)) {
      orig_coef <- coefs[i]
      parsed <- .parse_term(raw_nms[i])

      if (calc_mode == "equation") {
        display_nm <- .rebuild_term_name(parsed, predictor_names)
        term_str <- .fmt_coef(orig_coef, display_nm, first_term = first)
      } else if (calc_mode == "first_derivative") {
        result <- .differentiate_term(orig_coef, parsed, calc_var_idx)
        if (is.null(result)) next
        display_nm <- .rebuild_term_name(result$parsed, predictor_names)
        term_str <- .fmt_coef(result$coef, display_nm, first_term = first)
      } else if (calc_mode == "second_derivative") {
        result <- .differentiate2_term(orig_coef, parsed, calc_var_idx)
        if (is.null(result)) next
        display_nm <- .rebuild_term_name(result$parsed, predictor_names)
        term_str <- .fmt_coef(result$coef, display_nm, first_term = first)
      } else if (calc_mode == "antiderivative") {
        result <- .antidifferentiate_term(orig_coef, parsed, calc_var_idx)
        if (is.null(result)) next
        display_nm <- .rebuild_term_name(result$parsed, predictor_names)
        term_str <- .fmt_coef(result$coef, display_nm, first_term = first)
      }

      if (!is.null(term_str)) {
        terms_out <- c(terms_out, term_str)
        first <- FALSE
      }
    }

    if (length(terms_out) == 0) {
      formulas[k] <- paste0(lhs_label, " = 0")
    } else {
      rhs <- paste0(terms_out, collapse = "")
      if (calc_mode == "antiderivative") rhs <- paste0(rhs, " + C")
      formulas[k] <- paste0(lhs_label, " = ", rhs)
    }
  }

  ## Partition boundary info
  bounds_info <- NULL
  if (show_bounds && K > 0) {
    if (q == 1 && !is.null(object$knots)) {
      ## 1D interval bounds
      knot_vals <- sort(c(object$knots))
      bounds_info <- matrix(NA, nrow = K + 1, ncol = 2)
      colnames(bounds_info) <- c("lower", "upper")
      rownames(bounds_info) <- paste0("Partition ", 1:(K + 1))

      all_x <- unlist(lapply(object$X, function(Xk) {
        if (nrow(Xk) > 0 && length(object$power1_cols) > 0) {
          Xk[, object$power1_cols[1]]
        } else if (nrow(Xk) > 0 && length(object$nonspline_cols) > 0) {
          Xk[, object$nonspline_cols[1]]
        } else {
          NULL
        }
      }))
      if (length(all_x) > 0) {
        x_min <- min(all_x, na.rm = TRUE)
        x_max <- max(all_x, na.rm = TRUE)
      } else {
        x_min <- -Inf
        x_max <- Inf
      }

      bounds_info[1, 1] <- x_min
      bounds_info[1, 2] <- knot_vals[1]
      if (K >= 2) {
        for (kk in 2:K) {
          bounds_info[kk, 1] <- knot_vals[kk - 1]
          bounds_info[kk, 2] <- knot_vals[kk]
        }
      }
      bounds_info[K + 1, 1] <- knot_vals[K]
      bounds_info[K + 1, 2] <- x_max

    } else if (!is.null(object$make_partition_list)) {
      ## Multi-D bounding boxes from cluster centers
      centers <- object$make_partition_list$centers

      if (!is.null(centers)) {
        n_centers <- nrow(centers)
        n_dims <- ncol(centers)

        ## Only show spline predictor dimensions, not indicator columns
        n_spline <- length(object$power1_cols)
        if (n_spline > 0 && n_spline < n_dims) {
          show_dims <- seq_len(n_spline)
        } else {
          show_dims <- seq_len(n_dims)
        }

        center_colnames <- colnames(centers)
        if (is.null(center_colnames)) {
          if (!is.null(predictor_names) && n_dims <= length(predictor_names)) {
            center_colnames <- predictor_names[1:n_dims]
          } else {
            center_colnames <- paste0("x", 1:n_dims)
          }
        }

        ## Training predictor ranges per dimension
        training_ranges <- matrix(NA, 2, n_dims)
        rownames(training_ranges) <- c("min", "max")
        for (jj in seq_len(n_dims)) {
          col_vals <- unlist(lapply(object$X, function(Xk) {
            if (nrow(Xk) > 0 && jj <= length(all_linear_cols)) {
              Xk[, all_linear_cols[jj]]
            } else {
              NULL
            }
          }))
          if (length(col_vals) > 0) {
            training_ranges[1, jj] <- min(col_vals, na.rm = TRUE)
            training_ranges[2, jj] <- max(col_vals, na.rm = TRUE)
          }
        }

        ## Per-partition bounding boxes via midpoints between sorted centers
        bounds_lower <- matrix(NA, n_centers, n_dims)
        bounds_upper <- matrix(NA, n_centers, n_dims)
        colnames(bounds_lower) <- center_colnames
        colnames(bounds_upper) <- center_colnames

        for (jj in seq_len(n_dims)) {
          center_vals_j <- centers[, jj]
          ord_j <- order(center_vals_j)
          sorted_j <- center_vals_j[ord_j]
          if (length(sorted_j) > 1) {
            midpts_j <- (sorted_j[-length(sorted_j)] + sorted_j[-1]) / 2
          } else {
            midpts_j <- numeric(0)
          }
          for (rank in seq_along(ord_j)) {
            part_idx <- ord_j[rank]
            lo <- if (rank == 1 && !is.na(training_ranges[1, jj])) {
              training_ranges[1, jj]
            } else if (rank == 1) {
              -Inf
            } else {
              midpts_j[rank - 1]
            }
            hi <- if (rank == length(ord_j) && !is.na(training_ranges[2, jj])) {
              training_ranges[2, jj]
            } else if (rank == length(ord_j)) {
              Inf
            } else {
              midpts_j[rank]
            }
            bounds_lower[part_idx, jj] <- lo
            bounds_upper[part_idx, jj] <- hi
          }
        }

        bounds_info <- list(
          lower = bounds_lower[, show_dims, drop = FALSE],
          upper = bounds_upper[, show_dims, drop = FALSE],
          centers = centers,
          colnames = center_colnames[show_dims]
        )
      }
    }
  }

  ## Print
  cat("\n")
  div_major <- paste(rep("=", 72), collapse = "")
  div_minor <- paste(rep("-", 72), collapse = "")

  cat(div_major, "\n", sep = "")
  if (calc_mode == "equation") {
    cat("Fitted Equations: lgspline Model\n")
  } else if (calc_mode == "first_derivative") {
    cat("First Derivative w.r.t. ", calc_var_name,
        ": lgspline Model\n", sep = "")
  } else if (calc_mode == "second_derivative") {
    cat("Second Derivative w.r.t. ", calc_var_name,
        ": lgspline Model\n", sep = "")
  } else if (calc_mode == "antiderivative") {
    cat("Antiderivative w.r.t. ", calc_var_name,
        ": lgspline Model\n", sep = "")
  }
  cat(div_major, "\n", sep = "")
  cat("Family:", fam$family, "  Link:", fam$link, "\n")
  cat("Partitions:", K + 1, "  Basis functions per partition:", nc, "\n")
  cat(div_minor, "\n", sep = "")

  for (k in 1:(K + 1)) {
    cat("\n")
    cat("Partition", k, "\n")

    if (show_bounds && !is.null(bounds_info)) {
      if (q == 1 && is.matrix(bounds_info) && ncol(bounds_info) == 2) {
        lb <- bounds_info[k, 1]
        ub <- bounds_info[k, 2]
        lb_str <- if (is.finite(lb)) {
          format(round(lb, digits), nsmall = digits)
        } else { "-Inf" }
        ub_str <- if (is.finite(ub)) {
          format(round(ub, digits), nsmall = digits)
        } else { "+Inf" }
        lb_bracket <- if (k == 1) "[" else "("
        ub_bracket <- if (k == K + 1) "]" else ")"
        cat("  Bounds: ", predictor_names[1], " in ",
            lb_bracket, lb_str, ", ", ub_str, ub_bracket,
            "\n", sep = "")

      } else if (is.list(bounds_info) && !is.null(bounds_info$lower)) {
        n_dims <- ncol(bounds_info$lower)
        cnames <- bounds_info$colnames
        cat("  Bounds:\n")
        for (jj in seq_len(n_dims)) {
          lo <- bounds_info$lower[k, jj]
          hi <- bounds_info$upper[k, jj]
          lo_str <- if (is.finite(lo)) {
            format(round(lo, min(digits, 2)), nsmall = min(digits, 2))
          } else { "-Inf" }
          hi_str <- if (is.finite(hi)) {
            format(round(hi, min(digits, 2)), nsmall = min(digits, 2))
          } else { "+Inf" }
          lb_bracket <- if (lo == bounds_info$lower[
            which.min(bounds_info$lower[, jj]), jj]) "[" else "("
          ub_bracket <- if (hi == bounds_info$upper[
            which.max(bounds_info$upper[, jj]), jj]) "]" else ")"
          cat("    ", cnames[jj], " in ",
              lb_bracket, lo_str, ", ", hi_str, ub_bracket,
              "\n", sep = "")
        }
      }
    }

    cat("  ", formulas[k], "\n", sep = "")
  }

  cat("\n")
  cat(div_minor, "\n", sep = "")
  if (!is.null(object$sigmasq_tilde) && is.finite(object$sigmasq_tilde)) {
    cat("Dispersion (sigma^2):",
        format(object$sigmasq_tilde, digits = digits), "\n")
  }
  if (!is.null(object$trace_XUGX) && is.finite(object$trace_XUGX)) {
    cat("Effective df:",
        format(object$N - object$trace_XUGX, digits = 2), "\n")
  }
  cat(div_major, "\n", sep = "")
  cat("\n")

  invisible(list(
    formulas = formulas,
    bounds   = bounds_info,
    link     = fam$link,
    mode     = calc_mode,
    variable = calc_var_name
  ))
}


#' @rdname equation
#' @method print equation
#' @export
print.equation <- function(x, ...) {
  cat("\nFitted Equations")
  if (!is.null(x$mode) && x$mode != "equation") {
    mode_label <- switch(
      x$mode,
      "first_derivative"  = "First Derivative",
      "second_derivative" = "Second Derivative",
      "antiderivative"    = "Antiderivative"
    )
    cat(" (", mode_label, " w.r.t. ", x$variable, ")", sep = "")
  }
  cat(":\n")
  for (k in seq_along(x$formulas)) {
    cat("  Partition", k, ":", x$formulas[k], "\n")
  }
  if (!is.null(x$bounds)) {
    if (is.matrix(x$bounds)) {
      cat("\nPartition Bounds:\n")
      print(x$bounds)
    } else if (is.list(x$bounds) && !is.null(x$bounds$lower)) {
      cat("\nPartition Bounds (lower):\n")
      print(x$bounds$lower)
      cat("\nPartition Bounds (upper):\n")
      print(x$bounds$upper)
    }
  }
  invisible(x)
}

#' Generic for Numerical Integration
#'
#' @description
#' S3 generic that dispatches to \code{integrate.lgspline} for fitted
#' \code{lgspline} objects and falls back to \code{\link[stats]{integrate}}
#' for ordinary functions.
#'
#' @name integrate
#' @param f A fitted model object or a function.
#' @param ... Arguments passed to methods.
#'
#' @export
integrate <- function(f, ...) UseMethod("integrate")


#' Default method for the integrate generic
#'
#' @rdname integrate
#' @method integrate default
#' @export
integrate.default <- function(f, ...){
  stats::integrate(f, ...)
}





