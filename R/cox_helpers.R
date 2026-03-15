#' Cox Proportional Hazards Helpers for lgspline
#'
#' Functions for fitting Cox proportional hazards regression models within the
#' lgspline framework. Analogous to the Weibull AFT helpers, these provide
#' the partial log-likelihood, score, information, and all interface functions
#' needed by lgspline's unconstrained fitting, penalty tuning, and inference
#' machinery.
#'
#'
#' @name cox_helpers
NULL


#' Compute Cox Partial Log-Likelihood
#'
#' @description
#' Evaluates the Cox partial log-likelihood for a given coefficient vector,
#' using the Breslow approximation for tied event times.
#'
#' Observations must be sorted in ascending order of survival time before
#' calling this function. The internal helpers handle sorting automatically;
#' this function is exposed for diagnostics and testing.
#'
#' @param eta Numeric vector of linear predictors \eqn{\mathbf{X}\boldsymbol{\beta}},
#'   length N, sorted by ascending event time.
#' @param status Integer/logical vector of event indicators (1 = event,
#'   0 = censored), same length and order as \code{eta}.
#' @param y Optional numeric vector of observed event/censor times, same length
#'   and order as \code{eta}. When supplied, tied event times are handled by the
#'   Breslow approximation using a common risk-set denominator within each tied
#'   event-time block. When omitted, the function assumes there are no ties (or
#'   that ties have already been expanded appropriately).
#' @param weights Optional numeric vector of observation weights (default 1).
#'
#' @return Scalar partial log-likelihood value.
#'
#' @details
#' The partial log-likelihood (Breslow) is
#' \deqn{\ell(\boldsymbol{\beta}) =
#'   \sum_{g} \Bigl[\sum_{i \in D_g} w_i \eta_i -
#'   d_g^{(w)} \log\Bigl(\sum_{j \in R_g} w_j \exp(\eta_j)\Bigr)\Bigr]}
#' where \eqn{D_g} is the event set at tied event time \eqn{t_g},
#' \eqn{R_g = \{j : t_j \ge t_g\}} is the corresponding risk set, and
#' \eqn{d_g^{(w)} = \sum_{i \in D_g} w_i}.
#'
#' @examples
#' set.seed(1234)
#' eta <- rnorm(50)
#' status <- rbinom(50, 1, 0.6)
#' y <- rexp(50)
#' loglik_cox(eta, status, y)
#'
#' @export
loglik_cox <- function(eta, status, y = NULL, weights = 1) {

  N <- length(eta)
  if(length(weights) == 1) weights <- rep(weights, N)

  ## Weighted hazard contributions
  haz <- weights * exp(eta)

  ## No event-time vector supplied: fall back to no-ties convention
  if(is.null(y)) {
    rsk <- rev(cumsum(rev(haz)))
    return(sum(status * weights * (eta - log(rsk))))
  }

  ## Breslow ties by unique event times
  event_times <- sort(unique(y[status == 1]))
  if(length(event_times) == 0) {
    return(0)
  }

  ll <- 0
  for(tt in event_times) {
    died_g <- which(y == tt & status == 1)
    surv_g <- which(y >= tt)
    sum_w <- sum(weights[died_g])
    ll <- ll + sum(weights[died_g] * eta[died_g]) -
      sum_w * log(sum(haz[surv_g]))
  }

  ll
}


#' Compute Cox Partial Log-Likelihood Score Vector
#'
#' @description
#' Gradient of the Cox partial log-likelihood with respect to
#' \eqn{\boldsymbol{\beta}}. Data must be sorted by ascending event time.
#'
#' @param X Design matrix (N x p), sorted by ascending event time.
#' @param eta Linear predictor vector \eqn{\mathbf{X}\boldsymbol{\beta}}.
#' @param status Event indicator (1 = event, 0 = censored).
#' @param y Optional numeric vector of observed event/censor times, same length
#'   and order as \code{eta}. When supplied, tied event times are handled with
#'   the Breslow approximation. When omitted, the function assumes there are no
#'   ties.
#' @param weights Observation weights (default 1).
#'
#' @return Numeric column vector of length p (gradient).
#'
#' @details
#' Under the Breslow approximation, the score is
#' \deqn{\mathbf{u} =
#'   \sum_g \Bigl[\sum_{i \in D_g} w_i \mathbf{x}_i -
#'   d_g^{(w)} \frac{\sum_{j \in R_g} w_j e^{\eta_j}\mathbf{x}_j}
#'   {\sum_{j \in R_g} w_j e^{\eta_j}}\Bigr]}
#' where \eqn{D_g} is the tied event set at time \eqn{t_g},
#' \eqn{R_g = \{j : t_j \ge t_g\}} is the risk set, and
#' \eqn{d_g^{(w)} = \sum_{i \in D_g} w_i}.
#'
#' @keywords internal
#' @export
score_cox <- function(X, eta, status, y = NULL, weights = 1) {

  N <- length(eta)
  p <- ncol(X)
  if(length(weights) == 1) weights <- rep(weights, N)

  haz <- weights * exp(eta)

  ## No event-time vector supplied: fall back to no-ties convention
  if(is.null(y)) {
    rsk <- rev(cumsum(rev(haz)))
    d_over_rsk <- (weights * status) / rsk
    cum_d_over_rsk <- cumsum(d_over_rsk)
    Pd <- haz * cum_d_over_rsk
    return(crossprod(X, cbind(weights * status - Pd)))
  }

  event_times <- sort(unique(y[status == 1]))
  if(length(event_times) == 0) {
    return(matrix(0, p, 1))
  }

  sc <- matrix(0, p, 1)
  for(tt in event_times) {
    died_g <- which(y == tt & status == 1)
    surv_g <- which(y >= tt)
    sum_w <- sum(weights[died_g])
    score0 <- sum(haz[surv_g])
    score1 <- crossprod(X[surv_g, , drop = FALSE], cbind(haz[surv_g]))
    sc <- sc + crossprod(X[died_g, , drop = FALSE], cbind(weights[died_g])) -
      sum_w * (score1 / score0)
  }

  sc
}


#' Compute Cox Observed Information Matrix
#'
#' @description
#' Negative Hessian of the Cox partial log-likelihood under the Breslow
#' approximation for tied event times. Data must be sorted by ascending event
#' time.
#'
#' @param X Design matrix (N x p), sorted by ascending event time.
#' @param eta Linear predictor vector.
#' @param status Event indicator (1 = event, 0 = censored).
#' @param y Optional numeric vector of observed event/censor times, same length
#'   and order as \code{eta}. When supplied, tied event times are handled with
#'   the Breslow approximation. When omitted, the function assumes there are no
#'   ties.
#' @param weights Observation weights (default 1).
#'
#' @return Symmetric p x p observed information matrix.
#'
#' @details
#' Under the Breslow approximation, the observed information is
#' \deqn{\sum_g d_g^{(w)}
#'   \Bigl[\frac{S_{2g}}{S_{0g}} -
#'   \Bigl(\frac{S_{1g}}{S_{0g}}\Bigr)
#'   \Bigl(\frac{S_{1g}}{S_{0g}}\Bigr)^{\top}\Bigr]}
#' where
#' \eqn{S_{0g} = \sum_{j \in R_g} w_j e^{\eta_j}},
#' \eqn{S_{1g} = \sum_{j \in R_g} w_j e^{\eta_j}\mathbf{x}_j},
#' and \eqn{S_{2g} = \sum_{j \in R_g} w_j e^{\eta_j}\mathbf{x}_j\mathbf{x}_j^{\top}}.
#'
#' @keywords internal
#' @export
info_cox <- function(X, eta, status, y = NULL, weights = 1) {

  N <- length(eta)
  p <- ncol(X)
  if(length(weights) == 1) weights <- rep(weights, N)

  haz <- weights * exp(eta)

  ## No event-time vector supplied: fall back to no-ties convention
  if(is.null(y)) {
    info <- matrix(0, p, p)
    s0 <- sum(haz)
    s1 <- crossprod(X, cbind(haz))
    s2 <- crossprod(X, haz * X)

    for(i in 1:N) {
      if(status[i] == 1 && s0 > 0) {
        xbar <- s1 / s0
        info <- info + weights[i] * (s2 / s0 - tcrossprod(xbar))
      }
      xi <- X[i, ]
      s0 <- s0 - haz[i]
      s1 <- s1 - haz[i] * xi
      s2 <- s2 - haz[i] * tcrossprod(xi)
    }
    return(info)
  }

  event_times <- sort(unique(y[status == 1]))
  if(length(event_times) == 0) {
    return(matrix(0, p, p))
  }

  info <- matrix(0, p, p)
  for(tt in event_times) {
    died_g <- which(y == tt & status == 1)
    surv_g <- which(y >= tt)
    sum_w <- sum(weights[died_g])

    s0 <- sum(haz[surv_g])
    s1 <- crossprod(X[surv_g, , drop = FALSE],
                    cbind(haz[surv_g]))
    s2 <- crossprod(X[surv_g, , drop = FALSE],
                    haz[surv_g] * X[surv_g, , drop = FALSE])

    xbar <- s1 / s0
    info <- info + sum_w * (s2 / s0 - tcrossprod(xbar))
  }

  info
}


#' Cox Proportional Hazards Family for lgspline
#'
#' @description
#' Creates a family-like object for Cox PH models. The link function is
#' \code{log} (the linear predictor is log-relative-hazard), but unlike
#' standard GLM families there is no dispersion parameter and no closed-form
#' mean-variance relationship.
#'
#' The family provides \code{$loglik} and \code{$aic} methods compatible
#' with \code{\link{logLik.lgspline}}.
#'
#' @return A list with family components used by \code{lgspline}.
#'
#' @details
#' Cox PH is semiparametric: the baseline hazard is unspecified. The
#' partial log-likelihood depends only on the order of event times and the
#' linear predictor \eqn{\eta = \mathbf{X}\boldsymbol{\beta}}. Consequently:
#' \itemize{
#'   \item No dispersion parameter is estimated (\code{sigmasq_tilde} is
#'         fixed at 1).
#'   \item \code{dev.resids} returns martingale-style residuals for GCV
#'         tuning compatibility.
#'   \item The response variable is survival time (positive), and the link
#'         is \code{log}.
#' }
#'
#' @examples
#' fam <- cox_family()
#' fam$family
#' fam$link
#'
#' @export
cox_family <- function() list(
  family = "cox",
  link = "log",
  linkfun = log,
  linkinv = exp,

  ## Variance function: Cox PH has no natural mean-variance relationship
  #  in the GLM sense. Return 1 so that default glm_weight_function
  #  produces unit weights (overridden by cox_glm_weight_function anyway).
  variance = function(mu) rep(1, length(mu)),

  ## dev.resids for GCV tuning compatibility.
  #  Uses -2 * partial-log-likelihood contribution per observation as a
  #  proxy. This is NOT the standard martingale residual, but serves the
  #  same role as family$dev.resids in the penalty tuning loop.
  dev.resids = function(y, mu, wt) {
    ## Fallback: squared residual on log scale
    wt * (log(y) - log(mu))^2
  },

  ## Custom deviance used internally by tune_Lambda for GCV.
  #  Receives the full vectors in partition order and the status vector
  #  via .... Returns per-observation deviance contributions.
  custom_dev.resids = function(y, mu, order_indices, family,
                               observation_weights, status) {

    ## Sort by event time
    ord <- order(y)
    y_s <- y[ord]
    mu_s <- mu[ord]
    st_s <- status[order_indices][ord]
    wt_s <- observation_weights[ord]

    eta_s <- log(mu_s)
    haz_s <- wt_s * exp(eta_s)

    pll_i_s <- rep(0, length(y_s))
    event_times <- sort(unique(y_s[st_s == 1]))
    if(length(event_times) > 0) {
      for(tt in event_times) {
        died_g <- which(y_s == tt & st_s == 1)
        surv_g <- which(y_s >= tt)
        pll_i_s[died_g] <- wt_s[died_g] * (eta_s[died_g] -
                           log(sum(haz_s[surv_g])))
      }
    }

    pll_i <- rep(0, length(y))
    pll_i[ord] <- pll_i_s

    ## Return -2 * pll_i as deviance contribution
    -2 * pll_i
  },

  ## Log-likelihood extraction from fitted model
  loglik = function(model_fit) {
    y_fit   <- model_fit$y
    mu_fit  <- model_fit$ytilde
    wt_fit  <- if(!is.null(model_fit$weights)) model_fit$weights
    else rep(1, length(y_fit))

    ## Retrieve status
    status <- NULL
    if(!is.null(model_fit$status)) {
      status <- model_fit$status
    } else if(!is.null(model_fit$.fit_call_args$status)) {
      status <- model_fit$.fit_call_args$status
    } else {
      warning("Censoring status not found; assuming all events observed.")
      status <- rep(1, length(y_fit))
    }

    ## Sort by event time
    ord <- order(y_fit)
    y_s   <- y_fit[ord]
    eta_s <- log(mu_fit[ord])
    st_s  <- status[ord]
    wt_s  <- wt_fit[ord]

    loglik_cox(eta_s, st_s, y_s, wt_s)
  },

  ## AIC = -2*loglik + 2*edf (no +1 for scale, unlike Weibull)
  aic = function(model_fit) {
    ll <- cox_family()$loglik(model_fit)
    edf <- model_fit$trace_XUGX
    if(is.null(edf)) {
      warning("Cannot compute AIC: effective degrees of freedom not found.")
      return(NA)
    }
    -2 * ll + 2 * edf
  }
)


#' Cox PH GLM Weight Function
#'
#' @description
#' Computes working weights for the Cox PH information matrix, used by
#' lgspline when updating \eqn{\mathbf{G}} after obtaining constrained
#' estimates. The weights are a diagonal approximation built from the Breslow
#' tied-event information contributions.
#'
#' @param mu Predicted values (exp(eta), i.e., relative hazard).
#' @param y Observed survival times.
#' @param order_indices Observation indices in partition order.
#' @param family Cox family object (unused, for interface compatibility).
#' @param dispersion Dispersion parameter (fixed at 1 for Cox PH).
#' @param observation_weights Observation weights.
#' @param status Event indicator (1 = event, 0 = censored).
#'
#' @return Numeric vector of working weights, length N.
#'
#' @details
#' For a tied event-time block \eqn{g}, the diagonal approximation uses
#' \deqn{W_{jj}^{(g)} = d_g^{(w)} \frac{h_j}{S_g}
#'   \Bigl(1 - \frac{h_j}{S_g}\Bigr), \qquad j \in R_g}
#' where \eqn{h_j = w_j \exp(\eta_j)},
#' \eqn{S_g = \sum_{k \in R_g} h_k}, and
#' \eqn{d_g^{(w)} = \sum_{i \in D_g} w_i}.
#'
#' When the natural weights are degenerate (all zero or non-finite), falls
#' back to a vector of ones.
#'
#' @examples
#' ## Used internally by lgspline; see cox_helpers examples below.
#'
#' @export
cox_glm_weight_function <- function(mu,
                                    y,
                                    order_indices,
                                    family,
                                    dispersion,
                                    observation_weights,
                                    status) {

  N <- length(mu)
  if(length(observation_weights) == 1) {
    observation_weights <- rep(observation_weights, N)
  }

  ## Sort by event time
  ord <- order(y)
  y_s   <- y[ord]
  mu_s  <- mu[ord]
  st_s  <- status[order_indices][ord]
  wt_s  <- observation_weights[ord]

  haz_s <- wt_s * mu_s
  W_diag_s <- rep(0, N)

  event_times <- sort(unique(y_s[st_s == 1]))
  if(length(event_times) > 0) {
    for(tt in event_times) {
      died_g <- which(y_s == tt & st_s == 1)
      surv_g <- which(y_s >= tt)
      sum_w <- sum(wt_s[died_g])
      s0 <- sum(haz_s[surv_g])
      frac <- haz_s[surv_g] / s0
      W_diag_s[surv_g] <- W_diag_s[surv_g] + sum_w * frac * (1 - frac)
    }
  }

  ## Unsort back to original order
  W_diag <- rep(0, N)
  W_diag[ord] <- W_diag_s

  ## Fallback for degenerate cases
  if(any(!is.finite(W_diag)) || all(W_diag == 0)) {
    return(rep(1, N))
  }

  ## Floor at machine epsilon to avoid zero weights
  W_diag <- pmax(W_diag, .Machine$double.eps)
  return(W_diag)
}


#' Cox PH Dispersion Function
#'
#' @description
#' Returns 1 unconditionally. Cox PH has no dispersion parameter; this
#' function exists solely for interface compatibility with lgspline's
#' \code{dispersion_function} argument.
#'
#' @param mu Predicted values.
#' @param y Observed survival times.
#' @param order_indices Observation indices.
#' @param family Family object.
#' @param observation_weights Observation weights.
#' @param VhalfInv Inverse square root of correlation matrix.
#' @param ... Additional arguments (including \code{status}).
#'
#' @return Scalar 1.
#'
#' @export
cox_dispersion_function <- function(mu, y, order_indices, family,
                                    observation_weights, VhalfInv, ...) {
  1
}


#' Cox PH Score Function for Quadratic Programming and Blockfit
#'
#' @description
#' Computes the score (gradient of partial log-likelihood) in the format
#' expected by lgspline's \code{qp_score_function} interface. The block-
#' diagonal design matrix \eqn{\mathbf{X}} and response \eqn{\mathbf{y}} are
#' in partition order; this function internally sorts by event time, computes
#' the Cox score using the Breslow approximation for tied event times, and
#' returns the result in the original partition order.
#'
#' @param X Block-diagonal design matrix (N x P).
#' @param y Response vector (survival times, N x 1).
#' @param mu Predicted values (N x 1), same order as X and y.
#' @param order_list List of observation indices per partition.
#' @param dispersion Dispersion (fixed at 1).
#' @param VhalfInv Inverse square root correlation matrix (NULL for
#'   independent observations).
#' @param observation_weights Observation weights.
#' @param status Event indicator (1 = event, 0 = censored).
#'
#' @return Numeric column vector of length P (gradient w.r.t. coefficients).
#'
#' @export
cox_qp_score_function <- function(X, y, mu, order_list, dispersion,
                                  VhalfInv, observation_weights, status) {

  N <- length(y)
  if(length(observation_weights) == 1) {
    observation_weights <- rep(observation_weights, N)
  }
  order_indices <- unlist(order_list)

  ## Sort by event time
  ord <- order(y)
  X_s <- X[ord, , drop = FALSE]
  y_s <- y[ord]
  mu_s <- mu[ord]
  st_s <- status[order_indices][ord]
  wt_s <- observation_weights[ord]

  eta_s <- log(mu_s)
  score_cox(X_s, eta_s, st_s, y_s, wt_s)
}


#' Cox PH Schur Correction
#'
#' @description
#' Returns zero corrections for all partitions. Cox PH has no nuisance
#' dispersion parameter, so no Schur complement correction to the
#' information matrix is needed. This function exists for interface
#' compatibility with lgspline's \code{schur_correction_function}.
#'
#' @param X List of partition design matrices.
#' @param y List of partition response vectors.
#' @param B List of partition coefficient vectors.
#' @param dispersion Dispersion (fixed at 1).
#' @param order_list List of observation indices per partition.
#' @param K Number of knots.
#' @param family Family object.
#' @param observation_weights Observation weights.
#' @param ... Additional arguments.
#'
#' @return List of K+1 zeros.
#'
#' @export
cox_schur_correction <- function(X, y, B, dispersion, order_list, K,
                                 family, observation_weights, ...) {
  lapply(1:(K + 1), function(k) 0)
}


#' Unconstrained Cox PH Estimation for lgspline
#'
#' @description
#' Estimates penalized Cox PH coefficients for a single partition (or the
#' full model when K = 0) using damped Newton-Raphson on the penalized
#' partial log-likelihood with the Breslow approximation for tied event times.
#'
#' This function is passed to lgspline's \code{unconstrained_fit_fxn}
#' argument. It receives the partition-level design matrix and response,
#' sorts internally by event time, and returns the coefficient vector.
#'
#' @param X Design matrix (N_k x p) for partition k.
#' @param y Survival times for partition k.
#' @param LambdaHalf Square root of penalty matrix (unused here, retained
#'   for interface compatibility).
#' @param Lambda Penalty matrix.
#' @param keep_weighted_Lambda Logical; if TRUE, return hot-start estimates
#'   without refinement (not recommended for Cox).
#' @param family Cox family object.
#' @param tol Convergence tolerance.
#' @param K Number of knots.
#' @param parallel Logical for parallel processing.
#' @param cl Cluster object.
#' @param chunk_size, num_chunks, rem_chunks Parallelism parameters.
#' @param order_indices Observation indices mapping partition to full data.
#' @param weights Observation weights.
#' @param status Event indicator (1 = event, 0 = censored), full-data length.
#'
#' @return Numeric column vector of penalized partial-likelihood coefficient
#'   estimates.
#'
#' @details
#' The penalized partial log-likelihood is
#' \deqn{\ell_p(\boldsymbol{\beta}) = \ell(\boldsymbol{\beta})
#'   - \tfrac{1}{2}\boldsymbol{\beta}^{\top}
#'   \boldsymbol{\Lambda}\boldsymbol{\beta}}
#' Newton-Raphson updates use the Cox score and observed information
#' (computed via \code{score_cox} and \code{info_cox}), plus the penalty
#' term. The step is damped: the step size is halved until the penalized
#' log-likelihood improves.
#'
#' @examples
#' ## Used internally by lgspline; see the full-model example below.
#'
#' @keywords internal
#' @export
unconstrained_fit_cox <- function(X, y, LambdaHalf, Lambda,
                                  keep_weighted_Lambda, family,
                                  tol = 1e-8, K, parallel, cl,
                                  chunk_size, num_chunks, rem_chunks,
                                  order_indices, weights, status) {

  ## Handle empty partitions
  if(nrow(X) == 0) {
    return(cbind(rep(0, ncol(X))))
  }

  N <- nrow(X)
  p <- ncol(X)

  ## Observation weights
  if(any(!is.null(weights))) {
    weights <- c(weights)
  } else {
    weights <- rep(1, N)
  }

  ## Extract status for this partition
  st <- status[order_indices]

  ## Sort by event time within this partition
  ord <- order(y)
  X_s <- X[ord, , drop = FALSE]
  y_s <- y[ord]
  st_s <- st[ord]
  wt_s <- weights[ord]

  ## Initialize at zero (log-relative-hazard = 0 is the null model)
  init <- rep(0, p)

  ## Damped Newton-Raphson on penalized partial log-likelihood
  beta <- damped_newton_r(
    init,
    ## Penalized partial log-likelihood
    function(par) {
      eta <- c(X_s %**% cbind(par))
      loglik_cox(eta, st_s, y_s, wt_s) -
        0.5 * c(crossprod(par, Lambda) %**% cbind(par))
    },
    ## Score (gradient)
    function(par) {
      eta <- c(X_s %**% cbind(par))
      sc <- score_cox(X_s, eta, st_s, y_s, wt_s)
      c(sc - Lambda %**% cbind(par))
    },
    ## Information (negative Hessian)
    function(par) {
      eta <- c(X_s %**% cbind(par))
      info_cox(X_s, eta, st_s, y_s, wt_s) + Lambda
    },
    tol
  )

  cbind(beta)
}


#' Fit Cox Proportional Hazards Model via lgspline
#'
#' @description
#' Convenience wrapper that calls \code{lgspline} with the correct
#' family, weight, dispersion, score, and unconstrained-fit functions for
#' Cox proportional hazards regression. All standard lgspline arguments
#' (knots, penalties, constraints, parallel, etc.) are passed through.
#'
#' @param formula Formula specifying the model. The response should be
#'   survival time; \code{status} is passed separately.
#' @param data Data frame.
#' @param status Integer vector of event indicators (1 = event,
#'   0 = censored), same length as the number of rows in \code{data}.
#' @param ... Additional arguments passed to \code{lgspline}.
#'
#' @return An object of class \code{"lgspline"}.
#'
#' @details
#' Internally sets:
#' \itemize{
#'   \item \code{family = cox_family()}
#'   \item \code{unconstrained_fit_fxn = unconstrained_fit_cox}
#'   \item \code{glm_weight_function = cox_glm_weight_function}
#'   \item \code{qp_score_function = cox_qp_score_function}
#'   \item \code{dispersion_function = cox_dispersion_function}
#'   \item \code{schur_correction_function = cox_schur_correction}
#'   \item \code{need_dispersion_for_estimation = FALSE}
#'   \item \code{estimate_dispersion = FALSE}
#'   \item \code{standardize_response = FALSE}
#' }
#'
#' A formula interface is needed, e.g. \code{lgspline_cox(t, y, ...)} won't work,
#' unlike for ordinary \code{\link{lgspline}}.
#'
#' @examples
#' \donttest{
#' ## Cox PH with a nonlinear age effect on lung cancer survival
#' if(requireNamespace("survival", quietly = TRUE)) {
#'   library(survival)
#'   set.seed(1234)
#'   lung <- na.omit(lung[, c("time", "status", "age")])
#'   lung$age_std <- std(lung$age)
#'
#'   ## survival codes status as 1 = censored, 2 = dead
#'   event <- as.integer(lung$status == 2)
#'
#'   ## Spline on age
#'   fit <- lgspline_cox(
#'     time ~ spl(age_std),
#'     data = lung,
#'     status = event,
#'     K = 1
#'   )
#'   print(summary(fit))
#'   plot(fit,
#'        show_formulas = TRUE,
#'        custom_response_lab = 'HR',
#'        custom_predictor_lab = 'Standardized Age',
#'        ylim = c(0, 5))
#'
#' }
#' }
#'
#' @export
lgspline_cox <- function(formula, data, status, ...) {
  lgspline(
    formula = formula,
    data = data,
    family = cox_family(),
    unconstrained_fit_fxn = unconstrained_fit_cox,
    glm_weight_function = cox_glm_weight_function,
    qp_score_function = cox_qp_score_function,
    dispersion_function = cox_dispersion_function,
    schur_correction_function = cox_schur_correction,
    need_dispersion_for_estimation = FALSE,
    estimate_dispersion = FALSE,
    standardize_response = FALSE,
    status = status,
    ...
  )
}
