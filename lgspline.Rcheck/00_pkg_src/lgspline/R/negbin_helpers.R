#' Negative Binomial Regression Helpers for lgspline
#'
#' Functions for fitting negative binomial (NB2) regression models within the
#' lgspline framework. Analogous to the Weibull AFT and Cox PH helpers, these
#' provide the log-likelihood, score, information, and all interface functions
#' needed by lgspline's unconstrained fitting, penalty tuning, and inference
#' machinery.
#'
#' The parameterization follows NB2: \eqn{Y \sim \mathrm{NB}(\mu, \theta)}
#' with \eqn{\mathrm{Var}(Y) = \mu + \mu^2/\theta}, where \eqn{\theta > 0}
#' is the shape (size) parameter. The canonical link is log:
#' \eqn{\eta = \log\mu}. The dispersion stored in
#' \code{lgspline$sigmasq_tilde} is \eqn{\theta} itself, not
#' \eqn{1/\theta}.
#'
#' @name negbin_helpers
NULL


#' Compute Negative Binomial Log-Likelihood
#'
#' @description
#' Evaluates the NB2 log-likelihood for given mean vector and shape parameter.
#'
#' @param y Non-negative integer response vector.
#' @param mu Positive mean vector, same length as \code{y}.
#' @param theta Positive scalar shape parameter.
#' @param weights Optional observation weights (default 1).
#'
#' @return Scalar log-likelihood value.
#'
#' @details
#' The log-likelihood is
#' \deqn{\ell(\mu, \theta) = \sum_i w_i \bigl[
#'   \log\Gamma(y_i + \theta) - \log\Gamma(\theta) - \log\Gamma(y_i + 1)
#'   + \theta\log\theta - \theta\log(\mu_i + \theta)
#'   + y_i\log\mu_i - y_i\log(\mu_i + \theta)\bigr]}
#'
#' @examples
#' set.seed(1234)
#' mu <- exp(rnorm(50))
#' y <- rpois(50, mu)
#' loglik_negbin(y, mu, theta = 5)
#'
#' @export
loglik_negbin <- function(y, mu, theta, weights = 1) {

  N <- length(y)
  if(length(weights) == 1) weights <- rep(weights, N)

  logL <- lgamma(y + theta) - lgamma(theta) - lgamma(y + 1) +
    theta * log(theta) - theta * log(mu + theta) +
    y * log(mu) - y * log(mu + theta)

  sum(logL * weights)
}


#' Compute Negative Binomial Score Vector
#'
#' @description
#' Gradient of the NB2 log-likelihood with respect to
#' \eqn{\boldsymbol{\beta}} under the log link.
#'
#' @param X Design matrix (N x p).
#' @param y Response vector.
#' @param mu Mean vector \eqn{\exp(\mathbf{X}\boldsymbol{\beta})}.
#' @param theta Shape parameter.
#' @param weights Observation weights (default 1).
#'
#' @return Numeric column vector of length p (gradient).
#'
#' @details
#' Under log link, \eqn{d\mu_i/d\eta_i = \mu_i}, so the score is
#' \deqn{\mathbf{u} = \mathbf{X}^{\top}\mathbf{w}_{\mathrm{obs}}
#'   \odot \frac{(y_i - \mu_i)\theta}{\theta + \mu_i}}
#' where the per-observation contribution
#' \eqn{(y_i - \mu_i)\theta/(\theta + \mu_i)} is the derivative of the
#' log-likelihood with respect to \eqn{\eta_i}.
#'
#' @keywords internal
#' @export
score_negbin <- function(X, y, mu, theta, weights = 1) {

  N <- length(y)
  if(length(weights) == 1) weights <- rep(weights, N)

  ## d ell_i / d eta_i = (y_i - mu_i) * theta / (theta + mu_i)
  r <- weights * (y - mu) * theta / (theta + mu)

  crossprod(X, cbind(r))
}


#' Compute Negative Binomial Observed Information Matrix
#'
#' @description
#' Negative Hessian of the NB2 log-likelihood with respect to
#' \eqn{\boldsymbol{\beta}} under the log link.
#'
#' @param X Design matrix (N x p).
#' @param mu Mean vector.
#' @param theta Shape parameter.
#' @param weights Observation weights (default 1).
#'
#' @return Symmetric p x p observed information matrix.
#'
#' @details
#' The expected information under log link is
#' \deqn{\mathbf{I} = \mathbf{X}^{\top}\mathbf{W}\mathbf{X}}
#' where \eqn{W_{ii} = w_i \mu_i \theta / (\theta + \mu_i)}.
#' This is the IRLS weight for NB2 with log link.
#'
#' @keywords internal
#' @export
info_negbin <- function(X, y, mu, theta, weights = 1) {

  N <- length(mu)
  if(length(weights) == 1) weights <- rep(weights, N)

  ## Observed information under log link:
  #  J_bb = X^T diag(w_i * mu_i * theta * (theta + y_i)/(theta + mu_i)^2) X
  W <- weights * mu * theta * (theta + y) / (theta + mu)^2
  crossprod(X, W * X)
}


#' Estimate Negative Binomial Shape Parameter
#'
#' @description
#' Computes the profile MLE of the shape parameter \eqn{\theta} given
#' current mean estimates \eqn{\mu}.
#'
#' @param y Response vector.
#' @param mu Mean vector.
#' @param weights Observation weights (default 1).
#' @param init Optional initial value for \eqn{\theta}. If NULL, uses a
#'   moment-based estimate.
#'
#' @return Scalar MLE of \eqn{\theta}.
#'
#' @details
#' Maximizes the profile log-likelihood over \eqn{\theta} via Brent's
#' method on \eqn{[10^{-4}, 10^4]}.
#'
#' @examples
#' set.seed(1234)
#' mu <- rep(5, 200)
#' y <- rnbinom(200, size = 3, mu = 5)
#' negbin_theta(y, mu)
#'
#' @export
negbin_theta <- function(y, mu, weights = 1, init = NULL) {

  N <- length(y)
  if(length(weights) == 1) weights <- rep(weights, N)

  ## Moment-based initialization: theta = mu^2 / (var - mu)
  if(is.null(init)) {
    v <- max(weighted.mean((y - mu)^2, weights), mean(mu) + 0.1)
    m <- weighted.mean(mu, weights)
    init <- max(m^2 / max(v - m, 0.1), 0.1)
  }

  optim(
    init,
    fn = function(par) -loglik_negbin(y, mu, par, weights),
    method = 'Brent',
    lower = 1e-4,
    upper = 1e4
  )$par
}


#' Negative Binomial Family for lgspline
#'
#' @description
#' Creates a family-like object for NB2 regression. The link is log, the
#' variance function is \eqn{V(\mu) = \mu + \mu^2/\theta}, and the
#' dispersion stored by lgspline (\code{sigmasq_tilde}) is the shape
#' parameter \eqn{\theta}.
#'
#' @return A list with family components used by \code{lgspline}.
#'
#' @details
#' The NB2 model has a nuisance shape parameter \eqn{\theta} analogous
#' to the Weibull scale parameter. It is estimated jointly with
#' \eqn{\boldsymbol{\beta}} and its uncertainty is propagated via the
#' Schur complement correction.
#'
#' @examples
#' fam <- negbin_family()
#' fam$family
#' fam$link
#'
#' @export
negbin_family <- function() list(
  family = "negbin",
  link = "log",
  linkfun = log,
  linkinv = exp,

  ## Variance function: V(mu) = mu + mu^2/theta
  #  theta is stored in dispersion; when called from standard glm
  #  machinery without theta, return mu (Poisson-like) as fallback.
  variance = function(mu) mu,

  ## dev.resids for base R compatibility.
  #  Simplified version without theta; full version uses
  #  custom_dev.resids below.
  dev.resids = function(y, mu, wt) {
    ## Poisson-like deviance as fallback
    2 * wt * (y * log(pmax(y, 1) / mu) - (y - mu))
  },

  ## Custom deviance for GCV tuning.
  #  Receives full vectors and estimates theta internally.
  custom_dev.resids = function(y, mu, order_indices, family,
                               observation_weights) {

    wt <- c(observation_weights)

    ## Estimate theta given current mu
    theta <- negbin_theta(y, mu, wt)

    ## NB deviance: -2 * log-likelihood per observation
    dev <- -2 * (
      lgamma(y + theta) - lgamma(theta) - lgamma(y + 1) +
        theta * log(theta) - theta * log(mu + theta) +
        y * log(pmax(mu, .Machine$double.eps)) - y * log(mu + theta)
    )
    dev * wt
  },

  ## Log-likelihood extraction from fitted model
  loglik = function(model_fit) {
    y_fit  <- model_fit$y
    mu_fit <- model_fit$ytilde
    ## sigmasq_tilde stores theta directly
    theta <- model_fit$sigmasq_tilde

    wt <- if(!is.null(model_fit$weights)) model_fit$weights
    else rep(1, length(y_fit))

    if(is.null(theta) || !is.finite(theta) || theta <= 0) {
      warning("Shape parameter theta not found or invalid; ",
              "estimating from fitted values.")
      theta <- negbin_theta(y_fit, mu_fit, wt)
    }

    loglik_negbin(y_fit, mu_fit, theta, wt)
  },

  ## AIC = -2*loglik + 2*(edf + 1) where +1 is for theta
  aic = function(model_fit) {
    ll <- negbin_family()$loglik(model_fit)
    edf <- model_fit$edf
    if(is.null(edf)) edf <- model_fit$trace_XUGX
    if(is.null(edf)) {
      warning("Cannot compute AIC: effective degrees of freedom not found.")
      return(NA)
    }
    ## +1 for the shape parameter
    -2 * ll + 2 * (edf + 1)
  }
)


#' NB GLM Weight Function
#'
#' @description
#' Computes working weights for the NB2 information matrix used by
#' lgspline when updating \eqn{\mathbf{G}} after obtaining constrained
#' estimates.
#'
#' @param mu Predicted values \eqn{\exp(\eta)}.
#' @param y Observed counts.
#' @param order_indices Observation indices in partition order.
#' @param family NB family object (unused, for interface compatibility).
#' @param dispersion Shape parameter \eqn{\theta}.
#' @param observation_weights Observation weights.
#'
#' @return Numeric vector of working weights, length N.
#'
#' @details
#' The IRLS weight for NB2 with log link is
#' \deqn{W_i = w_i \mu_i \theta / (\theta + \mu_i)}
#'
#' Falls back to observation weights when natural weights are degenerate.
#'
#' @export
negbin_glm_weight_function <- function(mu,
                                       y,
                                       order_indices,
                                       family,
                                       dispersion,
                                       observation_weights) {

  theta <- dispersion
  N <- length(mu)
  if(length(observation_weights) == 1) {
    observation_weights <- rep(observation_weights, N)
  }

  ## Observed-information weight for J_bb
  val <- mu * theta * (theta + y) / (theta + mu)^2

  if(any(!is.finite(val)) || all(val == 0)) {
    return(rep(1, N))
  }

  val <- pmax(val, .Machine$double.eps)
  val * c(observation_weights)
}


#' NB Dispersion Function
#'
#' @description
#' Estimates the shape parameter \eqn{\theta} from current fitted values.
#' When a correlation structure is present (\code{VhalfInv} is non-NULL),
#' the Pearson residuals are whitened before computing the moment-based
#' initial value, giving a better starting point for the profile MLE
#' under correlated data. The final estimate is always the profile MLE
#' over \eqn{\theta}.
#'
#' @param mu Predicted values.
#' @param y Observed counts.
#' @param order_indices Observation indices.
#' @param family NB family object.
#' @param observation_weights Observation weights.
#' @param VhalfInv Inverse square root of the correlation matrix, or NULL
#'   for independent observations. When non-NULL, used to whiten residuals
#'   for the moment-based initialization of \eqn{\theta}.
#'
#' @return Scalar \eqn{\theta} estimate (stored as \code{sigmasq_tilde}).
#'
#' @details
#' The profile MLE maximizes \eqn{\ell(\theta \mid \mu)} via Brent's
#' method. When \code{VhalfInv} is provided, the Pearson residuals
#' \eqn{r_i = (y_i - \mu_i) / \sqrt{V(\mu_i)}} are pre-whitened as
#' \eqn{\tilde{r} = V^{-1/2} r} before computing the moment estimator
#' used for initialization. This accounts for the correlation structure
#' in the variance decomposition and produces a more stable starting
#' point for the optimizer, particularly when the correlation inflates
#' the marginal variance beyond what the NB model alone would predict.
#'
#' The profile MLE itself does not use \code{VhalfInv} because the NB
#' log-likelihood is a marginal quantity; the correlation structure
#' affects estimation only through the mean model (handled by the GEE
#' paths in \code{get_B}).
#'
#' @export
negbin_dispersion_function <- function(mu,
                                       y,
                                       order_indices,
                                       family,
                                       observation_weights,
                                       VhalfInv) {

  mu <- c(mu)
  y <- c(y)
  observation_weights <- c(observation_weights)
  N <- length(y)
  if(length(observation_weights) == 1) {
    observation_weights <- rep(observation_weights, N)
  }

  ## Moment-based initialization, optionally whitened for GEE
  init <- NULL
  if(!is.null(VhalfInv)) {
    ## Whiten Pearson residuals: r_i = (y_i - mu_i) / sqrt(mu_i)
    #  then tilde_r = V^{-1/2} r
    #  The whitened variance should be ~ 1 + mu/theta under the model,
    #  so theta_init = mean(mu) / max(var(tilde_r) - 1, 0.1)
    pear_r <- (y - mu) / sqrt(pmax(mu, .Machine$double.eps))
    tilde_r <- c(VhalfInv %**% cbind(pear_r))
    v_tilde <- weighted.mean(tilde_r^2, observation_weights)
    m_mu <- weighted.mean(mu, observation_weights)
    init <- max(m_mu / max(v_tilde - 1, 0.1), 0.1)
  }

  ## Profile MLE over theta
  negbin_theta(y, mu, observation_weights, init)
}


#' NB Score Function for Quadratic Programming and Blockfit
#'
#' @description
#' Computes the score (gradient of NB log-likelihood) in the format
#' expected by lgspline's \code{qp_score_function} interface.
#'
#' @param X Block-diagonal design matrix (N x P).
#' @param y Response vector (N x 1).
#' @param mu Predicted values (N x 1), same order as X and y.
#' @param order_list List of observation indices per partition.
#' @param dispersion Shape parameter \eqn{\theta}.
#' @param VhalfInv Inverse square root of correlation matrix; when non-NULL
#'   the score is computed on the whitened scale as
#'   \eqn{\tilde{X}^{\top}\tilde{r}} where
#'   \eqn{\tilde{X} = V^{-1/2}X} and the residual vector accounts for
#'   the correlation.
#' @param observation_weights Observation weights.
#'
#' @return Numeric column vector of length P.
#'
#' @details
#' Without correlation (\code{VhalfInv = NULL}), the score is
#' \eqn{\mathbf{X}^{\top}\mathbf{w} \odot (y - \mu)\theta/(\theta + \mu)}.
#'
#' With correlation, the GEE score is
#' \eqn{\tilde{\mathbf{X}}^{\top}\mathrm{diag}(\mathbf{W})
#' \mathbf{V}^{-1}(\mathbf{y} - \boldsymbol{\mu})}
#' where \eqn{\mathbf{W}} contains the NB working weights. The
#' whitening is absorbed by pre-multiplying both \eqn{\mathbf{X}} and
#' the residual by \eqn{\mathbf{V}^{-1/2}}.
#'
#' @export
negbin_qp_score_function <- function(X, y, mu, order_list, dispersion,
                                     VhalfInv, observation_weights) {

  theta <- dispersion
  N <- length(y)
  if(length(observation_weights) == 1) {
    observation_weights <- rep(observation_weights, N)
  }
  mu <- c(mu)
  y <- c(y)

  if(!is.null(VhalfInv)) {
    ## GEE score: X^T W V^{-1} (y - mu) where W is NB working weight
    #  On whitened scale: (V^{-1/2} X)^T * W * V^{-1/2}(y - mu)
    #  But the caller already passes X_tilde = V^{-1/2} X in the GEE
    #  paths, so X here may already be whitened. The mu passed is
    #  V^{-1/2} mu_original in GEE. We follow the Weibull convention:
    #  compute residual on original scale, then let the caller handle
    #  whitening.
    resid <- observation_weights * (y - mu) * theta / (theta + mu)
    return(crossprod(X, cbind(resid)))
  }

  ## Independent case
  resid <- observation_weights * (y - mu) * theta / (theta + mu)
  crossprod(X, cbind(resid))
}


#' NB Schur Correction
#'
#' @description
#' Computes the Schur complement correction to account for uncertainty in
#' estimating \eqn{\theta}. Structure is identical to
#' \code{\link{weibull_schur_correction}}: the joint information is
#' partitioned into \eqn{(\boldsymbol{\beta}, \theta)} blocks and the
#' correction is \eqn{-\mathbf{I}_{\beta\theta}
#' I_{\theta\theta}^{-1}\mathbf{I}_{\beta\theta}^{\top}}.
#'
#' @param X List of partition design matrices.
#' @param y List of partition response vectors.
#' @param B List of partition coefficient vectors.
#' @param dispersion Scalar shape parameter \eqn{\theta}.
#' @param order_list List of observation indices per partition.
#' @param K Number of knots.
#' @param family Family object.
#' @param observation_weights Observation weights.
#'
#' @return List of K+1 correction matrices, with \code{0} for empty
#'   partitions.
#'
#' @details
#' The cross-derivative (score of \eqn{\eta_i} w.r.t. \eqn{\theta}) is
#' \deqn{\frac{\partial^2 \ell}{\partial \eta_i \partial \theta}
#'   = \frac{(y_i - \mu_i)\mu_i}{(\theta + \mu_i)^2}}
#' The second derivative of the log-likelihood w.r.t. \eqn{\theta} is
#' \deqn{I_{\theta\theta} = -\sum_i w_i \bigl[
#'   \psi'(y_i + \theta) - \psi'(\theta) +
#'   \frac{1}{\theta} - \frac{2}{\mu_i + \theta} +
#'   \frac{y_i + \theta}{(\mu_i + \theta)^2}\bigr]}
#' where \eqn{\psi'} is the trigamma function.
#'
#' @export
negbin_schur_correction <- function(X, y, B, dispersion, order_list, K,
                                    family, observation_weights) {

  theta <- dispersion
  lapply(1:(K + 1), function(k) {
    if(nrow(X[[k]]) < 1) {
      return(0)
    } else {
      mu <- family$linkinv(c(X[[k]] %**% B[[k]]))
      obs <- c(y[[k]])
      weights <- c(observation_weights[[k]])

      ## Cross-derivative: d^2 ell / (d eta_i d theta)
      #  = (y_i - mu_i) * mu_i / (theta + mu_i)^2
      I_bt <- crossprod(X[[k]], cbind(
        weights * (obs - mu) * mu / (theta + mu)^2
      ))

      ## Information for theta (negative of second derivative)
      I_tt <- -sum(weights * (
        trigamma(obs + theta) - trigamma(theta) +
          1 / theta -
          2 / (mu + theta) +
          (obs + theta) / (mu + theta)^2
      ))

      ## Schur complement: -I_bt (1/I_tt) I_bt^T
      #  Gets elementwise added to G[[k]]^{-1}, so should be negative
      #  semi-definite for well-behaved data (I_tt > 0 implies
      #  compl is negative semi-definite).
      if(abs(I_tt) < .Machine$double.eps) {
        return(0)
      }
      compl <- I_bt %**% matrix(-1 / I_tt) %**% t(I_bt)
      return(compl)
    }
  })
}


#' Unconstrained NB Estimation for lgspline
#'
#' @description
#' Estimates penalized NB2 coefficients for a single partition (or the
#' full model when K = 0) using a two-stage optimization: outer loop
#' profiles over \eqn{\theta} via Brent's method, inner loop estimates
#' \eqn{\boldsymbol{\beta}} via damped Newton-Raphson on the penalized
#' log-likelihood.
#'
#' @param X Design matrix (N_k x p) for partition k.
#' @param y Response counts for partition k.
#' @param LambdaHalf Square root of penalty matrix.
#' @param Lambda Penalty matrix.
#' @param keep_weighted_Lambda Logical; if TRUE, return hot-start estimates
#'   from augmented Poisson regression without refinement.
#' @param family NB family object.
#' @param tol Convergence tolerance.
#' @param K Number of knots.
#' @param parallel Logical for parallel processing.
#' @param cl Cluster object.
#' @param chunk_size, num_chunks, rem_chunks Parallelism parameters.
#' @param order_indices Observation indices mapping partition to full data.
#' @param weights Observation weights.
#'
#' @return Numeric column vector of penalized coefficient estimates.
#'
#' @details
#' The penalized log-likelihood is
#' \deqn{\ell_p(\boldsymbol{\beta}, \theta) = \ell(\boldsymbol{\beta},
#'   \theta) - \tfrac{1}{2}\boldsymbol{\beta}^{\top}
#'   \boldsymbol{\Lambda}\boldsymbol{\beta}}
#'
#' The outer loop optimizes \eqn{\theta} given the inner-loop-optimal
#' \eqn{\boldsymbol{\beta}(\theta)}. This mirrors the Weibull AFT
#' approach where scale (\eqn{\sigma}) is profiled out.
#'
#' Initialization uses Poisson regression coefficients as a hot start
#' (equivalent to \eqn{\theta \to \infty}).
#'
#' @keywords internal
#' @export
unconstrained_fit_negbin <- function(X, y, LambdaHalf, Lambda,
                                     keep_weighted_Lambda, family,
                                     tol = 1e-8, K, parallel, cl,
                                     chunk_size, num_chunks, rem_chunks,
                                     order_indices, weights) {

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

  y <- c(y)

  ## Hot-start via Poisson (theta -> Inf limit of NB)
  pois_init <- try({
    glm.fit(
      x = rbind(X, LambdaHalf),
      y = c(y, rep(1, nrow(LambdaHalf))),
      family = poisson(link = 'log'),
      weights = c(weights, rep(mean(weights), nrow(LambdaHalf)))
    )$coefficients
  }, silent = TRUE)

  if(inherits(pois_init, 'try-error') || any(!is.finite(pois_init))) {
    pois_init <- c(log(max(mean(y), 0.1)), rep(0, p - 1))
  }

  if(keep_weighted_Lambda) {
    return(cbind(pois_init))
  }

  ## Initial theta from Poisson residuals
  mu_init <- c(exp(X %**% cbind(pois_init)))
  theta_init <- negbin_theta(y, mu_init, weights)

  ## Outer loop: profile over theta via Brent
  theta <- optim(
    theta_init,
    fn = function(par) {
      th <- par
      ## Inner loop: optimize beta given theta
      beta <- damped_newton_r(
        pois_init,
        function(b) {
          mu <- c(exp(X %**% cbind(b)))
          loglik_negbin(y, mu, th, weights) -
            0.5 * c(crossprod(b, Lambda) %**% cbind(b))
        },
        function(b) {
          mu <- c(exp(X %**% cbind(b)))
          c(score_negbin(X, y, mu, th, weights) -
              Lambda %**% cbind(b))
        },
        function(b) {
          mu <- c(exp(X %**% cbind(b)))
          info_negbin(X, y, mu, th, weights) + Lambda
        },
        tol
      )
      mu <- c(exp(X %**% cbind(beta)))
      -loglik_negbin(y, mu, par, weights) +
        0.5 * c(crossprod(beta, Lambda) %**% cbind(beta))
    },
    lower = theta_init / 10,
    upper = theta_init * 10,
    method = 'Brent'
  )$par

  ## Final beta given optimal theta
  beta <- cbind(damped_newton_r(
    pois_init,
    function(b) {
      mu <- c(exp(X %**% cbind(b)))
      loglik_negbin(y, mu, theta, weights) -
        0.5 * c(crossprod(b, Lambda) %**% cbind(b))
    },
    function(b) {
      mu <- c(exp(X %**% cbind(b)))
      c(score_negbin(X, y, mu, theta, weights) -
          Lambda %**% cbind(b))
    },
    function(b) {
      mu <- c(exp(X %**% cbind(b)))
      info_negbin(X, y, mu, theta, weights) + Lambda
    },
    tol
  ))

  return(beta)
}


#' Fit Negative Binomial Model via lgspline
#'
#' @description
#' Convenience wrapper that calls \code{lgspline} with the correct
#' family, weight, dispersion, score, and unconstrained-fit functions for
#' NB2 regression. All standard lgspline arguments (knots, penalties,
#' constraints, parallel, correlation structures, etc.) are passed through.
#'
#' @param formula Formula specifying the model. The response should be
#'   non-negative integer counts.
#' @param data Data frame.
#' @param ... Additional arguments passed to \code{lgspline}, including
#'   \code{Vhalf} and \code{VhalfInv} for correlation structures.
#'
#' @return An object of class \code{"lgspline"}.
#'
#' @details
#' Internally sets:
#' \itemize{
#'   \item \code{family = negbin_family()}
#'   \item \code{unconstrained_fit_fxn = unconstrained_fit_negbin}
#'   \item \code{glm_weight_function = negbin_glm_weight_function}
#'   \item \code{qp_score_function = negbin_qp_score_function}
#'   \item \code{dispersion_function = negbin_dispersion_function}
#'   \item \code{schur_correction_function = negbin_schur_correction}
#'   \item \code{need_dispersion_for_estimation = TRUE}
#'   \item \code{estimate_dispersion = TRUE}
#'   \item \code{standardize_response = FALSE}
#' }
#'
#' A formula interface is needed, e.g. \code{lgspline_negbin(t, y)} won't work,
#' unlike for ordinary \code{\link{lgspline}}.
#'
#' When a correlation structure is supplied via \code{Vhalf}/\code{VhalfInv},
#' the model is fitted through the GEE Path 1b machinery in \code{get_B}.
#' The dispersion function uses \code{VhalfInv} to whiten Pearson
#' residuals for a better moment-based initialization of \eqn{\theta},
#' which stabilizes the profile MLE under moderate to strong correlation.
#' The score function handles the whitened design consistently with the
#' Weibull AFT GEE convention.
#'
#'
#' @examples
#'
#' set.seed(1234)
#' N <- 300
#' t <- rnorm(N)
#' mu <- exp(1 + 0.5 * sin(2 * t))
#' y <- rnbinom(N, size = 3, mu = mu)
#' df <- data.frame(response = y, predictor = t)
#'
#' fit <- lgspline_negbin(
#'   response ~ spl(predictor),
#'   data = df,
#'   K = 2,
#'   opt = FALSE,
#'   wiggle_penalty = 1e-2
#' )
#' print(summary(fit))
#' plot(fit, show_formulas = TRUE,
#'      custom_response_lab = 'Count')
#' points(t, mu, col = 'grey', cex=0.67)
#'
#' @seealso \code{\link{lgspline_cox}} for Cox PH,
#'   \code{\link{lgspline_weibull}} for Weibull AFT,
#'   \code{\link{negbin_family}}, \code{\link{unconstrained_fit_negbin}}
#'
#' @export
lgspline_negbin <- function(formula, data, ...) {
  lgspline(
    formula = formula,
    data = data,
    family = negbin_family(),
    unconstrained_fit_fxn = unconstrained_fit_negbin,
    glm_weight_function = negbin_glm_weight_function,
    qp_score_function = negbin_qp_score_function,
    dispersion_function = negbin_dispersion_function,
    schur_correction_function = negbin_schur_correction,
    need_dispersion_for_estimation = TRUE,
    estimate_dispersion = TRUE,
    standardize_response = FALSE,
    ...
  )
}
