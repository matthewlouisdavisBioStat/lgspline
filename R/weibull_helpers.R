#' Compute Log-Likelihood for Weibull Accelerated Failure Time Model
#'
#' @description
#' Calculates the log-likelihood for a Weibull accelerated failure time (AFT)
#' survival model, supporting right-censored survival data.
#'
#' @param log_y Numeric vector of logarithmic response/survival times
#' @param log_mu Numeric vector of logarithmic predicted survival times
#' @param status Numeric vector of censoring indicators
#'   (1 = event, 0 = censored)
#'   Indicates whether an event of interest occurred (1) or the observation was
#'   right-censored (0). In survival analysis, right-censoring occurs when the
#'   full survival time is unknown, typically because the study ended or the
#'   subject was lost to follow-up before the event of interest occurred.
#' @param scale Numeric scalar representing the Weibull scale parameter
#'   (sigma), equivalent to \code{survreg$scale}. This is the square root of
#'   the dispersion stored in \code{lgspline$sigmasq_tilde}.
#' @param weights Optional numeric vector of observation weights (default = 1)
#'
#' @return
#' A numeric scalar representing the weighted total log-likelihood of the model
#'
#' @details
#' The function computes log-likelihood contributions for a Weibull AFT model,
#' explicitly accounting for right-censored observations. It supports optional
#' observation weighting to accommodate complex sampling designs.
#'
#' This both provides a tool for actually fitting Weibull AFT models, and
#' boilerplate code for users who wish to incorporate Lagrangian multiplier
#' smoothing splines into their own custom models.
#'
#' @examples
#'
#' ## Minimal example of fitting a Weibull Accelerated Failure Time model
#' # Simulating survival data with right-censoring
#' set.seed(1234)
#' x1 <- rnorm(1000)
#' x2 <- rbinom(1000, 1, 0.5)
#' yraw <- rexp(exp(0.01*x1 + 0.01*x2))
#' # status: 1 = event occurred, 0 = right-censored
#' status <- rbinom(1000, 1, 0.25)
#' yobs <- ifelse(status, runif(length(yraw), 0, yraw), yraw)
#' df <- data.frame(
#'   y = yobs,
#'   x1 = x1,
#'   x2 = x2
#' )
#'
#' ## Fit model using lgspline with Weibull AFT specifics
#' model_fit <- lgspline(y ~ spl(x1) + x2,
#'                       df,
#'                       unconstrained_fit_fxn = unconstrained_fit_weibull,
#'                       family = weibull_family(),
#'                       need_dispersion_for_estimation = TRUE,
#'                       dispersion_function = weibull_dispersion_function,
#'                       glm_weight_function = weibull_glm_weight_function,
#'                       schur_correction_function = weibull_schur_correction,
#'                       status = status,
#'                       opt = FALSE,
#'                       K = 1)
#'
#' loglik_weibull(log(model_fit$y), log(model_fit$ytilde), status,
#'   sqrt(model_fit$sigmasq_tilde))
#'
#' @export
loglik_weibull <- function(log_y, log_mu, status, scale, weights = 1) {

  ## Log-likelihood contributions
  ## scale = sigma (Weibull scale), NOT sigma^2 (dispersion)
  z <- (log_y - log_mu) / scale
  logL <- status * (-log(scale) +
                      z -
                      log_y) -
    exp(z)

  ## Return sum of log likelihood
  return(sum(logL * weights))
}


#' Compute Gradient of Log-Likelihood of Weibull Accelerated Failure Model
#'
#' @description
#' Calculates the gradient of log-likelihood for a Weibull accelerated failure
#' time (AFT) survival model, supporting right-censored survival data.
#'
#' @param X Design matrix
#' @param y Response vector
#' @param mu Predicted mean vector
#' @param order_list List of observation indices per partition
#' @param dispersion Dispersion parameter (sigma^2 = scale^2). The lgspline
#'   framework stores and passes dispersion (sigma^2); the Weibull scale
#'   (sigma) matching \code{survreg$scale} is \code{sqrt(dispersion)}.
#' @param VhalfInv Inverse square root of correlation matrix; unused here and
#'   retained for interface compatibility.
#' @param observation_weights Observation weights
#' @param status Censoring indicator (1 = event, 0 = censored)
#'
#' @return
#' Numeric column vector representing the gradient with respect to
#' coefficients.
#'
#' @details
#' Needed if using "blockfit", correlation structures, or quadratic programming
#' with Weibull AFT models.
#'
#' The gradient is computed on a scale that omits the 1/sigma prefactor.
#' Specifically, the true score is (1/sigma) * X^T diag(w) (exp(z) - status),
#' but both this function and the corresponding information matrix used
#' internally omit 1/sigma and 1/sigma^2 respectively, so the Newton-Raphson
#' step G*u remains correct. This matches the convention in
#' \code{\link{unconstrained_fit_weibull}}.
#'
#' @examples
#'
#' set.seed(1234)
#' t1 <- rnorm(1000)
#' t2 <- rbinom(1000, 1, 0.5)
#' yraw <- rexp(exp(0.01*t1 + 0.01*t2))
#' status <- rbinom(1000, 1, 0.25)
#' yobs <- ifelse(status, runif(length(yraw), 0, yraw), yraw)
#' df <- data.frame(
#'   y = yobs,
#'   t1 = t1,
#'   t2 = t2
#' )
#'
#' ## Example using blockfit for t2 as a linear term - output does not look
#' # different, but internal methods used for fitting change
#' model_fit <- lgspline(y ~ spl(t1) + t2,
#'                       df,
#'                       unconstrained_fit_fxn = unconstrained_fit_weibull,
#'                       family = weibull_family(),
#'                       need_dispersion_for_estimation = TRUE,
#'                       qp_score_function = weibull_qp_score_function,
#'                       dispersion_function = weibull_dispersion_function,
#'                       glm_weight_function = weibull_glm_weight_function,
#'                       schur_correction_function = weibull_schur_correction,
#'                       K = 1,
#'                       blockfit = TRUE,
#'                       opt = FALSE,
#'                       status = status,
#'                       verbose = TRUE)
#'
#' print(summary(model_fit))
#'
#' @export
weibull_qp_score_function = function(X,
                                     y,
                                     mu,
                                     order_list,
                                     dispersion,
                                     VhalfInv,
                                     observation_weights,
                                     status){
  ## dispersion = sigma^2; scale = sigma = sqrt(dispersion)
  scale <- sqrt(dispersion)
  order_indices <- unlist(order_list)
  ## Score of beta (scaled by sigma^2 for numerical convenience).
  # True score: (1/sigma) * X^T diag(w) (exp(z) - status)
  #  Scaled score (x sigma^2): sigma * X^T diag(w) (exp(z) - status)
  #  This matches the convention in unconstrained_fit_weibull, where
  #  the information is also unscaled (missing 1/sigma^2), so the
  #  Newton-Raphson step G*u remains correct.
  crossprod(X , cbind(
    observation_weights *
      (exp((log(y) - log(mu)) / scale) - status[order_indices]) *
      scale)
  )
}

#' Correction for the Variance-Covariance Matrix for Uncertainty in Scale
#'
#' @description
#' Computes the Schur complement \eqn{\textbf{S}} such that
#' \eqn{\textbf{G}^* = (\textbf{G}^{-1} + \textbf{S})^{-1}} properly
#' accounts for uncertainty in estimating the Weibull scale parameter when
#' estimating the variance-covariance matrix. Otherwise, the
#' variance-covariance matrix is optimistic and assumes the scale is known,
#' when it was in fact estimated. Note that the parameterization adds the
#' output of this function elementwise (not subtract) so for most cases, the
#' output of this function will be negative or a negative
#' definite/semi-definite matrix.
#'
#' @param X Block-diagonal matrices of spline expansions
#' @param y Block-vector of response
#' @param B Block-vector of coefficient estimates
#' @param dispersion Scalar, estimate of dispersion (sigma^2 = scale^2). The
#'   lgspline framework stores and passes dispersion (sigma^2); the Weibull
#'   scale (sigma) matching \code{survreg$scale} is \code{sqrt(dispersion)}.
#' @param order_list List of partition orders
#' @param K Number of partitions minus 1 (\eqn{K})
#' @param family Distribution family
#' @param observation_weights Optional observation weights (default = 1)
#' @param status Censoring indicator (1 = event, 0 = censored)
#'   Indicates whether an event of interest occurred (1) or the observation was
#'   right-censored (0). In survival analysis, right-censoring occurs when the
#'   full survival time is unknown, typically because the study ended or the
#'   subject was lost to follow-up before the event of interest occurred.
#'
#' @return
#' List of blockwise Schur-complement corrections \eqn{\textbf{S}_k} to be
#' elementwise added to each block of the information matrix before inversion,
#' with \code{0} returned for empty partitions.
#'
#' @details
#' Adjusts the variance-covariance matrix unscaled for coefficients to account
#' for uncertainty in estimating the Weibull scale parameter, that otherwise
#' would be lost if simply using
#' \eqn{\textbf{G}=(\textbf{X}^{T}\textbf{W}\textbf{X} + \textbf{L})^{-1}}.
#' This is accomplished using a correction based on the Schur complement so we
#' avoid having to construct the entire variance-covariance matrix, or
#' modifying the procedure for \code{\link{lgspline}} substantially.
#' For any model with nuisance parameters that must have uncertainty accounted
#' for, this tool will be helpful.
#'
#' This both provides a tool for actually fitting Weibull accelerated failure
#' time (AFT) models, and boilerplate code for users who wish to incorporate
#' Lagrangian multiplier smoothing splines into their own custom models.
#'
#' @examples
#'
#' ## Minimal example of fitting a Weibull Accelerated Failure Time model
#' # Simulating survival data with right-censoring
#' set.seed(1234)
#' t1 <- rnorm(1000)
#' t2 <- rbinom(1000, 1, 0.5)
#' yraw <- rexp(exp(0.01*t1 + 0.01*t2))
#' # status: 1 = event occurred, 0 = right-censored
#' status <- rbinom(1000, 1, 0.25)
#' yobs <- ifelse(status, runif(length(yraw), 0, yraw), yraw)
#' df <- data.frame(
#'   y = yobs,
#'   t1 = t1,
#'   t2 = t2
#' )
#'
#' ## Fit model using lgspline with Weibull Schur correction
#' model_fit <- lgspline(y ~ spl(t1) + t2,
#'                       df,
#'                       unconstrained_fit_fxn = unconstrained_fit_weibull,
#'                       family = weibull_family(),
#'                       need_dispersion_for_estimation = TRUE,
#'                       dispersion_function = weibull_dispersion_function,
#'                       glm_weight_function = weibull_glm_weight_function,
#'                       schur_correction_function = weibull_schur_correction,
#'                       status = status,
#'                       opt = FALSE,
#'                       K = 1)
#'
#' print(summary(model_fit))
#'
#' ## Fit model using lgspline without Weibull Schur correction
#' naive_fit <- lgspline(y ~ spl(t1) + t2,
#'                       df,
#'                       unconstrained_fit_fxn = unconstrained_fit_weibull,
#'                       family = weibull_family(),
#'                       need_dispersion_for_estimation = TRUE,
#'                       dispersion_function = weibull_dispersion_function,
#'                       glm_weight_function = weibull_glm_weight_function,
#'                       status = status,
#'                       opt = FALSE,
#'                       K = 1)
#'
#' print(summary(naive_fit))
#'
#' @export
weibull_schur_correction <- function(X,
                                     y,
                                     B,
                                     dispersion,
                                     order_list,
                                     K,
                                     family,
                                     observation_weights,
                                     status){
  ## dispersion = sigma^2; scale = sigma = sqrt(dispersion)
  scale <- sqrt(dispersion)
  lapply(1:(K+1), function(k){
    if(nrow(X[[k]]) < 1){
      return(0)
    } else {
      mu <- family$linkinv(c(X[[k]] %**% B[[k]]))
      s <- status[order_list[[k]]]
      obs <- y[[k]]
      z <- (log(obs) - log(mu)) / scale
      exp_z <- exp(z)
      zexp_z <- z * exp_z
      weights <- c(observation_weights[[k]])

      ## Correction via Schur complement
      # Extract true conditional variance-covariance of beta coefficients
      # conditional upon estimate of scale (sigma).
      # I = ( I_bb I_bs^{T} )
      #     ( I_bs I_ss     )
      # for b = beta, s = scale (sigma)
      # Note that: I_bb = invert(G[[k]]) is incorrect, I_bb is part of
      # Fisher info
      I_bs <- t(X[[k]]) %**% cbind(weights * zexp_z * scale)
      I_ss <- -sum(
        weights * (
          (s + 2*s*z + zexp_z + exp_z * z^2)
        )
      )
      # compl gets elementwise added to G[[k]] for all k = 1...K+1
      compl <- I_bs %**% matrix(-1/I_ss) %**% t(I_bs)
      # Schur complement correction to pass on to compute_G_eigen()
      return(compl)
    }
  })
}

#' @rdname weibull_schur_correction
#' @export
weibull_shur_correction <- weibull_schur_correction

#' Estimate Scale for Weibull Accelerated Failure Time Model
#'
#' @description
#' Computes maximum log-likelihood scale estimate (sigma) for a Weibull
#' accelerated failure time (AFT) survival model.
#'
#' This both provides a tool for actually fitting Weibull AFT Models, and
#' boilerplate code for users who wish to incorporate Lagrangian multiplier
#' smoothing splines into their own custom models.
#'
#' @param log_y Logarithm of response/survival times
#' @param log_mu Logarithm of predicted survival times
#' @param status Censoring indicator (1 = event, 0 = censored)
#'   Indicates whether an event of interest occurred (1) or the observation was
#'   right-censored (0). In survival analysis, right-censoring occurs when the
#'   full survival time is unknown, typically because the study ended or the
#'   subject was lost to follow-up before the event of interest occurred.
#' @param weights Optional observation weights (default = 1)
#'
#' @return
#' Scalar representing the estimated Weibull scale (sigma), equivalent to
#' \code{survreg$scale}. The dispersion (as stored in
#' \code{lgspline$sigmasq_tilde}) is sigma^2.
#'
#' @details
#' Calculates maximum log-likelihood estimate of scale (sigma) for Weibull AFT
#' model accounting for right-censored observations using Brent's method for
#' optimization.
#'
#' @examples
#'
#' ## Simulate exponential data with censoring
#' set.seed(1234)
#' mu <- 2  # mean of exponential distribution
#' n <- 500
#' y <- rexp(n, rate = 1/mu)
#'
#' ## Introduce censoring (25% of observations)
#' status <- rbinom(n, 1, 0.75)
#' y_obs <- ifelse(status, y, NA)
#'
#' ## Compute scale estimate
#' scale_est <- weibull_scale(
#'   log_y = log(y_obs[!is.na(y_obs)]),
#'   log_mu = log(mu),
#'   status = status[!is.na(y_obs)]
#' )
#'
#' print(scale_est)
#'
#'
#' @export
weibull_scale <- function(log_y, log_mu, status, weights = 1){
  optim(
    1,
    fn = function(par){
      ## par is scale (sigma), passed directly to loglik_weibull
      -loglik_weibull(log_y, log_mu, status, par, weights)
    },
    method = 'Brent',
    lower = 1e-64,
    upper = 100
  )$par
}

## weibull_family, with $aic, $loglik, $dev.resids methods
#  for logLik.lgspline compatibility

#' Weibull Family for Survival Model Specification
#'
#' @description
#' Creates a family-like object for Weibull accelerated failure time (AFT)
#' models, including custom log-likelihood, AIC, and deviance helpers.
#'
#' This both provides a tool for actually fitting Weibull AFT Models, and
#' boilerplate code for users who wish to incorporate Lagrangian multiplier
#' smoothing splines into their own custom models.
#'
#' @return
#' A family-like list containing the link functions and Weibull-specific
#' methods used by \code{lgspline}.
#'
#' @details
#' Provides a comprehensive family specification for Weibull AFT models,
#' including family name, link function, inverse link function, custom loss
#' function for model tuning, and methods for AIC and log-likelihood
#' computation compatible with \code{logLik.lgspline}.
#'
#' Supports right-censored survival data with flexible parameter estimation.
#'
#' Note on scale vs. dispersion: throughout this package, the lgspline object
#' stores \code{sigmasq_tilde} which equals sigma^2 (dispersion), where
#' sigma is the Weibull scale parameter matching \code{survreg$scale}.
#' Functions that accept a \code{dispersion} argument receive sigma^2;
#' functions that accept a \code{scale} argument receive sigma.
#'
#' @examples
#'
#' ## Simulate survival data with covariates
#' set.seed(1234)
#' n <- 1000
#' t1 <- rnorm(n)
#' t2 <- rbinom(n, 1, 0.5)
#'
#' ## Generate survival times with Weibull-like structure
#' lambda <- exp(0.5 * t1 + 0.3 * t2)
#' yraw <- rexp(n, rate = 1/lambda)
#'
#' ## Introduce right-censoring
#' status <- rbinom(n, 1, 0.75)
#' y <- ifelse(status, yraw, runif(length(yraw), 0, yraw))
#'
#' ## Prepare data
#' df <- data.frame(y = y, t1 = t1, t2 = t2, status = status)
#'
#' ## Fit model using custom Weibull family
#' model_fit <- lgspline(y ~ spl(t1) + t2,
#'                       df,
#'                       unconstrained_fit_fxn = unconstrained_fit_weibull,
#'                       family = weibull_family(),
#'                       need_dispersion_for_estimation = TRUE,
#'                       dispersion_function = weibull_dispersion_function,
#'                       glm_weight_function = weibull_glm_weight_function,
#'                       schur_correction_function = weibull_schur_correction,
#'                       status = status,
#'                       opt = FALSE,
#'                       K = 1)
#'
#' summary(model_fit)
#'
#' ## Log-likelihood now works via logLik.lgspline:
#' # logLik(model_fit)
#'
#' @export
weibull_family <- function() list(
  family = "weibull",
  link = "log",
  linkfun = log,
  linkinv = exp,

  ## Added dev.resids for base R
  #  compatibility when logLik.lgspline checks for family$dev.resids
  dev.resids = function(y, mu, wt) {
    ## Deviance residuals for Weibull AFT: -2 * log-likelihood
    #  This is a simplified version; full version uses status via
    #  custom_dev.resids below
    wt * (y/mu - log(y/mu) - 1) * 2
  },

  ## Custom loss used in place of MSE for computing GCV
  custom_dev.resids =
    function(y, mu, order_indices, family, observation_weights, status){
      log_mu <- log(mu)
      log_y <- log(y)
      status <- status[order_indices]

      ## Initialize scale (sigma) using null model
      init_scale <-
        weibull_scale(log_y, mean(log_y), status,
                      observation_weights)
      ## Find scale (sigma) given current mu
      scale <- optim(
        init_scale,
        fn = function(par){
          ## par is scale (sigma)
          -loglik_weibull(log_y, log_mu, status, par, observation_weights)
        },
        lower = init_scale/10,
        upper = init_scale*10,
        method = 'Brent'
      )$par

      ## -2 * log-likelihood
      dev <- -2*(
        status * (-log(scale) +
                    (1/scale - 1)*log_y -
                    log_mu/scale) -
          (exp((log_y - log_mu)/scale))
      )
      return(dev * observation_weights)
    },

  ## Added $loglik method for
  #  logLik.lgspline compatibility
  #  Accepts a fitted lgspline model object and returns the scalar
  #  log-likelihood value
  loglik = function(model_fit) {
    ## Extract necessary components from the fitted model
    log_y <- log(model_fit$y)
    log_mu <- log(model_fit$ytilde)
    ## sigmasq_tilde is dispersion (sigma^2); take sqrt to get scale (sigma)
    scale <- sqrt(model_fit$sigmasq_tilde)

    ## Status must be stored in model_fit$extra or model_fit$call_args
    #  Try common locations
    status <- NULL
    if(!is.null(model_fit$status)) {
      status <- model_fit$status
    } else if(!is.null(model_fit$extra_args$status)) {
      status <- model_fit$extra_args$status
    } else {
      ## If status not found, assume all events observed
      warning("Censoring status not found in model object; ",
              "assuming all events observed (status = 1).")
      status <- rep(1, length(log_y))
    }

    weights <- if(!is.null(model_fit$weights)){
      model_fit$weights
    } else if(!is.null(model_fit$observation_weights)){
      model_fit$observation_weights
    } else {
      rep(1, length(log_y))
    }

    loglik_weibull(log_y, log_mu, status, scale, weights)
  },

  ## Added $aic method for
  #  logLik.lgspline compatibility
  #  Returns AIC = -2*loglik + 2*edf
  #  edf = effective degrees of freedom (trace of hat matrix + 1 for scale)
  aic = function(model_fit) {
    ll <- weibull_family()$loglik(model_fit)
    ## edf: trace of smoother matrix + 1 for scale parameter
    edf <- model_fit$edf
    if(is.null(edf)){
      ## Fallback: use trace_XUGX if available
      edf <- model_fit$trace_XUGX
    }
    if(is.null(edf)){
      warning("Cannot compute AIC: effective degrees of freedom not found.")
      return(NA)
    }
    ## +1 for the scale parameter
    -2 * ll + 2 * (edf + 1)
  }
)


#' Estimate Weibull Dispersion for Accelerated Failure Time Model
#'
#' @description
#' Computes the dispersion parameter (sigma^2 = scale^2) for a Weibull
#' accelerated failure time (AFT) model, supporting right-censored survival
#' data. The returned value is sigma^2, where sigma is the Weibull scale
#' parameter matching \code{survreg$scale}.
#'
#' This both provides a tool for actually fitting Weibull AFT Models, and
#' boilerplate code for users who wish to incorporate Lagrangian multiplier
#' smoothing splines into their own custom models.
#'
#' @param mu Predicted survival times
#' @param y Observed response/survival times
#' @param order_indices Indices to align status with response
#' @param family Weibull AFT model family specification; unused here and
#'   retained for interface compatibility.
#' @param observation_weights Optional observation weights
#' @param VhalfInv Inverse square root of the correlation matrix; unused here
#'   and retained for interface compatibility.
#' @param status Censoring indicator (1 = event, 0 = censored)
#'   Indicates whether an event of interest occurred (1) or the observation was
#'   right-censored (0). In survival analysis, right-censoring occurs when the
#'   full survival time is unknown, typically because the study ended or the
#'   subject was lost to follow-up before the event of interest occurred.
#'
#' @return
#' Dispersion estimate (sigma^2) for the Weibull AFT model, i.e., the squared
#' scale parameter. The Weibull scale (sigma) matching \code{survreg$scale} is
#' \code{sqrt()} of this value.
#'
#' @seealso \code{\link{weibull_scale}} for the underlying scale estimation
#'   function
#'
#' @examples
#'
#' ## Simulate survival data with covariates
#' set.seed(1234)
#' n <- 1000
#' t1 <- rnorm(n)
#' t2 <- rbinom(n, 1, 0.5)
#'
#' ## Generate survival times with Weibull-like structure
#' lambda <- exp(0.5 * t1 + 0.3 * t2)
#' yraw <- rexp(n, rate = 1/lambda)
#'
#' ## Introduce right-censoring
#' status <- rbinom(n, 1, 0.75)
#' y <- ifelse(status, yraw, runif(length(yraw), 0, yraw))
#'
#' ## Example of using dispersion function
#' mu <- mean(y)
#' order_indices <- seq_along(y)
#' weights <- rep(1, n)
#'
#' ## Estimate dispersion (= scale^2 = sigma^2)
#' dispersion_est <- weibull_dispersion_function(
#'   mu = mu,
#'   y = y,
#'   order_indices = order_indices,
#'   family = weibull_family(),
#'   observation_weights = weights,
#'   VhalfInv = NULL,
#'   status = status
#' )
#'
#' print(dispersion_est)          # sigma^2
#' print(sqrt(dispersion_est))    # sigma (comparable to survreg$scale)
#'
#' @export
weibull_dispersion_function <- function(mu,
                                        y,
                                        order_indices,
                                        family,
                                        observation_weights,
                                        VhalfInv,
                                        status){

  ## Maximizes log-likelihood of right-censored data
  log_mu <- log(mu)
  log_y <- log(y)
  observation_weights <- c(observation_weights)
  status <- status[order_indices]

  ## Initialize scale (sigma) using null model
  init_scale <-
    weibull_scale(log_y,
                  mean(log_y),
                  status,
                  observation_weights)
  ## Find scale (sigma) given current mu
  scale <- optim(
    init_scale,
    fn = function(par){
      ## par is scale (sigma)
      -loglik_weibull(log_y,
                      log_mu,
                      status,
                      par,
                      observation_weights)
    },
    lower = init_scale/10, # Increased from 5 to 10
    upper = init_scale*10,
    method = 'Brent'
  )$par

  ## Return dispersion = sigma^2 = scale^2
  return(scale^2)
}

#' Weibull GLM Weight Function for Constructing Information Matrix
#'
#' @description
#' Computes the working weights used in the Weibull AFT information matrix,
#' including the observation-weight contribution returned on the vector scale.
#'
#' @param mu Predicted survival times
#' @param y Observed response/survival times
#' @param order_indices Order of observations when partitioned to match
#'   "status" to "response"
#' @param family Weibull AFT family; unused here and retained for interface
#'   compatibility.
#' @param dispersion Estimated dispersion parameter (sigma^2 = scale^2). The
#'   lgspline framework stores and passes dispersion (sigma^2); the Weibull
#'   scale (sigma) matching \code{survreg$scale} is \code{sqrt(dispersion)}.
#' @param observation_weights Weights of observations submitted to function
#' @param status Censoring indicator (1 = event, 0 = censored)
#'   Indicates whether an event of interest occurred (1) or the observation was
#'   right-censored (0). In survival analysis, right-censoring occurs when the
#'   full survival time is unknown, typically because the study ended or the
#'   subject was lost to follow-up before the event of interest occurred.
#'
#' @return
#' Numeric vector of working weights for the information matrix, including
#' observation weights when finite and a fallback of 1s when the natural
#' Weibull weights are non-finite.
#'
#' @details
#' This function generates weights used in constructing the information matrix
#' after unconstrained estimates have been found. Specifically, it is used in
#' the construction of the \eqn{\textbf{U}} and \eqn{\textbf{G}} matrices
#' following initial unconstrained parameter estimation.
#'
#' These weights are analogous to the variance terms in generalized linear
#' models (GLMs). Like logistic regression uses \eqn{\mu(1-\mu)}, Poisson
#' regression uses \eqn{e^{\mu}}, and Linear regression uses constant weights,
#' Weibull AFT models use \eqn{\exp((\log y - \log \mu)/\sigma)} where
#' \eqn{\sigma = \sqrt{\text{dispersion}}} is the scale parameter.
#'
#' @examples
#'
#' ## Demonstration of glm weight function in constrained model estimation
#' set.seed(1234)
#' n <- 1000
#' t1 <- rnorm(n)
#' t2 <- rbinom(n, 1, 0.5)
#'
#' ## Generate survival times
#' lambda <- exp(0.5 * t1 + 0.3 * t2)
#' yraw <- rexp(n, rate = 1/lambda)
#'
#' ## Introduce right-censoring
#' status <- rbinom(n, 1, 0.75)
#' y <- ifelse(status, yraw, runif(length(yraw), 0, yraw))
#'
#' ## Fit model demonstrating use of custom glm weight function
#' model_fit <- lgspline(y ~ spl(t1) + t2,
#'                       data.frame(y = y, t1 = t1, t2 = t2),
#'                       unconstrained_fit_fxn = unconstrained_fit_weibull,
#'                       family = weibull_family(),
#'                       need_dispersion_for_estimation = TRUE,
#'                       dispersion_function = weibull_dispersion_function,
#'                       glm_weight_function = weibull_glm_weight_function,
#'                       schur_correction_function = weibull_schur_correction,
#'                       status = status,
#'                       opt = FALSE,
#'                       K = 1)
#'
#' print(summary(model_fit))
#'
#'
#' @export
weibull_glm_weight_function <- function(mu,
                                        y,
                                        order_indices,
                                        family,
                                        dispersion,
                                        observation_weights,
                                        status){
  ## dispersion = sigma^2; scale = sigma = sqrt(dispersion)
  scale <- sqrt(dispersion)
  val <- exp((log(y) - log(mu)) / scale)
  if(any(!is.finite(val))){
    return(rep(1, length(val)))
  }
  newval <- val * c(observation_weights)
  return(newval)
}

#' Unconstrained Weibull Accelerated Failure Time Model Estimation
#'
#' @description
#' Estimates parameters for an unconstrained Weibull accelerated failure time
#' (AFT) model supporting right-censored survival data.
#'
#' This both provides a tool for actually fitting Weibull AFT Models, and
#' boilerplate code for users who wish to incorporate Lagrangian multiplier
#' smoothing splines into their own custom models.
#'
#' @param X Design matrix of predictors
#' @param y Survival/response times
#' @param LambdaHalf Square root of penalty matrix
#'   (\eqn{\boldsymbol{\Lambda}^{1/2}})
#' @param Lambda Penalty matrix (\eqn{\boldsymbol{\Lambda}})
#' @param keep_weighted_Lambda Flag to retain weighted penalties
#' @param family Distribution family specification
#' @param tol Convergence tolerance (default 1e-8)
#' @param K Number of partitions minus one (\eqn{K})
#' @param parallel Flag for parallel processing
#' @param cl Cluster object for parallel computation
#' @param chunk_size Processing chunk size
#' @param num_chunks Number of computational chunks
#' @param rem_chunks Remaining chunks
#' @param order_indices Observation ordering indices
#' @param weights Optional observation weights
#' @param status Censoring status indicator (1 = event, 0 = censored)
#'   Indicates whether an event of interest occurred (1) or the observation was
#'   right-censored (0). In survival analysis, right-censoring occurs when the
#'   full survival time is unknown, typically because the study ended or the
#'   subject was lost to follow-up before the event of interest occurred.
#'
#' @return
#' Numeric column vector of unconstrained coefficient estimates for the Weibull
#' AFT model.
#'
#' @details
#' Estimation Approach:
#' The function employs a two-stage optimization strategy for fitting
#' accelerated failure time models via maximum likelihood:
#'
#' 1. Outer Loop: Estimate Scale Parameter (sigma) using Brent's method
#'
#' 2. Inner Loop: Estimate Regression Coefficients (given sigma) using
#'    damped Newton-Raphson.
#'
#' Note: the score and information inside the Newton-Raphson are both scaled
#' by \eqn{\sigma^2} (i.e., both omit the \eqn{1/\sigma} and
#' \eqn{1/\sigma^2} prefactors respectively). Since the Newton-Raphson step is
#' \eqn{\mathbf{G}u = (\mathbf{X}^{\top}\mathbf{W}\mathbf{X})^{-1}
#' \mathbf{X}^{\top}\mathbf{v}}, the \eqn{\sigma^2} factors cancel and the
#' step remains correct.
#'
#' @examples
#'
#' ## Simulate survival data with covariates
#' set.seed(1234)
#' n <- 1000
#' t1 <- rnorm(n)
#' t2 <- rbinom(n, 1, 0.5)
#'
#' ## Generate survival times with Weibull-like structure
#' lambda <- exp(0.5 * t1 + 0.3 * t2)
#' yraw <- rexp(n, rate = 1/lambda)
#'
#' ## Introduce right-censoring
#' status <- rbinom(n, 1, 0.75)
#' y <- ifelse(status, yraw, runif(length(yraw), 0, yraw))
#' df <- data.frame(y = y, t1 = t1, t2 = t2)
#'
#' ## Fit model using lgspline with Weibull AFT unconstrained estimation
#' model_fit <- lgspline(y ~ spl(t1) + t2,
#'                       df,
#'                       unconstrained_fit_fxn = unconstrained_fit_weibull,
#'                       family = weibull_family(),
#'                       need_dispersion_for_estimation = TRUE,
#'                       dispersion_function = weibull_dispersion_function,
#'                       glm_weight_function = weibull_glm_weight_function,
#'                       schur_correction_function = weibull_schur_correction,
#'                       status = status,
#'                       opt = FALSE,
#'                       K = 1)
#'
#' ## Print model summary
#' summary(model_fit)
#'
#' @keywords internal
#' @export
unconstrained_fit_weibull <- function(X,
                                      y,
                                      LambdaHalf,
                                      Lambda,
                                      keep_weighted_Lambda,
                                      family,
                                      tol = 1e-8,
                                      K,
                                      parallel,
                                      cl,
                                      chunk_size,
                                      num_chunks,
                                      rem_chunks,
                                      order_indices,
                                      weights,
                                      status # passed through for censored Weibull fits
) {

  ## Weight if non-null
  if(any(!is.null(weights))){
    if(length(weights) > 0){
      weights <- c(weights)
    }
  } else {
    weights <- rep(1, length(y))
  }

  log_y <- log(y)

  ## Initialize scale (sigma)
  init_scale <- weibull_scale(log_y,
                              mean(log_y),
                              status[order_indices],
                              weights)

  ## First, use outer-loop to optimize scale (sigma)
  #  Then given scale, optimize beta.
  #  Note: the score and information inside damped_newton_r are both
  #  scaled by sigma^2 (omitting 1/sigma and 1/sigma^2 prefactors),
  #  so the Newton-Raphson step remains correct.
  scale <- optim(init_scale,
                 fn = function(par){
                   ## par is scale (sigma)
                   scale <- par
                   beta <- cbind(damped_newton_r(
                     c(mean(log_y), rep(0, ncol(X)-1)),
                     function(par){
                       beta <- cbind(par)
                       log_mu <- c(X %**% beta)
                       loglik_weibull(log_y,
                                      log_mu,
                                      status[order_indices],
                                      scale,
                                      weights) +
                         -0.5*c(crossprod(beta, Lambda) %**% beta)
                     },
                     function(par){
                       beta <- cbind(par)
                       eta <- c(X %**% beta)
                       z <- (log_y - eta) / scale
                       zeta <- exp(z)
                       ## Scaled score (x sigma^2): true score is
                       #  (1/sigma) * X^T diag(w)(exp(z) - status),
                       #  multiplied through by sigma^2 gives
                       #  sigma * X^T diag(w)(exp(z) - status)
                       grad_beta_sc <- crossprod(
                         X,
                         (weights*(zeta - status[order_indices])) * scale
                       ) +
                         Lambda %**% beta
                       cbind(grad_beta_sc)
                     },
                     function(par){
                       beta <- cbind(par)
                       eta <- X %**% beta
                       z <- (log_y - eta) / scale
                       zeta <- c(exp(z))
                       ## Scaled information (x sigma^2): true info is
                       #  (1/sigma^2) * X^T diag(w*exp(z)) X,
                       #  multiplied through by sigma^2 gives
                       #  X^T diag(w*exp(z)) X
                       info_sc <- crossprod(X, weights * zeta * X) + Lambda
                       info_sc
                     },
                     tol
                   ))
                   ## Return negated penalized log-likelihood
                   #  as the objective to minimize
                   log_mu <- X %**% beta
                   -loglik_weibull(log_y,
                                   log_mu,
                                   status[order_indices],
                                   par,
                                   weights) +
                     0.5*c(crossprod(beta, Lambda) %**% beta)
                 },
                 lower = init_scale/10, # Increased to 10
                 upper = init_scale*10, # (was 5)
                 method = 'Brent')$par

  ## Now optimize beta, given optimal scale (sigma)
  beta <- cbind(damped_newton_r(
    c(mean(log_y), rep(0, ncol(X)-1)),
    function(par){
      beta <- cbind(par)
      log_mu <- c(X %**% beta)
      loglik_weibull(log_y,
                     log_mu,
                     status[order_indices],
                     scale,
                     weights) +
        -0.5*c(crossprod(beta, Lambda) %**% beta)
    },
    function(par){
      beta <- cbind(par)
      eta <- c(X %**% beta)
      z <- (log_y - eta) / scale
      zeta <- exp(z)
      ## Scaled score (x sigma^2); see note above
      grad_beta_sc <- crossprod(
        X,
        (weights*(zeta - status[order_indices])) * scale
      ) +
        Lambda %**% beta
      cbind(grad_beta_sc)
    },
    function(par){
      beta <- cbind(par)
      eta <- X %**% beta
      z <- (log_y - eta) / scale
      zeta <- c(exp(z))
      ## Scaled information (x sigma^2); see note above
      info_sc <- crossprod(X, weights * zeta * X) + Lambda
      info_sc
    },
    tol
  ))

  return(beta)
}


#' Fit Weibull Accelerated Failure Time Model via lgspline
#'
#' @description
#' Convenience wrapper that calls \code{lgspline} with the correct
#' family, weight, dispersion, score, and unconstrained-fit functions for
#' Weibull accelerated failure time regression. All standard lgspline
#' arguments (knots, penalties, constraints, parallel, etc.) are passed
#' through.
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
#'   \item \code{family = weibull_family()}
#'   \item \code{unconstrained_fit_fxn = unconstrained_fit_weibull}
#'   \item \code{glm_weight_function = weibull_glm_weight_function}
#'   \item \code{qp_score_function = weibull_qp_score_function}
#'   \item \code{dispersion_function = weibull_dispersion_function}
#'   \item \code{schur_correction_function = weibull_schur_correction}
#'   \item \code{need_dispersion_for_estimation = TRUE}
#'   \item \code{estimate_dispersion = TRUE}
#'   \item \code{standardize_response = FALSE}
#' }
#'
#' A formula interface is needed, e.g. \code{lgspline_weibull(t, y, ...)} won't work,
#' unlike for ordinary \code{\link{lgspline}}.
#'
#' @examples
#'
#' ## Weibull AFT with a nonlinear age effect on lung cancer survival
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
#'   fit <- lgspline_weibull(
#'     time ~ spl(age_std),
#'     data = lung,
#'     status = event,
#'     K = 1,
#'     opt = FALSE,
#'     wiggle_penalty = 1e-4,
#'     flat_ridge_penalty = 1
#'   )
#'   print(summary(fit))
#'   plot(fit,
#'        show_formulas = TRUE,
#'        custom_response_lab = 'Survival Time',
#'        custom_predictor_lab = 'Standardized Age')
#' }
#'
#' @seealso \code{\link{lgspline_cox}} for Cox proportional hazards,
#'   \code{\link{weibull_family}}, \code{\link{unconstrained_fit_weibull}}
#'
#' @export
lgspline_weibull <- function(formula, data, status, ...) {
  lgspline(
    formula = formula,
    data = data,
    family = weibull_family(),
    unconstrained_fit_fxn = unconstrained_fit_weibull,
    glm_weight_function = weibull_glm_weight_function,
    qp_score_function = weibull_qp_score_function,
    dispersion_function = weibull_dispersion_function,
    schur_correction_function = weibull_schur_correction,
    need_dispersion_for_estimation = TRUE,
    estimate_dispersion = TRUE,
    standardize_response = FALSE,
    status = status,
    ...
  )
}
