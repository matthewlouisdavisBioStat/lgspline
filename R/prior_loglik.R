#' Log-Prior Distribution Evaluation for lgspline Models
#'
#' @description
#' Evaluates the log-prior on the spline coefficients conditional on the
#' dispersion and penalty matrices.
#'
#' @details
#' Returns the quadratic form of \eqn{\beta^{T}\Lambda\beta} evaluated at the
#' tuned or fixed penalties, scaled by negative one-half inverse dispersion.
#'
#' Assuming fixed penalties, the prior on \eqn{\beta} is taken to be
#' \deqn{\beta | \sigma^2 \sim \mathcal{N}(\textbf{0}, \sigma^2\Lambda^{-1})}
#' so that, up to a normalizing constant \eqn{C} with respect to \eqn{\beta},
#' \deqn{\implies \log P(\beta|\sigma^2) = C-\frac{1}{2\sigma^2}\beta^{T}\Lambda\beta}
#'
#' The value of \eqn{C} is included when \code{include_constant = TRUE}, and
#' omitted when \code{FALSE}.
#'
#' This is useful for computing joint penalized log-likelihoods and related
#' MAP-style diagnostics for a fitted \code{lgspline} object.
#'
#' @param model_fit An lgspline model object.
#' @param B_predict Optional list of coefficient vectors at which to evaluate
#'   the prior. Default NULL uses the fitted coefficients.
#' @param sigmasq_predict Numeric scalar dispersion parameter. If NULL,
#'   \code{model_fit$sigmasq_tilde} is used. Legacy \code{sigmasq} calls
#'   are still accepted.
#' @param include_constant Logical; if TRUE (default), include the
#'   multivariate normal normalizing constant.
#' @param ... Optional legacy arguments.
#'
#' @return A numeric scalar representing the prior log-likelihood.
#'
#' @examples
#'
#' \donttest{
#' ## Data
#' t <- sort(runif(100, -5, 5))
#' y <- sin(t) - 0.1*t^2 + rnorm(100)
#'
#' ## Model keeping penalties fixed
#' model_fit <- lgspline(t, y, opt = FALSE)
#'
#' ## Full joint log-likelihood, conditional upon known sigma^2 = 1
#' jntloglik <- sum(dnorm(model_fit$y,
#'                     model_fit$ytilde,
#'                     1,
#'                     log = TRUE)) +
#'           prior_loglik(model_fit, sigmasq_predict = 1)
#' print(jntloglik)
#' }
#'
#' @seealso \code{\link{lgspline}}
#'
#' @export
prior_loglik <- function(model_fit,
                         B_predict = NULL,
                         sigmasq_predict = NULL,
                         include_constant = TRUE,
                         ...){
  dots <- list(...)
  if(!is.null(B_predict) &&
     is.null(sigmasq_predict) &&
     is.numeric(B_predict) &&
     length(B_predict) == 1){

    sigmasq_predict <- B_predict

    B_predict <- NULL

  }
  if(is.null(sigmasq_predict) && !is.null(dots$sigmasq)){
    sigmasq_predict <- dots$sigmasq
  }
  if(is.null(sigmasq_predict)){
    sigmasq_predict <- model_fit$sigmasq_tilde
  }
  if(length(sigmasq_predict) != 1 ||
     !is.finite(sigmasq_predict) ||
     sigmasq_predict <= 0){
    stop('\n\t sigmasq_predict must be a positive finite scalar\n')
  }

  sigmasq_predict <- as.numeric(sigmasq_predict)

  if(is.null(B_predict)){
    B_raw_predict <- model_fit$B_raw
  } else {
    B_raw_predict <- lapply(B_predict, function(b){
      b_raw <- model_fit$forwtransform_coefficients(b)
      b_raw[1] <- b_raw[1] - model_fit$mean_y
      b_raw / model_fit$sd_y
    })
    names(B_raw_predict) <- names(B_predict)
  }

  ## If no partition-specific penalties were stored, treat them as zero.
  if(length(model_fit$penalties$L_partition_list) == 0){
    model_fit$penalties$L_partition_list <- lapply(
      1:(model_fit$K + 1), function(k) 0
    )
  }

  running_sum <- 0
  total_constant <- 0

  ## Sum the quadratic penalty contribution partition by partition.
  for(k in 1:(model_fit$K + 1)){
    beta_k <- B_raw_predict[[k]]
    Lambda_total <- model_fit$penalties$Lambda +
      model_fit$penalties$L_partition_list[[k]]
    p_k <- length(beta_k)
    quad_form <- as.numeric(
      t(beta_k) %**% Lambda_total %**% beta_k
    )
    running_sum <- running_sum -
      0.5 * quad_form / sigmasq_predict

    ## Add the multivariate-normal normalizing constant when requested.
    #  For beta | sigma^2 ~ N(0, sigma^2 Lambda^{-1}),
    #  the constant contribution is
    #    -0.5 * (p log(2 pi) + p log(sigma^2) - log |Lambda|).
    if(include_constant){
      logdet <- as.numeric(
        determinant(Lambda_total, logarithm = TRUE)$modulus
      )
      total_constant <- total_constant -
        0.5 * (
          p_k * log(2 * pi) +
            p_k * log(sigmasq_predict) -
            logdet
        )
    }
  }
  return(as.numeric(running_sum + total_constant))

}



