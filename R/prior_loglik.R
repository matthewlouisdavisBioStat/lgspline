#' Log-Prior Distribution Evaluation for lgspline Models
#'
#' @description
#' Evaluates the log-prior distribution on beta coefficients conditional upon dispersion and penalaties,
#'
#' @details
#' Returns the quadratic form of B^T(Lambda)B evaluated at the
#' tuned or fixed penalties, scaled by negative one-half inverse dispersion.
#'
#'
#' Assuming fixed penalties, the prior distribution of \eqn{\beta} is given as follows:
#'
#' \deqn{\beta | \sigma^2 \sim \mathcal{N}(\textbf{0}, \frac{1}{\sigma^2}\Lambda)}
#'
#' The log-likelihood obtained from this can be shown to be equivalent to the following,
#' with \eqn{C} a normalizing constant with respect to \eqn{\beta}.
#'
#' The value of \eqn{C} is returned when \code{include_constant=TRUE}, and ignored when \code{FALSE}.
#'
#' \deqn{\implies \log P(\beta|\sigma^2) = C-\frac{1}{2\sigma^2}\beta^{T}\Lambda\beta}
#'
#' This is useful for computing joint log-likelihoods and performing valid
#' likelihood ratio tests between nested lgspline models.
#'
#'
#' @param model_fit An lgspline model object.
#' @param sigmasq A scalar numeric representing the dispersion parameter. If NULL, model_fit$sigmasq_tilde is used.
#' @param include_constant Logical; if TRUE (default), include the multivariate normal normalizing constant.
#'
#' @return A numeric scalar representing the prior log-likelihood.
#'
#' @examples
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
#'           prior_loglik(model_fit, sigmasq = 1)
#' print(jntloglik)
#'
#' @seealso \code{\link{lgspline}}
#'
#' @export
prior_loglik <- function(model_fit,
                         sigmasq = NULL,
                         include_constant = TRUE){

  if(is.null(sigmasq)){
    sigmasq <- model_fit$sigmasq_tilde
  }

  if(length(model_fit$penalties$L_partition_list) == 0){
    model_fit$penalties$L_partition_list <- lapply(
      1:(model_fit$K + 1), function(k) 0
    )
  }

  running_sum <- 0
  total_constant <- 0

  for(k in 1:(model_fit$K + 1)){

    beta_k <- model_fit$B_raw[[k]]
    Lambda_total <- model_fit$penalties$Lambda +
      model_fit$penalties$L_partition_list[[k]]

    p_k <- length(beta_k)

    quad_form <- as.numeric(
      t(beta_k) %**% Lambda_total %**% beta_k
    )

    running_sum <- running_sum -
      0.5 * quad_form / sigmasq

    if(include_constant){

      logdet <- as.numeric(
        determinant(Lambda_total, logarithm = TRUE)$modulus
      )

      total_constant <- total_constant -
        0.5 * (
          p_k * log(2 * pi) -
            p_k * log(sigmasq) +
            logdet
        )
    }
  }

  return(as.numeric(running_sum + total_constant))
}
