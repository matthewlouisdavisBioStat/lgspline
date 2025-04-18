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
#' with \eqn{C} a constant with respect to \eqn{\beta}.
#'
#' \deqn{\implies \log P(\beta|\sigma^2) = C-\frac{1}{2\sigma^2}\beta^{T}\Lambda\beta}
#'
#' This is useful for computing joint log-likelihoods and performing valid
#' likelihood ratio tests between nested lgspline models.
#'
#'
#' @param model_fit An lgspline model object
#' @param sigmasq A scalar numeric representing the dispersion parameter. By default it is NULL, and the sigmasq_tilde associated with model_fit will be used. Otherwise, custom values can be supplied.
#'
#' @return A numeric scalar for the prior-loglikelihood (the penalty on beta coefficients
#' actually computed)
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
prior_loglik <- function(model_fit, sigmasq = NULL){
      if(is.null(sigmasq)){
        sigmasq <- model_fit$sigmasq_tilde
      }
      if(length(model_fit$penalties$L_partition_list) == 0){
        model_fit$penalties$L_partition_list <- lapply(
          1:(model_fit$K + 1), function(k)0
        )
      }
      running_sum <- -0.5 *
      (t(model_fit$B_raw[[1]]) %**%
        (model_fit$penalties$Lambda +
         model_fit$penalties$L_partition_list[[1]]) %**%
        model_fit$B_raw[[1]]) / sigmasq
      for(k in 2:(model_fit$K + 1)){
        running_sum <- running_sum +
          -0.5 *
         (t(model_fit$B_raw[[k]]) %**%
           (model_fit$penalties$Lambda +
            model_fit$penalties$L_partition_list[[k]]) %**%
            model_fit$B_raw[[k]]) / sigmasq
      }
  return(running_sum)
}
