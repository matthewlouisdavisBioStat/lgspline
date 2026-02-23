#' Compute Leave-One-Out Cross-Validated Predictions for
#' Gaussian Response/Identity Link under Constraint.
#'
#' @description
#' Computes the leave-one-out cross-validated predictions from a model
#' fit, assuming Gaussian-distributed response with identity link.
#'
#' The LOO closed-formula for observation \eqn{i} is \eqn{\hat{y}_{(-i)} = y_i -
#' \frac{1}{1 - H_{ii}}(y_i - \hat{y}_i)} where
#' \eqn{\mathbf{H}} is the effective hat matrix under
#' smoothing constraints, adjusted for weights and correlation structure if
#' present.
#'
#' Observations with leverage at or above \code{leverage_threshold} have
#' the denominator \eqn{1 - H_{ii}} near zero, making the shortcut
#' numerically unreliable. A warning is
#' issued to consider re-fitting with greater penalties or fewer knots if
#' many observations are flagged.
#'
#' @param model_fit A fitted lgspline model object.
#' @param leverage_threshold Numeric scalar in (0, 1). Observations with
#'   \eqn{H_{ii} \geq} \code{leverage_threshold} are treated as
#'   high-leverage and their LOO predictions are set to \code{NA}.
#'   Default \code{100}.
#'
#' @return A vector of leave-one-out cross-validated predictions
#'
#' @examples
#'
#' ## Basic usage with Gaussian response, computing PRESS
#' set.seed(1234)
#' t <- rnorm(50)
#' y <- sin(t) + rnorm(50, 0, .25)
#' model_fit <- lgspline(t, y)
#' loo <- leave_one_out(model_fit)
#' press <- mean((y - loo)^2, na.rm = TRUE)
#'
#' plot(loo, y,
#'    main = "LOO Cross-Validation Prediction vs. Observed Response",
#'    xlab = 'Prediction', ylab = 'Response')
#' abline(0, 1)
#'
#' @export
leave_one_out <- function(model_fit,
                          leverage_threshold = 100){

  ## Full N x P design matrix on the standardized scale
  X_block <- collapse_block_diagonal(
    lapply(model_fit$X, model_fit$std_X)
  )

  if(!is.null(model_fit$VhalfInv)){

    ## [Change 2026-02-17] Under correlated response the hat diagonal must
    #  come from G_correct = (X'V^{-1}X + Lambda)^{-1}, not the
    #  block-diagonal G. Uses the same G_correct path as varcovmat and
    #  trace_XUGX for consistency, which is the P x P dense matrix
    #  rather usual block-diagonal.
    X_ordered <- X_block[unlist(model_fit$og_order), , drop = FALSE]
    VinvhalfX <- t(t(model_fit$VhalfInv %**% X_ordered) *
                     model_fit$weights)

    has_part_pen <-
      length(model_fit$penalties$L_partition_list) == (model_fit$K + 1)
    Lambda_full <- collapse_block_diagonal(
      lapply(1:(model_fit$K + 1), function(k){
        if(has_part_pen){
          model_fit$penalties$Lambda +
            model_fit$penalties$L_partition_list[[k]]
        } else {
          model_fit$penalties$Lambda
        }
      })
    )

    Ghalf_correct <- matinvsqrt(gramMatrix(VinvhalfXU) + Lambda_full)

    ## hat diagonal: row-wise squared norms of V^{-1/2} X G_correct^{1/2}
    const <- sqrt(norm(Ghalf_correct, "2"))
    VinvhalfXUGhalf <- VinvhalfX %**% (model_fit$U %**% (Ghalf_correct / const))
    diag_hat <- rowSums(VinvhalfXUGhalf * VinvhalfXUGhalf) * const^2

    ## Reorder back to partition order used by y and ytilde
    diag_hat <- diag_hat[unlist(model_fit$order_list)]

  } else {

    ## [Change 2026-02-17] Fixed pre-existing bug: matmult_U called with
    #  Ghalf (not G) so that diag(X U G X') is computed rather than
    #  diag(X U G^2 X').
    UGhalf    <- matmult_U(model_fit$U,
                           model_fit$Ghalf,
                           model_fit$p,
                           model_fit$K)
    const  <- sqrt(norm(UGhalf, "2"))
    X_block <- t(t(X_block) * model_fit$weights)
    diag_hat <- rowSums((X_block %**% (UGhalf / const)) * X_block) * const
    diag_hat <- diag_hat[unlist(model_fit$og_order)]
  }

  ## Flag high-leverage observations where 1/(1 - H_ii) is unreliable
  n_high <- sum(diag_hat >= leverage_threshold, na.rm = TRUE)
  if(n_high > 0){
    warning(
      n_high, " observation(s) have leverage >= ", leverage_threshold,
      "; LOO predictions for these are unreliable.",
      "Consider re-fitting with greater penalties or fewer knots.",
      call. = FALSE
    )
  }
  denom <- 1 - diag_hat

  leave_one_out <- model_fit$y - (model_fit$y - model_fit$ytilde) / denom

  return(leave_one_out)
}
