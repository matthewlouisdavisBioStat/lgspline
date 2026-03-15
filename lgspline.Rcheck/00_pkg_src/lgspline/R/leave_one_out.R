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
#' Observations with leverage at or above \code{leverage_threshold} are flagged
#' in a warning, since extreme hat values can make the shortcut numerically
#' unreliable. The default \code{leverage_threshold = 100} is intentionally
#' permissive, so users who want diagnostic warnings for large \eqn{H_{ii}}
#' should set a smaller threshold explicitly.
#'
#' @param model_fit A fitted lgspline model object.
#' @param leverage_threshold Numeric scalar. Observations with
#'   \eqn{H_{ii} \geq} \code{leverage_threshold} are treated as
#'   high-leverage for the warning below. Default \code{100}.
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

    ## Under correlation, compute the hat diagonal from the same dense
    #  GLS system used by varcovmat and trace_XUGX.
    X_ordered <- X_block[unlist(model_fit$og_order), , drop = FALSE]
    VinvhalfX <- t(t(model_fit$VhalfInv %**% X_ordered) *
                     sqrt(model_fit$weights))

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

    Ghalf_correct <- matinvsqrt(gramMatrix(VinvhalfX) + Lambda_full)

    ## hat diagonal: diag(V^{-1/2} X U G X^T V^{-1/2})
    #  = rowSums((V^{-1/2} X U G^{1/2}) * (V^{-1/2} X G^{1/2}))
    const <- max(sqrt(norm(Ghalf_correct, "2")), sqrt(.Machine$double.eps))
    VinvhalfXUGhalf <-
      VinvhalfX %**% (model_fit$U %**% (Ghalf_correct / const))
    VinvhalfXGhalf <- VinvhalfX %**% (Ghalf_correct / const)
    diag_hat <- rowSums(VinvhalfXUGhalf * VinvhalfXGhalf) * const^2

  } else {

    ## Without correlation, use the block-diagonal shortcut for the
    #  hat diagonal.
    UGhalf <- matmult_U(model_fit$U,
                        model_fit$Ghalf,
                        model_fit$p,
                        model_fit$K)
    Ghalf_block <- collapse_block_diagonal(model_fit$Ghalf)
    const <- max(sqrt(norm(UGhalf, "2")), sqrt(.Machine$double.eps))
    X_block <- t(t(X_block) * sqrt(model_fit$weights))
    XUGhalf <- X_block %**% (UGhalf / const)
    XGhalf <- X_block %**% (Ghalf_block / const)
    diag_hat <- rowSums(XUGhalf * XGhalf) * const^2
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
  ## Apply the standard leave-one-out shortcut once the hat diagonal is in hand.
  denom <- 1 - diag_hat

  leave_one_out <- model_fit$y - (model_fit$y - model_fit$ytilde) / denom

  return(leave_one_out)
}
