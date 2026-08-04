#' Compute LOO GLM square-root weights
#'
#' @param model_fit A fitted lgspline model object.
#' @param order_indices Integer observation order used by the working system.
#'
#' @return Numeric vector of square-root working weights.
#' @noRd
.loo_glm_sqrt_weights <- function(model_fit, order_indices){
  obs_wt <- c(model_fit$weights)[order_indices]
  glm_weight_function <- NULL
  if(!is.null(model_fit$.fit_call_args$glm_weight_function)){
    glm_weight_function <- model_fit$.fit_call_args$glm_weight_function
  }
  if(is.null(glm_weight_function)){
    glm_weight_function <- default_glm_weight_function
  }

  extra <- list()
  if("status" %in% names(formals(glm_weight_function))){
    status <- model_fit$status
    if(is.null(status)) status <- model_fit$.fit_call_args$status
    if(!is.null(status)) extra$status <- status
  }

  args <- c(list(
    mu = model_fit$ytilde[order_indices],
    y = model_fit$y[order_indices],
    order_indices = order_indices,
    family = model_fit$family,
    dispersion = model_fit$sigmasq_tilde,
    observation_weights = obs_wt,
    glm_weight_function = glm_weight_function
  ), extra)

  out <- tryCatch(
    do.call(.gee_glm_working_components, args),
    error = function(e) NULL
  )
  if(is.null(out)) return(sqrt(obs_wt))
  out$sqrtW
}


.loo_hat_diag_lgspline <- function(model_fit){
  X_block <- collapse_block_diagonal(
    lapply(model_fit$X, model_fit$std_X)
  )

  if(!is.null(model_fit$VhalfInv)){
    X_ordered <- X_block[unlist(model_fit$og_order), , drop = FALSE]
    sqrtW <- .loo_glm_sqrt_weights(model_fit, seq_len(model_fit$N))
    VinvhalfX <- model_fit$VhalfInv %**% (X_ordered * c(sqrtW))

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

    Ghalf_correct <- matinvsqrt(crossprod(VinvhalfX) + Lambda_full)
    const <- max(sqrt(norm(Ghalf_correct, "2")), sqrt(.Machine$double.eps))
    VinvhalfXUGhalf <-
      VinvhalfX %**% (model_fit$U %**% (Ghalf_correct / const))
    VinvhalfXGhalf <- VinvhalfX %**% (Ghalf_correct / const)
    diag_hat <- rowSums(VinvhalfXUGhalf * VinvhalfXGhalf) * const^2
  } else {
    UGhalf <- matmult_U(model_fit$U,
                        model_fit$Ghalf,
                        model_fit$p,
                        model_fit$K)
    Ghalf_block <- collapse_block_diagonal(model_fit$Ghalf)
    const <- max(sqrt(norm(UGhalf, "2")), sqrt(.Machine$double.eps))
    sqrtW <- .loo_glm_sqrt_weights(model_fit, unlist(model_fit$order_list))
    X_block <- X_block * c(sqrtW)
    XUGhalf <- X_block %**% (UGhalf / const)
    XGhalf <- X_block %**% (Ghalf_block / const)
    diag_hat <- rowSums(XUGhalf * XGhalf) * const^2
    diag_hat <- diag_hat[unlist(model_fit$og_order)]
  }
  diag_hat
}


#' Compute Leave-One-Out Cross-Validated Predictions
#'
#' @description
#' Computes leave-one-out cross-validated predictions from a fitted
#' \code{lgspline} object. For \code{additive_lgspline} fits, the calculation
#' constructs the full appended additive design matrix, combines term penalties
#' and active constraints, and computes the hat diagonal from that single dense
#' system.
#'
#' The LOO shortcut is \eqn{\hat{y}_{(-i)} = y_i -
#' (y_i - \hat{y}_i)/(1 - H_{ii})}, where \eqn{\mathbf{H}} is the effective hat
#' matrix adjusted for weights and correlation structure when present.
#'
#' Observations with leverage at or above \code{leverage_threshold} are flagged
#' in a warning, since extreme hat values can make the calculation numerically
#' unreliable.
#'
#' @param model_fit A fitted \code{lgspline} or \code{additive_lgspline} object.
#' @param leverage_threshold Numeric scalar. Observations with
#'   \eqn{H_{ii} \geq} \code{leverage_threshold} are treated as high leverage.
#'   Default \code{100}.
#'
#' @return A vector of leave-one-out cross-validated predictions.
#'
#' @examples
#' set.seed(1234)
#' t <- rnorm(50)
#' y <- sin(t) + rnorm(50, 0, .25)
#' model_fit <- lgspline(t, y)
#' loo <- leave_one_out(model_fit)
#' press <- mean((y - loo)^2, na.rm = TRUE)
#'
#' plot(loo, y,
#'   main = "LOO Cross-Validation Prediction vs. Observed Response",
#'   xlab = "Prediction", ylab = "Response")
#' abline(0, 1)
#'
#' @references
#' Tarpey, T. (2000). A note on the prediction sum of squares statistic for
#' restricted least squares. \emph{The American Statistician}, 54(2), 116--118.
#' \doi{10.2307/2686028}
#'
#' @export
leave_one_out <- function(model_fit,
                          leverage_threshold = 100){

  if(inherits(model_fit, "additive_lgspline")){
    system <- .additive_combined_system(
      model_fit,
      VhalfInv = model_fit$VhalfInv,
      sigmasq_tilde = model_fit$sigmasq_tilde,
      use_glm_weights = TRUE
    )
    const <- max(sqrt(norm(system$Ghalf, "2")), sqrt(.Machine$double.eps))
    XUGhalf <- system$X_weighted %**% (system$U %**%
                                         (system$Ghalf / const))
    XGhalf <- system$X_weighted %**% (system$Ghalf / const)
    diag_hat <- rowSums(XUGhalf * XGhalf) * const^2

    n_high <- sum(diag_hat >= leverage_threshold, na.rm = TRUE)
    if(n_high > 0){
      warning(
        n_high, " observation(s) have leverage >= ", leverage_threshold,
        "; LOO predictions for these are unreliable.",
        " Consider re-fitting with greater penalties or fewer knots.",
        call. = FALSE
      )
    }
    denom <- 1 - diag_hat
    return(model_fit$y - (model_fit$y - model_fit$ytilde) / denom)
  }

  ## Full N x P design matrix on the standardized scale
  X_block <- collapse_block_diagonal(
    lapply(model_fit$X, model_fit$std_X)
  )

  if(!is.null(model_fit$VhalfInv)){

    ## Under correlation, compute the hat diagonal from the same dense
    #  GLS system used by varcovmat and trace_XUGX.
    X_ordered <- X_block[unlist(model_fit$og_order), , drop = FALSE]
    sqrtW <- .loo_glm_sqrt_weights(model_fit, seq_len(model_fit$N))
    VinvhalfX <- model_fit$VhalfInv %**% (X_ordered * c(sqrtW))

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

    Ghalf_correct <- matinvsqrt(crossprod(VinvhalfX) + Lambda_full)

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
    sqrtW <- .loo_glm_sqrt_weights(model_fit, unlist(model_fit$order_list))
    X_block <- X_block * c(sqrtW)
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
