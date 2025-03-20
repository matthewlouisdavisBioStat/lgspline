#' Compute Leave-One-Out Cross-Validated predictions for
#' Gaussian Response/Identity Link under Constraint.
#'
#' @description
#' Computes the leave-one-out cross-validated predictions from a model fit,
#' assuming Gaussian-distributed response with identity link.
#'
#' @param model_fit A fitted Lagrangian smoothing spline model
#'
#' @return A vector of leave-one-out cross-validated predictions
#'
#' @examples
#' \donttest{
#' ## Basic usage with Gaussian response, computing PRESS
#' set.seed(1234)
#' x <- rnorm(50)
#' y <- sin(x) + rnorm(50, 0, .25)
#' fit <- lgspline(x, y)
#' loo <- leave_one_out(fit)
#' press <- mean(loo^2)
#'
#' plot(loo, y,
#'    main = "Leave-One-Out Cross-Validation Prediction vs. Observed Response",
#'    xlab = 'Prediction', ylab = 'Response')
#' abline(0, 1)
#'
#' }
#'
#' @export
leave_one_out <- function(model_fit){

  ## Collapse X into a single block-diagonal matrix, unscaled
  X_block <- collapse_block_diagonal(
    lapply(model_fit$X, model_fit$std_X)
  )

  ## Multiply by VhalfInv if present
  if(!is.null(model_fit$VhalfInv)){
    X_block <- X_block[unlist(model_fit$og_order),]
    X_block <- model_fit$VhalfInv %**% X_block
    X_block <- X_block[unlist(model_fit$order_list),]
  }

  ## UG efficiently multiplied together
  UG <-  matmult_U(model_fit$U, model_fit$G, model_fit$p,  model_fit$K)

  ## The expensive operation
  const <- sqrt(norm(UG, '2')) # for computational stability
  diag_XUGX <- rowSums((X_block %**% (UG / const)) * X_block) * const

  ## Order it correctly
  diag_XUGX <- diag_XUGX[unlist(model_fit$og_order)]

  ## LOO predictions
  leave_one_out <-
    model_fit$y -
    1/(1 - diag_XUGX) *
    (model_fit$y - model_fit$ytilde)

  return(leave_one_out)
}
