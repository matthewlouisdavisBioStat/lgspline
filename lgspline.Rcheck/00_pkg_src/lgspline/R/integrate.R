#' Definite Integral of a Fitted lgspline
#'
#' @description
#' Given a fitted \code{lgspline} object, computes the definite integral of
#' the fitted surface over a rectangular domain using Gauss--Legendre
#' quadrature on \code{predict()}.
#'
#' @section Method:
#' The integration domain is discretised into a tensor-product grid of
#' Gauss--Legendre quadrature nodes.  Predicted values at each node come
#' from the model's \code{predict()} method, which correctly handles
#' partition assignment and piecewise polynomial evaluation.  The integral
#' is the weighted sum
#' \deqn{\int_{a_1}^{b_1} \cdots \int_{a_d}^{b_d}
#'       \hat{f}(\mathbf{t})\,\mathrm{d}t_1 \cdots \mathrm{d}t_d
#'       \;\approx\;
#'       \sum_{i=1}^{M} w_i\,\hat{f}(\mathbf{t}_i)}
#' where \eqn{M = n_{\mathrm{quad}}^d} and each weight incorporates the
#' Jacobian \eqn{(b_j - a_j)/2} for the affine map from \eqn{[-1, 1]} to
#' \eqn{[a_j, b_j]}.  Nodes and weights on \eqn{[-1, 1]} are computed via
#' the Golub--Welsch algorithm (eigenvalues of the symmetric tridiagonal
#' Jacobi matrix).
#'
#' For smooth polynomials, 30--50 nodes per dimension is typically
#' sufficient; highly partitioned models (large \eqn{K}) may benefit from
#' more.  Total evaluation points scale as \eqn{n_{\mathrm{quad}}^d}, so
#' problems with \eqn{d \ge 4} may require reducing \code{n_quad}.
#'
#' @section Integration scale:
#' By default (\code{link_scale = FALSE}), integration is on the response
#' scale \eqn{\mu = g^{-1}(\eta)}.  Setting \code{link_scale = TRUE}
#' integrates the linear predictor \eqn{\eta = \mathbf{x}^{\top}\boldsymbol{\beta}}
#' directly, which is useful when the quantity of interest is the area
#' under the link-transformed surface rather than the response.  For the
#' identity link the two coincide.
#'
#' @param f A fitted \code{lgspline} object.
#' @param lower Numeric vector of lower bounds, one per integration variable.
#'   Scalar values are recycled.
#' @param upper Numeric vector of upper bounds, one per integration variable.
#'   Scalar values are recycled.
#' @param vars Default: NULL.  Character or integer vector identifying which
#'   predictor(s) to integrate over.  When NULL all numeric predictors are
#'   integrated simultaneously.
#' @param initial_values Default: NULL.  Numeric vector of length \eqn{q}
#'   supplying fixed values for predictors not among the integration
#'   variables.  When NULL, non-integration predictors are held at the
#'   midpoint of their training range.
#' @param B_predict Default: NULL.  Optional list of coefficient vectors,
#'   one per partition.  When NULL the fitted coefficients are used.
#' @param link_scale Default: FALSE.  Logical; when TRUE the integral is
#'   computed on the link (linear predictor) scale \eqn{\eta} rather than
#'   the response scale \eqn{\mu}.
#' @param n_quad Default: 50.  Number of Gauss--Legendre nodes per
#'   integration dimension.
#' @param \ldots Additional arguments (currently unused; present for S3
#'   method compatibility).
#'
#' @return A numeric scalar: the estimated definite integral.
#'
#' @examples
#' ## 1-D: integral of fitted sin(t) over [-pi, pi] should be near 0
#' set.seed(1234)
#' t <- seq(-pi, pi, length.out = 1000)
#' y <- sin(t) + rnorm(length(t), 0, 0.01)
#' fit <- lgspline(t, y, K = 4, opt = FALSE)
#' integrate(fit, lower = -pi, upper = pi)
#'
#' ## Base R integrate still works as expected
#' integrate(sin, lower = -pi, upper = pi)
#'
#' ## 2-D: volume under fitted volcano surface
#' data(volcano)
#' vlong <- cbind(
#'   rep(seq_len(nrow(volcano)), ncol(volcano)),
#'   rep(seq_len(ncol(volcano)), each = nrow(volcano)),
#'   as.vector(volcano)
#' )
#' colnames(vlong) <- c("Length", "Width", "Height")
#' fit_v <- lgspline(vlong[, 1:2], vlong[, 3], K = 18,
#'                   include_quadratic_interactions = TRUE, opt = FALSE)
#' integrate(fit_v, lower = c(1, 1), upper = c(87, 61))
#'
#' @export
integrate.lgspline <- function(
    f,
    lower,
    upper,
    vars           = NULL,
    initial_values = NULL,
    B_predict      = NULL,
    link_scale     = FALSE,
    n_quad         = 50L,
    ...
){

  model_fit <- f

  ## Pull structural info and link function from the fitted model
  if(is.null(B_predict)) B_predict <- model_fit$B

  q_total <- model_fit$q
  link_name <- tryCatch(model_fit$family$link, error = function(e) "identity")
  if(is.null(link_name)) link_name <- "identity"

  ## Original predictor names (if the formula interface was used)
  og_cols <- tryCatch(
    get("og_cols", envir = environment(model_fit$predict)),
    error = function(e) NULL
  )

  ## Resolve integration variables from character names or integer indices
  #  When vars is NULL, all numeric predictors are integrated simultaneously.
  if(is.null(vars)){
    vars_idx <- model_fit$numerics
    if(length(vars_idx) == 0L) vars_idx <- seq_len(q_total)
  } else if(is.character(vars)){
    if(!is.null(og_cols) && all(vars %in% og_cols)){
      vars_idx <- match(vars, og_cols)
    } else {
      stop(
        "\n\tCharacter vars could not be resolved. Use integer predictor ",
        "indices unless original predictor names are available.\n"
      )
    }
  } else {
    vars_idx <- as.integer(vars)
  }

  n_integ <- length(vars_idx)
  if(n_integ == 0L) stop("\n\tNo integration variables identified.\n")
  if(any(is.na(vars_idx)) || any(vars_idx < 1L) || any(vars_idx > q_total)){
    stop("\n\tvars must be valid predictor indices between 1 and ", q_total, ".\n")
  }

  ## Validate and recycle bounds to match the number of integration variables
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  if(length(lower) == 1L) lower <- rep(lower, n_integ)
  if(length(upper) == 1L) upper <- rep(upper, n_integ)
  if(length(lower) != n_integ || length(upper) != n_integ){
    stop(
      "\n\tlower and upper must each have length equal to the number of ",
      "integration variables (", n_integ, ").\n"
    )
  }

  ## Golub--Welsch algorithm for Gauss--Legendre nodes and weights on [-1, 1]
  #  Eigenvalues of the symmetric tridiagonal Jacobi matrix give the nodes;
  #  squared first-row eigenvector entries (times 2) give the weights.
  .gauss_legendre <- function(n){
    i <- seq_len(n - 1L)
    b <- i / sqrt(4 * i^2 - 1)
    J <- matrix(0, n, n)
    for(k in seq_along(b)){
      J[k, k + 1L] <- b[k]
      J[k + 1L, k] <- b[k]
    }
    eig <- eigen(J, symmetric = TRUE)
    ord <- order(eig$values)
    list(nodes = eig$values[ord], weights = (2 * eig$vectors[1, ]^2)[ord])
  }

  ## Reconstruct raw training predictors from the stored design matrices
  #  Extracts the linear column for each predictor from each partition's
  #  X matrix and reassembles into an N x q matrix of raw predictor values.
  .reconstruct_training_predictors <- function(){
    X_parts <- model_fit$X
    ords <- model_fit$order_list
    nm <- model_fit$raw_expansion_names
    pred_mat <- matrix(0, nrow = model_fit$N, ncol = q_total)
    for(j in seq_len(q_total)){
      ci <- which(nm == paste0("_", j, "_"))
      if(length(ci) == 0L) next
      for(k in seq_along(X_parts)){
        if(is.null(X_parts[[k]]) || nrow(X_parts[[k]]) == 0L) next
        pred_mat[as.integer(unlist(ords[[k]])), j] <- X_parts[[k]][, ci[1L]]
      }
    }
    pred_mat
  }

  gl <- .gauss_legendre(as.integer(n_quad))

  ## Affine map from [-1, 1] to [lower[vi], upper[vi]] for each dimension
  mapped_nodes <- vector("list", n_integ)
  mapped_weights <- vector("list", n_integ)
  for(vi in seq_len(n_integ)){
    hw <- (upper[vi] - lower[vi]) / 2
    mp <- (upper[vi] + lower[vi]) / 2
    mapped_nodes[[vi]] <- hw * gl$nodes + mp
    mapped_weights[[vi]] <- gl$weights * hw
  }

  ## Tensor-product grid and corresponding product weights
  grid <- expand.grid(mapped_nodes)
  product_weights <- apply(expand.grid(mapped_weights), 1L, prod)
  n_grid <- nrow(grid)

  ## Fixed values for non-integration predictors
  #  When initial_values is NULL, each predictor is held at the midpoint
  #  of its training range. Otherwise the user-supplied values are used.
  if(is.null(initial_values)){
    training_preds <- .reconstruct_training_predictors()
    base_vec <- rep(0, q_total)
    for(j in seq_len(q_total)){
      xj <- training_preds[, j]
      if(length(unique(xj)) > 1L){
        base_vec[j] <- (min(xj, na.rm = TRUE) + max(xj, na.rm = TRUE)) / 2
      } else {
        base_vec[j] <- xj[1L]
      }
    }
  } else {
    base_vec <- as.numeric(initial_values)
    if(length(base_vec) != q_total){
      stop("\n\tinitial_values must have length ", q_total, ".\n")
    }
  }

  ## Build the evaluation matrix: base values everywhere, overwritten at
  #  integration-variable columns with the quadrature node coordinates
  pred_mat <- matrix(rep(base_vec, each = n_grid), nrow = n_grid)
  for(vi in seq_len(n_integ)) pred_mat[, vars_idx[vi]] <- grid[, vi]

  ## Evaluate the fitted surface and optionally transform to link scale
  predicted <- model_fit$predict(new_predictors = pred_mat,
                                 B_predict = B_predict)
  if(link_scale && link_name != "identity"){
    predicted <- model_fit$family$linkfun(predicted)
  }

  ## Weighted sum gives the quadrature estimate of the definite integral
  sum(as.numeric(predicted) * product_weights)
}
