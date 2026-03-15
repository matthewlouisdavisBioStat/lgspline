#' Compute Deviance for Non-GEE Models
#'
#' Evaluates the mean deviance (or mean squared error as fallback) at
#' the current fitted values.  Used inside the damped Newton-Raphson outer loop of
#' \code{\link{blockfit_solve}} for convergence monitoring.
#'
#' @param y_vec Numeric vector of observed responses.
#' @param mu_vec Numeric vector of fitted means.
#' @param obs_wt Numeric vector of observation weights.
#' @param ord Integer vector of observation indices (passed to
#'   custom deviance residual functions).
#' @param fam GLM family object.
#' @param ... Additional arguments forwarded to
#'   \code{fam$custom_dev.resids}.
#'
#' @return Scalar; the mean deviance contribution.
#'
#' @keywords internal
.bf_deviance <- function(y_vec, mu_vec, obs_wt, ord, fam, ...){
  if(!is.null(fam$custom_dev.resids)){
    mean(fam$custom_dev.resids(y_vec, mu_vec, ord,
                               fam, obs_wt, ...))
  } else if(!is.null(fam$dev.resids)){
    mean(obs_wt * fam$dev.resids(y_vec, cbind(mu_vec), wt = 1))
  } else {
    mean(obs_wt * (y_vec - mu_vec)^2)
  }
}


#' Compute Whitened Deviance for GEE Convergence Monitoring
#'
#' Evaluates deviance in the whitened (decorrelated) space when the
#' family supplies \code{fam$custom_dev.resids}.  In that case the raw
#' deviance residuals are divided by \eqn{\sqrt{W}} and then
#' pre-multiplied by \eqn{\mathbf{V}^{-1/2}} before squaring and
#' averaging.  If only \code{fam$dev.resids} is available, the function
#' falls back to the usual mean deviance (and otherwise mean squared
#' error).  Used in the GEE refinement loop (Case c) of
#' \code{\link{blockfit_solve}}.
#'
#' @inheritParams .bf_deviance
#' @param W_vec Numeric vector of damped Newton-Raphson weights at current iterate.
#' @param VhInv \eqn{\mathbf{V}^{-1/2}} matrix (permuted to
#'   partition ordering).
#'
#' @return Scalar; the mean whitened deviance.
#'
#' @keywords internal
.bf_gee_deviance <- function(y_vec, mu_vec, W_vec, ord,
                             obs_wt, VhInv, fam, ...){
  if(!is.null(fam$custom_dev.resids)){
    raw <- fam$custom_dev.resids(y_vec, mu_vec, ord,
                                 fam, obs_wt, ...)
    W_safe <- pmax(W_vec, sqrt(.Machine$double.eps))
    mean((VhInv %**%
            cbind(sign(raw) * sqrt(abs(raw)) / sqrt(c(W_safe))))^2)
  } else if(!is.null(fam$dev.resids)){
    mean(obs_wt * fam$dev.resids(y_vec, cbind(mu_vec), wt = 1))
  } else {
    mean(obs_wt * (y_vec - mu_vec)^2)
  }
}


#' Numerical Derivative of the Inverse Link Function
#'
#' Returns \eqn{\partial\mu / \partial\eta} for the supplied family
#' object, using the family's own \code{mu.eta} method when available
#' and falling back to analytic forms for common links or central
#' differences otherwise.
#'
#' @param eta Numeric vector of linear predictor values.
#' @param fam GLM family object.
#'
#' @return Numeric vector of the same length as \code{eta}.
#'
#' @keywords internal
.bf_get_mu_eta <- function(eta, fam){
  if(!is.null(fam$mu.eta)) return(fam$mu.eta(eta))
  mu <- fam$linkinv(eta)
  switch(fam$link,
         'identity' = rep(1, length(eta)),
         'log'      = mu,
         'logit'    = mu * (1 - mu),
         'inverse'  = -mu^2,
         'probit'   = dnorm(eta),
         'sqrt'     = 0.5 / sqrt(pmax(mu, .Machine$double.eps)),
         'cloglog'  = exp(eta - exp(eta)),
         {
           eps <- sqrt(.Machine$double.eps)
           (fam$linkinv(eta + eps) -
               fam$linkinv(eta - eps)) / (2 * eps)
         })
}


#' Construct Full Block-Diagonal Penalty Matrix
#'
#' Thin wrapper around \code{.solver_build_lambda_block()}.
#' The block/backfitting solver keeps its original helper name so the
#' surrounding code stays unchanged, while the shared implementation now
#' lives in \code{solver_utils.R}.
#'
#' @param Lambda \eqn{p_expansions \times p_expansions} shared penalty matrix.
#' @param K Integer; number of interior knots.
#' @param unique_penalty_per_partition Logical.
#' @param L_partition_list List of \eqn{K+1} partition-specific
#'   penalty matrices (or zeros).
#'
#' @return A \eqn{(p_expansions \cdot (K+1)) \times (p_expansions \cdot (K+1))} matrix.
#'
#' @keywords internal
.bf_make_Lambda_block <- function(Lambda, K,
                                  unique_penalty_per_partition,
                                  L_partition_list){
  .solver_build_lambda_block(
    Lambda = Lambda,
    K = K,
    unique_penalty_per_partition = unique_penalty_per_partition,
    L_partition_list = L_partition_list
  )
}


#' Split Design and Penalty into Spline and Flat Components
#'
#' Partitions the per-partition design matrices, penalty matrices, and
#' constraint matrix into "spline" (columns receiving partition-specific
#' coefficients) and "flat" (columns receiving a single shared
#' coefficient across all partitions) subsets.
#'
#' @param X List of \eqn{K+1} design matrices (\eqn{N_k \times p_expansions}).
#' @param flat_cols Integer vector of flat column indices.
#' @param p_expansions Integer; total columns per partition.
#' @param K Integer; number of interior knots.
#' @param Lambda \eqn{p_expansions \times p_expansions} shared penalty matrix.
#' @param L_partition_list List of partition-specific penalty matrices.
#' @param A Full \eqn{P \times R_constraints} constraint matrix.
#' @param constraint_values List of constraint right-hand sides.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{X_spline}{List of \eqn{K+1} matrices, spline columns only.}
#'   \item{X_flat}{List of \eqn{K+1} matrices, flat columns only.}
#'   \item{Lambda_spline}{\eqn{nc_s \times nc_s} penalty submatrix.}
#'   \item{Lambda_flat}{\eqn{nc_f \times nc_f} penalty submatrix.}
#'   \item{L_part_spline}{List of partition-specific penalty submatrices
#'     for the spline columns.}
#'   \item{A_spline}{Spline-only constraint matrix (columns pruned
#'     and rank-reduced via pivoted QR).}
#'   \item{nca_spline}{Integer; number of columns in
#'     \code{A_spline}.}
#'   \item{constraint_values_spline}{List of constraint RHS vectors
#'     restricted to spline rows.}
#'   \item{spline_cols}{Integer vector of spline column indices.}
#'   \item{nc_spline}{Integer; number of spline columns.}
#'   \item{nc_flat}{Integer; number of flat columns.}
#'   \item{A_flat}{Flat-only constraint matrix (rows for flat coefficients,
#'     columns pruned and rank-reduced). Used to enforce mixed constraints
#'     on the flat update step.}
#'   \item{nca_flat}{Integer; number of columns in \code{A_flat}.}
#'   \item{A_full}{The original full constraint matrix \code{A}, retained
#'     for mixed-constraint enforcement.}
#'   \item{flat_rows_all}{Integer vector of flat-coefficient row indices
#'     in the full P-space.}
#'   \item{spline_rows_all}{Integer vector of spline-coefficient row indices
#'     in the full P-space.}
#'   \item{has_mixed_constraints}{Logical; TRUE if any constraint column
#'     in A has nonzero entries on both spline and flat rows.}
#'   \item{mixed_constraint_cols}{Integer vector of column indices in the
#'     original A that are mixed (touch both spline and flat rows).}
#' }
#'
#' @keywords internal
.bf_split_components <- function(X, flat_cols, p_expansions, K, Lambda,
                                 L_partition_list, A,
                                 constraint_values){

  spline_cols <- setdiff(1:p_expansions, flat_cols)
  nc_spline <- length(spline_cols)
  nc_flat   <- length(flat_cols)

  X_spline <- lapply(X, function(Xk) Xk[, spline_cols, drop = FALSE])
  X_flat   <- lapply(X, function(Xk) Xk[, flat_cols,   drop = FALSE])

  Lambda_spline <- Lambda[spline_cols, spline_cols]
  Lambda_flat   <- Lambda[flat_cols,   flat_cols]

  L_part_spline <- lapply(L_partition_list, function(Lk){
    if(is.numeric(Lk) && length(Lk) == 1 && Lk == 0) return(0)
    Lk[spline_cols, spline_cols]
  })

  ## Extract spline-only and flat-only rows from full constraint matrix A
  flat_rows_all   <- unlist(lapply(0:K, function(k) k * p_expansions + flat_cols))
  spline_rows_all <- setdiff(1:(p_expansions * (K + 1)), flat_rows_all)
  A_spline <- A[spline_rows_all, , drop = FALSE]
  A_flat   <- A[flat_rows_all,   , drop = FALSE]

  ## Detect mixed constraints: columns of A that have nonzero entries
  #  on BOTH spline and flat rows
  spline_nz <- colSums(abs(A_spline)) > sqrt(.Machine$double.eps)
  flat_nz   <- colSums(abs(A_flat))   > sqrt(.Machine$double.eps)
  mixed_constraint_cols <- which(spline_nz & flat_nz)
  has_mixed_constraints <- length(mixed_constraint_cols) > 0

  ## Process A_spline: keep only columns with nonzero spline entries
  keep_cols <- which(colSums(abs(A_spline)) > sqrt(.Machine$double.eps))
  if(length(keep_cols) > 0){
    A_spline <- A_spline[, keep_cols, drop = FALSE]
    qr_As <- qr(A_spline)
    if(qr_As$rank < ncol(A_spline)){
      A_spline <- qr.Q(qr_As)[, 1:qr_As$rank, drop = FALSE]
    }
  } else {
    A_spline <- cbind(rep(0, nc_spline * (K + 1)))
  }
  nca_spline <- ncol(A_spline)

  ## Process A_flat: keep only columns with nonzero flat entries
  keep_cols_flat <- which(colSums(abs(A_flat)) > sqrt(.Machine$double.eps))
  if(length(keep_cols_flat) > 0){
    A_flat_kept <- A_flat[, keep_cols_flat, drop = FALSE]
    qr_Af <- qr(A_flat_kept)
    if(qr_Af$rank < ncol(A_flat_kept)){
      A_flat_kept <- qr.Q(qr_Af)[, 1:qr_Af$rank, drop = FALSE]
    }
    A_flat <- A_flat_kept
  } else {
    A_flat <- cbind(rep(0, length(flat_rows_all)))
  }
  nca_flat <- ncol(A_flat)

  ## Constraint values for spline-only subproblem
  if(length(constraint_values) > 0){
    constraint_values_spline <- lapply(constraint_values, function(cv){
      cv[spline_cols, , drop = FALSE]
    })
  } else {
    constraint_values_spline <- constraint_values
  }

  list(X_spline = X_spline,
       X_flat   = X_flat,
       Lambda_spline = Lambda_spline,
       Lambda_flat   = Lambda_flat,
       L_part_spline = L_part_spline,
       A_spline = A_spline,
       nca_spline = nca_spline,
       constraint_values_spline = constraint_values_spline,
       spline_cols = spline_cols,
       nc_spline = nc_spline,
       nc_flat   = nc_flat,
       A_flat = A_flat,
       nca_flat = nca_flat,
       A_full = A,
       flat_rows_all = flat_rows_all,
       spline_rows_all = spline_rows_all,
       has_mixed_constraints = has_mixed_constraints,
       mixed_constraint_cols = mixed_constraint_cols)
}


#' Constrained Flat Update Given Current Spline Solution
#'
#' Solves for flat coefficients subject to the residual equality
#' constraint imposed by the full constraint system, conditional on
#' the current spline block solution.
#'
#' When the full constraint system is \eqn{\mathbf{A}^\top \boldsymbol{\beta} = \mathbf{c}},
#' and we partition \eqn{\boldsymbol{\beta}} into spline and flat blocks, the flat
#' update must satisfy:
#' \deqn{
#' \mathbf{A}_{\mathrm{flat}}^\top \boldsymbol{\beta}_{\mathrm{flat}}
#' = \mathbf{c} - \mathbf{A}_{\mathrm{spline}}^\top \boldsymbol{\beta}_{\mathrm{spline}}
#' }
#'
#' This function solves the penalized least-squares problem for flat coefficients
#' subject to this residual equality constraint using a Lagrangian approach.
#'
#' @param XfWXf_pen Penalized flat Gram matrix
#'   (\eqn{nc_f \times nc_f}).
#' @param Xfr Right-hand side cross-product for flat update
#'   (\eqn{nc_f \times 1}).
#' @param A_full Full constraint matrix (\eqn{P \times R}).
#' @param constraint_values List of constraint RHS vectors.
#' @param beta_spline List of \eqn{K+1} current spline coefficient vectors.
#' @param flat_rows_all Integer vector of flat row indices in the full P-space.
#' @param spline_rows_all Integer vector of spline row indices in the full P-space.
#' @param spline_cols Integer vector of spline column indices within each partition.
#' @param flat_cols Integer vector of flat column indices within each partition.
#' @param p_expansions Integer; columns per partition.
#' @param K Integer; number of interior knots.
#' @param nc_flat Integer; number of flat columns.
#'
#' @return Numeric vector of flat coefficients (\eqn{nc_f \times 1}).
#'
#' @keywords internal
.bf_constrained_flat_update <- function(XfWXf_pen, Xfr,
                                        A_full, constraint_values,
                                        beta_spline,
                                        flat_rows_all, spline_rows_all,
                                        spline_cols, flat_cols,
                                        p_expansions, K, nc_flat){

  ## Build the full spline coefficient vector in full P-space ordering
  beta_spline_full <- rep(0, p_expansions * (K + 1))
  for(k in 1:(K + 1)){
    rows_k <- (k - 1) * p_expansions + spline_cols
    beta_spline_full[rows_k] <- beta_spline[[k]]
  }

  ## The full constraint is: A^T beta = c
  #  Partition: A_s^T beta_s + A_f^T beta_f = c
  #  So the flat constraint is: A_f^T beta_f = c - A_s^T beta_s
  #
  #  But flat coefficients are shared across partitions, so beta_f in
  #  the full P-space has the same nc_flat values replicated K+1 times.
  #  We need to express the constraint in terms of the unique nc_flat
  #  flat coefficients.

  ## Compute the spline contribution to each constraint
  spline_contribution <- crossprod(A_full, beta_spline_full)  # R x 1

  ## Get the constraint target
  if(length(constraint_values) > 0 && any(unlist(constraint_values) != 0)){
    cv_vec <- Reduce("rbind", constraint_values)
    constraint_target <- crossprod(A_full, cv_vec)  # R x 1
  } else {
    constraint_target <- rep(0, ncol(A_full))
  }

  ## Residual that the flat block must satisfy
  flat_target <- constraint_target - spline_contribution  # R x 1

  ## Build the flat constraint matrix in terms of unique flat coefficients.
  #  For each constraint column j of A, the flat contribution is:
  #    sum over k=0..K of A[flat_rows_k, j]^T * beta_flat_unique
  #  Since beta_flat is replicated, this equals:
  #    (sum over k of A[flat_rows_k, j])^T * beta_flat_unique
  A_flat_unique <- matrix(0, nrow = nc_flat, ncol = ncol(A_full))
  for(k in 1:(K + 1)){
    rows_k <- (k - 1) * p_expansions + flat_cols
    A_flat_unique <- A_flat_unique + A_full[rows_k, , drop = FALSE]
  }
  ## A_flat_unique is nc_flat x R

  ## Only keep constraints that actually involve flat coefficients
  active_cols <- which(colSums(abs(A_flat_unique)) > sqrt(.Machine$double.eps))

  if(length(active_cols) == 0){
    ## No constraints touch flat coefficients; unconstrained update
    return(c(invert(XfWXf_pen) %**% Xfr))
  }

  A_fc <- A_flat_unique[, active_cols, drop = FALSE]  # nc_flat x n_active
  b_fc <- flat_target[active_cols]                     # n_active x 1

  ## Solve the constrained problem via KKT / Lagrangian:
  #    minimize  0.5 * beta_f^T X'X beta_f - beta_f^T g
  #    subject to  A_fc^T beta_f = b_fc
  #
  #  KKT system:
  #    [ X'X    A_fc ] [ beta_f ] = [ g    ]
  #    [ A_fc^T  0 ] [ lambda ] = [ b_fc ]

  H <- XfWXf_pen
  g <- Xfr

  n_constr <- ncol(A_fc)

  ## Build and solve KKT system
  KKT <- rbind(
    cbind(H, A_fc),
    cbind(t(A_fc), matrix(0, n_constr, n_constr))
  )
  rhs <- rbind(g, cbind(b_fc))

  ## Use a robust solve
  kkt_sol <- try(solve(KKT, rhs), silent = TRUE)

  if(inherits(kkt_sol, "try-error")){
    ## If KKT is singular, try pseudoinverse approach
    #  This can happen if constraints are redundant
    kkt_sol <- try({
      ## Unconstrained solution
      beta_unc <- c(invert(H) %**% g)
      ## Project onto constraint
      residual <- b_fc - crossprod(A_fc, beta_unc)
      ## Correction
      HinvA <- invert(H) %**% A_fc
      correction <- HinvA %**% solve(crossprod(A_fc, HinvA), residual)
      cbind(beta_unc + c(correction))
    }, silent = TRUE)

    if(inherits(kkt_sol, "try-error")){
      ## Last resort: unconstrained
      return(c(invert(H) %**% g))
    }
    return(c(kkt_sol[1:nc_flat]))
  }

  c(kkt_sol[1:nc_flat])
}


#' Lagrangian Projection for the Spline-Only Subproblem
#'
#' Projects the adjusted cross-product \eqn{\mathbf{X}_s^{\top}
#' \tilde{\mathbf{y}}} through \eqn{\mathbf{G}_s^{1/2}} and
#' removes the component along \eqn{\mathbf{G}_s^{1/2}\mathbf{A}_s}
#' to enforce the smoothness constraints in the spline block.
#' When \code{cv_spline} contains nonzero values (e.g.\ from
#' \code{no_intercept}), the corresponding contribution is added
#' back after projection.
#'
#' @param Xy_adj List of \eqn{K+1} adjusted cross-product vectors
#'   (\eqn{nc_s \times 1}).
#' @param Ghalf_cur List of \eqn{K+1} \eqn{\mathbf{G}_s^{1/2}}
#'   matrices (each \eqn{nc_s \times nc_s}).
#' @param GhalfA_cur \eqn{(nc_s \cdot (K+1)) \times nca_s} matrix
#'   formed by row-stacking \eqn{\mathbf{G}_s^{1/2,(k)}
#'   \mathbf{A}_s^{(k)}}.
#' @param cv_spline List of constraint right-hand-side vectors for
#'   the spline subproblem (may be empty).
#' @param K Integer; number of interior knots.
#' @param nc_spline Integer; number of spline columns.
#' @param A_spline Spline-only constraint matrix.
#'
#' @return List of \eqn{K+1} coefficient vectors
#'   (\eqn{nc_s \times 1}).
#'
#' @keywords internal
.bf_lagrangian_project <- function(Xy_adj, Ghalf_cur, GhalfA_cur,
                                   cv_spline, K, nc_spline,
                                   A_spline){
  GhalfXy <- cbind(unlist(lapply(1:(K + 1), function(k){
    Ghalf_cur[[k]] %**% Xy_adj[[k]]
  })))

  sc <- 1 / sqrt(K + 1)
  resids_star <- .lm.fit(GhalfA_cur * sc,
                         GhalfXy * sc)$residuals / sc

  if(length(cv_spline) > 0){
    if(any(unlist(cv_spline) != 0)){
      cv_vec <- Reduce("rbind", cv_spline)
      preds_star <- GhalfA_cur %**%
        (invert(crossprod(GhalfA_cur) * sc) %**%
           (crossprod(A_spline, cv_vec) * sc))
      resids_star <- resids_star + c(preds_star)
    }
  }

  lapply(1:(K + 1), function(k){
    rows <- ((k - 1) * nc_spline + 1):(k * nc_spline)
    Ghalf_cur[[k]] %**% cbind(resids_star[rows])
  })
}


#' Gaussian Identity + GEE: Closed-Form Solve (Case a)
#'
#' When the response is Gaussian with identity link and a working
#' correlation structure is present, the whitened system
#' \eqn{\mathbf{V}^{-1/2}\mathbf{X}} is not block-diagonal, so
#' partition-wise backfitting does not apply.
#' This function performs a single closed-form solve that replicates
#' \code{get_B} Path 1a but in the backfitting context.
#'
#' @param X List of \eqn{K+1} unwhitened design matrices.
#' @param y List of \eqn{K+1} response vectors.
#' @param K Integer; number of interior knots.
#' @param p_expansions Integer; columns per partition.
#' @param order_list List of observation-index vectors by partition.
#' @param VhalfInv \eqn{\mathbf{V}^{-1/2}} matrix in the original
#'   observation ordering.
#' @param Lambda Shared penalty matrix.
#' @param L_partition_list Partition-specific penalty matrices.
#' @param unique_penalty_per_partition Logical.
#' @param A Full constraint matrix.
#' @param constraint_values List of constraint RHS vectors.
#' @param spline_cols,flat_cols Integer vectors of column indices.
#'
#' @return A named list:
#' \describe{
#'   \item{beta_spline}{List of \eqn{K+1} spline coefficient vectors.}
#'   \item{beta_flat}{Numeric vector of flat coefficients.}
#'   \item{result}{List of \eqn{K+1} full per-partition coefficient
#'     vectors.}
#' }
#'
#' @keywords internal
.bf_case_gauss_gee <- function(X, y, K, p_expansions, order_list, VhalfInv,
                               Lambda, L_partition_list,
                               unique_penalty_per_partition,
                               A, constraint_values,
                               spline_cols, flat_cols,
                               observation_weights){

  perm_bf <- unlist(order_list)
  VhalfInv_bf <- VhalfInv[perm_bf, perm_bf]
  obs_wt_bf <- c(unlist(observation_weights))
  if(length(obs_wt_bf) == 0L){
    obs_wt_bf <- rep(1, nrow(VhalfInv_bf))
  } else if(length(obs_wt_bf) == 1L){
    obs_wt_bf <- rep(obs_wt_bf, nrow(VhalfInv_bf))
  }

  X_block_bf <- collapse_block_diagonal(X)
  y_block_bf <- cbind(unlist(y))

  X_tilde_bf <- VhalfInv_bf %**% X_block_bf
  y_tilde_bf <- VhalfInv_bf %**% y_block_bf

  Gram_bf <- crossprod(X_tilde_bf, obs_wt_bf * X_tilde_bf)
  Lambda_block_bf <- .bf_make_Lambda_block(Lambda, K,
                                           unique_penalty_per_partition,
                                           L_partition_list)
  G_bf <- invert(Gram_bf + Lambda_block_bf)
  G_bf_half <- matsqrt(G_bf)
  Xy_bf <- crossprod(X_tilde_bf, obs_wt_bf * y_tilde_bf)

  ## Lagrangian projection in full P-space
  y_star <- G_bf_half %**% Xy_bf
  X_star <- G_bf_half %**% A

  sc <- 1 / sqrt(K + 1)
  resids_bf <- .lm.fit(X_star * sc, y_star * sc)$residuals / sc

  if(length(constraint_values) > 0){
    if(any(unlist(constraint_values) != 0)){
      cv_vec <- Reduce("rbind", constraint_values)
      preds_bf <- X_star %**%
        (invert(crossprod(X_star) * sc) %**%
           (crossprod(A, cv_vec) * sc))
      resids_bf <- resids_bf + c(preds_bf)
    }
  }

  beta_block <- G_bf_half %**% cbind(resids_bf)

  result <- lapply(1:(K + 1), function(k){
    cbind(beta_block[(k - 1) * p_expansions + 1:p_expansions])
  })

  beta_spline <- lapply(result, function(b) b[spline_cols, , drop = FALSE])
  beta_flat   <- result[[1]][flat_cols]

  list(beta_spline = beta_spline,
       beta_flat   = beta_flat,
       result      = result)
}


#' Gaussian Identity Backfitting Without Correlation (Case b)
#'
#' Standard block-coordinate descent on the convex quadratic
#' objective, alternating between the spline step (Lagrangian
#' projection) and the flat step (pooled penalized regression).
#'
#' @param X_spline,X_flat Lists of per-partition submatrices.
#' @param y List of response vectors.
#' @param K Integer.
#' @param Ghalf_sp List of \eqn{\mathbf{G}_s^{1/2}} matrices.
#' @param GhalfA_sp Pre-computed \eqn{\mathbf{G}_s^{1/2}
#'   \mathbf{A}_s}.
#' @param XfXf_pen_inv Penalised inverse of the pooled flat Gram
#'   matrix.
#' @param constraint_values_spline Spline-only constraint RHS.
#' @param nc_spline,nc_flat Integers.
#' @param A_spline Spline-only constraint matrix.
#' @param tol Convergence tolerance.
#' @param max_backfit_iter Maximum iterations.
#' @param verbose Logical.
#' @param split Output of \code{.bf_split_components}, needed for
#'   mixed-constraint enforcement on the flat step.
#' @param constraint_values Full constraint RHS list.
#' @param Lambda_flat Flat penalty submatrix.
#' @param p_expansions Integer; columns per partition.
#'
#' @return A named list with \code{beta_spline} and
#'   \code{beta_flat}.
#'
#' @keywords internal
.bf_case_gauss_no_corr <- function(X_spline, X_flat, y, K,
                                   Ghalf_sp, GhalfA_sp,
                                   XfXf_pen_inv,
                                   constraint_values_spline,
                                   nc_spline, nc_flat,
                                   A_spline,
                                   tol, max_backfit_iter,
                                   verbose,
                                   split = NULL,
                                   constraint_values = list(),
                                   Lambda_flat = NULL,
                                   p_expansions = NULL){

  beta_flat <- rep(0, nc_flat)
  beta_spline <- lapply(1:(K + 1), function(k) cbind(rep(0, nc_spline)))

  ## Determine whether we need constrained flat updates
  use_constrained_flat <- (!is.null(split) &&
                             split$has_mixed_constraints &&
                             !is.null(Lambda_flat) &&
                             !is.null(p_expansions))

  for(bf_iter in 1:max_backfit_iter){

    ## Spline step
    Xy_adj <- lapply(1:(K + 1), function(k){
      y_adj_k <- y[[k]] - X_flat[[k]] %**% cbind(beta_flat)
      crossprod(X_spline[[k]], y_adj_k)
    })
    beta_spline_new <- .bf_lagrangian_project(
      Xy_adj, Ghalf_sp, GhalfA_sp,
      constraint_values_spline, K, nc_spline, A_spline)

    ## Flat step
    Xfr <- Reduce("+", lapply(1:(K + 1), function(k){
      r_k <- y[[k]] - X_spline[[k]] %**% beta_spline_new[[k]]
      crossprod(X_flat[[k]], r_k)
    }))

    if(use_constrained_flat){
      ## Constrained flat update enforcing mixed equality constraints
      XfXf <- Reduce("+", lapply(1:(K + 1), function(k){
        crossprod(X_flat[[k]])
      }))
      XfXf_pen <- XfXf + Lambda_flat

      beta_flat_new <- .bf_constrained_flat_update(
        XfWXf_pen = XfXf_pen,
        Xfr = Xfr,
        A_full = split$A_full,
        constraint_values = constraint_values,
        beta_spline = beta_spline_new,
        flat_rows_all = split$flat_rows_all,
        spline_rows_all = split$spline_rows_all,
        spline_cols = split$spline_cols,
        flat_cols = setdiff(1:p_expansions, split$spline_cols),
        p_expansions = p_expansions,
        K = K,
        nc_flat = nc_flat)
    } else {
      beta_flat_new <- c(XfXf_pen_inv %**% Xfr)
    }

    ## Convergence check
    spline_diff <- max(abs(unlist(beta_spline_new) -
                             unlist(beta_spline)))
    flat_diff   <- max(abs(beta_flat_new - beta_flat))
    beta_spline <- beta_spline_new
    beta_flat   <- beta_flat_new

    if(verbose){
      cat('  Backfit iter', bf_iter,
          '| spline diff:', format(spline_diff, digits = 4),
          '| flat diff:',   format(flat_diff,   digits = 4), '\n')
    }

    if(max(spline_diff, flat_diff) < tol) break
  }

  list(beta_spline = beta_spline,
       beta_flat   = beta_flat)
}


#' GLM Backfitting Inner Loop (Weighted Working Response)
#'
#' Given damped Newton-Raphson weights and working response partitioned into lists,
#' runs the inner backfitting loop that alternates between a weighted
#' spline step and a weighted flat step until convergence.  Shared by
#' both Case c (GEE warm start) and Case d (GLM without GEE).
#'
#' @param X_spline,X_flat Lists of per-partition submatrices.
#' @param z_list List of \eqn{K+1} working-response vectors.
#' @param W_list List of \eqn{K+1} weight vectors.
#' @param K Integer.
#' @param Lambda_spline,Lambda_flat Penalty submatrices.
#' @param L_part_spline Partition-specific spline penalty list.
#' @param unique_penalty_per_partition Logical.
#' @param A_spline Spline-only constraint matrix.
#' @param constraint_values_spline Constraint RHS list.
#' @param nc_spline,nc_flat Integers.
#' @param beta_spline_init,beta_flat_init Initial coefficient values.
#' @param tol Convergence tolerance.
#' @param max_backfit_iter Maximum iterations.
#' @param parallel_eigen,cl,chunk_size,num_chunks,rem_chunks
#'   Parallel arguments forwarded to \code{compute_G_eigen}.
#' @param split Output of \code{.bf_split_components}, needed for
#'   mixed-constraint enforcement on the flat step.
#' @param constraint_values Full constraint RHS list.
#' @param p_expansions Integer; columns per partition.
#'
#' @return A named list with \code{beta_spline} and
#'   \code{beta_flat}.
#'
#' @keywords internal
.bf_inner_weighted <- function(X_spline, X_flat,
                               z_list, W_list, K,
                               Lambda_spline, Lambda_flat,
                               L_part_spline,
                               unique_penalty_per_partition,
                               A_spline, constraint_values_spline,
                               nc_spline, nc_flat,
                               beta_spline_init, beta_flat_init,
                               tol, max_backfit_iter,
                               parallel_eigen, cl,
                               chunk_size, num_chunks,
                               rem_chunks,
                               split = NULL,
                               constraint_values = list(),
                               p_expansions = NULL){

  beta_spline <- beta_spline_init
  beta_flat   <- beta_flat_init
  schur_zero <- lapply(1:(K + 1), function(k) 0)

  ## Determine whether we need constrained flat updates
  use_constrained_flat <- (!is.null(split) &&
                             split$has_mixed_constraints &&
                             !is.null(p_expansions))

  ## Weighted spline Gram and G
  X_gram_sp_w <- lapply(1:(K + 1), function(k){
    crossprod(X_spline[[k]] * sqrt(W_list[[k]]))
  })
  G_sp_w <- compute_G_eigen(X_gram_sp_w,
                            Lambda_spline, K,
                            parallel_eigen, cl,
                            chunk_size, num_chunks, rem_chunks,
                            gaussian(),
                            unique_penalty_per_partition,
                            L_part_spline,
                            keep_G = TRUE,
                            schur_zero)
  Ghalf_sp_w <- G_sp_w$Ghalf

  GhalfA_sp_w <- Reduce("rbind", lapply(1:(K + 1), function(k){
    rows <- ((k - 1) * nc_spline + 1):(k * nc_spline)
    Ghalf_sp_w[[k]] %**% A_spline[rows, , drop = FALSE]
  }))

  for(bf_iter in 1:max_backfit_iter){

    ## Spline step (weighted)
    Xy_adj <- lapply(1:(K + 1), function(k){
      z_adj_k <- z_list[[k]] - X_flat[[k]] %**% cbind(beta_flat)
      crossprod(X_spline[[k]] * W_list[[k]], cbind(z_adj_k))
    })
    beta_spline_new <- .bf_lagrangian_project(
      Xy_adj, Ghalf_sp_w, GhalfA_sp_w,
      constraint_values_spline, K, nc_spline, A_spline)

    ## Flat step (weighted, pooled)
    XfWXf <- Reduce("+", lapply(1:(K + 1), function(k){
      crossprod(X_flat[[k]] * sqrt(W_list[[k]]))
    }))
    XfWXf_pen <- XfWXf + Lambda_flat
    Xfr <- Reduce("+", lapply(1:(K + 1), function(k){
      r_k <- z_list[[k]] - X_spline[[k]] %**% beta_spline_new[[k]]
      crossprod(X_flat[[k]] * W_list[[k]], cbind(r_k))
    }))

    if(use_constrained_flat){
      ## Constrained flat update enforcing mixed equality constraints
      beta_flat_new <- .bf_constrained_flat_update(
        XfWXf_pen = XfWXf_pen,
        Xfr = Xfr,
        A_full = split$A_full,
        constraint_values = constraint_values,
        beta_spline = beta_spline_new,
        flat_rows_all = split$flat_rows_all,
        spline_rows_all = split$spline_rows_all,
        spline_cols = split$spline_cols,
        flat_cols = setdiff(1:p_expansions, split$spline_cols),
        p_expansions = p_expansions,
        K = K,
        nc_flat = nc_flat)
    } else {
      XfWXf_pen_inv <- invert(XfWXf_pen)
      beta_flat_new <- c(XfWXf_pen_inv %**% Xfr)
    }

    flat_diff <- max(abs(beta_flat_new - beta_flat))
    beta_spline <- beta_spline_new
    beta_flat   <- beta_flat_new

    if(flat_diff < tol) break
  }

  list(beta_spline = beta_spline,
       beta_flat   = beta_flat)
}


#' Damped SQP Refinement on the Full (Optionally Whitened) System
#'
#' Runs a damped sequential quadratic programming loop using
#' \code{quadprog::solve.QP} at each iteration, enforcing all
#' smoothness equalities, flat-equality constraints, and any
#' user-supplied inequality constraints.
#'
#' This function serves two roles within \code{\link{blockfit_solve}}:
#' \enumerate{
#'   \item \strong{GEE Stage 2} (Case c): operates on the whitened
#'     system \eqn{\tilde{\mathbf{X}} = \mathbf{V}^{-1/2}
#'     \mathbf{X}_{\mathrm{block}}}, initialized from the damped Newton-Raphson
#'     backfitting warm start.
#'   \item \strong{Non-GEE SQP refinement}: operates on the
#'     unwhitened block-diagonal design, applying inequality
#'     constraints after backfitting convergence.
#' }
#'
#' @param X_design \eqn{N \times P} design matrix (whitened for GEE,
#'   unwhitened otherwise).
#' @param y_design \eqn{N \times 1} response vector (whitened for
#'   GEE, unwhitened otherwise).
#' @param X_block_raw Unwhitened block-diagonal design matrix (used
#'   for computing the linear predictor on the original scale).
#'   For the non-GEE path this equals \code{X_design}.
#' @param beta_init \eqn{P \times 1} initial coefficient vector.
#' @param Lambda_block \eqn{P \times P} block-diagonal penalty.
#' @param qp_Amat_combined Combined constraint matrix (equalities
#'   then inequalities).
#' @param qp_bvec_combined Combined constraint RHS.
#' @param qp_meq_combined Number of leading equality constraints.
#' @param K,p_expansions Integers.
#' @param family GLM family object.
#' @param order_list Observation-index lists.
#' @param glm_weight_function,schur_correction_function,
#'   qp_score_function Functions.
#' @param need_dispersion_for_estimation Logical.
#' @param dispersion_function Function.
#' @param observation_weights List or vector of weights.
#' @param iterate Logical; if FALSE, takes at most two steps.
#' @param tol Convergence tolerance.
#' @param VhalfInv,VhalfInv_perm \eqn{\mathbf{V}^{-1/2}} matrices
#'   (original and permuted ordering).  NULL for non-GEE path.
#' @param is_gee Logical; TRUE when called from the GEE refinement
#'   path (affects deviance computation).
#' @param deviance_fun Function computing deviance; either
#'   \code{.bf_deviance} or \code{.bf_gee_deviance}.
#' @param X_partitions List of per-partition design matrices
#'   (unwhitened); needed for Schur correction in non-GEE path.
#' @param y_partitions List of per-partition response vectors
#'   (unwhitened); needed for Schur correction in non-GEE path.
#' @param verbose Logical.
#' @param ... Forwarded to weight, dispersion, and score functions.
#'
#' @return A named list with elements \code{beta_block} (final
#'   coefficient vector), \code{result} (per-partition coefficient
#'   list), \code{last_qp_sol} (last successful \code{solve.QP}
#'   output or NULL), \code{converged} (logical), and
#'   \code{final_deviance} (scalar).
#'
#' @keywords internal
.bf_sqp_loop <- function(X_design, y_design, X_block_raw,
                         beta_init, Lambda_block,
                         qp_Amat_combined, qp_bvec_combined,
                         qp_meq_combined,
                         K, p_expansions, family, order_list,
                         glm_weight_function,
                         schur_correction_function,
                         qp_score_function,
                         need_dispersion_for_estimation,
                         dispersion_function,
                         observation_weights,
                         iterate, tol,
                         VhalfInv, VhalfInv_perm,
                         is_gee, deviance_fun,
                         X_partitions, y_partitions,
                         verbose, ...){

  beta_block <- beta_init
  damp_cnt   <- 0
  master_cnt <- 0
  err <- Inf
  converged <- FALSE
  last_qp_sol <- NULL

  XB <- X_block_raw %**% beta_block
  y_block <- if(is_gee) cbind(unlist(y_partitions)) else y_design

  while(err > tol & damp_cnt < 10 & master_cnt < 100){
    master_cnt <- master_cnt + 1
    damp <- 2^(-(damp_cnt))

    ## Dispersion
    if(need_dispersion_for_estimation){
      dispersion_temp <- dispersion_function(
        mu = family$linkinv(XB),
        y = y_block,
        order_indices = unlist(order_list),
        family = family,
        observation_weights = unlist(observation_weights),
        VhalfInv = VhalfInv,
        ...
      )
    } else {
      dispersion_temp <- 1
    }

    ## Working weights
    W <- c(glm_weight_function(
      family$linkinv(XB), y_block,
      unlist(order_list), family, dispersion_temp,
      unlist(observation_weights), ...
    ))
    W <- pmax(W, sqrt(.Machine$double.eps))

    ## Schur correction
    result_temp <- lapply(1:(K + 1), function(k){
      cbind(beta_block[(k - 1) * p_expansions + 1:p_expansions])
    })

    if(is_gee){
      schur_correction <- schur_correction_function(
        list(X_block_raw), list(y_block),
        list(cbind(unlist(result_temp))),
        dispersion_temp, list(unlist(order_list)),
        0, family, unlist(observation_weights), ...
      )
    } else {
      schur_correction <- schur_correction_function(
        X_partitions, y_partitions,
        result_temp, dispersion_temp,
        order_list, K, family,
        observation_weights, ...
      )
    }

    if(!(any(unlist(schur_correction) != 0))){
      schur_correction_coll <- 0
    } else {
      schur_correction_coll <- collapse_block_diagonal(schur_correction)
    }

    ## Information matrix
    info <- crossprod(X_design, W * X_design) +
      Lambda_block + schur_correction_coll

    sc <- sqrt(mean(abs(info)))

    ## Score
    if(is_gee){
      qp_score <- qp_score_function(
        X_design, y_design,
        VhalfInv_perm %**% cbind(family$linkinv(c(XB))),
        unlist(order_list), dispersion_temp,
        VhalfInv_perm, unlist(observation_weights), ...
      )
    } else {
      qp_score <- qp_score_function(
        X_design, y_design,
        cbind(family$linkinv(XB)),
        unlist(order_list), dispersion_temp,
        NULL, unlist(observation_weights), ...
      )
    }

    ## Solve QP
    qp_sol <- try({quadprog::solve.QP(
      Dmat = info / sc,
      dvec = (qp_score -
                Lambda_block %**% beta_block +
                info %**% beta_block) / sc,
      Amat = qp_Amat_combined,
      bvec = qp_bvec_combined,
      meq  = qp_meq_combined
    )}, silent = TRUE)

    if(any(inherits(qp_sol, 'try-error'))){
      if(verbose){
        cat("    SQP iter", master_cnt,
            "- solve.QP failed, retaining current\n")
      }
      beta_new <- beta_block
    } else {
      last_qp_sol <- qp_sol
      beta_new    <- cbind(qp_sol$solution)
    }

    ## Step acceptance
    if(!iterate & master_cnt > 1){
      beta_block <- beta_new
      converged  <- TRUE
      damp_cnt   <- 11
      master_cnt <- 101
      err <- tol - 1
    } else {
      beta_new_damped <- (1 - damp) * beta_block + damp * beta_new
      XB <- X_block_raw %**% beta_new_damped

      if(need_dispersion_for_estimation){
        dispersion_temp <- dispersion_function(
          mu = family$linkinv(XB),
          y = y_block,
          order_indices = unlist(order_list),
          family = family,
          observation_weights = unlist(observation_weights),
          VhalfInv = VhalfInv,
          ...
        )
      } else {
        dispersion_temp <- 1
      }

      W <- c(glm_weight_function(
        family$linkinv(XB), y_block,
        unlist(order_list), family, dispersion_temp,
        unlist(observation_weights), ...
      ))
      W <- pmax(W, sqrt(.Machine$double.eps))

      if(is_gee){
        err_new <- deviance_fun(
          y_block, family$linkinv(c(XB)), W,
          unlist(order_list), unlist(observation_weights),
          VhalfInv_perm, family, ...
        )
      } else {
        err_new <- deviance_fun(
          y_block, family$linkinv(c(XB)),
          unlist(observation_weights),
          unlist(order_list), family, ...
        )
      }

      if(is.null(err_new) | is.na(err_new) | !is.finite(err_new)){
        damp_cnt <- damp_cnt + 1
      } else if(err_new <= err){
        prev_err <- err
        err <- err_new
        abs_diff <- max(abs(beta_new_damped - beta_block))
        beta_block <- beta_new_damped
        damp_cnt <- max(0, damp_cnt - 1)

        if((abs_diff < tol) &
           (prev_err - err < tol) &
           (master_cnt > 5)){
          converged  <- TRUE
          damp_cnt   <- 11
          master_cnt <- 101
          err <- tol - 1
        }
      } else {
        damp_cnt <- damp_cnt + 1
      }
    }

    if(verbose & (master_cnt <= 100)){
      cat("    SQP iter", master_cnt,
          "| deviance:", format(err, digits = 6),
          "| damp:", format(damp, digits = 4), "\n")
    }

    XB <- X_block_raw %**% beta_block
  }

  result <- lapply(1:(K + 1), function(k){
    cbind(beta_block[1:p_expansions + (k - 1) * p_expansions])
  })

  list(beta_block = beta_block,
       result     = result,
       last_qp_sol = last_qp_sol,
       converged  = converged,
       final_deviance = err)
}


#' Assemble qp_info from a solve.QP Solution
#'
#' Thin wrapper around \code{.solver_assemble_qp_info()}.
#' Keeping the local helper name avoids changing the public return shape
#' or downstream references inside \code{blockfit_solve()}.
#'
#' @param last_qp_sol Output of \code{quadprog::solve.QP}, or NULL.
#' @param beta_block Final coefficient vector.
#' @param qp_Amat_combined Combined constraint matrix.
#' @param qp_bvec_combined Combined constraint RHS.
#' @param qp_meq_combined Number of equalities.
#' @param converged Logical.
#' @param final_deviance Scalar.
#' @param info_matrix Information matrix (non-GEE path only; NULL
#'   otherwise).
#'
#' @return A list suitable for the \code{qp_info} slot, or NULL if
#'   \code{last_qp_sol} is NULL.
#'
#' @keywords internal
.bf_assemble_qp_info <- function(last_qp_sol, beta_block,
                                 qp_Amat_combined,
                                 qp_bvec_combined,
                                 qp_meq_combined,
                                 converged, final_deviance,
                                 info_matrix = NULL){
  .solver_assemble_qp_info(
    last_qp_sol = last_qp_sol,
    beta_block = beta_block,
    qp_Amat_combined = qp_Amat_combined,
    qp_bvec_combined = qp_bvec_combined,
    qp_meq_combined = qp_meq_combined,
    converged = converged,
    final_deviance = final_deviance,
    info_matrix = info_matrix
  )
}


#' GLM + GEE Two-Stage Estimation (Case c)
#'
#' Stage 1 runs damped Newton-Raphson backfitting on the unwhitened link-scale working
#' response to produce a warm start.  Stage 2 refines via damped SQP
#' on the full whitened system, replicating \code{get_B} Path 1b.
#'
#' @inheritParams blockfit_solve
#' @param split Output of \code{.bf_split_components}.
#' @param schur_zero List of zeros (one per partition).
#'
#' @return A named list with \code{result}, \code{beta_spline},
#'   \code{beta_flat}, and \code{qp_info}.
#'
#' @keywords internal
.bf_case_glm_gee <- function(X, y, K, p_expansions, flat_cols,
                             split, family, order_list,
                             glm_weight_function,
                             schur_correction_function,
                             qp_score_function,
                             need_dispersion_for_estimation,
                             dispersion_function,
                             observation_weights,
                             iterate, tol, max_backfit_iter,
                             parallel_eigen, cl,
                             chunk_size, num_chunks, rem_chunks,
                             schur_zero,
                             quadprog, qp_Amat, qp_bvec, qp_meq,
                             Lambda, L_partition_list,
                             unique_penalty_per_partition,
                             A, constraint_values,
                             Vhalf, VhalfInv,
                             verbose, ...){

  spline_cols <- split$spline_cols
  nc_spline <- split$nc_spline
  nc_flat   <- split$nc_flat

  beta_flat   <- rep(0, nc_flat)
  beta_spline <- lapply(1:(K + 1), function(k) cbind(rep(0, nc_spline)))

  dnr_err_ws <- Inf
  damp_cnt_ws  <- 0
  disp_ws      <- 1

  ## Stage 1: damped Newton-Raphson backfitting warm start
  for(dnr_iter_ws in 1:50){

    eta_ws <- unlist(lapply(1:(K + 1), function(k){
      c(split$X_spline[[k]] %**% beta_spline[[k]] +
          split$X_flat[[k]] %**% cbind(beta_flat))
    }))
    mu_ws  <- family$linkinv(eta_ws)
    mu_eta_ws <- .bf_get_mu_eta(eta_ws, family)
    mu_eta_ws <- ifelse(abs(mu_eta_ws) < sqrt(.Machine$double.eps),
                        sqrt(.Machine$double.eps), mu_eta_ws)

    y_vec_ws <- unlist(y)
    z_ws <- eta_ws + (y_vec_ws - mu_ws) / mu_eta_ws

    if(need_dispersion_for_estimation){
      disp_ws <- dispersion_function(
        mu = mu_ws, y = y_vec_ws,
        order_indices = unlist(order_list),
        family = family,
        observation_weights = unlist(observation_weights),
        VhalfInv = NULL, ...)
    }

    W_ws <- c(glm_weight_function(mu_ws, y_vec_ws,
                                  unlist(order_list),
                                  family, disp_ws,
                                  unlist(observation_weights), ...))
    W_ws <- pmax(W_ws, sqrt(.Machine$double.eps))

    ## Partition z and W
    idx_ws <- 0L
    z_list_ws <- W_list_ws <- vector("list", K + 1)
    for(k in 1:(K + 1)){
      nk <- length(y[[k]])
      rows_ws <- (idx_ws + 1):(idx_ws + nk)
      z_list_ws[[k]] <- z_ws[rows_ws]
      W_list_ws[[k]] <- W_ws[rows_ws]
      idx_ws <- idx_ws + nk
    }

    ## Inner weighted backfitting
    inner <- .bf_inner_weighted(
      split$X_spline, split$X_flat,
      z_list_ws, W_list_ws, K,
      split$Lambda_spline, split$Lambda_flat,
      split$L_part_spline,
      unique_penalty_per_partition,
      split$A_spline, split$constraint_values_spline,
      nc_spline, nc_flat,
      beta_spline, beta_flat,
      tol, max_backfit_iter,
      parallel_eigen, cl,
      chunk_size, num_chunks, rem_chunks,
      split = split,
      constraint_values = constraint_values,
      p_expansions = p_expansions)

    beta_spline <- inner$beta_spline
    beta_flat   <- inner$beta_flat

    ## Damped Newton-Raphson convergence
    eta_new_ws <- unlist(lapply(1:(K + 1), function(k){
      c(split$X_spline[[k]] %**% beta_spline[[k]] +
          split$X_flat[[k]] %**% cbind(beta_flat))
    }))
    mu_new_ws <- family$linkinv(eta_new_ws)
    err_new_ws <- .bf_deviance(y_vec_ws, mu_new_ws,
                               unlist(observation_weights),
                               unlist(order_list), family, ...)

    if(verbose){
      cat('  GEE warm-start DNR iter', dnr_iter_ws,
          '| deviance:', format(err_new_ws, digits = 6), '\n')
    }

    if(is.na(err_new_ws) | !is.finite(err_new_ws)){
      damp_cnt_ws <- damp_cnt_ws + 1
      if(damp_cnt_ws >= 10) break
    } else if(err_new_ws <= dnr_err_ws){
      prev_dnr_ws_err <- dnr_err_ws
      dnr_err_ws <- err_new_ws
      damp_cnt_ws <- 0
      if(abs(prev_dnr_ws_err - dnr_err_ws) < tol &
         dnr_iter_ws > 5) break
    } else {
      damp_cnt_ws <- damp_cnt_ws + 1
      if(damp_cnt_ws >= 10) break
    }

    if(!iterate & dnr_iter_ws >= 2) break
  }

  ## Assemble warm start into full block vector
  result_ws <- lapply(1:(K + 1), function(k){
    b <- rep(0, p_expansions)
    b[spline_cols] <- beta_spline[[k]]
    b[flat_cols]   <- beta_flat
    cbind(b)
  })
  beta_block <- cbind(unlist(result_ws))

  ## Stage 2: Damped SQP on full whitened system
  if(verbose) cat("  GEE full-system refinement (Path 1b)\n")

  perm <- unlist(order_list)
  VhalfInv_perm <- VhalfInv[perm, perm]

  X_block   <- collapse_block_diagonal(X)
  y_block   <- cbind(unlist(y))
  X_tilde   <- VhalfInv_perm %**% X_block
  y_tilde   <- VhalfInv_perm %**% y_block
  Lambda_block <- .bf_make_Lambda_block(Lambda, K,
                                        unique_penalty_per_partition,
                                        L_partition_list)

  ## Combined constraint matrix
  if(quadprog && !is.null(qp_Amat) && ncol(cbind(qp_Amat)) > 0){
    qp_Amat_combined <- cbind(A, qp_Amat)
    qp_bvec_combined <- c(rep(0, ncol(A)), qp_bvec)
    qp_meq_combined  <- ncol(A) + qp_meq
  } else {
    qp_Amat_combined <- A
    qp_bvec_combined <- rep(0, ncol(A))
    qp_meq_combined  <- ncol(A)
  }

  if(length(constraint_values) > 0){
    cv_vec <- Reduce("rbind", constraint_values)
    constr_rhs <- crossprod(A, cv_vec)
    if(length(constr_rhs) == ncol(A)){
      qp_bvec_combined[1:ncol(A)] <- c(constr_rhs)
    }
  }

  ## [Diagnostic 2026-03-05] Verify qp_Amat_combined dimensions and
  #  equality count for GEE + blockfit path.
  #
  #  Expected:
  #    qp_Amat_combined has nrow = P = p_expansions * (K+1)
  #    qp_meq_combined  = ncol(A) + qp_meq
  #    ncol(qp_Amat_combined) = ncol(A) + ncol(qp_Amat)  [if qp present]
  #                           = ncol(A)                    [if no ineq]
  #
  #  If qp_meq_combined > ncol(qp_Amat_combined), solve.QP will fail
  #  silently or produce garbage Lagrangians.
  #
  #  If ncol(A) doesn't match the number of smoothness constraints
  #  (should be ~K * p_expansions - K for typical models), something
  #  upstream changed A unexpectedly.
  if(verbose){
    cat("  GEE SQP setup:",
        "nrow(qp_Amat_combined) =", nrow(qp_Amat_combined),
        "| ncol(qp_Amat_combined) =", ncol(qp_Amat_combined),
        "| qp_meq_combined =", qp_meq_combined,
        "| ncol(A) =", ncol(A),
        "| qp_meq =", ifelse(is.null(qp_meq), "NULL", qp_meq),
        "| P =", p_expansions * (K + 1), "\n")
    if(qp_meq_combined > ncol(qp_Amat_combined)){
      warning("\n\t [BUG] qp_meq_combined (", qp_meq_combined,
              ") > ncol(qp_Amat_combined) (", ncol(qp_Amat_combined),
              "). solve.QP will treat inequality columns as equalities.\n")
    }
  }

  sqp <- .bf_sqp_loop(
    X_design = X_tilde,
    y_design = y_tilde,
    X_block_raw = X_block,
    beta_init = beta_block,
    Lambda_block = Lambda_block,
    qp_Amat_combined = qp_Amat_combined,
    qp_bvec_combined = qp_bvec_combined,
    qp_meq_combined  = qp_meq_combined,
    K = K, p_expansions = p_expansions, family = family,
    order_list = order_list,
    glm_weight_function = glm_weight_function,
    schur_correction_function = schur_correction_function,
    qp_score_function = qp_score_function,
    need_dispersion_for_estimation = need_dispersion_for_estimation,
    dispersion_function = dispersion_function,
    observation_weights = observation_weights,
    iterate = iterate, tol = tol,
    VhalfInv = VhalfInv, VhalfInv_perm = VhalfInv_perm,
    is_gee = TRUE, deviance_fun = .bf_gee_deviance,
    X_partitions = X, y_partitions = y,
    verbose = verbose, ...)

  qp_info <- .bf_assemble_qp_info(
    sqp$last_qp_sol, sqp$beta_block,
    qp_Amat_combined, qp_bvec_combined, qp_meq_combined,
    sqp$converged, sqp$final_deviance)

  ## [Diagnostic 2026-03-05] Verify Amat_active from GEE path.
  #  Should have ncol = qp_meq_combined + (number of binding inequalities).
  #  The number of binding inequalities should be small relative to total.
  #  If ncol(Amat_active) == ncol(qp_Amat_combined), all constraints
  #  are "active" which likely indicates the old iact bug resurfacing.
  if(verbose && !is.null(qp_info)){
    n_ineq_active <- ncol(qp_info$Amat_active) - qp_meq_combined
    cat("  GEE qp_info:",
        "ncol(Amat_active) =", ncol(qp_info$Amat_active),
        "| n_equalities =", qp_meq_combined,
        "| n_active_ineq =", n_ineq_active,
        "| lagrangian_nonzero =",
        sum(abs(qp_info$lagrangian) > sqrt(.Machine$double.eps)),
        "\n")
    if(n_ineq_active > 0.5 * (ncol(qp_Amat_combined) - qp_meq_combined)){
      warning("\n\t [SUSPECT] More than half of inequality constraints ",
              "are active (", n_ineq_active, " of ",
              ncol(qp_Amat_combined) - qp_meq_combined,
              "). Check .bf_assemble_qp_info logic.\n")
    }
  }

  list(result = sqp$result,
       beta_spline = lapply(sqp$result, function(b) b[spline_cols, , drop = FALSE]),
       beta_flat   = sqp$result[[1]][flat_cols],
       qp_info     = qp_info)
}


#' GLM Without GEE: Damped Newton-Raphson + Backfitting (Case d)
#'
#' Damped Newton-Raphson outer loop wrapping the weighted backfitting inner loop.
#' Each damped Newton-Raphson step computes working response and weights from the
#' current linear predictor, then calls \code{.bf_inner_weighted}.
#'
#' @inheritParams blockfit_solve
#' @param split Output of \code{.bf_split_components}.
#' @param schur_zero List of zeros.
#'
#' @return A named list with \code{beta_spline} and
#'   \code{beta_flat}.
#'
#' @keywords internal
.bf_case_glm_no_corr <- function(X, y, K, p_expansions,
                                 split, family, order_list,
                                 glm_weight_function,
                                 need_dispersion_for_estimation,
                                 dispersion_function,
                                 observation_weights,
                                 iterate, tol, max_backfit_iter,
                                 parallel_eigen, cl,
                                 chunk_size, num_chunks,
                                 rem_chunks, schur_zero,
                                 unique_penalty_per_partition,
                                 VhalfInv,
                                 verbose,
                                 constraint_values = list(),
                                 ...){

  nc_spline <- split$nc_spline
  nc_flat   <- split$nc_flat

  beta_flat   <- rep(0, nc_flat)
  beta_spline <- lapply(1:(K + 1), function(k) cbind(rep(0, nc_spline)))

  dnr_err  <- Inf
  damp_cnt   <- 0
  disp_temp  <- 1

  for(dnr_iter in 1:100){

    eta <- unlist(lapply(1:(K + 1), function(k){
      c(split$X_spline[[k]] %**% beta_spline[[k]] +
          split$X_flat[[k]] %**% cbind(beta_flat))
    }))
    mu  <- family$linkinv(eta)

    mu_eta <- .bf_get_mu_eta(eta, family)
    mu_eta <- ifelse(abs(mu_eta) < sqrt(.Machine$double.eps),
                     sqrt(.Machine$double.eps), mu_eta)

    y_vec <- unlist(y)
    z <- eta + (y_vec - mu) / mu_eta

    if(need_dispersion_for_estimation){
      disp_temp <- dispersion_function(
        mu = mu, y = y_vec,
        order_indices = unlist(order_list),
        family = family,
        observation_weights = unlist(observation_weights),
        VhalfInv = VhalfInv, ...)
    }

    W <- c(glm_weight_function(mu, y_vec, unlist(order_list),
                               family, disp_temp,
                               unlist(observation_weights), ...))
    W <- pmax(W, sqrt(.Machine$double.eps))

    ## Partition z and W
    idx <- 0L
    z_list <- W_list <- vector("list", K + 1)
    for(k in 1:(K + 1)){
      nk <- length(y[[k]])
      rows <- (idx + 1):(idx + nk)
      z_list[[k]] <- z[rows]
      W_list[[k]] <- W[rows]
      idx <- idx + nk
    }

    ## Inner weighted backfitting
    inner <- .bf_inner_weighted(
      split$X_spline, split$X_flat,
      z_list, W_list, K,
      split$Lambda_spline, split$Lambda_flat,
      split$L_part_spline,
      unique_penalty_per_partition,
      split$A_spline, split$constraint_values_spline,
      nc_spline, nc_flat,
      beta_spline, beta_flat,
      tol, max_backfit_iter,
      parallel_eigen, cl,
      chunk_size, num_chunks, rem_chunks,
      split = split,
      constraint_values = constraint_values,
      p_expansions = p_expansions)

    beta_spline <- inner$beta_spline
    beta_flat   <- inner$beta_flat

    ## Damped Newton-Raphson convergence
    eta_new <- unlist(lapply(1:(K + 1), function(k){
      c(split$X_spline[[k]] %**% beta_spline[[k]] +
          split$X_flat[[k]] %**% cbind(beta_flat))
    }))
    mu_new <- family$linkinv(eta_new)

    err_new <- .bf_deviance(y_vec, mu_new, unlist(observation_weights),
                            unlist(order_list), family, ...)

    if(verbose){
      cat('  Backfitting DNR iter', dnr_iter,
          '| deviance:', format(err_new, digits = 6), '\n')
    }

    if(is.na(err_new) | !is.finite(err_new)){
      damp_cnt <- damp_cnt + 1
      if(damp_cnt >= 10) break
    } else if(err_new <= dnr_err){
      prev_dnr_err <- dnr_err
      dnr_err <- err_new
      damp_cnt <- 0
      if(abs(prev_dnr_err - dnr_err) < tol & dnr_iter > 5) break
    } else {
      damp_cnt <- damp_cnt + 1
      if(damp_cnt >= 10) break
    }

    if(!iterate & dnr_iter >= 2) break
  }

  list(beta_spline = beta_spline,
       beta_flat   = beta_flat)
}


#' Backfitting Solver for Blockfit Models
#'
#' @description
#' Fits models with mixed spline and non-interactive linear ("flat") terms
#' using an iterative backfitting approach. Flat terms receive a single
#' shared coefficient across partitions (rather than \eqn{K+1}
#' partition-specific coefficients constrained to equality), reducing the
#' effective parameter count and improving efficiency when the number of
#' flat terms is large relative to spline terms.
#'
#' The backfitting loop alternates between:
#' \enumerate{
#'
#' \item \strong{Spline step}
#'
#' Fit spline terms on the response adjusted for the current flat
#' contribution using a constrained penalized least squares solve.
#'
#' If \eqn{\mathbf{y}_k} is the response vector for partition \eqn{k},
#' \eqn{\mathbf{Z}_k} the spline design matrix,
#' \eqn{\mathbf{X}_{\mathrm{flat}}^{(k)}} the flat design matrix,
#' and \eqn{\mathbf{v}} the flat coefficient vector, then
#'
#' \deqn{
#' \boldsymbol{\beta}_{\mathrm{spline}}^{(k)} =
#' \arg\min_{\boldsymbol{\beta}}
#' \left\|
#' \mathbf{y}_k
#' - \mathbf{Z}_k \boldsymbol{\beta}
#' - \mathbf{X}_{\mathrm{flat}}^{(k)} \mathbf{v}
#' \right\|^2
#' +
#' \boldsymbol{\beta}^{\top}
#' \mathbf{\Lambda}_{\mathrm{spline}}
#' \boldsymbol{\beta}
#' }
#'
#' subject to
#'
#' \deqn{
#' \mathbf{A}_{\mathrm{spline}}^{\top}
#' \boldsymbol{\beta}
#' =
#' \mathbf{c}_{\mathrm{spline}}.
#' }
#'
#' \item \strong{Flat step}
#'
#' Update flat coefficients via pooled penalized regression on residuals,
#' subject to the residual equality constraint from any mixed constraints:
#'
#' \deqn{
#' \mathbf{v}
#' =
#' \arg\min_{\mathbf{v}}
#' \sum_{k=0}^{K}
#' \left\|
#' \mathbf{y}_k
#' -
#' \mathbf{Z}_k
#' \boldsymbol{\beta}_{\mathrm{spline}}^{(k)}
#' -
#' \mathbf{X}_{\mathrm{flat}}^{(k)}
#' \mathbf{v}
#' \right\|^2
#' +
#' \mathbf{v}^{\top}
#' \mathbf{\Lambda}_{\mathrm{flat}}
#' \mathbf{v}
#' }
#'
#' subject to
#'
#' \deqn{
#' \tilde{\mathbf{A}}_{\mathrm{flat}}^{\top}
#' \mathbf{v}
#' =
#' \mathbf{c}
#' -
#' \mathbf{A}_{\mathrm{spline}}^{\top}
#' \boldsymbol{\beta}_{\mathrm{spline}}.
#' }
#' }
#'
#' Four estimation paths are selected automatically based on the model
#' configuration:
#' \describe{
#'   \item{Case (a)}{Gaussian identity + GEE: closed-form solve via
#'     \code{.bf_case_gauss_gee}.}
#'   \item{Case (b)}{Gaussian identity, no correlation: standard
#'     backfitting via \code{.bf_case_gauss_no_corr}.}
#'   \item{Case (c)}{GLM + GEE: two-stage (damped Newton-Raphson warm start then
#'     damped SQP) via \code{.bf_case_glm_gee}.}
#'   \item{Case (d)}{GLM without GEE: damped Newton-Raphson + backfitting via
#'     \code{.bf_case_glm_no_corr}.}
#' }
#'
#' When \code{quadprog = TRUE} without GEE, inequality constraints are
#' enforced after backfitting convergence.  The constraint handling
#' method is selected automatically by inspecting the sparsity pattern
#' of \code{qp_Amat}: if every column has nonzeros in only one
#' partition block (e.g.\ derivative sign or range constraints),
#' a partition-wise active-set method is used via
#' \code{.active_set_refine}, avoiding the dense \eqn{P \times P}
#' system.  If any column spans multiple partition blocks (e.g.\
#' cross-knot monotonicity), or if the active-set method does not
#' converge, the dense SQP fallback via \code{.bf_sqp_loop} is used.
#' GEE paths always use the dense system.
#'
#' After convergence, coefficients are reassembled into the standard
#' per-partition format (flat coefficients replicated across partitions)
#' for compatibility with downstream inference.
#'
#' @param X List of \eqn{K+1} design matrices
#'   (\eqn{N_k \times p_expansions}). Unwhitened even when GEE.
#' @param y List of \eqn{K+1} response vectors. Unwhitened even when GEE.
#' @param flat_cols Integer vector indicating flat columns of
#'   \eqn{\mathbf{X}^{(k)}}.
#' @param K Integer; number of interior knots.
#' @param p_expansions Integer; number of coefficients per partition.
#' @param Lambda \eqn{p_expansions \times p_expansions} penalty matrix.
#' @param L_partition_list List of partition-specific penalty matrices.
#' @param unique_penalty_per_partition Logical.
#' @param A Full \eqn{P \times R_constraints} constraint matrix.
#' @param R_constraints Integer; number of columns of \eqn{\mathbf{A}}.
#' @param constraint_values List of constraint right-hand sides.
#' @param X_gram List of Gram matrices
#'   \eqn{\mathbf{X}_k^{\top} \mathbf{X}_k}.
#' @param Ghalf_full,GhalfInv_full Lists of
#'   \eqn{\mathbf{G}^{1/2}} and \eqn{\mathbf{G}^{-1/2}} matrices.
#' @param family GLM family object.
#' @param order_list List of observation index vectors by partition.
#' @param glm_weight_function GLM weight function.
#' @param schur_correction_function Schur complement correction function.
#' @param need_dispersion_for_estimation Logical.
#' @param dispersion_function Dispersion estimation function.
#' @param observation_weights List of observation weights.
#' @param homogenous_weights Logical.
#' @param iterate Logical; if FALSE, single pass (no iteration).
#' @param tol Convergence tolerance.
#' @param parallel_eigen,cl,chunk_size,num_chunks,rem_chunks Parallel arguments.
#' @param return_G_getB Logical.
#' @param quadprog Logical; apply inequality constraint refinement if TRUE.
#' @param qp_Amat Inequality constraint matrix for
#'   \code{quadprog::solve.QP}.
#' @param qp_bvec Inequality constraint right-hand side.
#' @param qp_meq Number of leading equality constraints.
#' @param qp_score_function Score function for QP subproblem.
#' @param keep_weighted_Lambda Logical.
#' @param max_backfit_iter Integer.
#' @param Vhalf Square root of the working correlation matrix in the
#'   original observation ordering.
#' @param VhalfInv Inverse square root of the working correlation matrix
#'   in the original observation ordering.  When both are non-\code{NULL},
#'   the GEE paths are used.
#' @param include_warnings Logical.
#' @param verbose Logical.
#' @param ... Additional arguments passed to weight and dispersion functions.
#'
#' @return A list with elements:
#' \describe{
#'   \item{B}{
#'     List of \eqn{K+1} coefficient vectors
#'     (\eqn{p_expansions \times 1}), flat coefficients replicated across partitions.
#'   }
#'   \item{G_list}{
#'     List containing \eqn{\mathbf{G}},
#'     \eqn{\mathbf{G}^{1/2}}, and
#'     \eqn{\mathbf{G}^{-1/2}},
#'     each a list of \eqn{K+1}
#'     \eqn{p_expansions \times p_expansions} matrices,
#'     or NULL if \code{return_G_getB} is FALSE.
#'
#'     \eqn{\mathbf{G}} satisfies
#'     \eqn{\mathbf{G}
#'     =
#'     \mathbf{G}^{1/2}
#'     (\mathbf{G}^{1/2})^{\top}} exactly.
#'
#'     When a correlation structure is present, these matrices are obtained
#'     from the dense whitened information matrix and then split back into
#'     per-partition blocks.  When flat columns are present, the returned
#'     matrices also reflect the pooled flat-column information across
#'     partitions.
#'   }
#'   \item{qp_info}{
#'     NULL when no inequality-refinement metadata are produced.
#'     Otherwise a list of available QP or active-set diagnostics.
#'     Dense SQP solves include \code{solution}, \code{lagrangian},
#'     \code{active_constraints}, \code{iact}, \code{Amat_active},
#'     \code{bvec_active}, \code{meq_active}, \code{converged}, and
#'     \code{final_deviance}; non-GEE dense solves also carry
#'     \code{info_matrix}, \code{Amat_combined}, \code{bvec_combined},
#'     and \code{meq_combined}.  Partition-wise active-set refinement
#'     returns the corresponding active-constraint summary only.
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' ## Minimal verification example: Gaussian identity, no correlation,
#' ## 2 partitions, 3 spline columns + 1 flat column per partition.
#' ##
#' ## This confirms that blockfit_solve returns compatible output and
#' ## that flat coefficients are replicated across partitions.
#'
#' set.seed(1234)
#' n1 <- 50; n2 <- 50
#' p_expansions <- 4  # intercept + x + x^2 + z (flat)
#' K  <- 1  # 2 partitions
#'
#' X1 <- cbind(1, rnorm(n1), rnorm(n1)^2, rnorm(n1))
#' X2 <- cbind(1, rnorm(n2), rnorm(n2)^2, rnorm(n2))
#' beta_true <- c(1, 0.5, -0.3, 0.8)
#' y1 <- X1 %*% beta_true + rnorm(n1, 0, 0.5)
#' y2 <- X2 %*% beta_true + rnorm(n2, 0, 0.5)
#'
#' ## Constraint: spline coefficients equal across partitions
#' ## (columns 1:3 constrained, column 4 is flat)
#' A <- matrix(0, nrow = p_expansions * (K + 1), ncol = 3)
#' for(j in 1:3){
#'   A[j, j]      <-  1
#'   A[p_expansions + j, j]  <- -1
#' }
#' qr_A <- qr(A)
#' A <- qr.Q(qr_A)[, 1:qr_A$rank, drop = FALSE]
#'
#' Lambda <- diag(p_expansions) * 1e-4
#' L_partition_list <- list(0, 0)
#' X_gram <- list(crossprod(X1), crossprod(X2))
#'
#' G_list <- compute_G_eigen(
#'   X_gram, Lambda, K,
#'   parallel = FALSE, cl = NULL,
#'   chunk_size = NULL, num_chunks = NULL, rem_chunks = NULL,
#'   family = gaussian(),
#'   unique_penalty_per_partition = FALSE,
#'   L_partition_list = L_partition_list,
#'   keep_G = TRUE,
#'   schur_corrections = list(0, 0))
#'
#' result <- blockfit_solve(
#'   X = list(X1, X2),
#'   y = list(y1, y2),
#'   flat_cols = 4L,
#'   K = K, p_expansions = p_expansions,
#'   Lambda = Lambda,
#'   L_partition_list = L_partition_list,
#'   unique_penalty_per_partition = FALSE,
#'   A = A, R_constraints = ncol(A),
#'   constraint_values = list(),
#'   X_gram = X_gram,
#'   Ghalf_full = G_list$Ghalf,
#'   GhalfInv_full = G_list$GhalfInv,
#'   family = gaussian(),
#'   order_list = list(1:n1, (n1+1):(n1+n2)),
#'   glm_weight_function = function(mu, y, oi, fam, d, ow, ...) rep(1, length(y)),
#'   schur_correction_function = function(X, y, B, d, ol, K, fam, ow, ...) {
#'     lapply(1:(K + 1), function(k) 0)
#'   },
#'   need_dispersion_for_estimation = FALSE,
#'   dispersion_function = function(...) 1,
#'   observation_weights = list(rep(1, n1), rep(1, n2)),
#'   homogenous_weights = TRUE,
#'   iterate = TRUE, tol = 1e-8,
#'   parallel = FALSE, cl = NULL,
#'   chunk_size = NULL, num_chunks = NULL, rem_chunks = NULL,
#'   return_G_getB = TRUE,
#'   verbose = FALSE)
#'
#' ## Verify: flat coefficient (col 4) is identical across partitions
#' stopifnot(abs(result$B[[1]][4] - result$B[[2]][4]) < 1e-10)
#'
#' ## Verify: output structure
#' stopifnot(length(result$B) == K + 1)
#' stopifnot(length(result$B[[1]]) == p_expansions)
#' stopifnot(!is.null(result$G_list))
#' }
#'
#' @seealso \code{\link{lgspline}}, \code{\link{lgspline.fit}},
#'   \code{\link{get_B}}
#'
#' @keywords internal
#' @export
blockfit_solve <- function(
    X,
    y,
    flat_cols,
    K,
    p_expansions,
    Lambda,
    L_partition_list,
    unique_penalty_per_partition,
    A,
    R_constraints,
    constraint_values,
    X_gram,
    Ghalf_full,
    GhalfInv_full,
    family,
    order_list,
    glm_weight_function,
    schur_correction_function,
    need_dispersion_for_estimation,
    dispersion_function,
    observation_weights,
    homogenous_weights = TRUE,
    iterate,
    tol,
    parallel_eigen,
    cl,
    chunk_size,
    num_chunks,
    rem_chunks,
    return_G_getB,
    quadprog = FALSE,
    qp_Amat = NULL,
    qp_bvec = NULL,
    qp_meq = NULL,
    qp_score_function = NULL,
    keep_weighted_Lambda = FALSE,
    max_backfit_iter = 100,
    Vhalf = NULL,
    VhalfInv = NULL,
    include_warnings = TRUE,
    verbose = FALSE,
    ...
){

  ## Dimensions and case detection
  nc_flat <- length(flat_cols)
  spline_cols <- setdiff(1:p_expansions, flat_cols)
  nc_spline <- length(spline_cols)
  nr <- sum(sapply(y, length))
  is_gauss_id <- (paste0(family)[1] == 'gaussian' &
                    paste0(family)[2] == 'identity')

  has_corr <- !is.null(Vhalf) & !is.null(VhalfInv)
  needs_gee_glm <- has_corr & !is_gauss_id

  ## Auto-detect whether inequality constraints require global (dense)
  ## QP or can be handled partition-wise. GEE always forces global.
  if(quadprog){
    qp_global <- has_corr | .solver_detect_qp_global(qp_Amat, p_expansions, K)
  } else {
    qp_global <- FALSE
  }
  if(has_corr) qp_global <- TRUE

  ## Split design matrices, penalty, and constraints
  #  The spline-only objects are then reused across the case-specific solvers.
  split <- .bf_split_components(X, flat_cols, p_expansions, K, Lambda,
                                L_partition_list, A,
                                constraint_values)

  ## Compute spline-only G^{1/2}
  X_gram_spline <- lapply(1:(K + 1), function(k){
    crossprod(split$X_spline[[k]])
  })
  schur_zero <- lapply(1:(K + 1), function(k) 0)

  G_sp <- compute_G_eigen(X_gram_spline,
                          split$Lambda_spline,
                          K,
                          parallel_eigen,
                          cl,
                          chunk_size,
                          num_chunks,
                          rem_chunks,
                          gaussian(),
                          unique_penalty_per_partition,
                          split$L_part_spline,
                          keep_G = TRUE,
                          schur_zero)
  Ghalf_sp <- G_sp$Ghalf

  ## G^{1/2} A for Lagrangian projection
  GhalfA_sp <- Reduce("rbind", lapply(1:(K + 1), function(k){
    rows <- ((k - 1) * nc_spline + 1):(k * nc_spline)
    Ghalf_sp[[k]] %**% split$A_spline[rows, , drop = FALSE]
  }))

  ## Pooled flat Gram and penalized inverse
  XfXf <- Reduce("+", lapply(1:(K + 1), function(k){
    crossprod(split$X_flat[[k]])
  }))
  XfXf_pen_inv <- invert(XfXf + split$Lambda_flat)

  ## Initialize
  beta_flat <- rep(0, nc_flat)
  beta_spline <- lapply(1:(K + 1), function(k) cbind(rep(0, nc_spline)))
  qp_info <- NULL

  ## Dispatch to the appropriate blockfit solver path.
  #  Cases are ordered from most specialized to most general.

  ## Case (a): Gaussian identity + GEE
  if(is_gauss_id & has_corr){

    out_a <- .bf_case_gauss_gee(
      X, y, K, p_expansions, order_list, VhalfInv,
      Lambda, L_partition_list,
      unique_penalty_per_partition,
      A, constraint_values,
      spline_cols, flat_cols,
      observation_weights)
    beta_spline <- out_a$beta_spline
    beta_flat   <- out_a$beta_flat
    result      <- out_a$result

    ## Case (b): Gaussian identity, no correlation
  } else if(is_gauss_id & !has_corr){

    out_b <- .bf_case_gauss_no_corr(
      split$X_spline, split$X_flat, y, K,
      Ghalf_sp, GhalfA_sp, XfXf_pen_inv,
      split$constraint_values_spline,
      nc_spline, nc_flat, split$A_spline,
      tol, max_backfit_iter, verbose,
      split = split,
      constraint_values = constraint_values,
      Lambda_flat = split$Lambda_flat,
      p_expansions = p_expansions)

    beta_spline <- out_b$beta_spline
    beta_flat   <- out_b$beta_flat

    ## Case (c): GLM + GEE
  } else if(has_corr){

    out_c <- .bf_case_glm_gee(
      X, y, K, p_expansions, flat_cols, split,
      family, order_list,
      glm_weight_function, schur_correction_function,
      qp_score_function,
      need_dispersion_for_estimation, dispersion_function,
      observation_weights,
      iterate, tol, max_backfit_iter,
      parallel_eigen, cl,
      chunk_size, num_chunks, rem_chunks, schur_zero,
      quadprog, qp_Amat, qp_bvec, qp_meq,
      Lambda, L_partition_list,
      unique_penalty_per_partition,
      A, constraint_values,
      Vhalf, VhalfInv, verbose, ...)

    result      <- out_c$result
    beta_spline <- out_c$beta_spline
    beta_flat   <- out_c$beta_flat
    qp_info     <- out_c$qp_info

    ## Case (d): GLM without GEE
  } else {

    out_d <- .bf_case_glm_no_corr(
      X, y, K, p_expansions, split, family, order_list,
      glm_weight_function,
      need_dispersion_for_estimation, dispersion_function,
      observation_weights,
      iterate, tol, max_backfit_iter,
      parallel_eigen, cl,
      chunk_size, num_chunks, rem_chunks, schur_zero,
      unique_penalty_per_partition, VhalfInv,
      verbose,
      constraint_values = constraint_values,
      ...)

    beta_spline <- out_d$beta_spline
    beta_flat   <- out_d$beta_flat
  }

  ## Assemble full per-partition coefficients (Cases b and d)
  if(!(is_gauss_id & has_corr) & !needs_gee_glm){
    result <- lapply(1:(K + 1), function(k){
      b <- rep(0, p_expansions)
      b[spline_cols] <- beta_spline[[k]]
      b[flat_cols]   <- beta_flat
      cbind(b)
    })
  }

  ## Inequality constraint refinement (non-GEE cases only)
  if(quadprog & !needs_gee_glm & !(is_gauss_id & has_corr)){

    if(!qp_global){
      ## Partition-wise active-set refinement.
      # Reuse the full-dimensional Ghalf/GhalfInv from above.
      # For Gaussian identity, Xy_or_uncon = Xy (cross-products);
      # for GLM, use result as unconstrained proxy.
      if(verbose) cat("  QP refinement (partition-wise active-set)\n")

      ## Compute Xy = X_k^T y_adj_k per partition for Gaussian,
      # or use result as unconstrained proxy for GLM.
      if(is_gauss_id){
        Xy_for_as <- lapply(1:(K + 1), function(k){
          crossprod(X[[k]], y[[k]])
        })
        is_p3 <- FALSE
      } else {
        Xy_for_as <- result
        is_p3 <- TRUE
      }

      as_out <- .active_set_refine(
        result, X, y, K, p_expansions,
        A, R_constraints, constraint_values,
        Lambda, Ghalf_full, GhalfInv_full, family,
        qp_Amat, qp_bvec, qp_meq,
        Xy_or_uncon = Xy_for_as, is_path3 = is_p3,
        parallel_aga = FALSE, parallel_matmult = FALSE,
        cl = cl, chunk_size = chunk_size,
        num_chunks = num_chunks, rem_chunks = rem_chunks,
        tol = tol)

      if(as_out$converged){
        result <- as_out$result
        qp_info <- as_out$qp_info
      } else {
        ## Fallback to dense SQP.
        if(verbose) cat("  Active-set did not converge, falling back to dense SQP\n")

        X_block <- collapse_block_diagonal(X)
        beta_block <- cbind(unlist(result))
        Lambda_block <- .bf_make_Lambda_block(
          Lambda, K, unique_penalty_per_partition, L_partition_list)
        y_block <- cbind(unlist(y))

        qp_Amat_combined_std <- cbind(A, qp_Amat)
        qp_bvec_combined_std <- c(rep(0, ncol(A)), qp_bvec)
        qp_meq_combined_std  <- ncol(A) + qp_meq

        sqp <- .bf_sqp_loop(
          X_design = X_block,
          y_design = y_block,
          X_block_raw = X_block,
          beta_init = beta_block,
          Lambda_block = Lambda_block,
          qp_Amat_combined = qp_Amat_combined_std,
          qp_bvec_combined = qp_bvec_combined_std,
          qp_meq_combined  = qp_meq_combined_std,
          K = K, p_expansions = p_expansions, family = family,
          order_list = order_list,
          glm_weight_function = glm_weight_function,
          schur_correction_function = schur_correction_function,
          qp_score_function = qp_score_function,
          need_dispersion_for_estimation = need_dispersion_for_estimation,
          dispersion_function = dispersion_function,
          observation_weights = observation_weights,
          iterate = iterate, tol = tol,
          VhalfInv = VhalfInv, VhalfInv_perm = NULL,
          is_gee = FALSE, deviance_fun = .bf_deviance,
          X_partitions = X, y_partitions = y,
          verbose = verbose, ...)

        result <- sqp$result

        qp_info <- .bf_assemble_qp_info(
          sqp$last_qp_sol, cbind(unlist(sqp$result)),
          qp_Amat_combined_std, qp_bvec_combined_std,
          qp_meq_combined_std,
          sqp$converged, sqp$final_deviance,
          info_matrix = NULL)
      }

    } else {
      ## Dense SQP (original behavior, cross-partition constraints detected).
      if(verbose) cat("  QP refinement (SQP loop)\n")

      X_block <- collapse_block_diagonal(X)
      beta_block <- cbind(unlist(result))
      Lambda_block <- .bf_make_Lambda_block(
        Lambda, K, unique_penalty_per_partition, L_partition_list)
      y_block <- cbind(unlist(y))

      qp_Amat_combined_std <- cbind(A, qp_Amat)
      qp_bvec_combined_std <- c(rep(0, ncol(A)), qp_bvec)
      qp_meq_combined_std  <- ncol(A) + qp_meq

      sqp <- .bf_sqp_loop(
        X_design = X_block,
        y_design = y_block,
        X_block_raw = X_block,
        beta_init = beta_block,
        Lambda_block = Lambda_block,
        qp_Amat_combined = qp_Amat_combined_std,
        qp_bvec_combined = qp_bvec_combined_std,
        qp_meq_combined  = qp_meq_combined_std,
        K = K, p_expansions = p_expansions, family = family,
        order_list = order_list,
        glm_weight_function = glm_weight_function,
        schur_correction_function = schur_correction_function,
        qp_score_function = qp_score_function,
        need_dispersion_for_estimation = need_dispersion_for_estimation,
        dispersion_function = dispersion_function,
        observation_weights = observation_weights,
        iterate = iterate, tol = tol,
        VhalfInv = VhalfInv, VhalfInv_perm = NULL,
        is_gee = FALSE, deviance_fun = .bf_deviance,
        X_partitions = X, y_partitions = y,
        verbose = verbose, ...)

      result <- sqp$result

      qp_info <- .bf_assemble_qp_info(
        sqp$last_qp_sol, cbind(unlist(sqp$result)),
        qp_Amat_combined_std, qp_bvec_combined_std,
        qp_meq_combined_std,
        sqp$converged, sqp$final_deviance,
        info_matrix = NULL)
    }
  }

  ## If downstream inference was not requested, stop at the
  #  coefficient estimates and any qp metadata gathered above.
  if(!return_G_getB){
    return(list(B = result, G_list = NULL, qp_info = qp_info))
  }

  ## Recompute G at final estimates using the same path as get_B Path 3.
  #  Since point estimates coincide with get_B, this produces identical
  #  G/Ghalf/GhalfInv and therefore identical U, varcovmat, and
  #  generate_posterior draws.
  if(is_gauss_id){
    G_list_out <- compute_G_eigen(
      X_gram, Lambda, K,
      parallel_eigen, cl,
      chunk_size, num_chunks, rem_chunks,
      family,
      unique_penalty_per_partition,
      L_partition_list,
      keep_G = TRUE,
      lapply(1:(K + 1), function(k) 0)
    )
    G_list_out$G <- lapply(G_list_out$Ghalf, tcrossprod)

    ## [Debug 2026-03-05] Flat-column pooling for blockfit.
    #  The flat columns were estimated with pooled information across
    #  all K+1 partitions (single shared coefficient), but compute_G_eigen
    #  treats each partition independently. The per-partition G blocks
    #  therefore reflect only partition-k's X_flat_k^T X_flat_k, not the
    #  pooled sum. Rather than patching G/Ghalf here (which breaks the
    #  identity G = tcrossprod(Ghalf) when cross-terms are zeroed), the
    #  fix is applied in lgspline.fit by augmenting the constraint matrix
    #  A with flat-equality constraints before U is constructed. This
    #  ensures U projects out the flat-equality directions, correctly
    #  pooling flat-column information in varcovmat and generate_posterior.
    #
    #  If posterior draws or confint still show inflated variance for
    #  blockfit vs non-blockfit, check:
    #    1. Is A augmented with flat-equality constraints in lgspline.fit
    #       BEFORE U is constructed? (Search for "flat_eq_cols")
    #    2. Is return_list$A updated after augmentation?
    #    3. Are G_list$Ghalf[[k]] identical between blockfit and non-blockfit?
    #       They should be, since both call compute_G_eigen on the same
    #       X_gram and Lambda.
    #    4. Is model_fit$Ghalf populated? (requires return_Ghalf = TRUE)
    #       generate_posterior uses model_fit$Ghalf to build L_post.
    #    5. Compare diag(model_fit$varcovmat) between blockfit and non-blockfit.
    #       Spline-column variances should be similar; flat-column variances
    #       will differ slightly because blockfit pools via U projection
    #       while non-blockfit pools via the original A constraint.
    #
    #  Diagnostic snippet (run after fitting):
    #    cat("G[[1]] diag:", diag(G_list_out$G[[1]]), "\n")
    #    cat("Ghalf[[1]] diag:", diag(G_list_out$Ghalf[[1]]), "\n")
    #    cat("flat_cols:", flat_cols, "\n")
    #    cat("spline_cols:", spline_cols, "\n")

  } else {
    G_list_out <- .solver_recompute_G_at_estimate(
      X = X, y = y, result = result, K = K, Lambda = Lambda,
      family = family, order_list = order_list,
      glm_weight_function = glm_weight_function,
      schur_correction_function = schur_correction_function,
      need_dispersion_for_estimation = need_dispersion_for_estimation,
      dispersion_function = dispersion_function,
      observation_weights = observation_weights,
      VhalfInv = VhalfInv,
      parallel_eigen = parallel_eigen, parallel_matmult = FALSE,
      cl = cl, chunk_size = chunk_size,
      num_chunks = num_chunks, rem_chunks = rem_chunks,
      unique_penalty_per_partition = unique_penalty_per_partition,
      L_partition_list = L_partition_list, ...
    )

    ## [Debug 2026-03-05] Non-Gaussian blockfit G recomputation.
    #  .solver_recompute_G_at_estimate computes per-partition weighted Gram
    #  matrices using GLM working weights at the final estimates. For
    #  flat columns, this again reflects only per-partition information.
    #  The A augmentation in lgspline.fit handles pooling via U, same
    #  as the Gaussian case.
    #
    #  If non-Gaussian blockfit shows inflated variance:
    #    1. Verify .solver_recompute_G_at_estimate returns proper matrices
    #       (not NULL for GhalfInv -- it should exist for non-Gaussian).
    #    2. Check that the A augmentation in lgspline.fit fires for
    #       non-Gaussian families (use_blockfit && length(flat_cols) > 0).
    #    3. Compare G_list_out$G[[k]] diagonal between blockfit and
    #       non-blockfit for the spline columns -- these should match
    #       if working weights are identical (same point estimates).
  }

  if(verbose){
    cat("  blockfit_solve returning: K =", K,
        "| p =", p_expansions,
        "| flat_cols =", paste(flat_cols, collapse = ","),
        "| n_result =", length(result),
        "| has_G =", !is.null(G_list_out$G[[1]]),
        "| has_Ghalf =", !is.null(G_list_out$Ghalf[[1]]),
        "| has_GhalfInv =", !is.null(G_list_out$GhalfInv),
        "| has_qp_info =", !is.null(qp_info), "\n")
  }

  return(list(B = result, G_list = G_list_out, qp_info = qp_info))
}
