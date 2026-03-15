## solver_utils.R
#  Shared helpers used by get_B() and blockfit_solve().
#
## Main pieces
#   .solver_build_lambda_block()
#   .solver_detect_qp_global()
#   .solver_recompute_G_at_estimate()
#   .solver_assemble_qp_info()




#' Build Block-Diagonal Penalty Matrix
#'
#' Assembles the full \eqn{P \times P} block-diagonal penalty matrix
#' \eqn{\boldsymbol{\Lambda}} from a shared per-partition penalty
#' \code{Lambda} and optional partition-specific additive terms.
#' \eqn{P = p \times (K+1)} where \eqn{p} is the number of basis terms
#' per partition.
#'
#' When \code{unique_penalty_per_partition = TRUE}, the \eqn{k}-th
#' diagonal block is \code{Lambda + L_partition_list[[k]]}; otherwise
#' every block is \code{Lambda}.
#'
#' @param Lambda Shared \eqn{p \times p} penalty matrix.
#' @param K Integer; number of interior knots (\eqn{K+1} partitions).
#' @param unique_penalty_per_partition Logical; if \code{TRUE}, add the
#'   corresponding element of \code{L_partition_list} to each block.
#' @param L_partition_list List of \eqn{K+1} partition-specific
#'   \eqn{p \times p} penalty matrices.
#'
#' @return A \eqn{P \times P} block-diagonal matrix,
#'   \eqn{P = p \times (K+1)}.
#'
#' @keywords internal
.solver_build_lambda_block <- function(Lambda,
                                       K,
                                       unique_penalty_per_partition,
                                       L_partition_list) {
  ## Lay down one penalty block per partition, with zeros off the diagonal.
  Reduce("rbind", lapply(seq_len(K + 1L), function(k) {
    Reduce("cbind", lapply(seq_len(K + 1L), function(j) {
      if (k != j) return(0 * Lambda)
      if (unique_penalty_per_partition) Lambda + L_partition_list[[k]]
      else Lambda
    }))
  }))
}




#' Detect Whether Inequality Constraints Require a Dense Global QP
#'
#' Inspects the columns of \code{qp_Amat} to determine whether every
#' inequality constraint is confined to a single partition block
#' (block-separable) or whether any constraint couples coefficients
#' across partitions.
#'
#' Returns \code{FALSE} (partition-wise active-set is valid) when all
#' columns have nonzeros in at most one block.  Returns \code{TRUE}
#' (dense SQP required) when any column spans multiple blocks.
#'
#' @param qp_Amat Inequality constraint matrix
#'   (\eqn{P \times n_{\mathrm{ineq}}}), or \code{NULL}.
#' @param p_expansions Integer; number of basis terms per partition.
#' @param K Integer; number of interior knots.
#'
#' @return Logical scalar.
#'
#' @keywords internal
.solver_detect_qp_global <- function(qp_Amat, p_expansions, K) {
  if (is.null(qp_Amat)) return(FALSE)
  if (!is.matrix(qp_Amat) || ncol(qp_Amat) == 0L) return(FALSE)

  tol <- sqrt(.Machine$double.eps)

  ## Scan each constraint column and see how many partition blocks it touches.
  for (j in seq_len(ncol(qp_Amat))) {
    nz <- which(abs(qp_Amat[, j]) > tol)
    if (length(nz) == 0L) next
    blocks_hit <- unique(ceiling(nz / p_expansions))
    if (length(blocks_hit) > 1L) return(TRUE)
  }
  FALSE
}




#' Recompute G, Ghalf, and GhalfInv at a Supplied Coefficient Estimate
#'
#' Given current constrained coefficient estimates \code{result},
#' recomputes the penalized information matrix
#' \eqn{\mathbf{G}_k = (\mathbf{X}_k^{\top}\mathbf{W}_k
#' \mathbf{D}_k\mathbf{X}_k + \boldsymbol{\Lambda}_k)^{-1}} and its
#' matrix square roots for each partition.
#'
#' This is used after Newton-Raphson convergence in Path 3 of
#' \code{get_B()} and at the final return in \code{blockfit_solve()}
#' when \code{return_G_getB = TRUE}.  The implementation matches
#' \code{.recompute_G_at_estimate()} exactly; it exists here so both
#' solvers can call the same numerical core.
#'
#' @param X List of partition-specific design matrices
#'   \eqn{\mathbf{X}_k}.
#' @param y List of response vectors \eqn{\mathbf{y}_k} by partition.
#' @param result List of current coefficient column vectors
#'   \eqn{\tilde{\boldsymbol{\beta}}_k} by partition.
#' @param K Integer; number of interior knots.
#' @param Lambda Shared \eqn{p \times p} penalty matrix.
#' @param family GLM family object.
#' @param order_list List of index vectors mapping partition rows to
#'   original data ordering.
#' @param glm_weight_function Function computing GLM working weights
#'   \eqn{\mathbf{W}_k}.
#' @param schur_correction_function Function computing Schur corrections
#'   to the information matrix.
#' @param need_dispersion_for_estimation Logical; if \code{TRUE},
#'   estimate dispersion before computing weights.
#' @param dispersion_function Dispersion estimation function.
#' @param observation_weights List of observation weights
#'   \eqn{\mathbf{D}_k} by partition.
#' @param VhalfInv Inverse square root correlation matrix, or
#'   \code{NULL}.
#' @param parallel_eigen,parallel_matmult Logical flags for parallel
#'   computation.
#' @param cl Parallel cluster object.
#' @param chunk_size,num_chunks,rem_chunks Parallel distribution
#'   parameters.
#' @param unique_penalty_per_partition Logical.
#' @param L_partition_list List of partition-specific penalty matrices.
#' @param ... Passed to weight, correction, and dispersion functions.
#'
#' @return A list with components \code{G}, \code{Ghalf}, and
#'   \code{GhalfInv}, each a list of \eqn{K+1} matrices.
#'   \code{G[[k]]} is computed as \code{tcrossprod(Ghalf[[k]])} to
#'   guarantee exact symmetry.
#'
#' @keywords internal
.solver_recompute_G_at_estimate <- function(X,
                                             y,
                                             result,
                                             K,
                                             Lambda,
                                             family,
                                             order_list,
                                             glm_weight_function,
                                             schur_correction_function,
                                             need_dispersion_for_estimation,
                                             dispersion_function,
                                             observation_weights,
                                             VhalfInv,
                                             parallel_eigen,
                                             parallel_matmult,
                                             cl,
                                             chunk_size,
                                             num_chunks,
                                             rem_chunks,
                                             unique_penalty_per_partition,
                                             L_partition_list,
                                             ...) {

  ## Some families need a current dispersion estimate before their
  #  working weights can be recomputed.
  if (need_dispersion_for_estimation) {
    mu <- family$linkinv(
      unlist(
        matmult_block_diagonal(X, result, K, parallel_matmult,
                               cl, chunk_size, num_chunks, rem_chunks)))
    dispersion_temp <- dispersion_function(
      mu                  = mu,
      y                   = unlist(y),
      order_indices       = unlist(order_list),
      family              = family,
      observation_weights = unlist(observation_weights),
      VhalfInv            = VhalfInv,
      ...)
  } else {
    dispersion_temp <- 1
  }

  ## Form the weighted design X_k * sqrt(W_k) for each partition.
  #  The Gram matrix X_k^T W_k D_k X_k follows from crossprod of this.
  Xw <- lapply(seq_len(K + 1L), function(k) {
    if (nrow(X[[k]]) == 0L) return(X[[k]])
    var <- glm_weight_function(
      family$linkinv(X[[k]] %**% cbind(c(result[[k]]))),
      y[[k]], order_list[[k]], family, dispersion_temp,
      observation_weights[[k]], ...)
    if (length(var) == 1L) {
      if (c(var) == 0) return(X[[k]] * 0)
      return(X[[k]] * c(sqrt(var)))
    }
    X[[k]] * c(sqrt(var))
  })

  X_gram <- compute_gram_block_diagonal(Xw, parallel_matmult, cl,
                                        chunk_size, num_chunks, rem_chunks)

  ## Rebuild any Schur corrections at the current coefficient estimate
  #  before handing the weighted Gram matrices back to compute_G_eigen.
  schur_corrections <- schur_correction_function(
    X, y, result, dispersion_temp, order_list, K, family,
    observation_weights, ...)

  G_list <- compute_G_eigen(X_gram, Lambda, K, parallel_eigen, cl,
                            chunk_size, num_chunks, rem_chunks, family,
                            unique_penalty_per_partition, L_partition_list,
                            keep_G = TRUE, schur_corrections)

  ## Recompute G from Ghalf to avoid rounding asymmetry that can arise
  ## when compute_G_eigen returns G and Ghalf through separate paths.
  G_list$G <- lapply(G_list$Ghalf, tcrossprod)

  G_list
}




#' Assemble qp_info from a solve.QP Solution
#'
#' Packages the output of \code{quadprog::solve.QP} into the
#' \code{qp_info} list expected by downstream code (inference,
#' \code{generate_posterior}, and \code{varcovmat} construction).
#'
#' Active constraint columns are identified as all equality columns
#' (\code{1:qp_meq_combined}) plus any inequality column whose
#' Lagrange multiplier exceeds \code{sqrt(.Machine$double.eps)}.
#' This matches the convention used in \code{.qp_refine()}.
#'
#' @param last_qp_sol Output of \code{quadprog::solve.QP}, or
#'   \code{NULL}.  When \code{NULL} the function returns \code{NULL}.
#' @param beta_block Final \eqn{P \times 1} coefficient vector.
#' @param qp_Amat_combined Combined equality + inequality constraint
#'   matrix.
#' @param qp_bvec_combined Combined constraint right-hand side.
#' @param qp_meq_combined Integer; number of leading equality
#'   constraints.
#' @param converged Logical; whether the outer loop converged.
#' @param final_deviance Scalar deviance at convergence.
#' @param info_matrix Optional information matrix; included in the
#'   returned list when non-\code{NULL}.
#'
#' @return A list with elements \code{solution}, \code{lagrangian},
#'   \code{active_constraints}, \code{iact}, \code{Amat_active},
#'   \code{bvec_active}, \code{meq_active}, \code{converged}, and
#'   \code{final_deviance}, plus \code{info_matrix},
#'   \code{Amat_combined}, \code{bvec_combined}, and
#'   \code{meq_combined} when \code{info_matrix} is supplied.
#'   Returns \code{NULL} if \code{last_qp_sol} is \code{NULL}.
#'
#' @keywords internal
.solver_assemble_qp_info <- function(last_qp_sol,
                                      beta_block,
                                      qp_Amat_combined,
                                      qp_bvec_combined,
                                      qp_meq_combined,
                                      converged,
                                      final_deviance,
                                      info_matrix = NULL) {
  if (is.null(last_qp_sol)) return(NULL)

  lag_mult <- last_qp_sol$Lagrangian
  n_eq     <- qp_meq_combined
  n_con    <- ncol(qp_Amat_combined)

  ## Active inequality columns: those with non-negligible multipliers.
  if (n_eq < n_con) {
    ineq_active <- n_eq + which(
      abs(lag_mult[-(seq_len(n_eq))]) > sqrt(.Machine$double.eps))
  } else {
    ineq_active <- integer(0L)
  }
  active_cols <- sort(unique(c(seq_len(n_eq), ineq_active)))

  ## Package the solve.QP output in the form used downstream by inference
  #  and posterior drawing code.
  out <- list(
    solution           = c(beta_block),
    lagrangian         = lag_mult,
    active_constraints = active_cols,
    iact               = last_qp_sol$iact,
    Amat_active        = qp_Amat_combined[, active_cols, drop = FALSE],
    bvec_active        = qp_bvec_combined[active_cols],
    meq_active         = qp_meq_combined,
    converged          = converged,
    final_deviance     = final_deviance
  )

  ## Keep the full combined constraint objects when the caller also wants
  #  the information matrix back.
  if (!is.null(info_matrix)) {
    out$info_matrix   <- info_matrix
    out$Amat_combined <- qp_Amat_combined
    out$bvec_combined <- qp_bvec_combined
    out$meq_combined  <- qp_meq_combined
  }

  out
}
