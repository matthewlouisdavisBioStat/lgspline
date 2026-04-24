## get_B.R
#  Constrained coefficient estimation for the main fitting paths.
#
## Main pieces
#   get_B()
#   the three fitting paths
#   Lagrangian projection helpers
#   QP refinement helpers
#   Woodbury helpers for the GEE paths


## Shared-helper wrappers

#' Build Block-Diagonal Penalty Matrix
#'
#' Thin wrapper around \code{.solver_build_lambda_block()}.
#' Keeping the local helper name preserves existing call sites in
#' \code{get_B()} while centralizing the shared implementation in
#' \code{solver_utils.R}.
#'
#' @param Lambda Shared \eqn{p \times p} penalty matrix.
#' @param K Integer; number of interior knots (\eqn{K+1} partitions).
#' @param unique_penalty_per_partition Logical; if \code{TRUE}, add
#'   \code{L_partition_list[[k]]} to the \eqn{k}-th diagonal block.
#' @param L_partition_list List of partition-specific \eqn{p \times p}
#'   penalty matrices.
#'
#' @return A \eqn{P \times P} block-diagonal matrix where
#'   \eqn{P = p \times (K+1)}.
#'
#' @keywords internal
.build_lambda_block <- function(Lambda, K, unique_penalty_per_partition,
                                L_partition_list) {
  .solver_build_lambda_block(
    Lambda = Lambda,
    K = K,
    unique_penalty_per_partition = unique_penalty_per_partition,
    L_partition_list = L_partition_list
  )
}



#' Extract Per-Partition Diagonal Blocks from a Full Matrix
#'
#' Given a \eqn{P \times P} matrix (e.g., the full-system
#' \eqn{\mathbf{G} = (\mathbf{X}^{\top}\mathbf{V}^{-1}\mathbf{X} +
#' \boldsymbol{\Lambda})^{-1}}), extracts the \eqn{p \times p} diagonal
#' block for each of the \eqn{K+1} partitions.
#'
#' @param M A \eqn{P \times P} matrix where \eqn{P = p\_expansions \times (K+1)}.
#' @param p_expansions Integer; number of basis terms (columns) per partition.
#' @param K Integer; number of interior knots.
#'
#' @return A list of length \eqn{K+1}, each element a
#'   \eqn{p\_expansions \times p\_expansions} matrix corresponding to the
#'   diagonal block for partition \eqn{k}.
#'
#' @keywords internal
.extract_G_diagonal <- function(M, p_expansions, K) {
  lapply(1:(K + 1), function(k) {
    M[(k - 1) * p_expansions + 1:p_expansions,
      (k - 1) * p_expansions + 1:p_expansions]
  })
}


## Lagrangian projection helpers

#' Lagrangian Projection via OLS Reformulation
#'
#' Given unconstrained (or GEE-initialized) coefficient estimates
#' \eqn{\hat{\boldsymbol{\beta}}}, computes the constrained estimate
#' \eqn{\tilde{\boldsymbol{\beta}} = \mathbf{U}\hat{\boldsymbol{\beta}}}
#' where \eqn{\mathbf{U} = \mathbf{I} - \mathbf{G}\mathbf{A}
#' (\mathbf{A}^{\top}\mathbf{G}\mathbf{A})^{-1}\mathbf{A}^{\top}}.
#'
#' Rather than forming \eqn{\mathbf{U}} directly, the projection is
#' reformulated as a residual from an OLS problem (the
#' \eqn{\mathbf{G}^{1/2}\mathbf{r}^*} trick):
#' \deqn{\mathbf{y}^* = \mathbf{G}^{-1/2}\hat{\boldsymbol{\beta}}}
#' \deqn{\mathbf{X}^* = \mathbf{G}^{1/2}\mathbf{A}}
#' \deqn{\mathbf{r}^* = (\mathbf{I} - \mathbf{X}^*(\mathbf{X}^{*\top}
#'   \mathbf{X}^*)^{-1}\mathbf{X}^{*\top})\mathbf{y}^*}
#' \deqn{\tilde{\boldsymbol{\beta}} = \mathbf{G}^{1/2}\mathbf{r}^*}
#'
#' This avoids explicitly forming and inverting
#' \eqn{\mathbf{A}^{\top}\mathbf{G}\mathbf{A}}; the most expensive step is
#' the QR decomposition of the \eqn{R \times R} system inside
#' \code{.lm.fit}, which is far cheaper than the full \eqn{P \times P}
#' solve.
#'
#' @param GhalfXy Numeric column vector \eqn{\mathbf{y}^*} of length
#'   \eqn{P}; produced by multiplying \eqn{\mathbf{G}^{1/2}} (Path 2) or
#'   \eqn{\mathbf{G}^{-1/2}} (Path 3) into the unconstrained estimate.
#' @param Ghalf List of \eqn{\mathbf{G}^{1/2}_k} matrices (one per
#'   partition).
#' @param A Constraint matrix \eqn{\mathbf{A}} (\eqn{P \times R}).
#' @param K Integer; number of interior knots.
#' @param p_expansions Integer; number of basis terms per partition.
#' @param R_constraints Integer; number of columns of \eqn{\mathbf{A}}.
#' @param constraint_value_vectors List of constraint right-hand-side
#'   vectors encoding \eqn{\mathbf{A}^{\top}\boldsymbol{\beta} =
#'   \mathbf{c}}. When non-empty and nonzero, the particular solution is
#'   added.
#' @param family GLM family object (used for \code{linkinv(0)} fallback
#'   when \code{NA} values arise in \eqn{\mathbf{y}^*}).
#' @param parallel_aga,parallel_matmult Logical flags for parallel
#'   computation.
#' @param cl Parallel cluster object.
#' @param chunk_size,num_chunks,rem_chunks Parallel distribution parameters.
#'
#' @return A list of length \eqn{K+1}, each element a column vector of
#'   constrained coefficients \eqn{\tilde{\boldsymbol{\beta}}_k}.
#'
#' @keywords internal
.lagrangian_project <- function(GhalfXy, Ghalf, A, K, p_expansions,
                                R_constraints, constraint_value_vectors,
                                family, parallel_aga, parallel_matmult,
                                cl, chunk_size, num_chunks, rem_chunks) {

  ## Compute X* = G^{1/2} A by stacking partition blocks vertically.
  #  GAmult_wrapper handles the block-diagonal multiplication in parallel
  #  when parallel_aga = TRUE.
  GhalfA <- Reduce("rbind",
                   GAmult_wrapper(Ghalf,
                                  A,
                                  K,
                                  p_expansions,
                                  R_constraints,
                                  parallel_aga,
                                  cl,
                                  chunk_size,
                                  num_chunks,
                                  rem_chunks))

  ## Replace any NA values in y* with the link-transformed zero.
  #  NAs can arise from empty partitions (n_k = 0).
  GhalfXy <- ifelse(is.na(GhalfXy), family$linkinv(0), GhalfXy)

  ## OLS residuals: project y* onto the orthogonal complement of col(X*).
  #  Scaling by 1/sqrt(K+1) improves numerical stability of .lm.fit
  #  when the constraint matrix has many rows.
  col_scales <- .constraint_col_scales(GhalfA)
  GhalfA_scaled <- t(t(GhalfA) * col_scales)
  comp_stab_sc <- 1 / sqrt(K + 1)
  resids_star <- do.call('.lm.fit', list(
    x = GhalfA_scaled * comp_stab_sc,
    y = GhalfXy * comp_stab_sc
  ))$residuals / comp_stab_sc

  ## If constraint values are nonzero (A^T beta = j with j != 0),
  #  add the particular solution (I - U) b_0 where A^T b_0 = j.
  #  The full constrained solution is: U beta_hat + (I - U) b_0.
  if (length(constraint_value_vectors) > 0) {
    if (any(unlist(constraint_value_vectors) != 0)) {
      comp_stab_sc <- 1 / sqrt(K + 1)
      preds_star <- GhalfA %**%
        (invert(crossprod(GhalfA) * comp_stab_sc) %**%
           crossprod(A,
                     Reduce("rbind", constraint_value_vectors) *
                       comp_stab_sc))
      resids_star <- resids_star + c(preds_star)
    }
  }

  ## Back-transform: beta_constrained_k = G^{1/2}_k r*_k (per partition).
  if (parallel_matmult & !is.null(cl)) {
    ## Handle remainder partitions that don't fill a full chunk.
    if (rem_chunks > 0) {
      rem_indices <- num_chunks * chunk_size + 1:rem_chunks
      rem <- lapply(rem_indices, function(k) {
        Ghalf[[k]] %**% cbind(resids_star[(k - 1) * p_expansions +
                                            1:p_expansions])
      })
    } else {
      rem <- list()
    }
    ## Parallel back-transformation across partition chunks.
    result <- c(
      Reduce("c", parallel::parLapply(cl, 1:num_chunks, function(chunk) {
        start_idx <- (chunk - 1) * chunk_size + 1
        end_idx <- chunk * chunk_size
        lapply(start_idx:end_idx, function(k) {
          Ghalf[[k]] %**% cbind(resids_star[(k - 1) * p_expansions +
                                              1:p_expansions])
        })
      })),
      rem
    )
  } else {
    result <- lapply(1:(K + 1), function(k) {
      Ghalf[[k]] %**% cbind(resids_star[(k - 1) * p_expansions +
                                          1:p_expansions])
    })
  }

  result
}


## Woodbury Lagrangian projection

#' Woodbury-Corrected Lagrangian Projection
#'
#' Applies the constrained projection for the structured-correlation
#' Woodbury path. In supplement notation,
#' \eqn{\mathbf{G}^{1/2} = \mathbf{G}_{\mathrm{on}}^{1/2}\mathbf{F}^{1/2}},
#' where \eqn{\mathbf{G}_{\mathrm{off}}^{-1} = \mathbf{E}\mathbf{J}\mathbf{E}^{\top}}
#' and \eqn{\mathbf{F}^{1/2}} carries the low-rank correction.
#'
#' The routine follows the same four transformed-OLS steps as
#' \code{.lagrangian_project}: form \eqn{\mathbf{y}^{*}},
#' form \eqn{\mathbf{X}^{*}}, solve the OLS projection, then back-transform.
#' Only steps involving \eqn{\mathbf{G}^{1/2}} differ, and those are handled
#' as a block-diagonal multiplication through
#' \eqn{\mathbf{G}_{\mathrm{on}}^{1/2}} plus a rank-\eqn{r} update.
#'
#' Internally, \eqn{\mathbf{F}^{1/2}} is stored as
#' \eqn{\mathbf{I} - \mathbf{Q}\mathbf{C}\mathbf{Q}^{\top}}, where
#' \code{Q} is the basis corresponding to the supplement matrix
#' \eqn{\mathbf{Q}} and \code{C} is diagonal.
#'
#' @param GhalfXy_V Numeric column vector of length \eqn{P};
#'   the fully transformed right-hand side
#'   \eqn{\mathbf{G}^{1/2}\mathbf{X}^{\top}\mathbf{V}^{-1}\mathbf{y}}
#'   assembled through the Woodbury-corrected operator.
#' @param Ghalf_corrected List of \eqn{K+1} corrected
#'   \eqn{\mathbf{G}_{\mathrm{on},k}^{1/2}} matrices.
#' @param A Constraint matrix (\eqn{P \times R}).
#' @param K,p_expansions,R_constraints Integer dimensions.
#' @param constraint_value_vectors Constraint RHS list encoding
#'   \eqn{\mathbf{A}^{\top}\boldsymbol{\beta} = \mathbf{c}}.
#' @param family GLM family object.
#' @param wb_sqrt Output of \code{.woodbury_halfsqrt_components},
#'   containing the internal rank-\eqn{r} representation of
#'   \eqn{\mathbf{F}^{1/2}}: \code{Q}
#'   (basis matching supplement \eqn{\mathbf{Q}}), \code{C} and \code{C_inv} (diagonal correction
#'   coefficients), and \code{GhalfQ} (partition-wise copies of
#'   \eqn{\mathbf{G}_{\mathrm{on},k}^{1/2}\mathbf{Q}[rows_k,]}).
#' @param parallel_aga,parallel_matmult Logical flags.
#' @param cl,chunk_size,num_chunks,rem_chunks Parallel parameters.
#'
#' @return A list of \eqn{K+1} coefficient column vectors.
#'
#' @keywords internal
.lagrangian_project_woodbury <- function(GhalfXy_V, Ghalf_corrected,
                                         A, K, p_expansions,
                                         R_constraints,
                                         constraint_value_vectors,
                                         family,
                                         wb_sqrt,
                                         parallel_aga, parallel_matmult,
                                         cl, chunk_size, num_chunks,
                                         rem_chunks) {

  Q <- wb_sqrt$Q
  C <- wb_sqrt$C
  GhalfQ <- wb_sqrt$GhalfQ

  ## Step 1: transformed right-hand side
  #  The caller already formed y* = G^{1/2} X^T V^{-1} y through the
  #  same Woodbury-corrected operator used below for X*.
  y_star <- GhalfXy_V

  ## Replace NAs with link-transformed zero (same as .lagrangian_project).
  y_star <- ifelse(is.na(y_star), family$linkinv(0), y_star)

  ## Step 2: transformed constraint matrix
  #  Form X* = (G^{1/2})^T A with
  #  G^{1/2} = G_on^{1/2} F^{1/2}. For non-commuting factors,
  #  the transpose orientation is required.

  X_star <- .woodbury_transform_constraint_matrix_transpose(
    Ghalf_corrected = Ghalf_corrected,
    A = A,
    K = K,
    p_expansions = p_expansions,
    wb_sqrt = wb_sqrt,
    parallel_aga = parallel_aga,
    cl = cl,
    chunk_size = chunk_size,
    num_chunks = num_chunks,
    rem_chunks = rem_chunks
  )

  ## Step 3: transformed OLS residuals
  #  Same residual calculation as .lagrangian_project().
  col_scales <- .constraint_col_scales(X_star)
  X_star_scaled <- t(t(X_star) * col_scales)
  comp_stab_sc <- 1 / sqrt(K + 1)
  resids_star <- do.call('.lm.fit', list(
    x = X_star_scaled * comp_stab_sc,
    y = y_star * comp_stab_sc
  ))$residuals / comp_stab_sc

  ## Nonzero constraint values
  #  Same adjustment as .lagrangian_project() when A^T beta = j with j != 0.
  if (length(constraint_value_vectors) > 0) {
    if (any(unlist(constraint_value_vectors) != 0)) {
      comp_stab_sc <- 1 / sqrt(K + 1)
      preds_star <- X_star %**%
        (invert(crossprod(X_star) * comp_stab_sc) %**%
           crossprod(A,
                     Reduce("rbind", constraint_value_vectors) *
                       comp_stab_sc))
      resids_star <- resids_star + c(preds_star)
    }
  }

  ## Step 4: back-transform to beta_tilde
  #  Apply the same block-plus-rank-r operator used in Steps 1 and 2.
  v_r <- crossprod(Q, cbind(resids_star))  # r x 1
  C_v_r <- C %**% v_r                        # r x 1

  result <- lapply(1:(K + 1), function(k) {
    rows_k <- ((k - 1) * p_expansions + 1):(k * p_expansions)
    ## Block-diagonal part.
    block_part <- Ghalf_corrected[[k]] %**%
      cbind(resids_star[rows_k])
    ## Rank-r correction.
    correction_part <- GhalfQ[[k]] %**% C_v_r
    block_part - correction_part
  })

  result
}


## Full-system Lagrangian projection

#' Full-System Lagrangian Projection
#'
#' Applies the transformed OLS projection when \eqn{\mathbf{G}^{1/2}} is a
#' full \eqn{P \times P} matrix.  This is used as an exact active-set
#' refinement for correlated paths after the Woodbury solve has found the
#' current quadratic approximation.
#'
#' @keywords internal
.lagrangian_project_full <- function(GhalfXy, Ghalf_full, A, K,
                                     p_expansions, constraint_value_vectors,
                                     family) {

  y_star <- ifelse(is.na(GhalfXy), family$linkinv(0), GhalfXy)
  X_star <- Ghalf_full %**% A

  col_scales <- .constraint_col_scales(X_star)
  X_star_scaled <- t(t(X_star) * col_scales)
  comp_stab_sc <- 1 / sqrt(K + 1)
  resids_star <- do.call(".lm.fit", list(
    x = X_star_scaled * comp_stab_sc,
    y = y_star * comp_stab_sc
  ))$residuals / comp_stab_sc

  if (length(constraint_value_vectors) > 0) {
    if (any(unlist(constraint_value_vectors) != 0)) {
      cv_vec <- Reduce("rbind", constraint_value_vectors)
      preds_star <- X_star %**%
        (invert(crossprod(X_star) * comp_stab_sc) %**%
           crossprod(A, cv_vec * comp_stab_sc))
      resids_star <- resids_star + c(preds_star)
    }
  }

  beta_full <- Ghalf_full %**% cbind(resids_star)

  lapply(1:(K + 1), function(k) {
    rows_k <- ((k - 1) * p_expansions + 1):(k * p_expansions)
    beta_full[rows_k, , drop = FALSE]
  })
}


## Full-system KKT check

#' Check KKT Conditions for Full-System Active-Set Refinement
#'
#' @keywords internal
.check_kkt_full <- function(result, GhalfXy, Ghalf_full, A_aug,
                            n_eq_orig, qp_Amat, qp_bvec, active_ineq,
                            K, tol) {

  beta_full <- cbind(unlist(result))
  primal <- .solver_check_primal_feasibility(
    beta_full = beta_full,
    qp_Amat = qp_Amat,
    qp_bvec = qp_bvec,
    active_ineq = active_ineq,
    tol = tol
  )

  drop <- integer(0)
  multipliers <- numeric(0)

  if (length(active_ineq) > 0) {
    mult_out <- .solver_recover_active_multipliers(
      X_star = Ghalf_full %**% A_aug,
      y_star = GhalfXy,
      n_eq_orig = n_eq_orig,
      K = K,
      tol = tol,
      ineq_sign = -1
    )
    drop <- mult_out$drop
    multipliers <- mult_out$multipliers
  }

  list(feasible = primal$feasible,
       dual_feasible = length(drop) == 0,
       violated = primal$violated,
       drop = drop,
       multipliers = multipliers)
}


## Full-system active set

#' Full-System Active-Set Refinement for Correlated Paths
#'
#' Reuses the generic add/drop active-set controller but evaluates each
#' equality subproblem with the exact full \eqn{\mathbf{G}^{1/2}}.  This keeps
#' dense QP as a fallback while avoiding multiplier sign errors from using
#' only diagonal blocks of a full correlated square root.
#'
#' @keywords internal
.active_set_refine_full <- function(result, K, p_expansions,
                                    A, constraint_value_vectors,
                                    family,
                                    qp_Amat, qp_bvec, qp_meq,
                                    rhs_full, Ghalf_full,
                                    tol, max_as_iter = NULL,
                                    method = "active_set_full") {

  n_ineq <- if (is.null(qp_Amat)) 0L else ncol(qp_Amat)
  if (is.null(max_as_iter)) {
    max_as_iter <- max(200, min(1000, 10 * max(1, n_ineq)))
  }

  GhalfXy <- Ghalf_full %**% cbind(rhs_full)

  .solver_active_set(
    result = result,
    A = A,
    constraint_value_vectors = constraint_value_vectors,
    qp_Amat = qp_Amat,
    qp_bvec = qp_bvec,
    qp_meq = qp_meq,
    tol = tol,
    max_as_iter = max_as_iter,
    solve_subproblem = function(aug_state) {
      .lagrangian_project_full(
        GhalfXy = GhalfXy,
        Ghalf_full = Ghalf_full,
        A = aug_state$A_aug,
        K = K,
        p_expansions = p_expansions,
        constraint_value_vectors = aug_state$cv_for_proj,
        family = family
      )
    },
    kkt_subproblem = function(result_new, aug_state) {
      .check_kkt_full(
        result = result_new,
        GhalfXy = GhalfXy,
        Ghalf_full = Ghalf_full,
        A_aug = aug_state$A_aug,
        n_eq_orig = aug_state$n_eq_aug,
        qp_Amat = aug_state$qp_ineq_Amat,
        qp_bvec = aug_state$qp_ineq_bvec,
        active_ineq = match(aug_state$active_ineq_kept,
                            aug_state$ineq_cols),
        K = K,
        tol = tol
      )
    },
    method = method,
    debug_label = "full-as"
  )
}


#' Form the Woodbury-Corrected Constraint Design Matrix
#'
#' Constructs the transformed constraint matrix
#' \eqn{\mathbf{X}^{*} = \mathbf{G}^{1/2}\mathbf{A}} for the
#' Woodbury-accelerated correlated solver without forming the full dense
#' \eqn{\mathbf{G}^{1/2}} explicitly.
#'
#' @param Ghalf_corrected List of corrected
#'   \eqn{\mathbf{G}_{\mathrm{on},k}^{1/2}} matrices.
#' @param A Constraint matrix.
#' @param K,p_expansions Integer dimensions.
#' @param wb_sqrt Output of \code{.woodbury_halfsqrt_components()}.
#' @param parallel_aga Logical.
#' @param cl,chunk_size,num_chunks,rem_chunks Parallel parameters.
#'
#' @return Numeric matrix with \eqn{P} rows and \code{ncol(A)} columns.
#'
#' @keywords internal
.woodbury_transform_constraint_matrix <- function(Ghalf_corrected, A,
                                                  K, p_expansions,
                                                  wb_sqrt,
                                                  parallel_aga,
                                                  cl, chunk_size,
                                                  num_chunks,
                                                  rem_chunks) {
  GhalfA <- Reduce("rbind",
                   GAmult_wrapper(Ghalf_corrected, A, K, p_expansions,
                                  ncol(A), parallel_aga, cl,
                                  chunk_size, num_chunks, rem_chunks))

  QtA <- crossprod(wb_sqrt$Q, A)
  C_QtA <- wb_sqrt$C %**% QtA
  GhalfQ_full <- do.call(rbind, wb_sqrt$GhalfQ)
  correction_A <- GhalfQ_full %**% C_QtA

  GhalfA - correction_A
}


#' Form the Transpose Woodbury-Corrected Constraint Matrix
#'
#' Applies \eqn{(\mathbf{G}^{1/2})^{\top}\mathbf{A}} with
#' \eqn{\mathbf{G}^{1/2} = \mathbf{G}_{\mathrm{on}}^{1/2}\mathbf{F}^{1/2}}.
#' The transformed OLS system uses this orientation for \eqn{y^*} and
#' \eqn{X^*}, then maps residuals back with \eqn{\mathbf{G}^{1/2}}.
#'
#' @keywords internal
.woodbury_transform_constraint_matrix_transpose <- function(
    Ghalf_corrected, A, K, p_expansions, wb_sqrt,
    parallel_aga, cl, chunk_size, num_chunks, rem_chunks) {

  GhalfA <- Reduce("rbind",
                   GAmult_wrapper(Ghalf_corrected, A, K, p_expansions,
                                  ncol(A), parallel_aga, cl,
                                  chunk_size, num_chunks, rem_chunks))

  QtGhalfA <- crossprod(wb_sqrt$Q, GhalfA)
  C_QtGhalfA <- wb_sqrt$C %**% QtGhalfA

  GhalfA - wb_sqrt$Q %**% C_QtGhalfA
}


#' Compute Stable Column Rescaling for Constraint Matrices
#'
#' Returns multiplicative column scales that normalize nonzero columns to
#' unit Euclidean norm. This preserves column spaces and projections while
#' improving numerical conditioning for QR / OLS steps.
#'
#' @param X Numeric matrix.
#' @param eps Small positive threshold below which a column is treated as
#'   numerically zero.
#'
#' @return Numeric vector of length \code{ncol(X)} containing multiplicative
#'   scales.
#'
#' @keywords internal
.constraint_col_scales <- function(X,
                                   eps = sqrt(.Machine$double.eps)) {
  if (is.null(dim(X)) || ncol(X) == 0) return(numeric(0))
  col_norms <- sqrt(colSums(X^2))
  ifelse(is.finite(col_norms) & col_norms > eps, 1 / col_norms, 1)
}


#' Check KKT Conditions for Woodbury Partition-Wise Active-Set Method
#'
#' Woodbury analogue of \code{.check_kkt_partitionwise()}.  The primal
#' feasibility check is unchanged, while the dual multipliers are recovered
#' from the transformed OLS system used by
#' \code{.lagrangian_project_woodbury()}.
#'
#' @param result List of coefficient column vectors.
#' @param GhalfXy_V Numeric column vector already in the final
#'   \eqn{\mathbf{G}^{1/2}} coordinates used by the projection.
#' @param Ghalf_corrected List of corrected
#'   \eqn{\mathbf{G}_{\mathrm{on},k}^{1/2}} matrices.
#' @param A_aug Augmented constraint matrix.
#' @param n_eq_orig Integer; number of always-active equality constraints.
#' @param qp_Amat,qp_bvec Full inequality specification.
#' @param active_ineq Integer vector of active inequality indices.
#' @param K,p_expansions Integer dimensions.
#' @param family GLM family object.
#' @param wb_sqrt Output of \code{.woodbury_halfsqrt_components()}.
#' @param parallel_aga Logical.
#' @param cl,chunk_size,num_chunks,rem_chunks Parallel parameters.
#' @param tol Numeric feasibility tolerance.
#'
#' @return Same structure as \code{.check_kkt_partitionwise()}.
#'
#' @keywords internal
.check_kkt_partitionwise_woodbury <- function(result, GhalfXy_V,
                                              Ghalf_corrected,
                                              A_aug, n_eq_orig,
                                              qp_Amat, qp_bvec,
                                              active_ineq,
                                              K, p_expansions,
                                              family,
                                              wb_sqrt,
                                              parallel_aga,
                                              cl, chunk_size,
                                              num_chunks, rem_chunks,
                                              tol) {

  beta_full <- cbind(unlist(result))
  primal <- .solver_check_primal_feasibility(
    beta_full = beta_full,
    qp_Amat = qp_Amat,
    qp_bvec = qp_bvec,
    active_ineq = active_ineq,
    tol = tol
  )

  drop <- integer(0)
  multipliers <- numeric(0)

  if (length(active_ineq) > 0) {
    y_star <- GhalfXy_V
    y_star <- ifelse(is.na(y_star), family$linkinv(0), y_star)

    X_star <- .woodbury_transform_constraint_matrix_transpose(
      Ghalf_corrected = Ghalf_corrected,
      A = A_aug,
      K = K,
      p_expansions = p_expansions,
      wb_sqrt = wb_sqrt,
      parallel_aga = parallel_aga,
      cl = cl,
      chunk_size = chunk_size,
      num_chunks = num_chunks,
      rem_chunks = rem_chunks
    )

    mult_out <- .solver_recover_active_multipliers(
      X_star = X_star,
      y_star = y_star,
      n_eq_orig = n_eq_orig,
      K = K,
      tol = tol,
      ineq_sign = -1
    )
    drop <- mult_out$drop
    multipliers <- mult_out$multipliers
  }

  list(feasible = primal$feasible,
       dual_feasible = length(drop) == 0,
       violated = primal$violated,
       drop = drop,
       multipliers = multipliers)
}


#' Partition-Wise Active-Set Refinement for Woodbury GEE Paths
#'
#' Uses the same add/drop active-set logic as \code{.active_set_refine()},
#' but solves each equality-constrained subproblem through the
#' Woodbury-corrected Lagrangian projection rather than the plain
#' block-diagonal projection.
#'
#' @param result List of current coefficient column vectors.
#' @param K,p_expansions Integer dimensions.
#' @param A Original equality constraint matrix.
#' @param R_constraints Number of columns of \code{A}.
#' @param constraint_value_vectors Constraint RHS list.
#' @param family GLM family object.
#' @param qp_Amat,qp_bvec,qp_meq Inequality specification.
#' @param rhs_list List of partition-wise right-hand-side vectors on the
#'   original coefficient scale.  Gaussian GEE passes
#'   \eqn{\mathbf{X}_k^{\top}\mathbf{V}^{-1}\mathbf{y}}, while the
#'   non-Gaussian Woodbury path passes the partitioned score.
#' @param Ghalf_corrected List of corrected
#'   \eqn{\mathbf{G}_{\mathrm{on},k}^{1/2}} matrices.
#' @param wb_sqrt Output of \code{.woodbury_halfsqrt_components()}.
#' @param parallel_aga,parallel_matmult Logical flags.
#' @param cl,chunk_size,num_chunks,rem_chunks Parallel parameters.
#' @param tol Convergence tolerance.
#' @param max_as_iter Maximum active-set iterations (default 50).
#'
#' @return Same structure as \code{.active_set_refine()}.
#'
#' @keywords internal
.active_set_refine_woodbury <- function(result,
                                        K, p_expansions,
                                        A, R_constraints,
                                        constraint_value_vectors,
                                        family,
                                        qp_Amat, qp_bvec, qp_meq,
                                        rhs_list,
                                        Ghalf_corrected,
                                        wb_sqrt,
                                        parallel_aga, parallel_matmult,
                                        cl, chunk_size, num_chunks,
                                        rem_chunks, tol,
                                        max_as_iter = NULL) {
  n_ineq <- if (is.null(qp_Amat)) 0L else ncol(qp_Amat)
  rhs_full <- cbind(unlist(rhs_list))
  GhalfXy_V <- .woodbury_transform_constraint_matrix_transpose(
    Ghalf_corrected = Ghalf_corrected,
    A = rhs_full,
    K = K,
    p_expansions = p_expansions,
    wb_sqrt = wb_sqrt,
    parallel_aga = parallel_aga,
    cl = cl,
    chunk_size = chunk_size,
    num_chunks = num_chunks,
    rem_chunks = rem_chunks
  )

  ## Use the same larger budget as the full-system active-set path.
  if (is.null(max_as_iter)) {
    max_as_iter <- max(200, min(1000, 10 * max(1, n_ineq)))
  }

  .solver_active_set(
    result = result,
    A = A,
    constraint_value_vectors = constraint_value_vectors,
    qp_Amat = qp_Amat,
    qp_bvec = qp_bvec,
    qp_meq = qp_meq,
    tol = tol,
    max_as_iter = max_as_iter,
    solve_subproblem = function(aug_state) {
      .lagrangian_project_woodbury(
        GhalfXy_V = GhalfXy_V,
        Ghalf_corrected = Ghalf_corrected,
        A = aug_state$A_aug,
        K = K,
        p_expansions = p_expansions,
        R_constraints = ncol(aug_state$A_aug),
        constraint_value_vectors = aug_state$cv_for_proj,
        family = family,
        wb_sqrt = wb_sqrt,
        parallel_aga = parallel_aga,
        parallel_matmult = parallel_matmult,
        cl = cl,
        chunk_size = chunk_size,
        num_chunks = num_chunks,
        rem_chunks = rem_chunks
      )
    },
    kkt_subproblem = function(result_new, aug_state) {
      .check_kkt_partitionwise_woodbury(
        result = result_new,
        GhalfXy_V = GhalfXy_V,
        Ghalf_corrected = Ghalf_corrected,
        A_aug = aug_state$A_aug,
        n_eq_orig = aug_state$n_eq_aug,
        qp_Amat = aug_state$qp_ineq_Amat,
        qp_bvec = aug_state$qp_ineq_bvec,
        active_ineq = match(aug_state$active_ineq_kept, aug_state$ineq_cols),
        K = K,
        p_expansions = p_expansions,
        family = family,
        wb_sqrt = wb_sqrt,
        parallel_aga = parallel_aga,
        cl = cl,
        chunk_size = chunk_size,
        num_chunks = num_chunks,
        rem_chunks = rem_chunks,
        tol = tol
      )
    },
    method = "active_set_woodbury",
    debug_label = "woodbury-as"
  )
}


#' Attempt Woodbury Partition-Wise Active-Set Refinement
#'
#' Small wrapper around \code{.active_set_refine_woodbury()} used by both
#' \code{get_B()} and \code{blockfit_solve()}. When the refinement errors,
#' the caller can cleanly fall back to the dense QP path while optionally
#' surfacing the underlying failure.
#'
#' @param result Current coefficient list.
#' @param K,p_expansions Integer dimensions.
#' @param A Original equality constraint matrix.
#' @param R_constraints Number of columns of \code{A}.
#' @param constraint_value_vectors Constraint RHS list.
#' @param family GLM family object.
#' @param qp_Amat,qp_bvec,qp_meq Inequality specification.
#' @param rhs_list Partition-wise right-hand-side vectors for the current
#'   Woodbury-refined quadratic approximation.
#' @param Ghalf_corrected List of corrected
#'   \eqn{\mathbf{G}_{\mathrm{on},k}^{1/2}} matrices.
#' @param wb_sqrt Output of \code{.woodbury_halfsqrt_components()}.
#' @param parallel_aga,parallel_matmult Logical flags.
#' @param cl,chunk_size,num_chunks,rem_chunks Parallel parameters.
#' @param tol Numeric tolerance.
#' @param include_warnings Logical; emit warning if the active-set call
#'   errors and the caller will fall back.
#' @param warn_context Character scalar identifying the caller in warning text.
#'
#' @return Either a list with \code{result} and \code{qp_info}, or
#'   \code{NULL} to signal dense fallback.
#'
#' @keywords internal
.try_woodbury_active_set <- function(result,
                                     K, p_expansions,
                                     A, R_constraints,
                                     constraint_value_vectors,
                                     family,
                                     qp_Amat, qp_bvec, qp_meq,
                                     rhs_list,
                                     Ghalf_corrected,
                                     wb_sqrt,
                                     parallel_aga, parallel_matmult,
                                     cl, chunk_size, num_chunks,
                                     rem_chunks, tol,
                                     include_warnings = TRUE,
                                     warn_context = ".try_woodbury_active_set") {

  as_tol <- max(tol, 100 * sqrt(.Machine$double.eps))
  as_out <- try(.active_set_refine_woodbury(
    result = result,
    K = K,
    p_expansions = p_expansions,
    A = A,
    R_constraints = R_constraints,
    constraint_value_vectors = constraint_value_vectors,
    family = family,
    qp_Amat = qp_Amat,
    qp_bvec = qp_bvec,
    qp_meq = qp_meq,
    rhs_list = rhs_list,
    Ghalf_corrected = Ghalf_corrected,
    wb_sqrt = wb_sqrt,
    parallel_aga = parallel_aga,
    parallel_matmult = parallel_matmult,
    cl = cl,
    chunk_size = chunk_size,
    num_chunks = num_chunks,
    rem_chunks = rem_chunks,
    tol = as_tol
  ), silent = TRUE)

  if (inherits(as_out, "try-error")) {
    if (include_warnings) {
      err_msg <- if (!is.null(attr(as_out, "condition"))) {
        conditionMessage(attr(as_out, "condition"))
      } else {
        as.character(as_out)
      }
      warning(
        warn_context, " failed: ",
        err_msg,
        "\n  falling back to dense QP"
      )
    }
    return(NULL)
  }

  if (!isTRUE(as_out$converged)) {
    if (include_warnings) {
      warning(
        warn_context,
        " active-set did not converge; falling back to dense QP"
      )
    }
    return(NULL)
  }

  list(result = as_out$result, qp_info = as_out$qp_info)
}


## Shared G recomputation wrapper

#' Recompute G at Current Coefficient Estimates
#'
#' Thin wrapper around \code{.solver_recompute_G_at_estimate()}.
#' The numerical core is shared with \code{blockfit_solve()} so the
#' final \code{G}, \code{Ghalf}, and \code{GhalfInv} objects are
#' constructed through one code path.
#'
#' @inheritParams .solver_recompute_G_at_estimate
#'
#' @return A list with components \code{G}, \code{Ghalf}, and
#'   \code{GhalfInv}, each a list of \eqn{K+1} matrices.
#'
#' @keywords internal
.recompute_G_at_estimate <- function(X, y, result, K, Lambda, family,
                                     order_list, glm_weight_function,
                                     schur_correction_function,
                                     need_dispersion_for_estimation,
                                     dispersion_function,
                                     observation_weights, VhalfInv,
                                     parallel_eigen, parallel_matmult,
                                     cl, chunk_size, num_chunks, rem_chunks,
                                     unique_penalty_per_partition,
                                     L_partition_list, ...) {
  .solver_recompute_G_at_estimate(
    X = X,
    y = y,
    result = result,
    K = K,
    Lambda = Lambda,
    family = family,
    order_list = order_list,
    glm_weight_function = glm_weight_function,
    schur_correction_function = schur_correction_function,
    need_dispersion_for_estimation = need_dispersion_for_estimation,
    dispersion_function = dispersion_function,
    observation_weights = observation_weights,
    VhalfInv = VhalfInv,
    parallel_eigen = parallel_eigen,
    parallel_matmult = parallel_matmult,
    cl = cl,
    chunk_size = chunk_size,
    num_chunks = num_chunks,
    rem_chunks = rem_chunks,
    unique_penalty_per_partition = unique_penalty_per_partition,
    L_partition_list = L_partition_list,
    ...
  )
}


## QP helpers

#' Detect Whether Inequality Constraints Need the Global Bridge
#'
#' Thin wrapper around \code{.solver_detect_qp_global()}.
#' This preserves the original helper name inside \code{get_B()} while
#' using the shared detection logic defined in \code{solver_utils.R}.
#'
#' @param qp_Amat Inequality constraint matrix (\eqn{P \times n_{ineq}}),
#'   or \code{NULL}.
#' @param p_expansions Integer; number of basis terms per partition.
#' @param K Integer; number of interior knots.
#'
#' @return Logical; \code{TRUE} if the active-set loop must use the
#'   global equality bridge, \code{FALSE} if the equality re-solves can
#'   remain partition-wise.
#'
#' @keywords internal
.detect_qp_global <- function(qp_Amat, p_expansions, K) {
  .solver_detect_qp_global(
    qp_Amat = qp_Amat,
    p_expansions = p_expansions,
    K = K
  )
}


## Partition-wise KKT check

#' Check KKT Conditions for Partition-Wise Active-Set Method
#'
#' Given current constrained coefficient estimates and a set of active
#' inequality constraints (treated as equalities in the last Lagrangian
#' projection), checks primal feasibility of inactive constraints and
#' dual feasibility (non-negative multipliers) of active constraints.
#'
#' Multipliers for active inequality constraints are recovered from the
#' OLS fit used in the Lagrangian projection: the fitted coefficients
#' on \eqn{\mathbf{X}^* = \mathbf{G}^{1/2}\mathbf{A}_{\mathrm{aug}}}
#' give the Lagrangian multipliers (up to sign and scaling).
#'
#' @param result List of \eqn{K+1} coefficient column vectors.
#' @param Ghalf List of \eqn{\mathbf{G}^{1/2}_k} matrices.
#' @param GhalfInv List of \eqn{\mathbf{G}^{-1/2}_k} matrices.
#' @param Xy_or_uncon Either the list of cross-products
#'   \eqn{\mathbf{X}_k^{\top}\mathbf{y}_k} (Path 2) or unconstrained
#'   estimates (Path 3).
#' @param is_path3 Logical; if TRUE, \code{Xy_or_uncon} contains
#'   unconstrained estimates and GhalfInv is used for transformation.
#' @param A_aug Augmented constraint matrix (original A plus active
#'   inequality columns).
#' @param n_eq_orig Integer; number of original equality constraints
#'   (columns of A before augmentation).
#' @param qp_Amat Full inequality constraint matrix.
#' @param qp_bvec Full inequality constraint RHS.
#' @param active_ineq Integer vector; indices into columns of
#'   \code{qp_Amat} that are currently active.
#' @param K,p_expansions Integer dimensions.
#' @param family GLM family object.
#' @param parallel_matmult,parallel_aga Logical flags.
#' @param cl,chunk_size,num_chunks,rem_chunks Parallel parameters.
#' @param tol Numeric tolerance for feasibility and multiplier checks.
#'
#' @return A list with components:
#' \describe{
#'   \item{feasible}{Logical; TRUE if all inactive inequality constraints
#'     are satisfied within tolerance.}
#'   \item{dual_feasible}{Logical; TRUE if all active inequality
#'     multipliers are non-negative within tolerance.}
#'   \item{violated}{Integer vector; indices of violated inactive
#'     constraints (into \code{qp_Amat} columns).}
#'   \item{drop}{Integer vector; indices of active constraints with
#'     negative multipliers that should be dropped.}
#'   \item{multipliers}{Numeric vector of Lagrangian multipliers for
#'     the active inequality constraints.}
#' }
#'
#' @keywords internal
.check_kkt_partitionwise <- function(result, Ghalf, GhalfInv,
                                     Xy_or_uncon, is_path3,
                                     A_aug, n_eq_orig,
                                     qp_Amat, qp_bvec,
                                     active_ineq,
                                     K, p_expansions, family,
                                     parallel_matmult, parallel_aga,
                                     cl, chunk_size, num_chunks,
                                     rem_chunks, tol) {

  beta_full <- cbind(unlist(result))
  primal <- .solver_check_primal_feasibility(
    beta_full = beta_full,
    qp_Amat = qp_Amat,
    qp_bvec = qp_bvec,
    active_ineq = active_ineq,
    tol = tol
  )

  ## Compute Lagrangian multipliers for active inequality constraints.
  #  The multipliers come from the dual of the OLS problem in the
  #  Lagrangian projection.  Specifically, if y* = G^{-1/2} beta_hat
  #  and X* = G^{1/2} A_aug, then the OLS coefficients on X* give
  #  the multipliers (for the augmented system).
  drop <- integer(0)
  multipliers <- numeric(0)

  if (length(active_ineq) > 0) {
    ## Reconstruct y* from unconstrained or cross-product info.
    if (is_path3) {
      GhalfXy <- cbind(unlist(
        matmult_block_diagonal(GhalfInv, Xy_or_uncon, K,
                               parallel_matmult, cl,
                               chunk_size, num_chunks, rem_chunks)))
    } else {
      GhalfXy <- cbind(unlist(
        matmult_block_diagonal(Ghalf, Xy_or_uncon, K,
                               parallel_matmult, cl,
                               chunk_size, num_chunks, rem_chunks)))
    }
    GhalfXy <- ifelse(is.na(GhalfXy), family$linkinv(0), GhalfXy)

    ## X* = G^{1/2} A_aug (block-diagonal multiplication).
    GhalfA_aug <- Reduce("rbind",
                         GAmult_wrapper(Ghalf, A_aug, K, p_expansions,
                                        ncol(A_aug),
                                        parallel_aga, cl,
                                        chunk_size, num_chunks,
                                        rem_chunks))

    mult_out <- .solver_recover_active_multipliers(
      X_star = GhalfA_aug,
      y_star = GhalfXy,
      n_eq_orig = n_eq_orig,
      K = K,
      tol = tol,
      ineq_sign = -1
    )
    drop <- mult_out$drop
    multipliers <- mult_out$multipliers
  }

  list(feasible = primal$feasible,
       dual_feasible = length(drop) == 0,
       violated = primal$violated,
       drop = drop,
       multipliers = multipliers)
}


## Partition-wise active set

#' Partition-Wise Active-Set Refinement for Inequality Constraints
#'
#' Replaces the dense \code{.qp_refine} SQP loop with a partition-wise
#' active-set method that reuses the existing Lagrangian projection
#' machinery. The key insight is that an active inequality constraint
#' can be treated as an additional equality constraint and absorbed into
#' the augmented constraint matrix \eqn{\mathbf{A}_{\mathrm{aug}}},
#' after which the standard \eqn{\mathbf{G}^{1/2}\mathbf{r}^*} trick
#' applies without ever forming the full \eqn{P \times P} information
#' matrix.
#'
#' The active-set loop:
#' \enumerate{
#'   \item Solve the equality-constrained subproblem (Lagrangian projection)
#'     with current active set.
#'   \item Check primal feasibility: any violated inactive inequalities?
#'   \item Check dual feasibility: any negative multipliers on active
#'     inequalities?
#'   \item If both satisfied, return (KKT conditions met).
#'   \item Otherwise, add most-violated constraint or drop most-negative
#'     multiplier, and repeat.
#' }
#'
#' Falls back to \code{.qp_refine} if the active-set method does not
#' converge within \code{max_as_iter} iterations.
#'
#' @param result List of current coefficient column vectors by partition.
#' @param X,y Lists of partition-specific design matrices and responses.
#' @param K,p_expansions Integer dimensions.
#' @param A Original equality constraint matrix.
#' @param R_constraints Number of columns of \code{A}.
#' @param constraint_value_vectors Constraint RHS list.
#' @param Lambda Shared penalty matrix.
#' @param Ghalf,GhalfInv Lists of \eqn{\mathbf{G}^{1/2}_k} and
#'   \eqn{\mathbf{G}^{-1/2}_k} matrices.
#' @param family GLM family object.
#' @param qp_Amat Inequality constraint matrix.
#' @param qp_bvec Inequality constraint RHS.
#' @param qp_meq Number of equality constraints within \code{qp_Amat}
#'   (leading columns).
#' @param Xy_or_uncon Cross-products (Path 2) or unconstrained estimates
#'   (Path 3).
#' @param is_path3 Logical.
#' @param parallel_aga,parallel_matmult Logical flags.
#' @param cl,chunk_size,num_chunks,rem_chunks Parallel parameters.
#' @param tol Convergence tolerance.
#' @param max_as_iter Maximum active-set iterations (default 50).
#'
#' @return A list with components:
#' \describe{
#'   \item{result}{List of refined coefficient column vectors by partition.}
#'   \item{qp_info}{List with active constraint information, or NULL.}
#'   \item{converged}{Logical; TRUE if active-set method converged.}
#' }
#'
#' @keywords internal
.active_set_refine <- function(result, X, y, K, p_expansions,
                               A, R_constraints,
                               constraint_value_vectors,
                               Lambda, Ghalf, GhalfInv,
                               family,
                               qp_Amat, qp_bvec, qp_meq,
                               Xy_or_uncon, is_path3,
                               parallel_aga, parallel_matmult,
                               cl, chunk_size, num_chunks,
                               rem_chunks, tol,
                               max_as_iter = NULL) {
  n_ineq <- if (is.null(qp_Amat)) 0L else ncol(qp_Amat)

  ## Match the Woodbury active-set budget: many observation-local
  #  derivative constraints need more than 50 add/drop iterations.
  if (is.null(max_as_iter)) {
    max_as_iter <- max(200, min(1000, 10 * max(1, n_ineq)))
  }

  if (is_path3) {
    GhalfXy <- cbind(unlist(
      matmult_block_diagonal(GhalfInv, Xy_or_uncon, K,
                             parallel_matmult, cl,
                             chunk_size, num_chunks, rem_chunks)))
  } else {
    GhalfXy <- cbind(unlist(
      matmult_block_diagonal(Ghalf, Xy_or_uncon, K,
                             parallel_matmult, cl,
                             chunk_size, num_chunks, rem_chunks)))
  }

  .solver_active_set(
    result = result,
    A = A,
    constraint_value_vectors = constraint_value_vectors,
    qp_Amat = qp_Amat,
    qp_bvec = qp_bvec,
    qp_meq = qp_meq,
    tol = tol,
    max_as_iter = max_as_iter,
    solve_subproblem = function(aug_state) {
      .lagrangian_project(
        GhalfXy = GhalfXy,
        Ghalf = Ghalf,
        A = aug_state$A_aug,
        K = K,
        p_expansions = p_expansions,
        R_constraints = ncol(aug_state$A_aug),
        constraint_value_vectors = aug_state$cv_for_proj,
        family = family,
        parallel_aga = parallel_aga,
        parallel_matmult = parallel_matmult,
        cl = cl,
        chunk_size = chunk_size,
        num_chunks = num_chunks,
        rem_chunks = rem_chunks
      )
    },
    kkt_subproblem = function(result_new, aug_state) {
      .check_kkt_partitionwise(
        result = result_new,
        Ghalf = Ghalf,
        GhalfInv = GhalfInv,
        Xy_or_uncon = Xy_or_uncon,
        is_path3 = is_path3,
        A_aug = aug_state$A_aug,
        n_eq_orig = aug_state$n_eq_aug,
        qp_Amat = aug_state$qp_ineq_Amat,
        qp_bvec = aug_state$qp_ineq_bvec,
        active_ineq = match(aug_state$active_ineq_kept, aug_state$ineq_cols),
        K = K,
        p_expansions = p_expansions,
        family = family,
        parallel_matmult = parallel_matmult,
        parallel_aga = parallel_aga,
        cl = cl,
        chunk_size = chunk_size,
        num_chunks = num_chunks,
        rem_chunks = rem_chunks,
        tol = tol
      )
    },
    method = "active_set",
    debug_label = "dense-as"
  )
}


## Dense QP fallback

.solve_qp_stable <- function(Dmat, dvec, Amat, bvec, meq) {

  qp_result <- try(quadprog::solve.QP(
    Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = meq
  ), silent = TRUE)

  if (!inherits(qp_result, "try-error")) return(qp_result)

  ## Dense QP fallback: quadprog requires strictly positive-definite Dmat.
  Dsym <- 0.5 * (Dmat + t(Dmat))
  Dscale <- max(mean(abs(Dsym)), mean(abs(diag(Dsym))), 1)
  eig_min <- try(min(eigen(Dsym, symmetric = TRUE,
                           only.values = TRUE)$values), silent = TRUE)
  if (inherits(eig_min, "try-error") || !is.finite(eig_min)) {
    eig_min <- 0
  }
  ridge <- max(0, -eig_min) + sqrt(.Machine$double.eps) * Dscale
  Dstable <- Dsym + diag(ridge, nrow(Dsym))

  try(quadprog::solve.QP(
    Dmat = Dstable, dvec = dvec, Amat = Amat, bvec = bvec, meq = meq
  ), silent = TRUE)
}


.attach_qp_fallback_info <- function(out, reason, wb_decomp = NULL) {

  if (is.null(out$qp_info)) return(out)

  out$qp_info$fallback_reason <- reason
  if (!is.null(wb_decomp)) {
    out$qp_info$woodbury_reason <- wb_decomp$reason
    out$qp_info$woodbury_rank <- wb_decomp$r
    out$qp_info$woodbury_P <- wb_decomp$P
    out$qp_info$woodbury_rank_threshold <- wb_decomp$rank_threshold
    out$qp_info$woodbury_nnz_fraction <- wb_decomp$nnz_fraction
  }
  out
}


#' Quadratic Programming Refinement for Inequality Constraints
#'
#' After the Lagrangian projection handles smoothness equality constraints,
#' this function refines the estimate to satisfy additional inequality
#' constraints (monotonicity, derivative sign, range bounds, or user-supplied
#' constraints) via \code{quadprog::solve.QP}.
#'
#' The subproblem at each iteration is a second-order Taylor approximation
#' of the penalized log-likelihood around the current iterate
#' \eqn{\boldsymbol{\beta}^*}. Collecting terms, this yields:
#' \deqn{\tilde{\boldsymbol{\beta}} = \arg\min_{\boldsymbol{\beta}}
#'   \left\{-\mathbf{d}^{\top}\boldsymbol{\beta} + \frac{1}{2}
#'   \boldsymbol{\beta}^{\top}\mathbf{G}^{-1}\boldsymbol{\beta}\right\}
#'   \quad \text{s.t.} \quad
#'   \mathbf{A}^{\top}\boldsymbol{\beta} = \mathbf{0}, \quad
#'   \mathbf{C}^{\top}\boldsymbol{\beta} \succeq \mathbf{c}}
#' where \eqn{\mathbf{d}} is the score adjusted by the current iterate
#' and \eqn{\mathbf{c}} is the constraint value vector. Step acceptance
#' uses damped updates with deviance monitoring; see Nocedal and Wright
#' (2006) for the general SQP framework.
#'
#' @param result List of current coefficient column vectors by partition.
#' @param X List of partition-specific design matrices.
#' @param y List of response vectors by partition.
#' @param K Integer; number of interior knots.
#' @param p_expansions Integer; number of basis terms per partition.
#' @param A Equality constraint matrix \eqn{\mathbf{A}}.
#' @param Lambda Shared penalty matrix \eqn{\boldsymbol{\Lambda}}.
#' @param Lambda_block Full block-diagonal penalty matrix.
#' @param family GLM family object.
#' @param iterate Logical; if \code{TRUE}, iterate the SQP loop until
#'   convergence rather than taking a single step.
#' @param tol Convergence tolerance.
#' @param qp_Amat Inequality constraint matrix \eqn{\mathbf{C}} for
#'   \code{solve.QP}.
#' @param qp_bvec Inequality constraint vector \eqn{\mathbf{c}}.
#' @param qp_meq Number of equality constraints within \code{qp_Amat}.
#' @param qp_score_function Score function
#'   \eqn{\nabla_{\boldsymbol{\beta}}\ell(\boldsymbol{\beta}^*)} for the
#'   QP step.
#' @param order_list List of index vectors mapping partition rows to
#'   original data ordering.
#' @param glm_weight_function Function computing GLM working weights.
#' @param schur_correction_function Function computing Schur corrections.
#' @param need_dispersion_for_estimation Logical.
#' @param dispersion_function Dispersion estimation function.
#' @param observation_weights List of observation weights by partition.
#' @param VhalfInv Inverse square root of the working correlation matrix in
#'   the original observation ordering (or \code{NULL}); permuted internally
#'   via \code{order_list}.
#' @param ... Passed to weight, correction, dispersion, and score functions.
#'
#' @return A list with components:
#' \describe{
#'   \item{result}{List of refined coefficient column vectors by partition.}
#'   \item{qp_info}{List with QP solve metadata including Lagrangian
#'     multipliers and the active constraint matrix.}
#' }
#'
#' @keywords internal
.qp_refine <- function(result, X, y, K, p_expansions, A, Lambda,
                       Lambda_block, family, iterate, tol,
                       qp_Amat, qp_bvec, qp_meq,
                       qp_score_function,
                       order_list, glm_weight_function,
                       schur_correction_function,
                       need_dispersion_for_estimation,
                       dispersion_function, observation_weights,
                       VhalfInv, ...) {

  ## Assemble the full N x P block-diagonal design matrix by stacking
  #  partition blocks on the diagonal.
  X_block <- Reduce("rbind", lapply(1:(K + 1), function(k) {
    Reduce("cbind", lapply(1:(K + 1), function(j) {
      if (nrow(X[[k]]) == 0) {
        return(X[[k]])
      } else if (j == k) X[[k]] else 0 * X[[k]]
    }))
  }))
  beta_block <- cbind(unlist(result))
  y_block <- cbind(unlist(y))

  ## Damped SQP loop: iterate Newton-like QP steps with backtracking.
  damp_cnt <- 0
  master_cnt <- 0
  err <- Inf
  XB <- X_block %**% beta_block
  qp_info_out <- NULL

  while (err > tol & damp_cnt < 10 & master_cnt < 100) {
    master_cnt <- master_cnt + 1
    damp <- 2^(-(damp_cnt))

    if (need_dispersion_for_estimation) {
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
    W <- c(glm_weight_function(family$linkinv(XB),
                               y_block,
                               unlist(order_list),
                               family,
                               dispersion_temp,
                               unlist(observation_weights),
                               ...))
    W <- .solver_stabilize_working_weights(W)

    ## Schur correction in per-partition form, then collapse to full matrix.
    result <- lapply(1:(K + 1), function(k) {
      cbind(beta_block[(k - 1) * p_expansions + 1:p_expansions])
    })
    schur_correction <-
      schur_correction_function(
        X, y, result, dispersion_temp, order_list, K, family,
        observation_weights, ...
      )
    if (!(any(unlist(schur_correction) != 0))) {
      schur_correction <- 0
    } else {
      schur_correction <- collapse_block_diagonal(schur_correction)
    }

    ## Information matrix: X^T W X + Lambda + Schur correction.
    info <- crossprod(X_block, W * X_block) + Lambda_block + schur_correction

    sc <- sqrt(mean(abs(info)))

    ## QP step: equality (smoothness) and inequality constraints combined.
    #  The dvec encodes the score adjusted for the current iterate so that
    #  the solution is a full position estimate, not a step direction.
    qp_score <- qp_score_function(
      X_block, y_block,
      cbind(family$linkinv(XB)),
      unlist(order_list), dispersion_temp, NULL,
      unlist(observation_weights), ...
    )
    qp_result <- .solve_qp_stable(
      Dmat = info / sc,
      dvec = (qp_score -
                Lambda_block %**% beta_block +
                info %**% beta_block) / sc,
      Amat = cbind(A, qp_Amat),
      bvec = c(rep(0, ncol(A)), qp_bvec),
      meq = ncol(A) + qp_meq
    )

    if (any(inherits(qp_result, 'try-error'))) {
      beta_new <- 0 * beta_block
    } else {
      beta_new <- qp_result$solution
      ## Store most recent successful solve for multiplier extraction.
      qp_info_out <- list(
        lagrangian = qp_result$Lagrangian,
        Amat_active = cbind(A, qp_Amat)[,
                                        c(1:ncol(A),
                                          ncol(A) +
                                            which(abs(qp_result$Lagrangian[-(1:ncol(A))]) >
                                                    sqrt(.Machine$double.eps))),
                                        drop = FALSE]
      )
    }

    if (!iterate & master_cnt > 1) {
      ## Non-iterative mode: accept the single QP step and exit.
      beta_block <- beta_new
      damp_cnt <- 11
      master_cnt <- 101
      err <- tol - 1
    } else {
      ## Damped update: blend old and new iterates, check deviance.
      beta_new <- (1 - damp) * beta_block + damp * beta_new
      XB <- X_block %**% beta_new
      if (!is.null(family$custom_dev.resids) &
          is.null(family$dev.resids)) {
        err_new <- mean(
          family$custom_dev.resids(y_block,
                                   family$linkinv(c(XB)),
                                   unlist(order_list),
                                   family,
                                   unlist(observation_weights),
                                   ...))
      } else if (is.null(family$dev.resids)) {
        err_new <- mean(unlist(observation_weights) *
                          (y_block - family$linkinv(XB))^2)
      } else {
        err_new <-
          mean(unlist(observation_weights) *
                 family$dev.resids(y_block,
                                   family$linkinv(XB),
                                   wt = 1))
      }

      if (is.null(err_new) | is.na(err_new) | !is.finite(err_new)) {
        ## Non-finite deviance: increase damping and retry.
        damp_cnt <- damp_cnt + 1
      } else if (err_new <= err) {
        ## Deviance decreased: accept step and reset damping.
        prev_err <- err
        err <- err_new
        abs_diff <- max(abs(beta_new - beta_block))
        beta_block <- beta_new
        damp_cnt <- 0

        ## Early exit when both coefficient change and deviance decrease
        #  are below tolerance (after a burn-in of 10 iterations).
        if ((abs_diff < tol) &
            (prev_err - err < tol) &
            (master_cnt > 10)) {
          damp_cnt <- 11
          master_cnt <- 101
          err <- tol - 1
        }

      } else {
        ## Deviance increased: increase damping and retry.
        damp_cnt <- damp_cnt + 1
      }
    }

    ## Unpack full coefficient vector back into per-partition list.
    result <- lapply(1:(K + 1), function(k) {
      cbind(beta_block[1:p_expansions + (k - 1) * p_expansions])
    })
  }

  list(result = result, qp_info = qp_info_out)
}


## Woodbury decomposition

#' Woodbury Decomposition of Correlation Structure
#'
#' Builds the Woodbury factors for the correlated solver. The function
#' decomposes
#' \eqn{\mathbf{G}^{-1} = \mathbf{G}_{\mathrm{on}}^{-1} +
#' \mathbf{G}_{\mathrm{off}}^{-1}}
#' by absorbing the block-diagonal part of
#' \eqn{\mathbf{X}^{\top}\mathbf{V}^{-1}\mathbf{X} + \boldsymbol{\Lambda}}
#' into \eqn{\mathbf{G}_{\mathrm{on}}^{-1}}, then factoring the
#' cross-partition part as
#' \eqn{\mathbf{G}_{\mathrm{off}}^{-1} = \mathbf{E}\mathbf{J}\mathbf{E}^{\top}}.
#'
#' The output is exactly the information needed by
#' \code{.woodbury_halfsqrt_components()} and
#' \code{.lagrangian_project_woodbury()}:
#' corrected block-diagonal \eqn{\mathbf{G}_{\mathrm{on}}},
#' the low-rank factor \eqn{\mathbf{E}}, the sign matrix diagonal
#' \eqn{\mathbf{J}}, and the small inverse
#' \eqn{(\mathbf{J}^{-1} +
#' \mathbf{E}^{\top}\mathbf{G}_{\mathrm{on}}\mathbf{E})^{-1}}.
#'
#' When the cross-partition rank is zero or too large to be worthwhile,
#' the helper returns \code{use_woodbury = FALSE} so the caller can use
#' the dense correlated path instead.
#'
#' @param VhalfInv_perm \eqn{\mathbf{V}^{-1/2}} permuted to partition
#'   ordering.
#' @param X List of partition-specific design matrices.
#' @param K,p_expansions Integer dimensions.
#' @param Lambda,Lambda_block Shared and full block-diagonal penalty
#'   matrices.
#' @param unique_penalty_per_partition Logical.
#' @param L_partition_list Partition-specific penalty matrices.
#' @param order_list Partition-to-data index mapping.
#' @param parallel_eigen Logical; parallel eigendecomposition.
#' @param cl,chunk_size,num_chunks,rem_chunks Parallel parameters.
#' @param rank_threshold_fraction Numeric; Woodbury is used only when
#'   \eqn{r < P \times} this fraction (default 1/3).
#' @param family GLM family object.
#'
#' @return A list with components:
#' \describe{
#'   \item{use_woodbury}{Logical; FALSE triggers dense Path 1 fallback.}
#'   \item{Ghalf_corrected}{List of corrected \eqn{\mathbf{G}_{\mathrm{on},k}^{1/2}}.}
#'   \item{GhalfInv_corrected}{List of corrected \eqn{\mathbf{G}_{\mathrm{on},k}^{-1/2}}.}
#'   \item{G_corrected}{List of corrected \eqn{\mathbf{G}_{\mathrm{on},k}}.}
#'   \item{E}{Off-diagonal low-rank factor (\eqn{P \times r}) from the
#'     supplement factorization
#'     \eqn{\mathbf{G}_{\mathrm{off}}^{-1} = \mathbf{E}\mathbf{J}\mathbf{E}^{\top}}.}
#'   \item{J}{Length-\eqn{r} vector giving the diagonal of
#'     \eqn{\mathbf{J}}.}
#'   \item{r}{Integer effective rank.}
#'   \item{inner_inv}{Precomputed \eqn{r \times r} Woodbury inner inverse. In
#'     supplement notation this is
#'     \eqn{(\mathbf{J}^{-1} + \mathbf{E}^{\top}\mathbf{G}_{\mathrm{on}}
#'     \mathbf{E})^{-1}}.}
#'   \item{GL}{List of \eqn{K+1} matrices storing
#'     \eqn{\mathbf{G}_{\mathrm{on},k}\mathbf{E}[rows_k,]} in supplement notation,
#'     each \eqn{p \times r}.}
#' }
#'
#' @keywords internal
.woodbury_decompose_V <- function(VhalfInv_perm, X, K, p_expansions,
                                  Lambda, Lambda_block,
                                  unique_penalty_per_partition,
                                  L_partition_list,
                                  order_list,
                                  parallel_eigen, cl,
                                  chunk_size, num_chunks, rem_chunks,
                                  rank_threshold_fraction = 1/3,
                                  family) {

  P <- p_expansions * (K + 1)

  ## Step 1a: form V^{-1} and the correlation-induced perturbation
  ## V^{-1} - I in partition ordering.
  Vinv <- crossprod(VhalfInv_perm)
  N_obs <- nrow(Vinv)
  Delta_V <- Vinv - diag(N_obs)

  ## Sparsity check: if V^{-1} - I is too dense, Woodbury is not
  #  worthwhile (the off-diagonal rank will be high).
  nnz <- sum(abs(Delta_V) > sqrt(.Machine$double.eps))
  nnz_fraction <- nnz / (N_obs * N_obs)
  if (nnz_fraction > 0.5) {
    return(list(use_woodbury = FALSE,
                reason = "dense_correlation_perturbation",
                nnz_fraction = nnz_fraction,
                P = P))
  }

  ## Step 1b: assemble the cross-partition information correction
  ## X^T (V^{-1} - I) X.
  X_block <- collapse_block_diagonal(X)
  Delta_mat <- crossprod(X_block, Delta_V %**% X_block)

  ## Step 2: absorb the block-diagonal part into the corrected
  ## per-partition matrices defining G_on^{-1}.
  X_gram_corrected <- vector("list", K + 1)
  for (k in 1:(K + 1)) {
    rows_k <- ((k - 1) * p_expansions + 1):(k * p_expansions)
    Delta_diag_k <- Delta_mat[rows_k, rows_k]
    X_gram_corrected[[k]] <- crossprod(X[[k]]) + Delta_diag_k
  }

  ## Compute corrected G via eigendecomposition (reuses existing infra).
  #  Pass zero Schur corrections for the initial Gaussian decomposition.
  schur_zero <- lapply(1:(K + 1), function(k) {
    matrix(0, p_expansions, p_expansions)
  })
  G_list_c <- compute_G_eigen(X_gram_corrected, Lambda, K,
                              parallel_eigen, cl,
                              chunk_size, num_chunks, rem_chunks, family,
                              unique_penalty_per_partition,
                              L_partition_list,
                              keep_G = TRUE, schur_zero)

  ## Step 3: remove the absorbed diagonal blocks, leaving only the
  ## cross-partition remainder.
  Delta_offdiag <- Delta_mat
  for (k in 1:(K + 1)) {
    rows_k <- ((k - 1) * p_expansions + 1):(k * p_expansions)
    Delta_offdiag[rows_k, rows_k] <- 0
  }

  ## Step 4: eigentruncate the remainder and encode it as E J E^T.
  eig_off <- eigen(Delta_offdiag, symmetric = TRUE)
  r <- sum(abs(eig_off$values) >
             max(abs(eig_off$values)) * sqrt(.Machine$double.eps))

  if (r == 0) {
    ## No cross-partition coupling; corrected G is exact.
    return(list(
      use_woodbury = FALSE,
      reason = "no_cross_partition_coupling",
      r = r,
      P = P,
      rank_threshold = P * rank_threshold_fraction,
      nnz_fraction = nnz_fraction,
      Ghalf_corrected = G_list_c$Ghalf,
      GhalfInv_corrected = G_list_c$GhalfInv,
      G_corrected = G_list_c$G
    ))
  }

  if (r >= P * rank_threshold_fraction) {
    ## Low-rank approximation not worthwhile; fall back to dense.
    return(list(use_woodbury = FALSE,
                reason = "offdiagonal_rank_too_high",
                r = r,
                P = P,
                rank_threshold = P * rank_threshold_fraction,
                nnz_fraction = nnz_fraction))
  }

  ## Supplement notation: E J E^T for the retained off-diagonal remainder.
  E <- eig_off$vectors[, 1:r, drop = FALSE] %**%
    diag(sqrt(abs(eig_off$values[1:r])), r, r)
  J <- sign(eig_off$values[1:r])

  ## Precompute G_{on,k} E_k blockwise for the later Woodbury inverse and
  ## square-root steps.
  GL <- lapply(1:(K + 1), function(k) {
    rows_k <- ((k - 1) * p_expansions + 1):(k * p_expansions)
    G_list_c$G[[k]] %**% E[rows_k, , drop = FALSE]
  })

  ## Assemble E^T G_on E on the retained subspace.
  GL_full <- do.call(rbind, GL)
  LtGL <- crossprod(E, GL_full)

  ## Precompute (J^{-1} + E^T G_on E)^{-1}; this is the only dense inverse
  ## in the Woodbury correction.
  inner_inv <- invert(diag(1 / J, r, r) + LtGL)

  list(
    use_woodbury = TRUE,
    Ghalf_corrected = G_list_c$Ghalf,
    GhalfInv_corrected = G_list_c$GhalfInv,
    G_corrected = G_list_c$G,
    E = E,
    J = J,
    r = r,
    P = P,
    rank_threshold = P * rank_threshold_fraction,
    nnz_fraction = nnz_fraction,
    inner_inv = inner_inv,
    GL = GL
  )
}


## Weighted Woodbury update

#' Per-Iteration Weighted Woodbury Redecomposition
#'
#' Weighted analogue of \code{.woodbury_decompose_V()}. At each GLM
#' iteration, this helper rebuilds the current correlated information,
#' absorbs the block-diagonal part into
#' \eqn{\mathbf{G}_{\mathrm{on}}^{-1}}, and refactors the cross-partition
#' remainder as \eqn{\mathbf{E}\mathbf{J}\mathbf{E}^{\top}}.
#'
#' It is standalone in the sense that it returns the full Woodbury state
#' for the current working weights:
#' corrected block-diagonal factors, the low-rank basis, and the small
#' inverse \eqn{(\mathbf{J}^{-1} +
#' \mathbf{E}^{\top}\mathbf{G}_{\mathrm{on}}\mathbf{E})^{-1}}.
#'
#' @param Delta_V Fixed \eqn{N \times N} matrix
#'   \eqn{\mathbf{V}^{-1} - \mathbf{I}}.
#' @param X_block Full \eqn{N \times P} block-diagonal design matrix.
#' @param DV_X Precomputed \eqn{N \times P} product
#'   of \eqn{(\mathbf{V}^{-1} - \mathbf{I})\mathbf{X}}.
#' @param W Length-\eqn{N} vector of current GLM working weights.
#' @param X List of partition-specific design matrices.
#' @param K,p_expansions Integer dimensions.
#' @param X_gram_weighted List of \eqn{K+1} weighted Gram matrices
#'   \eqn{\mathbf{X}_k^{\top}\mathrm{diag}(\mathbf{W}_k)\mathbf{X}_k}.
#' @param Lambda Shared \eqn{p \times p} penalty matrix.
#' @param schur_corrections List of \eqn{K+1} Schur correction matrices.
#' @param unique_penalty_per_partition Logical.
#' @param L_partition_list Partition-specific penalty matrices.
#' @param parallel_eigen Logical; parallel eigendecomposition.
#' @param cl,chunk_size,num_chunks,rem_chunks Parallel parameters.
#' @param rank_threshold_fraction Numeric; Woodbury threshold (default 1/3).
#' @param family GLM family object.
#'
#' @return A list with the same structure as \code{.woodbury_decompose_V}:
#' \describe{
#'   \item{use_woodbury}{Logical; FALSE if rank exceeds threshold.}
#'   \item{Ghalf_corrected}{List of corrected \eqn{\mathbf{G}_{\mathrm{on},k}^{1/2}}.}
#'   \item{GhalfInv_corrected}{List of corrected \eqn{\mathbf{G}_{\mathrm{on},k}^{-1/2}}.}
#'   \item{G_corrected}{List of corrected \eqn{\mathbf{G}_{\mathrm{on},k}}.}
#'   \item{E}{Off-diagonal low-rank factor (\eqn{P \times r}), denoted by
#'     \eqn{\mathbf{E}} in the supplement.}
#'   \item{J}{Length-\eqn{r} sign vector, i.e.\ the diagonal of
#'     \eqn{\mathbf{J}}.}
#'   \item{r}{Integer effective rank.}
#'   \item{inner_inv}{Precomputed \eqn{r \times r} Woodbury inner inverse; in
#'     supplement notation,
#'     \eqn{(\mathbf{J}^{-1} + \mathbf{E}^{\top}\mathbf{G}_{\mathrm{on}}
#'     \mathbf{E})^{-1}}.}
#'   \item{GL}{List of per-partition matrices corresponding to
#'     \eqn{\mathbf{G}_{\mathrm{on},k}\mathbf{E}[rows_k,]} in supplement notation.}
#' }
#'
#' @keywords internal
.woodbury_redecompose_weighted <- function(Delta_V, X_block, DV_X, W,
                                           X, K, p_expansions,
                                           X_gram_weighted,
                                           Lambda, schur_corrections,
                                           unique_penalty_per_partition,
                                           L_partition_list,
                                           parallel_eigen, cl,
                                           chunk_size, num_chunks,
                                           rem_chunks,
                                           rank_threshold_fraction = 1/3,
                                           family) {

  P <- p_expansions * (K + 1)
  W <- .solver_stabilize_working_weights(W)

  ## Weighted analogue of the cross-partition information correction:
  ## X^T diag(W) (V^{-1} - I) X.
  ## Since DV_X = (V^{-1} - I) X is precomputed, this is:
  #    crossprod(X_block, W * DV_X)
  #  which costs O(N * P^2) -- cheaper than the P^3 eigen that the
  #  dense path requires.
  Delta_mat_W <- crossprod(X_block, c(W) * DV_X)
  if (any(!is.finite(Delta_mat_W))) {
    return(list(use_woodbury = FALSE))
  }
  Delta_mat_W <- 0.5 * (Delta_mat_W + t(Delta_mat_W))

  ## Absorb the weighted block-diagonal portion into the corrected
  ## per-partition information defining G_on^{-1} at the current iterate.
  X_gram_corrected <- vector("list", K + 1)
  for (k in 1:(K + 1)) {
    rows_k <- ((k - 1) * p_expansions + 1):(k * p_expansions)
    Delta_diag_W_k <- Delta_mat_W[rows_k, rows_k]
    X_gram_corrected[[k]] <- X_gram_weighted[[k]] + Delta_diag_W_k
  }

  ## Eigendecompose per partition via compute_G_eigen.
  #  This adds Lambda + Schur corrections and produces Ghalf, GhalfInv, G.
  G_list_c <- compute_G_eigen(X_gram_corrected, Lambda, K,
                              parallel_eigen, cl,
                              chunk_size, num_chunks, rem_chunks, family,
                              unique_penalty_per_partition,
                              L_partition_list,
                              keep_G = TRUE, schur_corrections)
  if (any(vapply(G_list_c$G, is.null, logical(1))) ||
      any(vapply(G_list_c$Ghalf, is.null, logical(1))) ||
      any(vapply(G_list_c$GhalfInv, is.null, logical(1)))) {
    return(list(use_woodbury = FALSE))
  }

  ## Remove the absorbed diagonal blocks, leaving the weighted
  ## cross-partition remainder only.
  Delta_offdiag_W <- Delta_mat_W
  for (k in 1:(K + 1)) {
    rows_k <- ((k - 1) * p_expansions + 1):(k * p_expansions)
    Delta_offdiag_W[rows_k, rows_k] <- 0
  }

  ## Eigentruncate the weighted remainder and encode it as E J E^T.
  eig_off <- eigen(Delta_offdiag_W, symmetric = TRUE)
  r <- sum(abs(eig_off$values) >
             max(abs(eig_off$values)) * sqrt(.Machine$double.eps))

  if (r == 0) {
    ## No cross-partition coupling at this W; corrected G is exact.
    return(list(
      use_woodbury = TRUE,
      Ghalf_corrected = G_list_c$Ghalf,
      GhalfInv_corrected = G_list_c$GhalfInv,
      G_corrected = G_list_c$G,
      E = matrix(0, P, 0),
      J = numeric(0),
      r = 0L,
      inner_inv = matrix(0, 0, 0),
      GL = lapply(1:(K + 1), function(k) matrix(0, p_expansions, 0))
    ))
  }

  if (r >= P * rank_threshold_fraction) {
    ## Rank too high at this iteration; signal fallback to dense.
    return(list(use_woodbury = FALSE))
  }

  ## Supplement notation: E J E^T for the weighted off-diagonal remainder.
  E <- eig_off$vectors[, 1:r, drop = FALSE] %**%
    diag(sqrt(abs(eig_off$values[1:r])), r, r)
  J <- sign(eig_off$values[1:r])

  ## Precompute G_{c,k} E_k blockwise for the current iterate.
  GL <- lapply(1:(K + 1), function(k) {
    rows_k <- ((k - 1) * p_expansions + 1):(k * p_expansions)
    G_list_c$G[[k]] %**% E[rows_k, , drop = FALSE]
  })

  ## Assemble E^T G_on E on the retained subspace.
  GL_full <- do.call(rbind, GL)
  LtGL <- crossprod(E, GL_full)

  ## Precompute (J^{-1} + E^T G_on E)^{-1} on the retained subspace.
  inner_inv <- invert(diag(1 / J, r, r) + LtGL)

  list(
    use_woodbury = TRUE,
    Ghalf_corrected = G_list_c$Ghalf,
    GhalfInv_corrected = G_list_c$GhalfInv,
    G_corrected = G_list_c$G,
    E = E,
    J = J,
    r = r,
    inner_inv = inner_inv,
    GL = GL
  )
}


## Woodbury half-square-root pieces

#' Compute Woodbury Half-Square-Root Components
#'
#' Converts the low-rank Woodbury decomposition into square-root objects
#' used directly by the solver. Starting from
#' \eqn{\mathbf{N} = \mathbf{G}_{\mathrm{on}}^{1/2}\mathbf{E}}, it builds
#' a basis for the supplement matrix \eqn{\mathbf{Q}} and stores
#' \eqn{\mathbf{F}^{1/2}} and \eqn{\mathbf{F}^{-1/2}} in the internal forms
#' \eqn{\mathbf{I} - \mathbf{Q}\mathbf{C}\mathbf{Q}^{\top}} and
#' \eqn{\mathbf{I} + \mathbf{Q}\mathbf{C}_{\mathrm{inv}}\mathbf{Q}^{\top}}.
#'
#' This helper is standalone: its output is the complete square-root state
#' consumed by \code{.lagrangian_project_woodbury()}.
#'
#' @param Ghalf_corrected List of corrected \eqn{\mathbf{G}_{\mathrm{on},k}^{1/2}}.
#' @param E Off-diagonal low-rank factor (\eqn{P \times r}), denoted
#'   \eqn{\mathbf{E}} in the supplement.
#' @param J Length-\eqn{r} sign vector corresponding to the
#'   diagonal of \eqn{\mathbf{J}}.
#' @param r Integer effective rank.
#' @param K,p_expansions Integer dimensions.
#'
#' @return A list with components:
#' \describe{
#'   \item{valid}{Logical; FALSE if \eqn{\mathbf{F}} is not positive
#'     definite (signals fallback to dense).}
#'   \item{Q}{\eqn{P \times r} basis corresponding to the supplement factor
#'     \eqn{\mathbf{Q}} in the QR decomposition of
#'     \eqn{\mathbf{N} = \mathbf{G}_{\mathrm{on}}^{1/2}\mathbf{E}}. This spans the same
#'     subspace as that QR factor.}
#'   \item{C}{\eqn{r \times r} diagonal matrix
#'     encoding the rank-\eqn{r} update in the internal form
#'     \eqn{\mathbf{F}^{1/2} = \mathbf{I}_P - \mathbf{Q}\mathbf{C}
#'     \mathbf{Q}^{\top}}.}
#'   \item{C_inv}{\eqn{r \times r} diagonal matrix
#'     encoding the inverse update in the internal form
#'     \eqn{\mathbf{F}^{-1/2} = \mathbf{I}_P +
#'     \mathbf{Q}\mathbf{C}_{\mathrm{inv}}\mathbf{Q}^{\top}}.}
#'   \item{GhalfQ}{List of \eqn{K+1} matrices, each \eqn{p \times r}:
#'     \eqn{\mathbf{G}_{\mathrm{on},k}^{1/2} \mathbf{Q}[rows_k,]}.}
#' }
#'
#' @keywords internal
.woodbury_halfsqrt_components <- function(Ghalf_corrected, E, J,
                                          r, K, p_expansions) {
  P <- p_expansions * (K + 1)

  ## Handle the degenerate case r = 0 (no off-diagonal coupling).
  #  Return an identity-like structure with trivial correction matrices.
  if (r == 0) {
    GhalfQ <- lapply(1:(K + 1), function(k) {
      matrix(0, p_expansions, 0)
    })
    return(list(
      valid = TRUE,
      Q = matrix(0, P, 0),
      C = matrix(0, 0, 0),
      C_inv = matrix(0, 0, 0),
      GhalfQ = GhalfQ
    ))
  }

  ## Form the supplement matrix N = G_on^{1/2} E.
  N_mat <- matrix(0, P, r)
  for (k in 1:(K + 1)) {
    rows_k <- ((k - 1) * p_expansions + 1):(k * p_expansions)
    N_mat[rows_k, ] <- Ghalf_corrected[[k]] %**% E[rows_k, , drop = FALSE]
  }

  ## Thin SVD of manuscript N. The resulting basis spans the same subspace as
  ## the manuscript QR factor Q. The square-root itself is computed from
  #  the small r x r F matrix, preserving mixed-sign Woodbury updates.
  svd_Q <- svd(N_mat, nu = r, nv = r)
  Q <- svd_Q$u[, 1:r, drop = FALSE]

  N_in_Q <- crossprod(Q, N_mat)
  NtN <- crossprod(N_mat)
  N_inner_inv <- invert(diag(1 / J, r, r) + NtN)
  F_sub <- diag(1, r, r) - N_in_Q %**% N_inner_inv %**% t(N_in_Q)
  F_sub <- 0.5 * (F_sub + t(F_sub))

  eig_F <- eigen(F_sub, symmetric = TRUE)
  if (any(eig_F$values <= sqrt(.Machine$double.eps))) {
    return(list(valid = FALSE))
  }

  F_sub_half <- eig_F$vectors %**%
    (t(eig_F$vectors) * sqrt(eig_F$values))
  F_sub_inv_half <- eig_F$vectors %**%
    (t(eig_F$vectors) * (1 / sqrt(eig_F$values)))

  C <- diag(1, r, r) - F_sub_half
  C_inv <- F_sub_inv_half - diag(1, r, r)

  ## Precompute G_{on,k}^{1/2} Q blockwise so the Woodbury-corrected
  ## projection can apply the rank-r update partition by partition.
  GhalfQ <- lapply(1:(K + 1), function(k) {
    rows_k <- ((k - 1) * p_expansions + 1):(k * p_expansions)
    Ghalf_corrected[[k]] %**% Q[rows_k, , drop = FALSE]
  })

  list(valid = TRUE,
       Q = Q,
       C = C,
       C_inv = C_inv,
       GhalfQ = GhalfQ)
}


## Partitioned GLS cross-products

#' Compute Partitioned GLS Cross-Products
#'
#' Computes \eqn{\mathbf{X}_k^{\top}\mathbf{V}^{-1}\mathbf{y}} for each
#' partition, accounting for cross-partition contributions from the
#' correlation structure.
#'
#' @param X List of partition-specific design matrices.
#' @param y List of response vectors by partition.
#' @param VhalfInv_perm \eqn{\mathbf{V}^{-1/2}} permuted to partition
#'   ordering.
#' @param K,p_expansions Integer dimensions.
#' @param order_list Partition-to-data index mapping (unused here but
#'   retained for interface consistency).
#'
#' @return A list of \eqn{K+1} column vectors, each
#'   \eqn{p\_expansions \times 1}.
#'
#' @keywords internal
.compute_Xy_V <- function(X, y, VhalfInv_perm, K, p_expansions,
                          order_list) {

  X_block <- collapse_block_diagonal(X)
  y_vec <- cbind(unlist(y))

  ## V^{-1} y in partition order (two-step multiply avoids forming V^{-1}).
  Vinv_y <- crossprod(VhalfInv_perm, VhalfInv_perm %**% y_vec)

  ## X^T V^{-1} y, split by partition.
  XtVinvy_full <- crossprod(X_block, Vinv_y)

  lapply(1:(K + 1), function(k) {
    rows_k <- ((k - 1) * p_expansions + 1):(k * p_expansions)
    XtVinvy_full[rows_k, , drop = FALSE]
  })
}


## Partitioned GEE score

#' Compute Partitioned GEE Score Vector
#'
#' Computes the GEE score vector \eqn{\mathbf{X}^{\top}\mathrm{diag}
#' (\mathbf{W})\mathbf{V}^{-1}(\mathbf{y} - \boldsymbol{\mu})} split
#' into per-partition pieces of dimension \eqn{p \times 1} each. Uses
#' the identity \eqn{\mathbf{V}^{-1} = \mathbf{I} + \boldsymbol{\Delta}_V}
#' where \eqn{\boldsymbol{\Delta}_V = \mathbf{V}^{-1} - \mathbf{I}} is
#' precomputed and fixed across iterations, avoiding explicit formation
#' of \eqn{\mathbf{V}^{-1}}.
#'
#' @param X List of partition-specific design matrices.
#' @param X_block Full \eqn{N \times P} block-diagonal design.
#' @param y List of response vectors by partition.
#' @param result List of current coefficient column vectors by partition.
#' @param K,p_expansions Integer dimensions.
#' @param family GLM family object.
#' @param W Length-\eqn{N} vector of GLM working weights.
#' @param Delta_V Fixed \eqn{N \times N} matrix
#'   \eqn{\mathbf{V}^{-1} - \mathbf{I}}.
#' @param observation_weights List of observation weights by partition.
#'
#' @return A list of \eqn{K+1} column vectors, each
#'   \eqn{p\_expansions \times 1}, representing the partition-wise
#'   components of the full GEE score.
#'
#' @keywords internal
.compute_score_V_partitioned <- function(X, X_block, y, result,
                                         K, p_expansions,
                                         family, W, Delta_V,
                                         observation_weights) {

  ## Current linear predictor and mean.
  beta_block <- cbind(unlist(result))
  mu <- family$linkinv(X_block %**% beta_block)

  ## Raw residual: y - mu.
  y_block <- cbind(unlist(y))
  residual <- y_block - mu

  ## V^{-1} (y - mu) = (I + Delta_V)(y - mu) = residual + Delta_V residual.
  #  This avoids forming V^{-1} explicitly.
  Vinv_resid <- residual + Delta_V %**% residual

  ## Full score: X^T W V^{-1} (y - mu).
  score_full <- crossprod(X_block, c(W) * Vinv_resid)

  ## Split by partition.
  lapply(1:(K + 1), function(k) {
    rows_k <- ((k - 1) * p_expansions + 1):(k * p_expansions)
    score_full[rows_k, , drop = FALSE]
  })
}


## Path 1a: Gaussian GEE

#' Path 1a: Gaussian Identity + GEE (Closed-Form Full-System Solve)
#'
#' Computes the constrained penalized GLS estimate for Gaussian response
#' with identity link when a correlation structure is present. Because
#' \eqn{\mathbf{V}^{-1/2}} couples all partitions, fitting must operate on
#' the full \eqn{P}-dimensional whitened system rather than partition-wise.
#'
#' The unconstrained GLS estimate is:
#' \deqn{\hat{\boldsymbol{\beta}} = \mathbf{G}(\mathbf{X}^{*\top}
#'   \mathbf{y}^{*}), \quad
#'   \mathbf{G} = (\mathbf{X}^{*\top}\mathbf{X}^{*} +
#'   \boldsymbol{\Lambda})^{-1}}
#' where \eqn{\mathbf{X}^{*} = \mathbf{V}^{-1/2}\mathbf{X}} and
#' \eqn{\mathbf{y}^{*} = \mathbf{V}^{-1/2}\mathbf{y}}. The constrained
#' estimate is then obtained via the \eqn{\mathbf{G}^{1/2}\mathbf{r}^*}
#' trick in the full \eqn{P}-dimensional space.
#'
#' @param X_block Full \eqn{N \times P} unwhitened block-diagonal design.
#' @param X_tilde Whitened design \eqn{\mathbf{V}^{-1/2}\mathbf{X}}.
#' @param y_tilde Whitened response \eqn{\mathbf{V}^{-1/2}\mathbf{y}}.
#' @param VhalfInv_perm \eqn{\mathbf{V}^{-1/2}} permuted to partition
#'   ordering.
#' @param Lambda_block Full \eqn{P \times P} block-diagonal penalty.
#' @param A Constraint matrix (\eqn{P \times R}).
#' @param K,p_expansions Integer dimensions.
#' @param constraint_value_vectors Constraint RHS list encoding
#'   \eqn{\mathbf{A}^{\top}\boldsymbol{\beta} = \mathbf{c}}.
#' @param family GLM family object.
#' @param return_G_getB Logical; return covariance components.
#' @param quadprog Logical; apply QP refinement.
#' @param qp_Amat,qp_bvec,qp_meq QP constraint specification.
#' @param qp_score_function Score function for QP step.
#' @param order_list,observation_weights Standard partition arguments.
#' @param ... Passed to score function.
#'
#' @return If \code{return_G_getB = TRUE}: list with \code{B},
#'   \code{G_list}, and \code{qp_info}. Otherwise: list of \code{B} and
#'   \code{qp_info}.
#'
#' @keywords internal
.get_B_gee_gaussian <- function(X_block, X_tilde, y_tilde, VhalfInv_perm,
                                Lambda_block, A, K, p_expansions,
                                constraint_value_vectors, family,
                                return_G_getB,
                                quadprog, qp_Amat, qp_bvec, qp_meq,
                                qp_score_function,
                                order_list, observation_weights, ...) {

  ## Observation weights enter the Gaussian GEE solve in the whitened system,
  #  exactly as they do in the iterative GLM GEE paths.
  obs_wt <- c(unlist(observation_weights))
  if (length(obs_wt) == 0L) {
    obs_wt <- rep(1, nrow(X_tilde))
  } else if (length(obs_wt) == 1L) {
    obs_wt <- rep(obs_wt, nrow(X_tilde))
  }

  ## Full P x P Gram matrix from the weighted whitened design:
  #  X_tilde^T D X_tilde = X^T V^{-1/2} D V^{-1/2} X.
  Gram_full <- crossprod(X_tilde, obs_wt * X_tilde)

  ## G = (X_tilde^T X_tilde + Lambda)^{-1}, a dense P x P matrix.
  G_full_inv <- Gram_full + Lambda_block
  G_full <- invert(G_full_inv)

  ## G^{1/2} via eigendecomposition of the full P x P G.
  eig_G <- eigen(G_full, symmetric = TRUE)
  vals_G <- eig_G$values
  vals_G[vals_G <= 0] <- 0
  G_full_half <- eig_G$vectors %**%
    (t(eig_G$vectors) * sqrt(vals_G))

  ## G^{-1/2} for per-partition extraction.
  inv_sqrt_vals_G <- ifelse(sqrt(pmax(eig_G$values, 0)) > 0,
                            1 / sqrt(pmax(eig_G$values, 0)), 0)
  G_full_half_inv <- eig_G$vectors %**%
    (t(eig_G$vectors) * inv_sqrt_vals_G)

  ## X_tilde^T D y_tilde: the P x 1 sufficient statistic for beta.
  Xy_tilde <- crossprod(X_tilde, obs_wt * y_tilde)

  ## Lagrangian projection in the full P-dimensional space.
  if (K == 0 | any(all.equal(unique(A), 0) == TRUE)) {
    A_proj <- cbind(rep(0, p_expansions * (K + 1)))
  } else {
    A_proj <- A
  }

  ## G^{1/2} r* trick.
  y_star <- G_full_half %**% Xy_tilde
  X_star <- G_full_half %**% A_proj

  comp_stab_sc <- 1 / sqrt(K + 1)
  resids_star <- .lm.fit(
    X_star * comp_stab_sc,
    y_star * comp_stab_sc
  )$residuals / comp_stab_sc

  ## Nonzero constraint values: add particular solution b_0.
  if (length(constraint_value_vectors) > 0) {
    if (any(unlist(constraint_value_vectors) != 0)) {
      preds_star <- X_star %**%
        (invert(crossprod(X_star) * comp_stab_sc) %**%
           crossprod(A_proj,
                     Reduce("rbind",
                            constraint_value_vectors) *
                       comp_stab_sc))
      resids_star <- resids_star + c(preds_star)
    }
  }

  ## beta_constrained = G^{1/2} r*.
  beta_block <- G_full_half %**% cbind(resids_star)

  ## Unpack into per-partition column vectors.
  result <- lapply(1:(K + 1), function(k) {
    cbind(beta_block[(k - 1) * p_expansions + 1:p_expansions])
  })

  qp_info <- NULL

  ## Optional QP refinement for inequality constraints.
  if (quadprog) {
    if (K == 0 | any(all.equal(unique(A), 0) == TRUE)) {
      A_qp <- cbind(rep(0, p_expansions * (K + 1)))
    } else {
      A_qp <- A
    }
    if (!is.null(qp_Amat)) {
      qp_Amat_c <- cbind(A_qp, qp_Amat)
    } else {
      qp_Amat_c <- A_qp
    }
    if (length(constraint_value_vectors) > 0) {
      constr_rhs <- Reduce('rbind', constraint_value_vectors)
      if (nrow(constr_rhs) < ncol(A_qp)) {
        constr_rhs <- c(rep(0, ncol(A_qp) - nrow(constr_rhs)),
                        constr_rhs)
      }
    } else {
      constr_rhs <- rep(0, ncol(A_qp))
    }
    if (!is.null(qp_bvec)) {
      qp_bvec_c <- c(constr_rhs, qp_bvec)
    } else {
      qp_bvec_c <- constr_rhs
    }
    if (!is.null(qp_meq)) {
      qp_meq_c <- ncol(A_qp) + qp_meq
    } else {
      qp_meq_c <- ncol(A_qp)
    }

    ## Information matrix for the QP: use full whitened Gram.
    info <- Gram_full + Lambda_block
    sc <- sqrt(mean(abs(info)))

    qp_score <- qp_score_function(
      X_tilde, y_tilde,
      VhalfInv_perm %**% cbind(family$linkinv(
        c(X_block %**% beta_block))),
      unlist(order_list), 1, VhalfInv_perm,
      unlist(observation_weights), ...
    )
    qp_result <- .solve_qp_stable(
      Dmat = info / sc,
      dvec = (qp_score -
                Lambda_block %**% beta_block +
                info %**% beta_block) / sc,
      Amat = qp_Amat_c,
      bvec = qp_bvec_c,
      meq = qp_meq_c
    )

    if (!any(inherits(qp_result, 'try-error'))) {
      active_ineq <- if (ncol(qp_Amat_c) > ncol(A_qp)) {
        which(abs(qp_result$Lagrangian[-(1:ncol(A_qp))]) >
                sqrt(.Machine$double.eps))
      } else {
        integer(0)
      }
      beta_block <- cbind(qp_result$solution)
      result <- lapply(1:(K + 1), function(k) {
        cbind(beta_block[(k - 1) * p_expansions + 1:p_expansions])
      })
      qp_info <- list(
        lagrangian = qp_result$Lagrangian,
        Amat_active = qp_Amat_c[,
                                c(1:ncol(A_qp),
                                  ncol(A_qp) + which(abs(qp_result$Lagrangian[-(1:ncol(A_qp))]) >
                                                       sqrt(.Machine$double.eps))),
                                drop = FALSE],
        active_ineq = active_ineq,
        converged = TRUE,
        method = "dense_qp_gee_gaussian"
      )
    }
  }

  ## Return with per-partition G components.
  if (return_G_getB) {
    G_diag <- .extract_G_diagonal(G_full, p_expansions, K)
    Ghalf_diag <- .extract_G_diagonal(G_full_half, p_expansions, K)
    GhalfInv_diag <- .extract_G_diagonal(G_full_half_inv, p_expansions, K)
    G_diag_from_half <- lapply(1:(K + 1), function(k) {
      tcrossprod(Ghalf_diag[[k]])
    })
    return(list(
      B = result,
      G_list = list(G = G_diag_from_half,
                    Ghalf = Ghalf_diag,
                    GhalfInv = GhalfInv_diag),
      qp_info = qp_info
    ))
  } else {
    return(list(B = result, qp_info = qp_info))
  }
}


#' Path 1a-Woodbury: Gaussian GEE with Woodbury Acceleration
#'
#' Gaussian correlated solve using the Woodbury factorization from the
#' supplement. This replaces the dense GEE path when
#' \eqn{\mathbf{G}_{\mathrm{off}}^{-1} = \mathbf{E}\mathbf{J}\mathbf{E}^{\top}}
#' has low rank and computes the constrained estimate through the same
#' transformed OLS projection as the independent case, with
#' \eqn{\mathbf{G}^{1/2} =
#' \mathbf{G}_{\mathrm{on}}^{1/2}\mathbf{F}^{1/2}}.
#'
#' This function is standalone for the Gaussian Woodbury path: it takes the
#' decomposition from \code{.woodbury_decompose_V()}, the square-root state
#' from \code{.woodbury_halfsqrt_components()}, performs the constrained
#' solve, and returns the same object structure as the dense Gaussian GEE
#' solver.
#'
#' Cost is \eqn{O(Kp^3 + Pr^2)} compared to \eqn{O(P^3)} for the dense
#' Path 1a.
#'
#' @param X,y Lists of partition-specific design matrices and responses.
#' @param K,p_expansions Integer dimensions.
#' @param VhalfInv_perm \eqn{\mathbf{V}^{-1/2}} permuted to partition
#'   ordering.
#' @param order_list Partition-to-data index mapping.
#' @param A Constraint matrix (\eqn{P \times R}).
#' @param R_constraints Number of columns of A.
#' @param constraint_value_vectors Constraint RHS list.
#' @param family GLM family object.
#' @param return_G_getB Logical; return covariance components.
#' @param quadprog Logical; apply QP refinement.
#' @param qp_Amat,qp_bvec,qp_meq QP constraint specification.
#' @param qp_score_function Score function for QP step.
#' @param observation_weights Observation weights.
#' @param wb_decomp Output of \code{.woodbury_decompose_V}.
#' @param wb_sqrt Output of \code{.woodbury_halfsqrt_components}.
#' @param Lambda_block Full block-diagonal penalty matrix.
#' @param parallel_aga,parallel_matmult Logical flags.
#' @param cl,chunk_size,num_chunks,rem_chunks Parallel parameters.
#' @param qp_global Logical; TRUE when inequality constraints couple
#'   partitions and therefore require the global equality bridge inside
#'   the active-set refinement.
#' @param tol Numeric tolerance used by the Woodbury active-set refinement.
#' @param ... Passed to sub-functions.
#'
#' @return Same structure as \code{.get_B_gee_gaussian}.
#'
#' @keywords internal
.get_B_gee_woodbury <- function(X, y, K, p_expansions,
                                VhalfInv_perm, order_list,
                                A, R_constraints,
                                constraint_value_vectors, family,
                                return_G_getB,
                                quadprog, qp_Amat, qp_bvec, qp_meq,
                                qp_score_function,
                                observation_weights,
                                wb_decomp, wb_sqrt,
                                Lambda_block,
                                parallel_aga, parallel_matmult,
                                cl, chunk_size, num_chunks, rem_chunks,
                                qp_global, tol, ...) {

  dots <- list(...)
  include_warnings <- TRUE
  if (!is.null(dots$include_warnings)) {
    include_warnings <- isTRUE(dots$include_warnings)
  }

  P <- p_expansions * (K + 1)

  ## Work with the Woodbury factor L = G_on^{1/2} F^{1/2}.
  X_block <- collapse_block_diagonal(X)
  y_block <- cbind(unlist(y))

  ## Observation weights for the whitened system.
  obs_wt <- c(unlist(observation_weights))
  if (length(obs_wt) == 0L) {
    obs_wt <- rep(1, nrow(VhalfInv_perm))
  } else if (length(obs_wt) == 1L) {
    obs_wt <- rep(obs_wt, nrow(VhalfInv_perm))
  }

  Xy_V_list <- .compute_Xy_V(X, y, VhalfInv_perm, K, p_expansions,
                             order_list)
  Xy_full <- cbind(unlist(Xy_V_list))

  ## y* = L^T X^T V^{-1} y.
  GhalfXy_V <- .woodbury_transform_constraint_matrix_transpose(
    Ghalf_corrected = wb_decomp$Ghalf_corrected,
    A = Xy_full,
    K = K,
    p_expansions = p_expansions,
    wb_sqrt = wb_sqrt,
    parallel_aga = parallel_aga,
    cl = cl,
    chunk_size = chunk_size,
    num_chunks = num_chunks,
    rem_chunks = rem_chunks
  )

  ## Constrained projection via Lagrangian.
  result <- .lagrangian_project_woodbury(
    GhalfXy_V, wb_decomp$Ghalf_corrected,
    A, K, p_expansions, R_constraints,
    constraint_value_vectors, family,
    wb_sqrt,
    parallel_aga, parallel_matmult,
    cl, chunk_size, num_chunks, rem_chunks)

  qp_info <- NULL

  ## QP refinement: try the Woodbury active-set first; dense paths
  #  remain fallbacks.
  if (quadprog) {

    as_out <- try(.active_set_refine_woodbury(
      result = result,
      K = K,
      p_expansions = p_expansions,
      A = A,
      R_constraints = R_constraints,
      constraint_value_vectors = constraint_value_vectors,
      family = family,
      qp_Amat = qp_Amat,
      qp_bvec = qp_bvec,
      qp_meq = qp_meq,
      rhs_list = Xy_V_list,
      Ghalf_corrected = wb_decomp$Ghalf_corrected,
      wb_sqrt = wb_sqrt,
      parallel_aga = parallel_aga,
      parallel_matmult = parallel_matmult,
      cl = cl,
      chunk_size = chunk_size,
      num_chunks = num_chunks,
      rem_chunks = rem_chunks,
      tol = max(tol, 100 * sqrt(.Machine$double.eps))
    ), silent = TRUE)

    if (!inherits(as_out, "try-error") && isTRUE(as_out$converged)) {
      result <- as_out$result
      qp_info <- as_out$qp_info
    } else {
      X_tilde <- VhalfInv_perm %**% X_block
      y_tilde <- VhalfInv_perm %**% y_block
      Gram_full <- crossprod(X_tilde, obs_wt * X_tilde)
      G_full_inv <- Gram_full + Lambda_block
      G_full <- invert(G_full_inv)
      eig_G <- eigen(G_full, symmetric = TRUE)
      vals_G <- pmax(eig_G$values, 0)
      G_full_half <- eig_G$vectors %**%
        (t(eig_G$vectors) * sqrt(vals_G))
      Xy_tilde <- crossprod(X_tilde, obs_wt * y_tilde)

      as_full <- try(.active_set_refine_full(
        result = result,
        K = K,
        p_expansions = p_expansions,
        A = A,
        constraint_value_vectors = constraint_value_vectors,
        family = family,
        qp_Amat = qp_Amat,
        qp_bvec = qp_bvec,
        qp_meq = qp_meq,
        rhs_full = Xy_tilde,
        Ghalf_full = G_full_half,
        tol = max(tol, 100 * sqrt(.Machine$double.eps)),
        method = "active_set_woodbury"
      ), silent = TRUE)

      if (!inherits(as_full, "try-error") && isTRUE(as_full$converged)) {
        result <- as_full$result
        qp_info <- as_full$qp_info
      } else {
      ## Fall back to dense solve.QP.
      beta_block <- cbind(unlist(result))

      if (K == 0 | any(all.equal(unique(A), 0) == TRUE)) {
        A_qp <- cbind(rep(0, p_expansions * (K + 1)))
      } else {
        A_qp <- A
      }
      if (!is.null(qp_Amat)) {
        qp_Amat_c <- cbind(A_qp, qp_Amat)
      } else {
        qp_Amat_c <- A_qp
      }
      constr_rhs <- if (length(constraint_value_vectors) > 0) {
        cr <- Reduce('rbind', constraint_value_vectors)
        if (nrow(cr) < ncol(A_qp)) c(rep(0, ncol(A_qp) - nrow(cr)), cr)
        else cr
      } else {
        rep(0, ncol(A_qp))
      }
      qp_bvec_c <- if (!is.null(qp_bvec)) c(constr_rhs, qp_bvec) else constr_rhs
      qp_meq_c <- ncol(A_qp) + if (!is.null(qp_meq)) qp_meq else 0

      info <- Gram_full + Lambda_block
      sc <- sqrt(mean(abs(info)))

      qp_score <- qp_score_function(
        X_tilde, y_tilde,
        VhalfInv_perm %**% cbind(family$linkinv(
          c(X_block %**% beta_block))),
        unlist(order_list), 1, VhalfInv_perm,
        unlist(observation_weights), ...
      )
      qp_result <- .solve_qp_stable(
        Dmat = info / sc,
        dvec = (qp_score -
                  Lambda_block %**% beta_block +
                  info %**% beta_block) / sc,
        Amat = qp_Amat_c,
        bvec = qp_bvec_c,
        meq = qp_meq_c
      )

      if (!any(inherits(qp_result, 'try-error'))) {
        active_ineq <- if (ncol(qp_Amat_c) > ncol(A_qp)) {
          which(abs(qp_result$Lagrangian[-(1:ncol(A_qp))]) >
                  sqrt(.Machine$double.eps))
        } else {
          integer(0)
        }
        beta_block <- cbind(qp_result$solution)
        result <- lapply(1:(K + 1), function(k) {
          cbind(beta_block[(k - 1) * p_expansions + 1:p_expansions])
        })
        qp_info <- list(
          lagrangian = qp_result$Lagrangian,
          Amat_active = qp_Amat_c[,
                                  c(1:ncol(A_qp),
                                    ncol(A_qp) + which(abs(qp_result$Lagrangian[-(1:ncol(A_qp))]) >
                                                         sqrt(.Machine$double.eps))),
                                  drop = FALSE],
          active_ineq = active_ineq,
          converged = TRUE,
          method = "dense_qp_gee_gaussian",
          fallback_reason = "woodbury_active_set_failed",
          woodbury_rank = wb_decomp$r,
          woodbury_P = wb_decomp$P,
          woodbury_rank_threshold = wb_decomp$rank_threshold,
          woodbury_nnz_fraction = wb_decomp$nnz_fraction
        )
      }
      }
    }
  }

  ## Return with per-partition G components.
  if (return_G_getB) {
    if (exists("G_full_half", inherits = FALSE)) {
      Ghalf_diag <- .extract_G_diagonal(G_full_half, p_expansions, K)
      inv_sqrt_vals_G2 <- ifelse(sqrt(pmax(eig_G$values, 0)) > 0,
                                 1 / sqrt(pmax(eig_G$values, 0)), 0)
      G_full_half_inv2 <- eig_G$vectors %**%
        (t(eig_G$vectors) * inv_sqrt_vals_G2)
      GhalfInv_diag <- .extract_G_diagonal(G_full_half_inv2, p_expansions, K)
      G_diag_from_half <- lapply(1:(K + 1), function(k) {
        tcrossprod(Ghalf_diag[[k]])
      })
    } else {
      Ghalf_diag <- wb_decomp$Ghalf_corrected
      GhalfInv_diag <- wb_decomp$GhalfInv_corrected
      G_diag_from_half <- lapply(Ghalf_diag, tcrossprod)
    }
    return(list(
      B = result,
      G_list = list(G = G_diag_from_half,
                    Ghalf = Ghalf_diag,
                    GhalfInv = GhalfInv_diag),
      qp_info = qp_info
    ))
  } else {
    return(list(B = result, qp_info = qp_info))
  }
}


## Path 1b: non-Gaussian GEE

#' Path 1b: Non-Gaussian GEE (Damped SQP with Full Whitened Design)
#'
#' Estimates constrained coefficients for non-Gaussian GLMs with correlation
#' structures using damped Sequential Quadratic Programming (SQP) in the
#' full \eqn{P}-dimensional whitened space. The whitened design
#' \eqn{\mathbf{X}^{*} = \mathbf{V}^{-1/2}\mathbf{X}} is required
#' because \eqn{\mathbf{V}^{-1/2}} couples all partitions.
#'
#' @inheritParams .get_B_gee_gaussian
#' @param y_block Unwhitened response vector \eqn{\mathbf{y}}.
#' @param iterate Logical; if \code{TRUE}, iterate until convergence.
#' @param tol Convergence tolerance.
#' @param glm_weight_function Function computing GLM working weights.
#' @param schur_correction_function Function computing Schur corrections.
#' @param need_dispersion_for_estimation Logical.
#' @param dispersion_function Dispersion estimation function.
#' @param VhalfInv Inverse square root correlation matrix.
#' @param ... Passed to weight, correction, dispersion, and score functions.
#'
#' @return Same structure as \code{.get_B_gee_gaussian}.
#'
#' @keywords internal
.get_B_gee_glm <- function(X_block, X_tilde, y_block, y_tilde,
                           VhalfInv_perm, Lambda_block, A, K,
                           p_expansions, constraint_value_vectors,
                           family, return_G_getB, iterate, tol,
                           qp_Amat, qp_bvec, qp_meq,
                           qp_score_function,
                           order_list, observation_weights,
                           glm_weight_function, schur_correction_function,
                           need_dispersion_for_estimation,
                           dispersion_function, VhalfInv, ...) {

  beta_block <- cbind(rep(0, ncol(X_block)))

  ## SQP iteration control.
  damp_cnt <- 0
  master_cnt <- 0
  err <- Inf

  ## Original-scale linear predictor (X_block is unwhitened).
  XB <- X_block %**% beta_block

  ## Combine smoothness equality and user inequality constraints.
  if (K == 0 | any(all.equal(unique(A), 0) == TRUE)) {
    A <- cbind(rep(0, p_expansions * (K + 1)))
  }
  if (!is.null(qp_Amat)) {
    qp_Amat <- cbind(A, qp_Amat)
  } else {
    qp_Amat <- A
  }
  if (length(constraint_value_vectors) > 0) {
    constr_A <- Reduce('rbind', constraint_value_vectors)
    if (nrow(constr_A) < ncol(A)) {
      constr_A <- c(rep(0, ncol(A) - nrow(constr_A)), constr_A)
    }
  } else {
    constr_A <- rep(0, ncol(A))
  }
  if (!is.null(qp_bvec)) {
    qp_bvec <- c(constr_A, qp_bvec)
  } else {
    qp_bvec <- constr_A
  }
  if (!is.null(qp_meq)) {
    qp_meq <- ncol(A) + qp_meq
  } else {
    qp_meq <- ncol(A)
  }

  qp_info <- NULL

  ## Damped SQP loop.
  while (err > tol & damp_cnt < 10 & master_cnt < 100) {
    master_cnt <- master_cnt + 1
    damp <- 2^(-(damp_cnt))

    if (need_dispersion_for_estimation) {
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

    W <- c(glm_weight_function(family$linkinv(XB),
                               y_block,
                               unlist(order_list),
                               family, dispersion_temp,
                               unlist(observation_weights), ...))
    W <- .solver_stabilize_working_weights(W)

    ## Schur correction using the original-scale design.
    result <- lapply(1:(K + 1), function(k) {
      cbind(beta_block[(k - 1) * p_expansions + 1:p_expansions])
    })
    schur_correction <-
      schur_correction_function(
        list(X_block), list(y_block),
        list(cbind(unlist(result))),
        dispersion_temp, list(unlist(order_list)),
        0, family, unlist(observation_weights), ...
      )
    if (!(any(unlist(schur_correction) != 0))) {
      schur_correction_collapsed <- 0
    } else {
      schur_correction_collapsed <-
        collapse_block_diagonal(schur_correction)
    }

    info <- crossprod(X_tilde, W * X_tilde) +
      Lambda_block +
      schur_correction_collapsed

    sc <- sqrt(mean(abs(info)))

    if (master_cnt == 1) {
      ## First iteration: constrained Newton step.
      infoinv_block <- sc * invert(sc * info)
      U <- (diag(1, nrow(info)) -
              tcrossprod(infoinv_block %**%
                           A %**%
                           invert(crossprod(A,
                                            infoinv_block %**% A)),
                         A))
      mu_vec <- cbind(family$linkinv(c(XB)))
      beta_new <- c(
        U %**% infoinv_block %**%
          crossprod(X_tilde,
                    y_tilde - VhalfInv_perm %**% mu_vec)
      )
      infoinv_block <- NULL
    } else {
      ## Subsequent iterations: full QP solve.
      qp_score <- qp_score_function(
        X_tilde, y_tilde,
        VhalfInv_perm %**% cbind(family$linkinv(XB)),
        unlist(order_list), dispersion_temp,
        VhalfInv_perm,
        unlist(observation_weights), ...
      )
      qp_result <- .solve_qp_stable(
        Dmat = info / sc,
        dvec = (qp_score -
                  Lambda_block %**% beta_block +
                  info %**% beta_block) / sc,
        Amat = qp_Amat,
        bvec = qp_bvec,
        meq = qp_meq
      )
      if (any(inherits(qp_result, 'try-error'))) {
        beta_new <- 0 * beta_block
      } else {
        active_ineq <- if (ncol(qp_Amat) > ncol(A)) {
          which(abs(qp_result$Lagrangian[-(1:ncol(A))]) >
                  sqrt(.Machine$double.eps))
        } else {
          integer(0)
        }
        beta_new <- qp_result$solution
        qp_info <- list(
          lagrangian = qp_result$Lagrangian,
          Amat_active = qp_Amat[,
                                c(1:ncol(A),
                                  ncol(A) + which(abs(qp_result$Lagrangian[-(1:ncol(A))]) >
                                                    sqrt(.Machine$double.eps))),
                                drop = FALSE],
          active_ineq = active_ineq,
          converged = TRUE,
          method = "dense_qp_gee_glm"
        )
      }
    }

    ## Non-iterative mode: accept after the first Newton step.
    if (!iterate & master_cnt > 2) {
      beta_block <- beta_new
      damp_cnt <- 11
      master_cnt <- 101
      err <- tol - 1
    } else {
      beta_new <- (1 - damp) * beta_block + damp * beta_new
      XB <- X_block %**% cbind(beta_new)

      if (need_dispersion_for_estimation) {
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
      W <- c(glm_weight_function(family$linkinv(XB),
                                 y_block,
                                 unlist(order_list),
                                 family, dispersion_temp,
                                 unlist(observation_weights), ...))
      W <- .solver_stabilize_working_weights(W)

      mu_new <- family$linkinv(c(XB))
      if (!is.null(family$custom_dev.resids)) {
        raw <- family$custom_dev.resids(
          y_block, mu_new, unlist(order_list),
          family, unlist(observation_weights), ...
        )
        W_safe <- pmax(W, sqrt(.Machine$double.eps))
        err_new <- mean((
          VhalfInv_perm %**%
            cbind(sign(raw) * sqrt(abs(raw)) / sqrt(c(W_safe)))
        )^2)
      } else if (is.null(family$dev.resids)) {
        err_new <- mean((unlist(observation_weights) *
                           (y_block - cbind(mu_new)))^2)
      } else {
        err_new <-
          mean(unlist(observation_weights) *
                 family$dev.resids(y_block, cbind(mu_new),
                                   wt = 1 / W))
      }

      if (is.null(err_new) | is.na(err_new) | !is.finite(err_new)) {
        damp_cnt <- damp_cnt + 1
      } else if (err_new <= err) {
        prev_err <- err
        err <- err_new
        abs_diff <- max(abs(beta_new - beta_block))
        beta_block <- beta_new
        damp_cnt <- 0
        if ((abs_diff < tol) & (prev_err - err < tol) &
            (master_cnt > 10)) {
          damp_cnt <- 11
          master_cnt <- 101
          err <- tol - 1
        }
      } else {
        damp_cnt <- damp_cnt + 1
      }
    }

    result <- lapply(1:(K + 1), function(k) {
      cbind(beta_block[(k - 1) * p_expansions + 1:p_expansions])
    })
  }

  ## Return with G components if requested.
  if (return_G_getB) {
    XB <- X_block %**% cbind(unlist(result))
    if (need_dispersion_for_estimation) {
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
    W <- c(glm_weight_function(family$linkinv(XB),
                               y_block,
                               unlist(order_list),
                               family, dispersion_temp,
                               unlist(observation_weights), ...))
    W <- .solver_stabilize_working_weights(W)

    schur_correction <-
      schur_correction_function(
        list(X_block), list(y_block),
        list(cbind(unlist(result))),
        dispersion_temp, list(unlist(order_list)),
        0, family, unlist(observation_weights), ...
      )
    if (!(any(unlist(schur_correction) != 0))) {
      schur_correction_collapsed <- 0
    } else {
      schur_correction_collapsed <-
        collapse_block_diagonal(schur_correction)
    }

    info <- crossprod(X_tilde, W * X_tilde) +
      Lambda_block +
      schur_correction_collapsed

    G <- lapply(1:(K + 1), function(k) NA)
    Ghalf <- lapply(1:(K + 1), function(k) NA)
    GhalfInv <- lapply(1:(K + 1), function(k) NA)
    for (k in 1:(K + 1)) {
      info_kk <- info[(k - 1) * p_expansions + 1:p_expansions,
                      (k - 1) * p_expansions + 1:p_expansions]
      eig <- eigen(info_kk, symmetric = TRUE)
      vals <- eig$values
      vals_safe <- pmax(vals, 0)
      vals_safe[vals <= sqrt(.Machine$double.eps)] <- 0
      sqrt_vals <- sqrt(vals_safe)
      Ghalf[[k]] <- eig$vectors %**%
        (t(eig$vectors) * sqrt_vals)
      inv_sqvals <- ifelse(sqrt_vals > 0, 1 / sqrt_vals, 0)
      GhalfInv[[k]] <- eig$vectors %**%
        (t(eig$vectors) * inv_sqvals)
      G[[k]] <- tcrossprod(Ghalf[[k]])
    }

    return(list(
      B = result,
      G_list = list(Ghalf = Ghalf,
                    GhalfInv = GhalfInv,
                    G = G),
      qp_info = qp_info
    ))
  } else {
    return(list(B = result, qp_info = qp_info))
  }
}


#' Path 1b-Woodbury: Non-Gaussian GEE with Woodbury Acceleration
#'
#' Non-Gaussian correlated solve using the Woodbury factorization. At each
#' damped Newton step, the current weighted information matrix is split into
#' a block-diagonal part defining \eqn{\mathbf{G}_{\mathrm{on}}} and a
#' low-rank cross-partition part.  The constrained step is then taken
#' through \code{.lagrangian_project_woodbury()}.
#'
#' @inheritParams .get_B_gee_glm
#' @param X,y Lists of partition-specific design matrices and responses.
#' @param Lambda,Lambda_block Shared and full penalty matrices.
#' @param unique_penalty_per_partition Logical.
#' @param L_partition_list Partition-specific penalty matrices.
#' @param parallel_eigen,parallel_aga,parallel_matmult Logical flags.
#' @param cl,chunk_size,num_chunks,rem_chunks Parallel parameters.
#' @param ... Passed to sub-functions.
#'
#' @return Same structure as \code{.get_B_gee_glm}.
#'
#' @keywords internal
.get_B_gee_glm_woodbury <- function(X, y, K, p_expansions,
                                    VhalfInv_perm, order_list,
                                    A, R_constraints,
                                    constraint_value_vectors, family,
                                    return_G_getB, iterate, tol,
                                    quadprog,
                                    qp_Amat, qp_bvec, qp_meq,
                                    qp_score_function,
                                    observation_weights,
                                    Lambda, Lambda_block,
                                    unique_penalty_per_partition,
                                    L_partition_list,
                                    wb_decomp_init, wb_sqrt_init,
                                    glm_weight_function,
                                    schur_correction_function,
                                    need_dispersion_for_estimation,
                                    dispersion_function, VhalfInv,
                                    parallel_eigen, parallel_aga,
                                    parallel_matmult,
                                    cl, chunk_size, num_chunks,
                                    rem_chunks, ...) {

  dots <- list(...)
  include_warnings <- TRUE
  if (!is.null(dots$include_warnings)) {
    include_warnings <- isTRUE(dots$include_warnings)
  }

  P <- p_expansions * (K + 1)

  ## Precomputation (done once, fixed across Newton iterations).
  X_block <- collapse_block_diagonal(X)
  y_block <- cbind(unlist(y))
  N_obs <- nrow(X_block)

  ## Delta_V = V^{-1} - I: fixed across iterations.
  Delta_V <- crossprod(VhalfInv_perm) - diag(N_obs)

  ## DV_X = (V^{-1} - I) X: fixed across iterations.
  DV_X <- Delta_V %**% X_block

  ## Partition sizes for splitting W.
  n_per_partition <- sapply(X, nrow)

  ## Initialize beta at zero.
  beta_block <- cbind(rep(0, P))

  ## Damped Newton control.
  damp_cnt <- 0
  master_cnt <- 0
  err <- Inf
  XB <- X_block %**% beta_block

  ## Track fallback state.
  fell_back_to_dense <- FALSE

  qp_info <- NULL

  ## Newton loop.
  while (err > tol & damp_cnt < 10 & master_cnt < 100) {
    master_cnt <- master_cnt + 1
    damp <- 2^(-(damp_cnt))

    if (need_dispersion_for_estimation) {
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

    W <- c(glm_weight_function(family$linkinv(XB),
                               y_block,
                               unlist(order_list),
                               family, dispersion_temp,
                               unlist(observation_weights), ...))
    W <- .solver_stabilize_working_weights(W)

    result <- lapply(1:(K + 1), function(k) {
      cbind(beta_block[(k - 1) * p_expansions + 1:p_expansions])
    })
    schur_corrections <- schur_correction_function(
      X, y, result, dispersion_temp, order_list, K, family,
      observation_weights, ...
    )

    W_split <- split(W, rep(1:(K + 1), n_per_partition))
    Xw <- lapply(1:(K + 1), function(k) {
      if (nrow(X[[k]]) == 0) return(X[[k]])
      X[[k]] * sqrt(W_split[[k]])
    })
    X_gram_weighted <- compute_gram_block_diagonal(
      Xw, parallel_matmult, cl, chunk_size, num_chunks, rem_chunks)

    wb <- .woodbury_redecompose_weighted(
      Delta_V, X_block, DV_X, W,
      X, K, p_expansions,
      X_gram_weighted,
      Lambda, schur_corrections,
      unique_penalty_per_partition, L_partition_list,
      parallel_eigen, cl,
      chunk_size, num_chunks, rem_chunks,
      rank_threshold_fraction = 1/2, family)

    if (!wb$use_woodbury) {
      fell_back_to_dense <- TRUE
      break
    }

    wb_sqrt <- .woodbury_halfsqrt_components(
    wb$Ghalf_corrected, wb$E, wb$J,
      wb$r, K, p_expansions)

    if (!wb_sqrt$valid) {
      fell_back_to_dense <- TRUE
      break
    }

    score_V <- .compute_score_V_partitioned(
      X, X_block, y, result, K, p_expansions,
      family, W, Delta_V, observation_weights)

    GhalfXy_V <- .woodbury_transform_constraint_matrix_transpose(
      Ghalf_corrected = wb$Ghalf_corrected,
      A = cbind(unlist(score_V)),
      K = K,
      p_expansions = p_expansions,
      wb_sqrt = wb_sqrt,
      parallel_aga = parallel_aga,
      cl = cl,
      chunk_size = chunk_size,
      num_chunks = num_chunks,
      rem_chunks = rem_chunks
    )

    result_new <- .lagrangian_project_woodbury(
      GhalfXy_V, wb$Ghalf_corrected,
      A, K, p_expansions, R_constraints,
      constraint_value_vectors, family,
      wb_sqrt,
      parallel_aga, parallel_matmult,
      cl, chunk_size, num_chunks, rem_chunks)

    beta_new <- cbind(unlist(result_new))

    if (!iterate & master_cnt > 1) {
      beta_block <- beta_new
      damp_cnt <- 11
      master_cnt <- 101
      err <- tol - 1
    } else {
      beta_new <- (1 - damp) * beta_block + damp * beta_new
      XB <- X_block %**% beta_new

      mu_new <- family$linkinv(c(XB))
      if (!is.null(family$custom_dev.resids) &
          is.null(family$dev.resids)) {
        err_new <- mean(
          family$custom_dev.resids(y_block,
                                   mu_new,
                                   unlist(order_list),
                                   family,
                                   unlist(observation_weights),
                                   ...))
      } else if (is.null(family$dev.resids)) {
        err_new <- mean(unlist(observation_weights) *
                          (y_block - mu_new)^2)
      } else {
        err_new <-
          mean(unlist(observation_weights) *
                 family$dev.resids(y_block,
                                   cbind(mu_new),
                                   wt = 1))
      }

      if (is.null(err_new) | is.na(err_new) | !is.finite(err_new)) {
        damp_cnt <- damp_cnt + 1
      } else if (err_new <= err) {
        prev_err <- err
        err <- err_new
        abs_diff <- max(abs(beta_new - beta_block))
        beta_block <- beta_new
        damp_cnt <- 0

        if ((abs_diff < tol) &
            (prev_err - err < tol) &
            (master_cnt > 10)) {
          damp_cnt <- 11
          master_cnt <- 101
          err <- tol - 1
        }

      } else {
        damp_cnt <- damp_cnt + 1
      }
    }

    result <- lapply(1:(K + 1), function(k) {
      cbind(beta_block[(k - 1) * p_expansions + 1:p_expansions])
    })
  }

  ## Fallback to dense Path 1b.
  if (fell_back_to_dense) {
    X_tilde <- VhalfInv_perm %**% X_block
    y_tilde <- VhalfInv_perm %**% y_block

    return(.get_B_gee_glm(
      X_block, X_tilde, y_block, y_tilde,
      VhalfInv_perm, Lambda_block, A, K,
      p_expansions, constraint_value_vectors,
      family, return_G_getB, iterate, tol,
      qp_Amat, qp_bvec, qp_meq,
      qp_score_function,
      order_list, observation_weights,
      glm_weight_function, schur_correction_function,
      need_dispersion_for_estimation,
      dispersion_function, VhalfInv, ...))
  }

  ## QP refinement after the Woodbury Newton-Raphson loop.
  #  Build the final quadratic exactly and use active-set equality
  #  re-solves before falling back to dense SQP.
  if (quadprog) {

    beta_block <- cbind(unlist(result))
    XB_final <- X_block %**% beta_block

    if (need_dispersion_for_estimation) {
      dispersion_final <- dispersion_function(
        mu = family$linkinv(XB_final),
        y = y_block,
        order_indices = unlist(order_list),
        family = family,
        observation_weights = unlist(observation_weights),
        VhalfInv = VhalfInv,
        ...
      )
    } else {
      dispersion_final <- 1
    }

    W_final <- c(glm_weight_function(
      family$linkinv(XB_final),
      y_block,
      unlist(order_list),
      family,
      dispersion_final,
      unlist(observation_weights),
      ...
    ))
    W_final <- .solver_stabilize_working_weights(W_final)

    ## Build the final Woodbury quadratic.
    X_tilde <- VhalfInv_perm %**% X_block
    y_tilde <- VhalfInv_perm %**% y_block

    schur_corrections_final <- schur_correction_function(
      X, y, result, dispersion_final, order_list, K, family,
      observation_weights, ...
    )
    schur_coll <- if (any(unlist(schur_corrections_final) != 0)) {
      collapse_block_diagonal(schur_corrections_final)
    } else {
      0
    }

    ## QP RHS: score - Lambda beta + info beta.
    #  The penalty stays in the quadratic operator; this leaves the non-penalty
    #  information contribution plus the current score.
    qp_score <- qp_score_function(
      X_tilde, y_tilde,
      VhalfInv_perm %**% cbind(family$linkinv(XB_final)),
      unlist(order_list), dispersion_final,
      VhalfInv_perm,
      unlist(observation_weights), ...
    )

    info_nopen_beta <- crossprod(
      X_tilde, W_final * (X_tilde %**% beta_block)
    )
    if (any(unlist(schur_corrections_final) != 0)) {
      info_nopen_beta <- info_nopen_beta + schur_coll %**% beta_block
    }
    dvec_full <- qp_score + info_nopen_beta

    W_split_final <- split(W_final, rep(1:(K + 1), n_per_partition))
    Xw_final <- lapply(1:(K + 1), function(k) {
      if (nrow(X[[k]]) == 0) return(X[[k]])
      X[[k]] * sqrt(W_split_final[[k]])
    })
    X_gram_weighted_final <- compute_gram_block_diagonal(
      Xw_final, parallel_matmult, cl, chunk_size, num_chunks, rem_chunks
    )
    wb_final <- .woodbury_redecompose_weighted(
      Delta_V, X_block, DV_X, W_final,
      X, K, p_expansions,
      X_gram_weighted_final,
      Lambda, schur_corrections_final,
      unique_penalty_per_partition, L_partition_list,
      parallel_eigen, cl,
      chunk_size, num_chunks, rem_chunks,
      rank_threshold_fraction = 1/2, family
    )

    as_out <- NULL
    if (isTRUE(wb_final$use_woodbury)) {
      wb_sqrt_final <- .woodbury_halfsqrt_components(
        wb_final$Ghalf_corrected, wb_final$E, wb_final$J,
        wb_final$r, K, p_expansions
      )
      if (isTRUE(wb_sqrt_final$valid)) {
        rhs_list <- lapply(1:(K + 1), function(k) {
          rows_k <- ((k - 1) * p_expansions + 1):(k * p_expansions)
          dvec_full[rows_k, , drop = FALSE]
        })
        as_out <- try(.active_set_refine_woodbury(
          result = result,
          K = K,
          p_expansions = p_expansions,
          A = A,
          R_constraints = R_constraints,
          constraint_value_vectors = constraint_value_vectors,
          family = family,
          qp_Amat = qp_Amat,
          qp_bvec = qp_bvec,
          qp_meq = qp_meq,
          rhs_list = rhs_list,
          Ghalf_corrected = wb_final$Ghalf_corrected,
          wb_sqrt = wb_sqrt_final,
          parallel_aga = parallel_aga,
          parallel_matmult = parallel_matmult,
          cl = cl,
          chunk_size = chunk_size,
          num_chunks = num_chunks,
          rem_chunks = rem_chunks,
          tol = max(tol, 100 * sqrt(.Machine$double.eps))
        ), silent = TRUE)
      }
    }

    if (!is.null(as_out) &&
        !inherits(as_out, "try-error") &&
        isTRUE(as_out$converged)) {
      result <- as_out$result
      qp_info <- as_out$qp_info
      wb <- wb_final
    } else {
      ## Full-system active-set fallback before dense SQP.
      info_full <- crossprod(X_tilde, W_final * X_tilde) +
        Lambda_block + schur_coll
      G_full <- invert(info_full)
      eig_G <- eigen(G_full, symmetric = TRUE)
      vals_G <- pmax(eig_G$values, 0)
      G_full_half <- eig_G$vectors %**%
        (t(eig_G$vectors) * sqrt(vals_G))
      inv_sqrt_vals <- ifelse(sqrt(vals_G) > 0, 1 / sqrt(vals_G), 0)
      G_full_half_inv <- eig_G$vectors %**%
        (t(eig_G$vectors) * inv_sqrt_vals)

      as_out <- try(.active_set_refine_full(
      result = result,
      K = K,
      p_expansions = p_expansions,
      A = A,
      constraint_value_vectors = constraint_value_vectors,
      family = family,
      qp_Amat = qp_Amat,
      qp_bvec = qp_bvec,
      qp_meq = qp_meq,
      rhs_full = dvec_full,
      Ghalf_full = G_full_half,
      tol = max(tol, 100 * sqrt(.Machine$double.eps)),
      method = "active_set_woodbury"
      ), silent = TRUE)

      if (!inherits(as_out, "try-error") && isTRUE(as_out$converged)) {
        result <- as_out$result
        qp_info <- as_out$qp_info
        Ghalf_list <- .extract_G_diagonal(G_full_half, p_expansions, K)
        GhalfInv_list <- .extract_G_diagonal(G_full_half_inv, p_expansions, K)
      } else {
        ## Fall back to dense SQP.
        dense_out <- .get_B_gee_glm(
          X_block, X_tilde, y_block, y_tilde,
          VhalfInv_perm, Lambda_block, A, K,
          p_expansions, constraint_value_vectors,
          family, return_G_getB, iterate, tol,
          qp_Amat, qp_bvec, qp_meq,
          qp_score_function,
          order_list, observation_weights,
          glm_weight_function, schur_correction_function,
          need_dispersion_for_estimation,
          dispersion_function, VhalfInv, ...)
        return(dense_out)
      }
    }
  }

  ## Return.
  if (return_G_getB) {
    if (!is.null(qp_info) && exists("Ghalf_list", inherits = FALSE)) {
      ## QP refinement built the full info; use its blocks.
      G_list <- list(
        G = lapply(Ghalf_list, tcrossprod),
        Ghalf = Ghalf_list,
        GhalfInv = GhalfInv_list
      )
    } else {
      ## No QP refinement; use Woodbury-corrected G.
      G_list <- list(
        G = lapply(wb$Ghalf_corrected, function(m) tcrossprod(m)),
        Ghalf = wb$Ghalf_corrected,
        GhalfInv = wb$GhalfInv_corrected
      )
    }
    return(list(B = result, G_list = G_list, qp_info = qp_info))
  } else {
    return(list(B = result, qp_info = qp_info))
  }
}


#' Path 2: Gaussian Identity Link, No Correlation
#'
#' For the canonical Gaussian case without correlation structures, the
#' unconstrained estimate has the closed form
#' \eqn{\hat{\boldsymbol{\beta}}_k = \mathbf{G}_k \mathbf{X}_k^{\top}
#' \mathbf{y}_k} and the constrained estimate follows by a single
#' Lagrangian projection. No iteration is needed.
#'
#' When \eqn{K = 0} and there are no inequality constraints or nonzero
#' constraint values, the function returns
#' \eqn{\hat{\boldsymbol{\beta}} = \mathbf{G}\mathbf{X}^{\top}\mathbf{y}}
#' directly without forming the full \eqn{P}-dimensional OLS system.
#'
#' @param Xy List of cross-products \eqn{\mathbf{X}_k^{\top}\mathbf{y}_k}
#'   by partition.
#' @param Ghalf,GhalfInv Lists of \eqn{\mathbf{G}^{1/2}_k} and
#'   \eqn{\mathbf{G}^{-1/2}_k} matrices by partition.
#' @param A Constraint matrix \eqn{\mathbf{A}}.
#' @param K,p_expansions,R_constraints Integer dimensions.
#' @param constraint_value_vectors Constraint RHS list encoding
#'   \eqn{\mathbf{A}^{\top}\boldsymbol{\beta} = \mathbf{c}}.
#' @param family GLM family object.
#' @param return_G_getB Logical; return covariance components.
#' @param quadprog Logical; apply QP refinement for inequality constraints.
#' @param qp_Amat,qp_bvec,qp_meq QP constraint specification.
#' @param qp_score_function Score function for QP step.
#' @param parallel_aga,parallel_matmult Logical flags.
#' @param cl,chunk_size,num_chunks,rem_chunks Parallel parameters.
#' @param X,y,Lambda,order_list,observation_weights Standard arguments
#'   passed through for optional QP refinement.
#' @param iterate,tol,glm_weight_function,schur_correction_function,
#'   need_dispersion_for_estimation,dispersion_function,
#'   unique_penalty_per_partition,L_partition_list,VhalfInv,
#'   homogenous_weights Standard arguments passed through.
#' @param ... Passed to sub-functions.
#'
#' @return Same structure as \code{get_B}.
#'
#' @keywords internal
.get_B_gaussian_nocorr <- function(Xy, Ghalf, GhalfInv, A, K,
                                   p_expansions, R_constraints,
                                   constraint_value_vectors, family,
                                   return_G_getB,
                                   quadprog, qp_Amat, qp_bvec, qp_meq,
                                   qp_score_function,
                                   parallel_aga, parallel_matmult,
                                   cl, chunk_size, num_chunks, rem_chunks,
                                   X, y, Lambda,
                                   order_list, observation_weights,
                                   iterate, tol,
                                   glm_weight_function,
                                   schur_correction_function,
                                   need_dispersion_for_estimation,
                                   dispersion_function,
                                   unique_penalty_per_partition,
                                   L_partition_list, VhalfInv,
                                   homogenous_weights, ...) {

  ## Fast path: K = 0, no inequality constraints, homogeneous RHS.
  #  Direct solution: beta = G X^T y (one partition, no projection needed).
  if (K == 0 & !quadprog & length(constraint_value_vectors) == 0) {
    G <- list(tcrossprod(Ghalf[[1]]))
    result <- list(G[[1]] %**% Xy[[1]])
    if (return_G_getB) {
      Ghalf_exact <- lapply(G, function(Gk) {
        matsqrt(Gk)
      })
      GhalfInv_exact <- lapply(G, function(Gk) {
        matinvsqrt(Gk)
      })
      return(list(
        B = result,
        G_list = list(G = G,
                      Ghalf = Ghalf_exact,
                      GhalfInv = GhalfInv_exact),
        qp_info = NULL
      ))
    } else {
      return(list(B = result, qp_info = NULL))
    }
  }

  ## General path: y* = G^{1/2} X^T y, then Lagrangian projection.
  GhalfXy <- cbind(
    unlist(
      matmult_block_diagonal(Ghalf, Xy, K, parallel_matmult, cl,
                             chunk_size, num_chunks, rem_chunks)))

  result <- .lagrangian_project(GhalfXy, Ghalf, A, K, p_expansions,
                                R_constraints, constraint_value_vectors,
                                family, parallel_aga, parallel_matmult,
                                cl, chunk_size, num_chunks, rem_chunks)

  qp_info <- NULL

  ## Optional QP refinement for inequality constraints.
  if (quadprog) {
    ## Active set first.
    #  If repeated equality re-solves fail, drop to the dense QP fallback.
    as_out <- try(.active_set_refine(
      result, X, y, K, p_expansions,
      A, R_constraints, constraint_value_vectors,
      Lambda, Ghalf, GhalfInv, family,
      qp_Amat, qp_bvec, qp_meq,
      Xy_or_uncon = Xy, is_path3 = FALSE,
      parallel_aga, parallel_matmult,
      cl, chunk_size, num_chunks, rem_chunks, tol), silent = TRUE)
    if (!inherits(as_out, "try-error") && as_out$converged) {
      result <- as_out$result
      qp_info <- as_out$qp_info
    } else {
      Lambda_block <- .build_lambda_block(Lambda, K,
                                          unique_penalty_per_partition,
                                          L_partition_list)
      qp_out <- .qp_refine(result, X, y, K, p_expansions, A, Lambda,
                           Lambda_block, family, iterate, tol,
                           qp_Amat, qp_bvec, qp_meq,
                           qp_score_function,
                           order_list, glm_weight_function,
                           schur_correction_function,
                           need_dispersion_for_estimation,
                           dispersion_function, observation_weights,
                           VhalfInv, ...)
      result <- qp_out$result
      qp_info <- qp_out$qp_info
    }
  }

  if (return_G_getB) {
    G_exact <- lapply(Ghalf, function(mat) tcrossprod(mat))
    return(list(
      B = result,
      G_list = list(G = G_exact,
                    Ghalf = Ghalf,
                    GhalfInv = GhalfInv),
      qp_info = qp_info
    ))
  } else {
    return(list(B = result, qp_info = qp_info))
  }
}


## Path 3: non-Gaussian without correlation

#' Path 3: Non-Gaussian GLM, No Correlation
#'
#' For GLMs with non-identity links or non-Gaussian families (without
#' correlation structures), unconstrained partition-wise estimates are
#' first obtained via Newton-Raphson (or OLS for Gaussian identity), then
#' the constrained estimate is computed by a single Lagrangian projection.
#' For non-canonical links, \eqn{\mathbf{G}} depends on the current
#' fitted values through the GLM working weights \eqn{\mathbf{W}}; the
#' projection is therefore iterated, updating \eqn{\mathbf{G}} at the
#' current constrained estimate until convergence.
#'
#' @param X,y Lists of partition-specific design matrices and responses.
#' @param X_gram List of Gram matrices.
#' @param Xy List of cross-products.
#' @param Lambda Shared penalty matrix.
#' @param Ghalf,GhalfInv Lists of matrix square roots by partition.
#' @param A Constraint matrix.
#' @param K,p_expansions,R_constraints Integer dimensions.
#' @param constraint_value_vectors Constraint RHS list encoding
#'   \eqn{\mathbf{A}^{\top}\boldsymbol{\beta} = \mathbf{c}}.
#' @param family GLM family object.
#' @param return_G_getB Logical; return covariance components.
#' @param iterate Logical; iterate for non-canonical links.
#' @param tol Convergence tolerance.
#' @param quadprog Logical; apply QP refinement.
#' @param qp_Amat,qp_bvec,qp_meq QP constraint specification.
#' @param qp_score_function Score function for QP step.
#' @param unconstrained_fit_fxn Function for partition-wise unconstrained
#'   estimation.
#' @param keep_weighted_Lambda,unique_penalty_per_partition Logical flags.
#' @param L_partition_list Partition-specific penalty matrices.
#' @param parallel_eigen,parallel_aga,parallel_matmult,
#'   parallel_unconstrained Logical flags.
#' @param cl,chunk_size,num_chunks,rem_chunks Parallel parameters.
#' @param order_list,observation_weights Standard partition arguments.
#' @param glm_weight_function,schur_correction_function,
#'   need_dispersion_for_estimation,dispersion_function GLM customization.
#' @param VhalfInv Inverse square root correlation (unused here).
#' @param homogenous_weights Logical.
#' @param ... Passed to fitting, weight, and score functions.
#'
#' @return Same structure as \code{get_B}.
#'
#' @keywords internal
.get_B_glm_nocorr <- function(X, y, X_gram, Xy, Lambda, Ghalf, GhalfInv,
                              A, K, p_expansions, R_constraints,
                              constraint_value_vectors, family,
                              return_G_getB, iterate, tol,
                              quadprog, qp_Amat, qp_bvec, qp_meq,
                              qp_score_function,
                              unconstrained_fit_fxn,
                              keep_weighted_Lambda,
                              unique_penalty_per_partition,
                              L_partition_list,
                              parallel_eigen, parallel_aga,
                              parallel_matmult, parallel_unconstrained,
                              cl, chunk_size, num_chunks, rem_chunks,
                              order_list, observation_weights,
                              glm_weight_function,
                              schur_correction_function,
                              need_dispersion_for_estimation,
                              dispersion_function, VhalfInv,
                              homogenous_weights, ...) {

  ## Lambda^{1/2} for augmented regression in unconstrained_fit_fxn.
  LambdaHalf <- matsqrt(Lambda)

  ## Obtain unconstrained estimates partition-wise.
  if (parallel_unconstrained) {
    unconstrained_estimate <-
      parallel::parLapply(cl, 1:(K + 1),
                          function(k,
                                   unique_penalty_per_partition,
                                   unconstrained_fit_fxn,
                                   keep_weighted_Lambda,
                                   family, tol, K,
                                   parallel_unconstrained, cl,
                                   chunk_size, num_chunks, rem_chunks,
                                   observation_weights) {
                            if (unique_penalty_per_partition) {
                              Lambda_temp <- Lambda + L_partition_list[[k]]
                              LambdaHalf_temp <- matsqrt(Lambda_temp)
                            } else {
                              LambdaHalf_temp <- LambdaHalf
                              Lambda_temp <- Lambda
                            }
                            cbind(c(unconstrained_fit_fxn(
                              X[[k]], y[[k]], LambdaHalf_temp, Lambda_temp,
                              keep_weighted_Lambda, family, tol, K,
                              parallel_unconstrained, cl,
                              chunk_size, num_chunks, rem_chunks,
                              order_list[[k]], observation_weights[[k]],
                              ...)))
                          },
                          unique_penalty_per_partition,
                          unconstrained_fit_fxn,
                          keep_weighted_Lambda,
                          family, tol, K,
                          parallel_unconstrained, cl,
                          chunk_size, num_chunks, rem_chunks,
                          observation_weights)
  } else {
    unconstrained_estimate <- lapply(1:(K + 1), function(k) {
      if (unique_penalty_per_partition) {
        Lambda_temp <- Lambda + L_partition_list[[k]]
        LambdaHalf_temp <- matsqrt(Lambda_temp)
      } else {
        LambdaHalf_temp <- LambdaHalf
        Lambda_temp <- Lambda
      }
      cbind(c(unconstrained_fit_fxn(
        X[[k]], y[[k]], LambdaHalf_temp, Lambda_temp,
        keep_weighted_Lambda, family, tol, K,
        parallel_unconstrained, cl,
        chunk_size, num_chunks, rem_chunks,
        order_list[[k]], observation_weights[[k]],
        ...)))
    })
  }

  ## Fast path: K = 0, no inequality constraints, homogeneous RHS.
  if (K == 0 & !quadprog & length(constraint_value_vectors) == 0) {
    if (return_G_getB) {
      G_list <- .recompute_G_at_estimate(
        X, y, unconstrained_estimate, K, Lambda, family,
        order_list, glm_weight_function, schur_correction_function,
        need_dispersion_for_estimation, dispersion_function,
        observation_weights, VhalfInv,
        parallel_eigen, parallel_matmult, cl,
        chunk_size, num_chunks, rem_chunks,
        unique_penalty_per_partition, L_partition_list, ...
      )
      return(list(B = unconstrained_estimate, G_list = G_list,
                  qp_info = NULL))
    } else {
      return(list(B = unconstrained_estimate, qp_info = NULL))
    }
  }

  ## Transform: y* = G^{-1/2} beta_hat.
  GhalfXy <- cbind(
    unlist(
      matmult_block_diagonal(
        GhalfInv, unconstrained_estimate, K, parallel_matmult,
        cl, chunk_size, num_chunks, rem_chunks)))

  ## Initial Lagrangian projection.
  result <- .lagrangian_project(GhalfXy, Ghalf, A, K, p_expansions,
                                R_constraints, constraint_value_vectors,
                                family, parallel_aga, parallel_matmult,
                                cl, chunk_size, num_chunks, rem_chunks)

  ## Iterative projection for non-canonical links
  if (iterate) {
    prevB_IRWLS <- NULL
    prev_diff_IRWLS <- Inf

    for (IRWLS_iter in 1:100) {
      if (!is.null(prevB_IRWLS)) {
        diff_IRWLS <- mean(
          unlist(
            sapply(1:(K + 1),
                   function(k) mean(
                     abs(result[[k]] - prevB_IRWLS[[k]])))))
        if (diff_IRWLS < tol) break
        if (prev_diff_IRWLS <= diff_IRWLS) {
          result <- prevB_IRWLS
          break
        }
        prev_diff_IRWLS <- diff_IRWLS
      }

      prevB_IRWLS <- result

      if (need_dispersion_for_estimation) {
        dispersion_temp <- dispersion_function(
          mu = family$linkinv(unlist(matmult_block_diagonal(
            X, result, K, parallel_matmult, cl,
            chunk_size, num_chunks, rem_chunks))),
          y = unlist(y),
          order_indices = unlist(order_list),
          family = family,
          observation_weights = unlist(observation_weights),
          VhalfInv = VhalfInv,
          ...
        )
      } else {
        dispersion_temp <- 1
      }

      Xw <- lapply(1:(K + 1), function(k) {
        if (nrow(X[[k]]) == 0) {
          return(X[[k]])
        }
        var <- glm_weight_function(
          family$linkinv(X[[k]] %**% cbind(c(result[[k]]))),
          y[[k]], order_list[[k]], family, dispersion_temp,
          observation_weights[[k]], ...
        )
        if (length(var) == 1) {
          if (c(var) == 0) {
            return(X[[k]] * 0)
          } else {
            return(X[[k]] * c(sqrt(var)))
          }
        } else {
          var <- c(sqrt(var))
        }
        X[[k]] * var
      })
      X_gram_IRWLS <- compute_gram_block_diagonal(
        Xw, parallel_matmult, cl, chunk_size, num_chunks, rem_chunks)

      schur_corrections_IRWLS <- schur_correction_function(
        X, y, result, dispersion_temp, order_list, K, family,
        observation_weights, ...
      )
      G_list_IRWLS <- compute_G_eigen(
        X_gram_IRWLS, Lambda, K, parallel_eigen, cl,
        chunk_size, num_chunks, rem_chunks, family,
        unique_penalty_per_partition, L_partition_list,
        keep_G = FALSE, schur_corrections_IRWLS)

      Ghalf <- G_list_IRWLS$Ghalf
      GhalfInv <- G_list_IRWLS$GhalfInv

      GhalfXy <- cbind(
        unlist(
          matmult_block_diagonal(
            GhalfInv, unconstrained_estimate, K, parallel_matmult,
            cl, chunk_size, num_chunks, rem_chunks)))

      result <- .lagrangian_project(GhalfXy, Ghalf, A, K, p_expansions,
                                    R_constraints, constraint_value_vectors,
                                    family, parallel_aga, parallel_matmult,
                                    cl, chunk_size, num_chunks, rem_chunks)
    }
  }

  qp_info <- NULL

  ## Optional QP refinement for inequality constraints.
  if (quadprog) {
    ## Active set first.
    #  If repeated equality re-solves fail, drop to the dense QP fallback.
    as_out <- try(.active_set_refine(
      result, X, y, K, p_expansions,
      A, R_constraints, constraint_value_vectors,
      Lambda, Ghalf, GhalfInv, family,
      qp_Amat, qp_bvec, qp_meq,
      Xy_or_uncon = unconstrained_estimate, is_path3 = TRUE,
      parallel_aga, parallel_matmult,
      cl, chunk_size, num_chunks, rem_chunks, tol), silent = TRUE)
    if (!inherits(as_out, "try-error") && as_out$converged) {
      result <- as_out$result
      qp_info <- as_out$qp_info
    } else {
      Lambda_block <- .build_lambda_block(Lambda, K,
                                          unique_penalty_per_partition,
                                          L_partition_list)
      qp_out <- .qp_refine(result, X, y, K, p_expansions, A, Lambda,
                           Lambda_block, family, iterate, tol,
                           qp_Amat, qp_bvec, qp_meq,
                           qp_score_function,
                           order_list, glm_weight_function,
                           schur_correction_function,
                           need_dispersion_for_estimation,
                           dispersion_function, observation_weights,
                           VhalfInv, ...)
      result <- qp_out$result
      qp_info <- qp_out$qp_info
    }
  }

  if (return_G_getB) {
    if (paste0(family)[1] == 'gaussian' &
        paste0(family)[2] == 'identity') {
      G_exact <- lapply(Ghalf, function(mat) tcrossprod(mat))
      return(list(
        B = result,
        G_list = list(G = G_exact,
                      Ghalf = Ghalf,
                      GhalfInv = GhalfInv),
        qp_info = qp_info
      ))
    }
    G_list <- .recompute_G_at_estimate(
      X, y, result, K, Lambda, family,
      order_list, glm_weight_function, schur_correction_function,
      need_dispersion_for_estimation, dispersion_function,
      observation_weights, VhalfInv,
      parallel_eigen, parallel_matmult, cl,
      chunk_size, num_chunks, rem_chunks,
      unique_penalty_per_partition, L_partition_list, ...
    )
    return(list(B = result, G_list = G_list, qp_info = qp_info))
  } else {
    return(list(B = result, qp_info = qp_info))
  }
}


## Main function: get_B

#' Compute Constrained GLM Coefficient Estimates via Lagrangian Multipliers
#'
#' @description
#' Core estimation function for Lagrangian multiplier smoothing splines.
#' Computes penalized coefficient estimates subject to smoothness constraints
#' (continuity and derivative matching at knots) and optional user-supplied
#' linear equality or inequality constraints. Dispatches to one of three
#' computational paths depending on the model structure:
#'
#' \bold{Path 1. GEE (correlation structure present):}
#' When both \code{Vhalf} and \code{VhalfInv} are provided, the full
#' \eqn{N \times P} whitened design
#' \eqn{\mathbf{X}^{*} = \mathbf{V}_{\mathrm{perm}}^{-1/2}\mathbf{X}} is used
#' after permuting observations to partition ordering via \code{order_list}.
#' For structured correlation, the solver tries the supplement's
#' decomposition
#' \eqn{\mathbf{G}^{-1} =
#' \mathbf{G}_{\mathrm{on}}^{-1} + \mathbf{G}_{\mathrm{off}}^{-1}}
#' with
#' \eqn{\mathbf{G}_{\mathrm{off}}^{-1} = \mathbf{E}\mathbf{J}\mathbf{E}^{\top}}.
#' This gives a Woodbury path with cost \eqn{O(Kp^3 + Pr^2)} when the
#' cross-partition rank \eqn{r} is small; otherwise the dense
#' \eqn{O(P^3)} path is used.
#' Two sub-paths:
#' \itemize{
#'   \item \bold{Path 1a.} Gaussian identity + GEE.
#'     See \code{\link{.get_B_gee_gaussian}} (dense) and
#'     \code{\link{.get_B_gee_woodbury}} (accelerated).
#'   \item \bold{Path 1b.} Non-Gaussian GEE.
#'     See \code{\link{.get_B_gee_glm}} (dense) and
#'     \code{\link{.get_B_gee_glm_woodbury}} (accelerated).
#' }
#'
#' \bold{Path 2. Gaussian identity link, no correlation:}
#' Unconstrained estimate \eqn{\hat{\boldsymbol{\beta}}_k =
#' \mathbf{G}_k\mathbf{X}_k^{\top}\mathbf{y}_k} per partition, then a
#' single Lagrangian projection. See \code{\link{.get_B_gaussian_nocorr}}.
#'
#' \bold{Path 3. Non-Gaussian GLM, no correlation:}
#' Partition-wise unconstrained estimates via Newton-Raphson, then
#' Lagrangian projection. See \code{\link{.get_B_glm_nocorr}}.
#'
#' \bold{Inequality constraint handling:}
#' When inequality constraints \eqn{\mathbf{C}^{\top}\boldsymbol{\beta}
#' \succeq \mathbf{c}} are present, the sparsity pattern of
#' \eqn{\mathbf{C}} is inspected automatically. If every constraint
#' column has nonzeros in only a single partition block (block-separable),
#' a partition-wise active-set method is used at cost
#' \eqn{O(Kp^3)} per iteration. If any constraint spans multiple
#' partition blocks (e.g.\ cross-knot monotonicity), the dense
#' \code{quadprog::solve.QP} SQP fallback is invoked.
#'
#' @param X List of length \eqn{K+1}; partition-specific design matrices.
#' @param X_gram List of Gram matrices by partition.
#' @param Lambda Combined penalty matrix.
#' @param keep_weighted_Lambda Logical.
#' @param unique_penalty_per_partition Logical.
#' @param L_partition_list List of partition-specific penalty matrices.
#' @param A Constraint matrix \eqn{\mathbf{A}} (\eqn{P \times R}).
#' @param Xy List of cross-products by partition.
#' @param y List of response vectors by partition.
#' @param K Integer; number of interior knots.
#' @param p_expansions Integer; basis terms per partition.
#' @param R_constraints Integer; columns of \eqn{\mathbf{A}}.
#' @param Ghalf List of \eqn{\mathbf{G}^{1/2}_k} matrices.
#' @param GhalfInv List of \eqn{\mathbf{G}^{-1/2}_k} matrices.
#' @param parallel_eigen,parallel_aga,parallel_matmult,parallel_unconstrained
#'   Logical flags.
#' @param cl Cluster object.
#' @param chunk_size,num_chunks,rem_chunks Parallel distribution parameters.
#' @param family GLM family object.
#' @param unconstrained_fit_fxn Partition-wise unconstrained estimator.
#' @param iterate Logical; iterate for non-canonical links.
#' @param qp_score_function Score function for QP steps.
#' @param quadprog Logical; use \code{quadprog::solve.QP} for inequality
#'   constraints.
#' @param qp_Amat,qp_bvec,qp_meq Inequality constraint specification
#'   \eqn{\mathbf{C}^{\top}\boldsymbol{\beta} \succeq \mathbf{c}}.
#' @param prevB,prevUnconB,iter_count,prev_diff Retired; ignored.
#' @param tol Convergence tolerance.
#' @param constraint_value_vectors List encoding nonzero RHS
#'   \eqn{\mathbf{A}^{\top}\boldsymbol{\beta} = \mathbf{c}}.
#' @param order_list Partition-to-data index mapping.
#' @param glm_weight_function Function computing working weights.
#' @param schur_correction_function Function computing Schur corrections.
#' @param need_dispersion_for_estimation Logical.
#' @param dispersion_function Dispersion estimator.
#' @param observation_weights List of observation weights by partition.
#' @param homogenous_weights Logical.
#' @param return_G_getB Logical; return \eqn{\mathbf{G}},
#'   \eqn{\mathbf{G}^{1/2}}, \eqn{\mathbf{G}^{-1/2}} per partition.
#' @param blockfit,just_linear_without_interactions Retired; retained for
#'   call-site compatibility.
#' @param Vhalf,VhalfInv Square root and inverse square root of the
#'   working correlation matrix in the original observation ordering.
#'   When both are non-\code{NULL}, dispatches to Path 1; \code{VhalfInv}
#'   is permuted internally to partition ordering via \code{order_list}.
#' @param ... Passed to fitting, weight, correction, and dispersion
#'   functions.
#'
#' @return
#' If \code{return_G_getB = FALSE}: a list with \code{B} (coefficient
#' column vectors by partition) and \code{qp_info}.
#'
#' If \code{return_G_getB = TRUE}: a list with elements:
#' \describe{
#'   \item{B}{List of constrained coefficient column vectors
#'     \eqn{\tilde{\boldsymbol{\beta}}_k} by partition.}
#'   \item{G_list}{List with \code{G}, \code{Ghalf}, \code{GhalfInv},
#'     each a list of \eqn{K+1} matrices.}
#'   \item{qp_info}{QP or active-set metadata, including Lagrangian
#'     multipliers and active-constraint information when available, or
#'     \code{NULL}.}
#' }
#'
#' @keywords internal
#' @export
get_B <- function(X,
                  X_gram,
                  Lambda,
                  keep_weighted_Lambda,
                  unique_penalty_per_partition,
                  L_partition_list,
                  A,
                  Xy,
                  y,
                  K,
                  p_expansions,
                  R_constraints,
                  Ghalf,
                  GhalfInv,
                  parallel_eigen,
                  parallel_aga,
                  parallel_matmult,
                  parallel_unconstrained,
                  cl,
                  chunk_size,
                  num_chunks,
                  rem_chunks,
                  family,
                  unconstrained_fit_fxn,
                  iterate,
                  qp_score_function,
                  quadprog,
                  qp_Amat,
                  qp_bvec,
                  qp_meq,
                  prevB = NULL,
                  prevUnconB = NULL,
                  iter_count = 0,
                  prev_diff = Inf,
                  tol,
                  constraint_value_vectors,
                  order_list,
                  glm_weight_function,
                  schur_correction_function,
                  need_dispersion_for_estimation,
                  dispersion_function,
                  observation_weights,
                  homogenous_weights,
                  return_G_getB,
                  blockfit,
                  just_linear_without_interactions,
                  Vhalf,
                  VhalfInv,
                  ...) {

  ## Path 1: GEE (correlation structure present)
  if (!is.null(Vhalf) & !is.null(VhalfInv)) {

    ## Permute correlation matrices so rows/columns align with partition
    #  ordering rather than the original data ordering.
    perm <- unlist(order_list)
    Vhalf_perm    <- Vhalf[perm, perm]
    VhalfInv_perm <- VhalfInv[perm, perm]

    ## Default to unit observation weights when none were supplied.
    if (is.null(observation_weights) |
        length(unlist(observation_weights)) == 0) {
      observation_weights <- 1
    }

    ## Assemble the full N x P block-diagonal design.
    X_block <- collapse_block_diagonal(X)
    y_block <- cbind(unlist(y))

    ## Whiten: X_tilde = V^{-1/2} X, y_tilde = V^{-1/2} y.
    X_tilde <- VhalfInv_perm %**% X_block
    y_tilde <- VhalfInv_perm %**% y_block

    ## Full P x P block-diagonal penalty matrix.
    Lambda_block <- .build_lambda_block(Lambda, K,
                                        unique_penalty_per_partition,
                                        L_partition_list)

    ## Split the GEE path into the Gaussian one-step solve versus
    #  the iterative non-Gaussian path.
    is_gauss_id <- (paste0(family)[1] == 'gaussian' &
                      paste0(family)[2] == 'identity')

    ## Woodbury gate
    #  If the off-block correlated piece is low rank, keep the solve on the
    #  reduced Woodbury path instead of forming the full dense system.
    wb_decomp <- .woodbury_decompose_V(
      VhalfInv_perm, X, K, p_expansions,
      Lambda, Lambda_block,
      unique_penalty_per_partition, L_partition_list,
      order_list,
      parallel_eigen, cl,
      chunk_size, num_chunks, rem_chunks,
      rank_threshold_fraction = 1/3,
      family)

    if (is_gauss_id) {
      obs_wt_full <- c(unlist(observation_weights))
      if (length(obs_wt_full) == 0L) {
        obs_wt_full <- rep(1, nrow(VhalfInv_perm))
      } else if (length(obs_wt_full) == 1L) {
        obs_wt_full <- rep(obs_wt_full, nrow(VhalfInv_perm))
      }
      ## The Gaussian Woodbury path currently assumes unit observation weights.
      can_use_gauss_woodbury <- all(abs(obs_wt_full - 1) <
                                      sqrt(.Machine$double.eps))

      if (wb_decomp$use_woodbury && can_use_gauss_woodbury) {
        ## Compute the low-rank F^{1/2} pieces once for this solve.
        wb_sqrt <- .woodbury_halfsqrt_components(
          wb_decomp$Ghalf_corrected, wb_decomp$E, wb_decomp$J,
          wb_decomp$r, K, p_expansions)

        if (wb_sqrt$valid) {
          ## Path 1a-Woodbury.
          out <- .get_B_gee_woodbury(
            X, y, K, p_expansions,
            VhalfInv_perm, order_list,
            A, R_constraints,
            constraint_value_vectors, family,
            return_G_getB,
            quadprog, qp_Amat, qp_bvec, qp_meq,
            qp_score_function,
            observation_weights,
            wb_decomp, wb_sqrt,
            Lambda_block,
            parallel_aga, parallel_matmult,
            cl, chunk_size, num_chunks, rem_chunks,
            qp_global = .detect_qp_global(qp_Amat, p_expansions, K),
            tol = tol, ...)
        } else {
          ## F failed the half-square-root step, so use the dense path.
          out <- .get_B_gee_gaussian(
            X_block, X_tilde, y_tilde, VhalfInv_perm,
            Lambda_block, A, K, p_expansions,
            constraint_value_vectors, family,
            return_G_getB,
            quadprog, qp_Amat, qp_bvec, qp_meq,
            qp_score_function,
            order_list, observation_weights, ...)
          out <- .attach_qp_fallback_info(
            out, "woodbury_square_root_not_positive_definite", wb_decomp)
        }
      } else {
        ## If the gate fails, stay on the dense Gaussian GEE path.
        out <- .get_B_gee_gaussian(
          X_block, X_tilde, y_tilde, VhalfInv_perm,
          Lambda_block, A, K, p_expansions,
          constraint_value_vectors, family,
          return_G_getB,
          quadprog, qp_Amat, qp_bvec, qp_meq,
          qp_score_function,
          order_list, observation_weights, ...)
        fallback_reason <- if (!can_use_gauss_woodbury) {
          "nonunit_observation_weights"
        } else {
          wb_decomp$reason
        }
        out <- .attach_qp_fallback_info(out, fallback_reason, wb_decomp)
      }
    } else {
      if (wb_decomp$use_woodbury) {
        wb_sqrt <- .woodbury_halfsqrt_components(
          wb_decomp$Ghalf_corrected, wb_decomp$E, wb_decomp$J,
          wb_decomp$r, K, p_expansions)

        if (wb_sqrt$valid) {
          ## Path 1b-Woodbury.
          out <- .get_B_gee_glm_woodbury(
            X, y, K, p_expansions,
            VhalfInv_perm, order_list,
            A, R_constraints,
            constraint_value_vectors, family,
            return_G_getB, iterate, tol,
            quadprog,
            qp_Amat, qp_bvec, qp_meq,
            qp_score_function,
            observation_weights,
            Lambda, Lambda_block,
            unique_penalty_per_partition, L_partition_list,
            wb_decomp, wb_sqrt,
            glm_weight_function, schur_correction_function,
            need_dispersion_for_estimation, dispersion_function,
            VhalfInv,
            parallel_eigen, parallel_aga, parallel_matmult,
            cl, chunk_size, num_chunks, rem_chunks, ...)
        } else {
          ## F failed the half-square-root step, so use the dense path.
          out <- .get_B_gee_glm(
            X_block, X_tilde, y_block, y_tilde,
            VhalfInv_perm, Lambda_block, A, K, p_expansions,
            constraint_value_vectors, family,
            return_G_getB, iterate, tol,
            qp_Amat, qp_bvec, qp_meq,
            qp_score_function,
            order_list, observation_weights,
            glm_weight_function, schur_correction_function,
            need_dispersion_for_estimation, dispersion_function,
            VhalfInv, ...)
          out <- .attach_qp_fallback_info(
            out, "woodbury_square_root_not_positive_definite", wb_decomp)
        }
      } else {
        ## If the gate fails, stay on the dense non-Gaussian GEE path.
        out <- .get_B_gee_glm(
          X_block, X_tilde, y_block, y_tilde,
          VhalfInv_perm, Lambda_block, A, K, p_expansions,
          constraint_value_vectors, family,
          return_G_getB, iterate, tol,
          qp_Amat, qp_bvec, qp_meq,
          qp_score_function,
          order_list, observation_weights,
          glm_weight_function, schur_correction_function,
          need_dispersion_for_estimation, dispersion_function,
          VhalfInv, ...)
        out <- .attach_qp_fallback_info(out, wb_decomp$reason, wb_decomp)
      }
    }

    ## Normalize return format.
    if (return_G_getB) {
      return(list(B = out$B, G_list = out$G_list, qp_info = out$qp_info))
    } else {
      return(out$B)
    }
  }

  ## Paths 2 and 3: no correlation structure

  ## Without correlation, the family decides between the Gaussian
  #  closed-form path and the iterative GLM path.
  is_gauss_id <- (paste0(family)[2] == 'identity') &
    (paste0(family)[1] == 'gaussian')

  if (is_gauss_id) {
    ## Path 2: Gaussian identity, no correlation.
    out <- .get_B_gaussian_nocorr(
      Xy, Ghalf, GhalfInv, A, K, p_expansions, R_constraints,
      constraint_value_vectors, family,
      return_G_getB,
      quadprog, qp_Amat, qp_bvec, qp_meq,
      qp_score_function,
      parallel_aga, parallel_matmult,
      cl, chunk_size, num_chunks, rem_chunks,
      X, y, Lambda,
      order_list, observation_weights,
      iterate, tol,
      glm_weight_function, schur_correction_function,
      need_dispersion_for_estimation, dispersion_function,
      unique_penalty_per_partition, L_partition_list,
      VhalfInv, homogenous_weights, ...
    )
  } else {
    ## Path 3: Non-Gaussian GLM, no correlation.
    out <- .get_B_glm_nocorr(
      X, y, X_gram, Xy, Lambda, Ghalf, GhalfInv,
      A, K, p_expansions, R_constraints,
      constraint_value_vectors, family,
      return_G_getB, iterate, tol,
      quadprog, qp_Amat, qp_bvec, qp_meq,
      qp_score_function,
      unconstrained_fit_fxn,
      keep_weighted_Lambda,
      unique_penalty_per_partition, L_partition_list,
      parallel_eigen, parallel_aga, parallel_matmult,
      parallel_unconstrained,
      cl, chunk_size, num_chunks, rem_chunks,
      order_list, observation_weights,
      glm_weight_function, schur_correction_function,
      need_dispersion_for_estimation, dispersion_function,
      VhalfInv, homogenous_weights, ...
    )
  }

  ## Normalize return format.
  if (return_G_getB) {
    return(list(B = out$B, G_list = out$G_list, qp_info = out$qp_info))
  } else {
    return(out$B)
  }
}



#' @title Verification Examples for get_B
#'
#' @description
#' Simple, self-contained examples that reviewers can run to verify that
#' \code{get_B} produces correct output. These exercise Path 2 (Gaussian,
#' no correlation) and Path 3 (binomial GLM).
#'
#' @examples
#' \dontrun{
#' ## Example 1: Path 2 - Gaussian identity, with knots
#' set.seed(42)
#' t <- runif(200, -5, 5)
#' y <- sin(t) + rnorm(200, 0, 0.5)
#' fit1 <- lgspline(t, y, K = 3, opt = FALSE, wiggle_penalty = 1e-4)
#' stopifnot(inherits(fit1, "lgspline"))
#' stopifnot(length(fit1$B) == 4)  # K+1 = 4 partitions
#' cat("Example 1 passed: Gaussian identity, K=3\n")
#'
#' preds1 <- predict(fit1, new_predictors = rnorm(10))
#' stopifnot(all(is.finite(preds1)))
#' cat("  Predictions finite: OK\n")
#'
#' ## Example 2: Path 2 - Gaussian identity, K=0 (no constraints)
#' fit2 <- lgspline(t, y, K = 0, opt = FALSE, wiggle_penalty = 1e-4)
#' stopifnot(inherits(fit2, "lgspline"))
#' stopifnot(length(fit2$B) == 1)
#' preds2 <- predict(fit2, new_predictors = rnorm(10))
#' stopifnot(all(is.finite(preds2)))
#' cat("Example 2 passed: Gaussian identity, K=0\n")
#'
#' ## Example 3: Path 3 - Binomial GLM
#' y_bin <- rbinom(200, 1, plogis(sin(t)))
#' fit3 <- lgspline(t, y_bin, K = 2, family = binomial(),
#'                  opt = FALSE, wiggle_penalty = 1e-3)
#' stopifnot(inherits(fit3, "lgspline"))
#' preds3 <- predict(fit3, new_predictors = rnorm(10))
#' stopifnot(all(preds3 >= 0 & preds3 <= 1))
#' cat("Example 3 passed: Binomial GLM, K=2\n")
#'
#' ## Example 4: Path 2 with QP constraints (monotonic increase)
#' t_sorted <- sort(runif(100, -3, 3))
#' y_mono <- t_sorted + 0.5 * sin(t_sorted) + rnorm(100, 0, 0.3)
#' fit4 <- lgspline(t_sorted, y_mono, K = 2,
#'                  qp_monotonic_increase = TRUE,
#'                  opt = FALSE, wiggle_penalty = 1e-4)
#' preds4 <- predict(fit4, new_predictors = cbind(sort(rnorm(50))))
#' stopifnot(all(diff(preds4) >= -sqrt(.Machine$double.eps)))
#' cat("Example 4 passed: Monotonic increase QP constraint\n")
#'
#' ## Example 5: Path 2 with range constraints
#' fit5 <- lgspline(t, y, K = 3,
#'                  qp_range_lower = -2, qp_range_upper = 2,
#'                  opt = FALSE, wiggle_penalty = 1e-4)
#' stopifnot(all(fit5$ytilde >= -2 - 0.01))
#' stopifnot(all(fit5$ytilde <= 2 + 0.01))
#' cat("Example 5 passed: Range constraints\n")
#'
#' ## Example 6: Multi-predictor
#' set.seed(123)
#' x1 <- rnorm(300)
#' x2 <- rnorm(300)
#' y6 <- sin(x1) + cos(x2) + rnorm(300, 0, 0.5)
#' fit6 <- lgspline(cbind(x1, x2), y6, K = 4,
#'                  opt = FALSE, wiggle_penalty = 1e-5)
#' stopifnot(inherits(fit6, "lgspline"))
#' stopifnot(fit6$q_predictors == 2)
#' cat("Example 6 passed: 2D predictor, K=4\n")
#'
#' ## Example 7: Coefficient consistency (determinism)
#' set.seed(999)
#' t7 <- runif(150, -4, 4)
#' y7 <- 2 * cos(t7) + rnorm(150, 0, 0.4)
#' fit7a <- lgspline(t7, y7, K = 2, opt = FALSE, wiggle_penalty = 1e-5)
#' set.seed(999)
#' fit7b <- lgspline(t7, y7, K = 2, opt = FALSE, wiggle_penalty = 1e-5)
#' max_diff <- max(abs(unlist(fit7a$B) - unlist(fit7b$B)))
#' stopifnot(max_diff < 1e-12)
#' cat("Example 7 passed: Deterministic coefficient reproduction\n")
#'
#' cat("\nAll verification examples passed.\n")
#' }
#'
#' @name get_B_verification_examples
NULL
