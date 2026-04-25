## solver_utils.R
#  Shared solver helpers used by get_B() and blockfit_solve().
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




#' Detect Whether Inequality Constraints Need the Global Bridge
#'
#' Inspects the columns of \code{qp_Amat} to determine whether every
#' inequality constraint is confined to a single partition block
#' (block-separable) or whether any constraint couples coefficients
#' across partitions.
#'
#' Returns \code{FALSE} when all columns have nonzeros in at most one
#' block, so the active-set loop can keep partition-wise equality
#' re-solves. Returns \code{TRUE} when any column spans multiple blocks,
#' so the active-set loop must use the global equality bridge instead.
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




#' Build Equality Right-Hand Side for Constraint Matrix
#'
#' Converts the list-form equality targets used throughout
#' \code{get_B()} / \code{blockfit_solve()} into the vector
#' \eqn{\mathbf{c}} satisfying \eqn{\mathbf{A}^{\top}\beta = \mathbf{c}}.
#'
#' @param A Equality constraint matrix.
#' @param constraint_value_vectors List of coefficient-space vectors
#'   encoding nonzero equality targets.
#'
#' @return Numeric vector of length \code{ncol(A)}.
#'
#' @keywords internal
.solver_build_constraint_rhs <- function(A, constraint_value_vectors) {
  if (is.null(A) || !is.matrix(A) || ncol(A) == 0L) return(numeric(0))
  if (length(constraint_value_vectors) == 0L) return(rep(0, ncol(A)))

  cv_vec <- Reduce("rbind", constraint_value_vectors)
  if (length(cv_vec) == 0L) return(rep(0, ncol(A)))

  c(crossprod(A, cv_vec))
}




#' Build Augmented Equality System for an Active-Set Iteration
#'
#' Forms
#' \eqn{\mathbf{A}_{aug} = [\mathbf{A} \mid \mathbf{C}_{meq} \mid
#' \mathbf{C}_{active}]},
#' removes linearly dependent columns while preserving the original
#' column identities, and assembles the corresponding right-hand side.
#'
#' @param A Base equality constraint matrix.
#' @param constraint_value_vectors Equality right-hand-side vectors for
#'   \code{A}.
#' @param qp_Amat Full QP constraint matrix.
#' @param qp_bvec Full QP right-hand side.
#' @param qp_meq Number of leading equality columns in \code{qp_Amat}.
#' @param active_ineq Integer vector of currently active inequality
#'   indices on the original \code{qp_Amat} indexing.
#' @param rank_tol Numeric QR tolerance used for rank filtering.
#'
#' @return List with the augmented constraint matrix, reduced active set,
#'   and bookkeeping needed by the shared active-set loop.
#'
#' @keywords internal
.solver_build_augmented_constraints <- function(A,
                                                constraint_value_vectors,
                                                qp_Amat,
                                                qp_bvec,
                                                qp_meq,
                                                active_ineq = integer(0),
                                                rank_tol = 1e-10,
                                                parallel_qr = FALSE,
                                                cl = NULL) {
  n_eq_orig <- if (is.null(A)) 0L else ncol(A)
  n_qp <- if (is.null(qp_Amat)) 0L else ncol(qp_Amat)

  if (n_qp == 0L) {
    qp_Amat <- matrix(0, nrow = nrow(A), ncol = 0)
    qp_bvec <- numeric(0)
    qp_meq <- 0L
  }

  if (is.null(qp_meq)) qp_meq <- 0L
  qp_meq <- max(0L, min(as.integer(qp_meq), n_qp))

  if (qp_meq > 0L) {
    permanent_eq_cols <- seq_len(qp_meq)
    ineq_cols <- if (qp_meq < n_qp) (qp_meq + 1L):n_qp else integer(0)
  } else {
    permanent_eq_cols <- integer(0)
    ineq_cols <- if (n_qp > 0L) seq_len(n_qp) else integer(0)
  }

  if (length(permanent_eq_cols) > 0L) {
    A_base <- cbind(A, qp_Amat[, permanent_eq_cols, drop = FALSE])
  } else {
    A_base <- A
  }
  n_eq_base <- ncol(A_base)

  active_ineq <- unique(as.integer(active_ineq))
  active_ineq <- active_ineq[active_ineq %in% ineq_cols]

  if (length(active_ineq) > 0L) {
    A_aug_full <- cbind(A_base, qp_Amat[, active_ineq, drop = FALSE])
  } else {
    A_aug_full <- A_base
  }

  keep_cols <- integer(0)
  if (!is.null(A_aug_full) && ncol(A_aug_full) > 0L) {
    for (j in seq_len(ncol(A_aug_full))) {
      cand_cols <- c(keep_cols, j)
      cand_A <- A_aug_full[, cand_cols, drop = FALSE]
      cand_scales <- .constraint_col_scales(cand_A)
      cand_A_scaled <- t(t(cand_A) * cand_scales)
      if (.constraint_rank(cand_A_scaled,
                           tol = rank_tol,
                           parallel_qr = parallel_qr,
                           cl = cl) > length(keep_cols)) {
        keep_cols <- cand_cols
      }
    }
  }

  A_aug <- A_aug_full[, keep_cols, drop = FALSE]
  n_eq_aug <- sum(keep_cols <= n_eq_base)
  active_ineq_kept <- if (length(active_ineq) > 0L) {
    kept_active_pos <- keep_cols[keep_cols > n_eq_base] - n_eq_base
    active_ineq[kept_active_pos]
  } else {
    integer(0)
  }

  rhs_orig <- .solver_build_constraint_rhs(A, constraint_value_vectors)
  combined_rhs_full <- rep(0, ncol(A_aug_full))
  if (length(rhs_orig) > 0L) {
    combined_rhs_full[seq_len(n_eq_orig)] <- rhs_orig
  }
  if (length(permanent_eq_cols) > 0L) {
    combined_rhs_full[n_eq_orig + seq_along(permanent_eq_cols)] <-
      qp_bvec[permanent_eq_cols]
  }
  if (length(active_ineq) > 0L) {
    combined_rhs_full[n_eq_base + seq_along(active_ineq)] <-
      qp_bvec[active_ineq]
  }
  combined_rhs <- combined_rhs_full[keep_cols]

  if (length(combined_rhs) > 0L && any(combined_rhs != 0)) {
    b0_aug <- A_aug %**%
      (invert(crossprod(A_aug)) %**% cbind(combined_rhs))
    cv_for_proj <- list(b0_aug)
  } else {
    cv_for_proj <- list()
  }

  list(
    A_base = A_base,
    A_aug = A_aug,
    n_eq_orig = n_eq_orig,
    n_eq_base = n_eq_base,
    n_eq_aug = n_eq_aug,
    keep_cols = keep_cols,
    ineq_cols = ineq_cols,
    qp_ineq_Amat = qp_Amat[, ineq_cols, drop = FALSE],
    qp_ineq_bvec = qp_bvec[ineq_cols],
    active_ineq_kept = active_ineq_kept,
    combined_rhs = combined_rhs,
    cv_for_proj = cv_for_proj
  )
}




#' Check Primal Feasibility of Inactive Inequality Constraints
#'
#' Computes the slack
#' \eqn{\mathbf{C}^{\top}\beta - \mathbf{c}}
#' and reports which inactive inequalities are violated.
#'
#' @param beta_full Full coefficient column vector.
#' @param qp_Amat Inequality constraint matrix.
#' @param qp_bvec Inequality right-hand side.
#' @param active_ineq Integer vector of active inequality indices, using
#'   the local indexing of \code{qp_Amat}.
#' @param tol Numeric tolerance.
#'
#' @return List with \code{slack}, \code{violated}, and \code{feasible}.
#'
#' @keywords internal
.solver_check_primal_feasibility <- function(beta_full,
                                             qp_Amat,
                                             qp_bvec,
                                             active_ineq = integer(0),
                                             tol = sqrt(.Machine$double.eps)) {
  if (is.null(qp_Amat) || !is.matrix(qp_Amat) || ncol(qp_Amat) == 0L) {
    return(list(
      slack = numeric(0),
      violated = integer(0),
      feasible = TRUE
    ))
  }

  slack <- crossprod(qp_Amat, beta_full) - cbind(qp_bvec)
  violated_mask <- c(slack) < -tol
  inactive_idx <- setdiff(seq_len(ncol(qp_Amat)), active_ineq)
  violated <- inactive_idx[violated_mask[inactive_idx]]

  list(
    slack = slack,
    violated = violated,
    feasible = length(violated) == 0L
  )
}




#' Recover Active-Constraint Multipliers from Transformed OLS
#'
#' Given the transformed regression pair
#' \eqn{(\mathbf{y}^{*}, \mathbf{X}^{*})}, extracts the implied KKT
#' multipliers for the appended active inequalities.
#'
#' @param X_star Transformed constraint matrix.
#' @param y_star Transformed right-hand side.
#' @param n_eq_orig Number of leading equality constraints in
#'   \code{X_star}.
#' @param K Integer; number of interior knots.
#' @param tol Numeric tolerance.
#' @param ineq_sign Numeric sign mapping from transformed-OLS
#'   coefficients to KKT multipliers. For
#'   \eqn{\mathbf{C}^{\top}\beta \ge c}, use \code{-1}.
#'
#' @return List with \code{multipliers}, \code{drop}, and the full
#'   transformed-OLS coefficients.
#'
#' @keywords internal
.solver_recover_active_multipliers <- function(X_star,
                                               y_star,
                                               n_eq_orig,
                                               K,
                                               tol = sqrt(.Machine$double.eps),
                                               ineq_sign = -1,
                                               parallel_qr = FALSE,
                                               cl = NULL) {
  if (is.null(X_star) || !is.matrix(X_star) || ncol(X_star) <= n_eq_orig) {
    return(list(
      multipliers = numeric(0),
      drop = integer(0),
      coefficients = numeric(0)
    ))
  }

  col_scales <- .constraint_col_scales(X_star)
  X_star_scaled <- t(t(X_star) * col_scales)
  comp_stab_sc <- 1 / sqrt(K + 1)
  ols_fit <- .parallel_qr_lm_fit(
    X = X_star_scaled * comp_stab_sc,
    y = y_star * comp_stab_sc,
    parallel_qr = parallel_qr,
    cl = cl
  )

  all_coef <- ols_fit$coefficients * col_scales
  ineq_coef <- all_coef[(n_eq_orig + 1L):length(all_coef)]
  ineq_multipliers <- ineq_sign * ineq_coef

  list(
    multipliers = ineq_multipliers,
    drop = which(ineq_multipliers < -tol),
    coefficients = all_coef
  )
}




#' Build qp_info for an Active-Set Solve
#'
#' Packages the final active working set using the same fields returned
#' by the legacy active-set helpers.
#'
#' @param A Base equality constraint matrix.
#' @param qp_Amat Full QP constraint matrix.
#' @param active_ineq Final active inequality indices on the original
#'   \code{qp_Amat} indexing.
#' @param lagrangian Numeric vector of active inequality multipliers.
#' @param method Character solver label.
#'
#' @return A \code{qp_info} list.
#'
#' @keywords internal
.solver_build_active_qp_info <- function(A,
                                         qp_Amat,
                                         active_ineq,
                                         lagrangian,
                                         method) {
  if (length(active_ineq) > 0L) {
    Amat_active <- cbind(A, qp_Amat[, active_ineq, drop = FALSE])
  } else {
    Amat_active <- A
  }

  list(
    lagrangian = lagrangian,
    Amat_active = Amat_active,
    active_ineq = active_ineq,
    converged = TRUE,
    method = method
  )
}




#' Re-Solve an Equality Subproblem for a Fixed Active Set
#'
#' Uses the generic augmented-constraint builder together with a supplied
#' equality-only solver to validate and package a fixed active set.
#'
#' @param result Current coefficient list.
#' @param A Base equality constraint matrix.
#' @param constraint_value_vectors Equality right-hand-side vectors.
#' @param qp_Amat Full QP constraint matrix.
#' @param qp_bvec Full QP right-hand side.
#' @param qp_meq Number of leading equality columns in \code{qp_Amat}.
#' @param active_ineq Candidate active inequality indices on the original
#'   \code{qp_Amat} indexing.
#' @param solve_subproblem Callback taking one augmented-constraint state
#'   and returning the equality-constrained coefficient list.
#' @param kkt_subproblem Callback taking \code{(result_new, aug_state)}
#'   and returning KKT diagnostics.
#' @param method Character solver label for \code{qp_info}.
#'
#' @return Either a converged \code{list(result, qp_info, converged)} or
#'   \code{NULL} if the supplied active set does not satisfy KKT.
#'
#' @keywords internal
.solver_resolve_active_set <- function(result,
                                       A,
                                       constraint_value_vectors,
                                       qp_Amat,
                                       qp_bvec,
                                       qp_meq,
                                       active_ineq,
                                       solve_subproblem,
                                       kkt_subproblem,
                                       method,
                                       parallel_qr = FALSE,
                                       cl = NULL) {
  if (length(active_ineq) == 0L) return(NULL)

  aug_state <- .solver_build_augmented_constraints(
    A = A,
    constraint_value_vectors = constraint_value_vectors,
    qp_Amat = qp_Amat,
    qp_bvec = qp_bvec,
    qp_meq = qp_meq,
    active_ineq = active_ineq,
    parallel_qr = parallel_qr,
    cl = cl
  )

  if (length(aug_state$active_ineq_kept) == 0L) return(NULL)

  result_new <- solve_subproblem(aug_state)
  kkt <- kkt_subproblem(result_new, aug_state)
  if (!isTRUE(kkt$feasible) || !isTRUE(kkt$dual_feasible)) {
    if (isTRUE(getOption("lgspline.debug_active_set", FALSE))) {
      message(
        "[seeded-active] active=", paste(aug_state$active_ineq_kept, collapse = ","),
        " feasible=", isTRUE(kkt$feasible),
        " dual_feasible=", isTRUE(kkt$dual_feasible),
        " violated=", length(kkt$violated),
        " drop=", length(kkt$drop)
      )
    }
    return(NULL)
  }

  list(
    result = result_new,
    qp_info = .solver_build_active_qp_info(
      A = A,
      qp_Amat = qp_Amat,
      active_ineq = aug_state$active_ineq_kept,
      lagrangian = if (length(kkt$multipliers) > 0L) kkt$multipliers else numeric(0),
      method = method
    ),
    converged = TRUE
  )
}




#' Stabilize GLM Working Weights Before Linear Algebra
#'
#' Keeps iteratively reweighted least-squares weights finite and bounded
#' away from zero so the weighted Gram and Woodbury updates stay
#' numerically well posed.
#'
#' @param W Numeric working-weight vector.
#' @param min_weight Lower bound applied after replacing non-finite
#'   values.
#' @param max_weight Upper bound for extremely large finite weights.
#'
#' @return Numeric vector of stabilized working weights.
#'
#' @keywords internal
.solver_stabilize_working_weights <- function(
    W,
    min_weight = sqrt(.Machine$double.eps),
    max_weight = 1e4) {

  W <- c(W)
  if (length(W) == 0L) return(W)

  finite_pos <- W[is.finite(W) & W > 0]
  fallback_hi <- if (length(finite_pos) > 0L) max(finite_pos) else min_weight
  fallback_hi <- min(max(fallback_hi, min_weight), max_weight)

  W[is.na(W) | W <= 0] <- min_weight
  W[is.infinite(W)] <- fallback_hi
  W[!is.finite(W)] <- fallback_hi

  pmin(pmax(W, min_weight), max_weight)
}




#' Generic Outer Active-Set Loop Around Equality-Constrained Solves
#'
#' Central add/drop controller shared by the dense transformed-OLS and
#' Woodbury-transformed solvers. The inner equality solve is provided by
#' \code{solve_subproblem}; this helper only manages the working set.
#'
#' @param result Current coefficient list.
#' @param A Base equality constraint matrix.
#' @param constraint_value_vectors Equality right-hand-side vectors.
#' @param qp_Amat Full QP constraint matrix.
#' @param qp_bvec Full QP right-hand side.
#' @param qp_meq Number of leading equality columns in \code{qp_Amat}.
#' @param tol Numeric feasibility tolerance.
#' @param max_as_iter Maximum active-set iterations.
#' @param solve_subproblem Callback taking one augmented-constraint state
#'   and returning the equality-constrained coefficient list.
#' @param kkt_subproblem Callback taking \code{(result_new, aug_state)}
#'   and returning a KKT-status list with components
#'   \code{feasible}, \code{dual_feasible}, \code{violated},
#'   \code{drop}, and \code{multipliers}.
#' @param method Character solver label for \code{qp_info}.
#' @param debug_label Character tag used in optional debug tracing.
#'
#' @return List with \code{result}, \code{qp_info}, and \code{converged}.
#'
#' @keywords internal
.solver_active_set <- function(result,
                               A,
                               constraint_value_vectors,
                               qp_Amat,
                               qp_bvec,
                               qp_meq,
                               tol,
                               max_as_iter,
                               solve_subproblem,
                               kkt_subproblem,
                               method,
                               debug_label = "active-set",
                               parallel_qr = FALSE,
                               cl = NULL) {
  if (is.null(qp_Amat) || !is.matrix(qp_Amat) || ncol(qp_Amat) == 0L) {
    return(list(result = result, qp_info = NULL, converged = TRUE))
  }

  active_ineq <- integer(0)
  active_ineq_kept <- integer(0)
  converged <- FALSE
  kkt <- list(multipliers = numeric(0))
  recent_drop <- integer(0)
  visited_states <- character(0)

  for (as_iter in seq_len(max_as_iter)) {
    ## Build the equality system for the current working set.
    aug_state <- .solver_build_augmented_constraints(
      A = A,
      constraint_value_vectors = constraint_value_vectors,
      qp_Amat = qp_Amat,
      qp_bvec = qp_bvec,
      qp_meq = qp_meq,
      active_ineq = active_ineq,
      parallel_qr = parallel_qr,
      cl = cl
    )

    active_ineq <- aug_state$active_ineq_kept
    active_ineq_kept <- aug_state$active_ineq_kept
    state_key <- paste(sort(active_ineq_kept), collapse = ",")
    visited_states <- unique(c(visited_states, state_key))

    ## Re-solve the equality-only problem, then check KKT.
    result_new <- solve_subproblem(aug_state)
    kkt <- kkt_subproblem(result_new, aug_state)

    if (isTRUE(getOption("lgspline.debug_active_set", FALSE))) {
      beta_new <- cbind(unlist(result_new))
      slack_all <- if (ncol(aug_state$qp_ineq_Amat) > 0L) {
        crossprod(aug_state$qp_ineq_Amat, beta_new) -
          cbind(aug_state$qp_ineq_bvec)
      } else {
        numeric(0)
      }
      message(
        "[", debug_label, "] iter=", as_iter,
        " active=", length(active_ineq),
        " min_slack=",
        if (length(slack_all) > 0L) signif(min(c(slack_all)), 6) else NA_character_,
        " violated=", length(kkt$violated),
        " dual_drop=", length(kkt$drop),
        " dual_min=",
        if (length(kkt$multipliers) > 0L) signif(min(kkt$multipliers), 6) else NA_character_
      )
    }

    if (isTRUE(kkt$feasible) && isTRUE(kkt$dual_feasible)) {
      result <- result_new
      converged <- TRUE
      break
    }

    if (!isTRUE(kkt$dual_feasible)) {
      ## Drop the most-negative multiplier unless that would cycle.
      drop_pos <- kkt$drop[which.min(kkt$multipliers[kkt$drop])]
      drop_idx <- active_ineq_kept[drop_pos]
      tentative <- active_ineq[!(active_ineq %in% drop_idx)]
      tentative_key <- paste(sort(tentative), collapse = ",")

      if (tentative_key %in% visited_states && length(kkt$drop) > 1L) {
        ## If one-at-a-time dropping would revisit a prior state,
        # drop the whole dual-infeasible set and move on.
        drop_all_idx <- active_ineq_kept[kkt$drop]
        active_ineq <- active_ineq[!(active_ineq %in% drop_all_idx)]
        recent_drop <- drop_all_idx
      } else {
        active_ineq <- tentative
        recent_drop <- drop_idx
      }
      result <- result_new
      next
    }

    if (!isTRUE(kkt$feasible)) {
      ## Add the most-violated inactive inequality, avoiding fresh cycles
      #  when there are ties or a recent drop.
      beta_new <- cbind(unlist(result_new))
      slack_viol <- crossprod(
        aug_state$qp_ineq_Amat[, kkt$violated, drop = FALSE],
        beta_new
      ) - cbind(aug_state$qp_ineq_bvec[kkt$violated])

      viol_order <- order(c(slack_viol))
      worst_pool <- aug_state$ineq_cols[kkt$violated[viol_order]]
      worst_idx <- worst_pool[1]
      if (length(worst_pool) > 1L || length(recent_drop) > 0L) {
        candidate_pool <- worst_pool
        if (length(recent_drop) > 0L) {
          fresh_pool <- candidate_pool[!(candidate_pool %in% recent_drop)]
          if (length(fresh_pool) > 0L) {
            candidate_pool <- c(fresh_pool, setdiff(candidate_pool, fresh_pool))
          }
        }

        found_unvisited <- FALSE
        for (cand_idx in candidate_pool) {
          cand_state <- .solver_build_augmented_constraints(
            A = A,
            constraint_value_vectors = constraint_value_vectors,
            qp_Amat = qp_Amat,
            qp_bvec = qp_bvec,
            qp_meq = qp_meq,
            active_ineq = c(cand_idx, active_ineq),
            parallel_qr = parallel_qr,
            cl = cl
          )
          cand_key <- paste(sort(cand_state$active_ineq_kept), collapse = ",")
          if (!(cand_key %in% visited_states)) {
            worst_idx <- cand_idx
            found_unvisited <- TRUE
            break
          }
        }
        ## Every candidate leads to a visited state; terminate to
        # avoid an infinite cycle.
        if (!found_unvisited) break
      }
      cand_active <- c(worst_idx, active_ineq)
      cand_A <- cbind(aug_state$A_base, qp_Amat[, cand_active, drop = FALSE])
      cand_scales <- .constraint_col_scales(cand_A)
      cand_A_scaled <- t(t(cand_A) * cand_scales)
      cand_rank <- .constraint_rank(cand_A_scaled,
                                    tol = 1e-10,
                                    parallel_qr = parallel_qr,
                                    cl = cl)

      if (!(worst_idx %in% active_ineq) &&
          cand_rank <= (aug_state$n_eq_base + length(active_ineq)) &&
          length(active_ineq_kept) > 0L &&
          length(kkt$multipliers) == length(active_ineq_kept)) {
        weakest_pos <- which.min(kkt$multipliers)
        weakest_idx <- active_ineq_kept[weakest_pos]
        active_ineq <- c(
          worst_idx,
          setdiff(active_ineq, c(worst_idx, weakest_idx))
        )
      } else {
        active_ineq <- c(worst_idx, setdiff(active_ineq, worst_idx))
      }
      recent_drop <- integer(0)
      result <- result_new
    }
  }

  qp_info <- NULL
  if (converged) {
    qp_info <- .solver_build_active_qp_info(
      A = A,
      qp_Amat = qp_Amat,
      active_ineq = active_ineq_kept,
      lagrangian = if (length(kkt$multipliers) > 0L) kkt$multipliers else numeric(0),
      method = method
    )
  }

  list(result = result, qp_info = qp_info, converged = converged)
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
