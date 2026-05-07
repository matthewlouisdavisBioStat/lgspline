## process_qp.R
## Quadratic programming constraint preparation for lgspline
##
## [Change 2026-03-06] Refactored from lgspline.fit into a standalone
## helper so that (a) the main fitting function is easier to review,
## (b) the QP setup can be tested independently, and (c) per-variable
## derivative sign constraints are supported.
##
## Convention notes (for reviewers):
## - First comment line: ##   continuation: #
## - Long delimiter lines (## ## ## ##) used sparingly
## - LaTeX notation in Roxygen uses \eqn{}


## Parse qp_observations into a uniform lookup map.
#   NULL           -> no filtering
#   atomic vector  -> one shared subset for every built-in constraint
#   named list     -> per-key subsets, keyed by "var:qp_<type>" or
#                     "qp_<type>"
.parse_qp_observations <- function(qp_observations) {

  if (is.null(qp_observations)) {
    return(list(mode = "none", global = NULL, per_key = list()))
  }

  if (is.atomic(qp_observations) && !is.list(qp_observations)) {
    obs <- try(as.integer(c(qp_observations)), silent = TRUE)
    if (inherits(obs, "try-error") || any(is.na(obs))) {
      stop("\n \t qp_observations must be coercible to an integer vector, ",
           "or a named list of observation subsets. \n")
    }
    return(list(mode = "global", global = obs, per_key = list()))
  }

  if (is.list(qp_observations)) {
    keys <- names(qp_observations)
    if (is.null(keys) || any(!nzchar(keys))) {
      stop("\n \t qp_observations list entries must be named using ",
           "'var:qp_<type>' or 'qp_<type>'. \n")
    }
    per_key <- lapply(qp_observations, function(v) {
      obs <- try(as.integer(c(v)), silent = TRUE)
      if (inherits(obs, "try-error") || any(is.na(obs))) {
        stop("\n \t qp_observations list values must be coercible to ",
             "integer vectors of observation indices. \n")
      }
      obs
    })
    names(per_key) <- keys
    return(list(mode = "keyed", global = NULL, per_key = per_key))
  }

  stop("\n \t qp_observations must be NULL, an atomic vector, or a ",
       "named list. \n")
}


## Look up observation subsets for one constraint key.
#   Preference: var-specific key, then bare key, then global subset.
.lookup_qp_obs <- function(obs_map, constraint_type, var_label = NULL) {
  if (obs_map$mode == "none") return(NULL)
  if (obs_map$mode == "global") return(obs_map$global)

  if (!is.null(var_label)) {
    key_var <- paste0(var_label, ":", constraint_type)
    if (key_var %in% names(obs_map$per_key)) {
      return(obs_map$per_key[[key_var]])
    }
  }

  if (constraint_type %in% names(obs_map$per_key)) {
    return(obs_map$per_key[[constraint_type]])
  }

  NULL
}


## Look up non-variable subsets and union matching keyed entries.
#   Used for range and monotonicity constraints, where the variable prefix
#   is only a label and does not change the constraint itself.
.lookup_qp_obs_nonvar <- function(obs_map, constraint_type) {
  if (obs_map$mode == "none") return(NULL)
  if (obs_map$mode == "global") return(obs_map$global)

  nm <- names(obs_map$per_key)
  suffix <- paste0(":", constraint_type)
  hits <- obs_map$per_key[nm == constraint_type | endsWith(nm, suffix)]

  if (length(hits) == 0) return(NULL)
  sort(unique(unlist(hits)))
}


## Subset the block-diagonal design by original observation index.
.subset_X_block <- function(X_block_full, obs_order, obs) {
  if (is.null(obs)) return(X_block_full)
  X_block_full[obs_order %in% obs, , drop = FALSE]
}


## Keep only pivot inequality columns within each partition block.
#   Equality columns (the first meq) are left untouched. Inequality columns
#   spanning multiple partitions, such as monotonicity constraints, are also
#   left untouched because their feasibility set is not partition-separable.
#   For the remaining partition-local inequalities, only exact same-direction
#   duplicates / positive scalar multiples are removed. A general QR basis is
#   not cone-preserving for inequalities and can change the feasible set.
.reduce_qp_constraint_block <- function(qp_Amat,
                                        qp_bvec,
                                        qp_meq,
                                        p_expansions,
                                        K,
                                        parallel_qr_qp = FALSE,
                                        cl = NULL,
                                        tol = 1e-10) {
  if (is.null(qp_Amat) || !is.matrix(qp_Amat) || ncol(qp_Amat) == 0L) {
    return(list(Amat = qp_Amat, bvec = qp_bvec, meq = qp_meq))
  }

  n_con <- ncol(qp_Amat)
  qp_meq <- max(0L, min(as.integer(qp_meq), n_con))
  if (qp_meq >= n_con) {
    return(list(Amat = qp_Amat, bvec = qp_bvec, meq = qp_meq))
  }

  ineq_cols <- seq.int(qp_meq + 1L, n_con)
  if (length(ineq_cols) == 0L) {
    return(list(Amat = qp_Amat, bvec = qp_bvec, meq = qp_meq))
  }

  eps <- sqrt(.Machine$double.eps)
  block_id <- rep(NA_integer_, length(ineq_cols))
  keep_ineq <- integer(0)

  for (j in seq_along(ineq_cols)) {
    col_idx <- ineq_cols[[j]]
    nz <- which(abs(qp_Amat[, col_idx]) > eps)

    if (length(nz) == 0L) {
      next
    }

    blocks_hit <- unique(ceiling(nz / p_expansions))
    if (length(blocks_hit) == 1L) {
      block_id[[j]] <- blocks_hit[[1]]
    } else {
      keep_ineq <- c(keep_ineq, col_idx)
    }
  }

  block_jobs <- lapply(sort(unique(stats::na.omit(block_id))), function(k) {
    cols_k <- ineq_cols[which(block_id == k)]
    row_k <- ((k - 1L) * p_expansions + 1L):(k * p_expansions)
    list(
      cols = cols_k,
      A_local = qp_Amat[row_k, cols_k, drop = FALSE],
      b_local = qp_bvec[cols_k]
    )
  })

  reduce_job <- function(job) {
    if (length(job$cols) <= 1L) return(job$cols)

    round_digits <- max(6L, min(14L, ceiling(-log10(tol))))
    keep_cols <- integer(0)
    seen_keys <- character(0)

    for (j in seq_along(job$cols)) {
      col_vec <- c(job$A_local[, j], job$b_local[[j]])
      col_norm <- sqrt(sum(col_vec^2))
      if (!is.finite(col_norm) || col_norm <= eps) {
        next
      }

      key <- paste(round(col_vec / col_norm, round_digits), collapse = ",")
      if (!(key %in% seen_keys)) {
        seen_keys <- c(seen_keys, key)
        keep_cols <- c(keep_cols, job$cols[[j]])
      }
    }

    sort(keep_cols)
  }

  if (isTRUE(parallel_qr_qp) &&
      !is.null(cl) &&
      inherits(cl, "cluster") &&
      length(block_jobs) > 1L) {
    kept_blocks <- unlist(parallel::parLapply(cl, block_jobs, reduce_job))
  } else {
    kept_blocks <- unlist(lapply(block_jobs, reduce_job))
  }

  keep_cols <- sort(unique(c(seq_len(qp_meq), keep_ineq, kept_blocks)))

  list(
    Amat = qp_Amat[, keep_cols, drop = FALSE],
    bvec = qp_bvec[keep_cols],
    meq = qp_meq
  )
}


## Resolve keyed subsets for one derivative variable.
#   The derivative list may use either the original variable names or the
#   internal "_j_" labels, so try both forms.
.get_obs_for_var <- function(obs_map, constraint_type, var_nm, og_cols) {
  if (is.null(obs_map) || obs_map$mode != "keyed") return(NULL)

  obs_var <- .lookup_qp_obs(obs_map, constraint_type, var_nm)
  if (!is.null(obs_var)) return(obs_var)

  if (is.null(og_cols)) return(NULL)

  if (grepl("^_[0-9]+_$", var_nm)) {
    j <- suppressWarnings(as.integer(gsub("_", "", var_nm)))
    if (!is.na(j) && j >= 1L && j <= length(og_cols)) {
      obs_var <- .lookup_qp_obs(obs_map, constraint_type, og_cols[j])
      if (!is.null(obs_var)) return(obs_var)
    }
  } else if (var_nm %in% og_cols) {
    j <- which(og_cols == var_nm)
    if (length(j) == 1L) {
      obs_var <- .lookup_qp_obs(
        obs_map,
        constraint_type,
        paste0("_", j, "_")
      )
      if (!is.null(obs_var)) return(obs_var)
    }
  }

  NULL
}


#' Build Derivative QP Constraints in Full \eqn{P}-Dimensional Space
#'
#' @description
#' Internal helper that constructs the \code{Amat} / \code{bvec} pair
#' enforcing derivative sign constraints at every row of a block-diagonal
#' design matrix.
#'
#' Given an \eqn{N_{\mathrm{sub}} \times P} block-diagonal design matrix
#' \code{X_block}, where \eqn{P = p \times (K+1)}, this function:
#' \enumerate{
#'   \item Recovers the per-partition \eqn{p}-column expansion matrix
#'     \code{C_qp} and records each row's partition assignment.
#'   \item Calls \code{make_derivative_matrix} on \code{C_qp} to obtain
#'     first or second derivative matrices with respect to each predictor.
#'   \item Optionally selects only derivatives for a subset of predictor
#'     variables (\code{target_vars}).
#'   \item Maps each derivative row into the full \eqn{P}-dimensional
#'     coefficient space, yielding one constraint column per (observation,
#'     variable) pair.
#' }
#'
#' The result is a constraint pair \eqn{\mathbf{A}^{\top}\boldsymbol{\beta}
#' \ge \mathbf{b}} (with \eqn{\mathbf{b} = \mathbf{0}}) suitable for
#' \code{\link[quadprog]{solve.QP}}.
#'
#' @param X_block Numeric matrix, \eqn{N_{\mathrm{sub}} \times P}.
#'   Block-diagonal design matrix for the QP-selected observations.
#' @param sign_mult Numeric scalar, \code{+1} for positive constraints
#'   or \code{-1} for negative constraints.
#' @param just_first Logical. If \code{TRUE}, constrain first derivatives;
#'   if \code{FALSE}, constrain second derivatives.
#' @param p_expansions Integer. Number of basis expansions per partition.
#' @param K Integer. Number of interior knots (partitions minus 1).
#' @param colnm_expansions Character vector of length \code{p_expansions}.
#'   Column names of the expansion matrix.
#' @param power1_cols Integer vector of linear-term column indices.
#' @param power2_cols Integer vector of quadratic-term column indices.
#' @param nonspline_cols Integer vector of non-spline linear column indices.
#' @param interaction_single_cols Integer vector of linear-by-linear
#'   interaction column indices.
#' @param interaction_quad_cols Integer vector of linear-by-quadratic
#'   interaction column indices.
#' @param triplet_cols Integer vector of three-way interaction column indices.
#' @param include_2way_interactions Logical switch forwarded to
#'   \code{make_derivative_matrix}.
#' @param include_3way_interactions Logical switch forwarded to
#'   \code{make_derivative_matrix}.
#' @param include_quadratic_interactions Logical switch forwarded to
#'   \code{make_derivative_matrix}.
#' @param expansion_scales Numeric vector of length \code{p_expansions - 1}.
#'   Passed as unit scales (\code{1 + 0 * expansion_scales}) so that the
#'   constraint matrix stays on the expansion-standardized coefficient scale.
#' @param target_vars Optional. Integer vector of predictor column indices
#'   or character vector of predictor names identifying which predictors
#'   to constrain. When \code{NULL} (default), all derivative variables
#'   are constrained.
#' @param og_cols Optional character vector of original predictor column
#'   names, used to resolve character \code{target_vars} to integer indices.
#' @param obs_map Optional parsed \code{qp_observations} map. In keyed mode,
#'   derivative rows are filtered per variable using
#'   \code{constraint_type} and \code{obs_order}.
#' @param constraint_type Optional character scalar naming the current
#'   derivative constraint, for example \code{"qp_negative_derivative"}.
#' @param obs_order Optional integer vector of original observation indices,
#'   one per row of \code{X_block}.
#'
#' @return A list with components:
#' \describe{
#'   \item{Amat}{\eqn{P \times M} constraint matrix, where \eqn{M} is the
#'     number of unique constraint columns (after deduplication).}
#'   \item{bvec}{Numeric vector of length \eqn{M}, all zeros.}
#'   \item{meq}{Integer, always \code{0} (inequality constraints).}
#' }
#'
#' @keywords internal
#' @seealso \code{\link{process_qp}}, \code{\link[quadprog]{solve.QP}}
.build_deriv_qp <- function(X_block,
                            sign_mult,
                            just_first,
                            p_expansions,
                            K,
                            colnm_expansions,
                            power1_cols,
                            power2_cols,
                            nonspline_cols,
                            interaction_single_cols,
                            interaction_quad_cols,
                            triplet_cols,
                            include_2way_interactions,
                            include_3way_interactions,
                            include_quadratic_interactions,
                            expansion_scales,
                            target_vars = NULL,
                            og_cols = NULL,
                            obs_map = NULL,
                            constraint_type = NULL,
                            obs_order = NULL) {

  N_sub <- nrow(X_block)
  if (N_sub == 0) {
    return(list(
      Amat = matrix(0, nrow = p_expansions * (K + 1), ncol = 0),
      bvec = numeric(0),
      meq = 0L
    ))
  }

  ## Recover the p-column expansion matrix from the block-diagonal design.
  # Each row of X_block belongs to exactly one partition k, and its
  # p_expansions nonzero entries sit in columns [(k-1)*p + 1] : [k*p].
  C_qp <- matrix(0, N_sub, p_expansions)
  colnames(C_qp) <- colnm_expansions
  partition_assignment <- integer(N_sub)
  for (i in seq_len(N_sub)) {
    for (k in 1:(K + 1)) {
      cols_k <- ((k - 1) * p_expansions + 1):(k * p_expansions)
      if (any(X_block[i, cols_k] != 0)) {
        C_qp[i, ] <- X_block[i, cols_k]
        partition_assignment[i] <- k
        break
      }
    }
  }

  ## Compute p-column derivative matrices on the standardized scale.
  # Unit scales ensure the constraint matrix matches the coefficient scale.
  derivs <- make_derivative_matrix(
    p_expansions, C_qp, power1_cols, power2_cols, nonspline_cols,
    interaction_single_cols, interaction_quad_cols, triplet_cols,
    K, include_2way_interactions, include_3way_interactions,
    include_quadratic_interactions, colnm_expansions,
    1 + 0 * expansion_scales,
    just_first_derivatives = just_first,
    just_spline_effects = FALSE
  )

  ## Select first or second derivative list
  if (just_first) {
    deriv_list <- derivs$first_derivatives
  } else {
    deriv_list <- derivs$second_derivatives
  }

  ## Resolve target_vars to integer indices into deriv_list.
  # deriv_list is named by predictor; names(deriv_list) holds the
  # expansion-name keys (e.g. "_1_", "_2_", or "Girth", "Height").
  if (!is.null(target_vars)) {
    deriv_names <- names(deriv_list)

    ## Build a permissive set of candidate labels so character targets
    # can match either the original predictor names (e.g. "Dose") or
    # the internal expansion names (e.g. "_1_").
    candidate_targets <- character(0)

    if (is.character(target_vars)) {
      candidate_targets <- unique(c(candidate_targets, target_vars))

      if (!is.null(og_cols)) {
        target_numeric <- unique(unlist(lapply(target_vars, function(v) {
          idx <- which(og_cols == v)
          if (length(idx) == 0) idx <- grep(v, og_cols, fixed = TRUE)
          idx
        })))

        if (length(target_numeric) > 0) {
          candidate_targets <- unique(c(
            candidate_targets,
            paste0("_", target_numeric, "_"),
            og_cols[target_numeric]
          ))
        }
      }

    } else if (is.numeric(target_vars)) {
      target_numeric <- unique(as.integer(target_vars))
      candidate_targets <- unique(c(
        candidate_targets,
        paste0("_", target_numeric, "_")
      ))

      if (!is.null(og_cols)) {
        target_numeric <- target_numeric[
          target_numeric >= 1 & target_numeric <= length(og_cols)
        ]
        if (length(target_numeric) > 0) {
          candidate_targets <- unique(c(
            candidate_targets,
            og_cols[target_numeric]
          ))
        }
      }

    } else {
      ## Fallback: treat as-is
      candidate_targets <- unique(c(candidate_targets, target_vars))
    }

    ## Keep only derivative components whose name matches one of the
    # candidate labels assembled above.
    keep_idx <- which(sapply(deriv_names, function(nm) {
      any(sapply(candidate_targets, function(tv) {
        if (is.na(tv) || !nzchar(tv)) return(FALSE)
        identical(nm, tv) || grepl(tv, nm, fixed = TRUE)
      }))
    }))

    if (length(keep_idx) == 0) {
      ## No matching variables: return empty constraint
      return(list(
        Amat = matrix(0, nrow = p_expansions * (K + 1), ncol = 0),
        bvec = numeric(0),
        meq = 0L
      ))
    }
    deriv_list <- deriv_list[keep_idx]
  }

  ## Map each row of each derivative matrix into P-dimensional space.
  # One constraint column per (observation, variable) pair.
  constraint_cols <- list()
  deriv_names_kept <- names(deriv_list)
  for (v in seq_along(deriv_list)) {
    deriv_mat <- deriv_list[[v]]
    var_nm <- deriv_names_kept[[v]]

    row_mask <- rep(TRUE, nrow(deriv_mat))
    if (!is.null(obs_map) &&
        !is.null(constraint_type) &&
        !is.null(obs_order)) {
      obs_var <- .get_obs_for_var(obs_map, constraint_type, var_nm, og_cols)
      if (!is.null(obs_var)) {
        row_mask <- obs_order %in% obs_var
      }
    }

    for (i in which(row_mask)) {
      k <- partition_assignment[i]
      col_P <- rep(0, p_expansions * (K + 1))
      col_P[((k - 1) * p_expansions + 1):(k * p_expansions)] <-
        sign_mult * deriv_mat[i, ]
      constraint_cols[[length(constraint_cols) + 1]] <- col_P
    }
  }

  if (length(constraint_cols) == 0) {
    return(list(
      Amat = matrix(0, nrow = p_expansions * (K + 1), ncol = 0),
      bvec = numeric(0),
      meq = 0L
    ))
  }

  qp_Amat <- do.call(cbind, constraint_cols)
  qp_Amat <- t(unique(t(qp_Amat)))
  qp_bvec <- rep(0, ncol(qp_Amat))

  list(Amat = qp_Amat, bvec = qp_bvec, meq = 0L)
}


#' Prepare Quadratic Programming Constraints for lgspline
#'
#' @description
#' Builds the linear inequality system used for shape-restricted fitting.
#' The helper collects built-in range, derivative-sign, second-derivative,
#' and monotonicity restrictions, together with any user-supplied custom
#' constraint functions, and returns the resulting
#' \eqn{\mathbf{C}^{\top}\boldsymbol{\beta} \ge \mathbf{c}} objects.
#'
#' This logic was refactored out of \code{\link{lgspline.fit}} so the
#' constraint construction can be reviewed and tested on its own. Existing
#' calls that pass \code{TRUE}/\code{FALSE} for derivative flags remain
#' backward-compatible.
#'
#' @section Per-Variable Derivative Constraints:
#' The arguments \code{qp_positive_derivative}, \code{qp_negative_derivative},
#' \code{qp_positive_2ndderivative}, and \code{qp_negative_2ndderivative}
#' now accept three forms:
#' \describe{
#'   \item{\code{FALSE}}{No constraint (default).}
#'   \item{\code{TRUE}}{Constrain \strong{all} predictor variables (backward
#'     compatible with previous behavior).}
#'   \item{Character or integer vector}{Constrain only the specified
#'     predictor variables. Character entries are resolved via
#'     \code{og_cols}; integer entries refer to column positions in the
#'     predictor matrix.}
#' }
#'
#' This allows, for example, enforcing a nonnegative first derivative for
#' \code{"Dose"} and a nonpositive first derivative for \code{"Time"}
#' simultaneously:
#' \preformatted{
#'   lgspline(...,
#'            qp_positive_derivative = "Dose",
#'            qp_negative_derivative = "Time")
#' }
#'
#' The arguments \code{qp_monotonic_increase} and
#' \code{qp_monotonic_decrease} remain \code{TRUE}/\code{FALSE} only,
#' because they constrain fitted values in observation order (not
#' per-variable).
#'
#' @param X List of per-partition design matrices.
#' @param K Integer. Number of interior knots.
#' @param p_expansions Integer. Number of basis expansions per partition.
#' @param order_list List of per-partition observation index vectors.
#' @param colnm_expansions Character vector of expansion column names.
#' @param expansion_scales Numeric vector of expansion standardization scales.
#' @param power1_cols Integer vector of linear-term column indices.
#' @param power2_cols Integer vector of quadratic-term column indices.
#' @param nonspline_cols Integer vector of non-spline linear column indices.
#' @param interaction_single_cols Integer vector of linear-by-linear
#'   interaction column indices.
#' @param interaction_quad_cols Integer vector of linear-by-quadratic
#'   interaction column indices.
#' @param triplet_cols Integer vector of three-way interaction column indices.
#' @param include_2way_interactions Logical switch controlling whether
#'   two-way interactions are included in derivative construction.
#' @param include_3way_interactions Logical switch controlling whether
#'   three-way interactions are included in derivative construction.
#' @param include_quadratic_interactions Logical switch controlling whether
#'   quadratic interactions are included in derivative construction.
#' @param family GLM family object.
#' @param mean_y,sd_y Numeric scalars for response standardization.
#' @param N_obs Integer. Total sample size.
#' @param qp_observations Optional observation subset. Supply either a single
#'   integer vector to thin every active built-in QP constraint the same way,
#'   or a named list keyed by \code{"var:qp_<type>"} (or bare
#'   \code{"qp_<type>"}) to give different constraints different subsets.
#'   Known types are \code{qp_range_lower}, \code{qp_range_upper},
#'   \code{qp_positive_derivative}, \code{qp_negative_derivative},
#'   \code{qp_positive_2ndderivative}, \code{qp_negative_2ndderivative},
#'   \code{qp_monotonic_increase}, and \code{qp_monotonic_decrease}. For
#'   range and monotonicity, use the bare keys such as
#'   \code{"qp_range_lower"} and \code{"qp_monotonic_increase"} because
#'   those constraints are variable-independent; prefixed versions are still
#'   accepted and unioned internally.
#'   Unknown keyed entries are ignored with a warning when
#'   \code{include_warnings = TRUE}.
#' @param qp_positive_derivative,qp_negative_derivative Logical scalar,
#'   character vector, or integer vector. See section
#'   \emph{Per-Variable Derivative Constraints}.
#' @param qp_positive_2ndderivative,qp_negative_2ndderivative Same as above
#'   but for second derivatives.
#' @param qp_monotonic_increase,qp_monotonic_decrease Logical.
#'   Constrain fitted values to be monotonic in observation order.
#' @param qp_range_upper,qp_range_lower Optional numeric upper/lower bounds
#'   for fitted values.
#' @param qp_Amat_fxn,qp_bvec_fxn,qp_meq_fxn Optional user-supplied
#'   constraint-generating functions.
#' @param qp_Amat,qp_bvec,qp_meq Optional pre-built QP objects. Their presence
#'   marks QP handling as active, but this helper does not append them to the
#'   built-in constraints it constructs; they are expected to be handled
#'   outside this constructor.
#' @param qr_pivot_inequality_constraints Logical. If \code{TRUE},
#'   partition-local inequality columns are thinned to QR pivot columns
#'   before the final \code{solve.QP} objects are stacked. Leading equality
#'   columns are left untouched; built-in range and monotonicity blocks are
#'   also left untouched; and columns spanning multiple partitions are left
#'   untouched.
#' @param parallel_qr_qp Logical. If \code{TRUE} and a cluster is supplied
#'   in \code{cl}, the partition-local QR pivot steps used by
#'   \code{qr_pivot_inequality_constraints} are parallelized across
#'   partitions.
#' @param all_derivatives_fxn Function to compute derivatives from expansion
#'   matrices (the \code{all_derivatives} closure from \code{lgspline.fit}).
#' @param og_cols Optional character vector of original predictor column names.
#' @param cl Optional parallel cluster object used only when
#'   \code{parallel_qr_qp = TRUE}.
#' @param include_warnings Logical. Whether to issue warnings.
#' @param ... Additional arguments forwarded to custom constraint functions.
#'
#' @return A list with components:
#' \describe{
#'   \item{qp_Amat}{\eqn{P \times M} combined constraint matrix.}
#'   \item{qp_bvec}{Numeric vector of length \eqn{M}.}
#'   \item{qp_meq}{Integer. Number of leading equality constraints.}
#'   \item{quadprog}{Logical. \code{TRUE} if any QP constraints are active.}
#' }
#'
#' @examples
#' \dontrun{
#' ## Standalone verification: simple 1-D monotonic increase
#' set.seed(1234)
#' t <- seq(-5, 5, length.out = 200)
#' y <- 3 * sin(t) + t + rnorm(200, 0, 0.5)
#'
#' ## Fit with positive first-derivative constraint on all variables
#' fit1 <- lgspline(t, y, K = 3,
#'                  qp_positive_derivative = TRUE)
#'
#' ## Verify: first derivative should be >= 0 everywhere
#' derivs1 <- predict(fit1, new_predictors = sort(t),
#'                    take_first_derivatives = TRUE)
#' stopifnot(all(derivs1$first_deriv >= -1e-8))
#'
#' ## Fit with monotonic increase (observation-order)
#' fit2 <- lgspline(t, y, K = 3,
#'                  qp_monotonic_increase = TRUE)
#' preds2 <- predict(fit2, new_predictors = sort(t))
#' stopifnot(all(diff(preds2) >= -1e-8))
#'
#' ## Per-variable constraints: 2-D example
#' t1 <- runif(500, -5, 5)
#' t2 <- runif(500, -5, 5)
#' y2 <- t1 + sin(t2) + rnorm(500, 0, 0.5)
#' dat2 <- data.frame(t1 = t1, t2 = t2, y = y2)
#'
#' ## Constrain t1 to have positive derivative, t2 to have negative
#' fit3 <- lgspline(y ~ spl(t1, t2), data = dat2, K = 2,
#'                  qp_positive_derivative = "t1",
#'                  qp_negative_derivative = "t2")
#'
#' ## Verify per-variable derivatives
#' newdat <- expand.grid(t1 = seq(-4, 4, length.out = 50),
#'                       t2 = seq(-4, 4, length.out = 50))
#' d3 <- predict(fit3, new_predictors = newdat,
#'               take_first_derivatives = TRUE)
#'
#' ## t1 derivative should be >= 0
#' stopifnot(all(unlist(d3$first_deriv[["_1_"]]) >= -1e-6))
#' ## t2 derivative should be <= 0
#' stopifnot(all(unlist(d3$first_deriv[["_2_"]]) <= 1e-6))
#'
#' ## Numeric column indices work identically
#' fit4 <- lgspline(y ~ spl(t1, t2), data = dat2, K = 2,
#'                  qp_positive_derivative = 1,
#'                  qp_negative_derivative = 2)
#'
#' ## Range + derivative constraints simultaneously
#' fit5 <- lgspline(t, y, K = 3,
#'                  qp_positive_derivative = TRUE,
#'                  qp_range_lower = -5,
#'                  qp_range_upper = 15)
#' preds5 <- predict(fit5)
#' stopifnot(all(preds5 >= -5 - 1e-6))
#' stopifnot(all(preds5 <= 15 + 1e-6))
#' }
#'
#' @export
process_qp <- function(X,
                       K,
                       p_expansions,
                       order_list,
                       colnm_expansions,
                       expansion_scales,
                       power1_cols,
                       power2_cols,
                       nonspline_cols,
                       interaction_single_cols,
                       interaction_quad_cols,
                       triplet_cols,
                       include_2way_interactions,
                       include_3way_interactions,
                       include_quadratic_interactions,
                       family,
                       mean_y,
                       sd_y,
                       N_obs,
                       qp_observations = NULL,
                       qp_positive_derivative = FALSE,
                       qp_negative_derivative = FALSE,
                       qp_positive_2ndderivative = FALSE,
                       qp_negative_2ndderivative = FALSE,
                       qp_monotonic_increase = FALSE,
                       qp_monotonic_decrease = FALSE,
                       qp_range_upper = NULL,
                       qp_range_lower = NULL,
                       qp_Amat_fxn = NULL,
                       qp_bvec_fxn = NULL,
                       qp_meq_fxn = NULL,
                       qp_Amat = NULL,
                       qp_bvec = NULL,
                       qp_meq = 0,
                       qr_pivot_inequality_constraints = FALSE,
                       parallel_qr_qp = FALSE,
                       all_derivatives_fxn = NULL,
                       og_cols = NULL,
                       cl = NULL,
                       include_warnings = TRUE,
                       ...) {

  ## Determine whether any QP constraints are active.
  # The derivative flags can now be TRUE, FALSE, or a vector of
  # variable names/indices, so we check with .is_active().
  .is_active <- function(x) {
    if (is.logical(x)) return(any(x))
    if (is.numeric(x) || is.character(x)) return(length(x) > 0)
    FALSE
  }

  quadprog <- (
    .is_active(qp_negative_derivative) |
    .is_active(qp_monotonic_decrease) |
    .is_active(qp_negative_2ndderivative) |
    .is_active(qp_positive_derivative) |
    .is_active(qp_monotonic_increase) |
    .is_active(qp_positive_2ndderivative) |
    any(!is.null(qp_range_upper)) |
    any(!is.null(qp_range_lower)) |
    (any(!is.null(qp_Amat_fxn)) &
     any(!is.null(qp_bvec_fxn)) &
     any(!is.null(qp_meq_fxn))) |
    (any(!is.null(qp_Amat)) &
     any(!is.null(qp_bvec)) &
     any(!is.null(qp_meq)))
  )

  if (!quadprog) {
    return(list(
      qp_Amat = NULL,
      qp_bvec = NULL,
      qp_meq = 0L,
      quadprog = FALSE
    ))
  }

  ## Accumulate constraint triples from each source
  qp_Amat_list <- list()
  qp_bvec_list <- list()
  qp_meq_list <- list()
  qp_type_list <- list()

  ## Parse qp_observations once so each constraint block can pull the rows
  #  it needs while keeping the legacy shared-subset path intact.
  obs_map <- .parse_qp_observations(qp_observations)

  ## Warn on unknown keyed entries so typos do not silently do nothing.
  if (obs_map$mode == "keyed" && include_warnings) {
    known_qp_types <- c(
      "qp_range_lower",
      "qp_range_upper",
      "qp_positive_derivative",
      "qp_negative_derivative",
      "qp_positive_2ndderivative",
      "qp_negative_2ndderivative",
      "qp_monotonic_increase",
      "qp_monotonic_decrease"
    )

    for (key in names(obs_map$per_key)) {
      key_parts <- strsplit(key, ":", fixed = TRUE)[[1]]
      key_type <- key_parts[[length(key_parts)]]
      if (!(key_type %in% known_qp_types)) {
        warning("\n\t qp_observations key '", key,
                "' is not a known constraint type and will be ignored.\n")
      }
    }
  }

  ## Build the block-diagonal design matrix for QP observations.
  # This is N_full x P where P = p_expansions * (K+1).
  X_block_full <- Reduce("rbind", lapply(1:(K + 1), function(k) {
    if (nrow(X[[k]]) == 0) return(X[[k]])
    Reduce("cbind", lapply(1:(K + 1), function(j) {
      if (j == k) X[[k]] else 0 * X[[k]]
    }))
  }))

  obs_order_full <- unlist(order_list)

  if (obs_map$mode == "global") {
    X_block <- .subset_X_block(X_block_full, obs_order_full, obs_map$global)
    obs_order <- obs_order_full[obs_order_full %in% obs_map$global]
  } else {
    X_block <- X_block_full
    obs_order <- obs_order_full
  }

  ## Resolve per-variable derivative flag to target_vars (NULL = all).
  # If TRUE, target_vars = NULL (constrain all); if FALSE, skip;
  # if character/numeric, pass through as target_vars.
  .resolve_target <- function(flag) {
    if (is.logical(flag)) {
      if (any(flag)) return(NULL)  # NULL = all variables
      else return(NA)              # NA = skip
    }
    ## Character or numeric: specific variables
    flag
  }

  ## Range constraints
  if (!any(is.null(qp_range_upper)) | !any(is.null(qp_range_lower))) {
    if (paste0(family)[2] != 'identity') {
      if (any(!is.null(qp_range_upper))) {
        qp_range_upper <- family$linkfun(qp_range_upper)
      }
      if (any(!is.null(qp_range_lower))) {
        qp_range_lower <- family$linkfun(qp_range_lower)
      }
      if (any(is.na(c(qp_range_upper, qp_range_lower))) |
          any(!is.finite(c(qp_range_upper, qp_range_lower)))) {
        stop('\n\tQuadratic programming upper/lower range constraints ',
             'are NA or not finite after link function transformation. ',
             'Supply bounds on raw-response scale, which are in-range of ',
             'link function transformations.\n')
      }
    }

    if (obs_map$mode == "keyed") {
      obs_lo <- .lookup_qp_obs_nonvar(obs_map, "qp_range_lower")
      obs_up <- .lookup_qp_obs_nonvar(obs_map, "qp_range_upper")
      X_block_lo <- .subset_X_block(X_block_full, obs_order_full, obs_lo)
      X_block_up <- .subset_X_block(X_block_full, obs_order_full, obs_up)

      if (!any(is.null(qp_range_lower)) && nrow(X_block_lo) > 0) {
        range_Amat <- t(X_block_lo)
        if (length(qp_range_lower) == 1) {
          range_Amat <- t(unique(t(range_Amat)))
          range_bvec <- rep(qp_range_lower, ncol(range_Amat))
        } else {
          range_bvec <- qp_range_lower
        }
        range_bvec <- (range_bvec - mean_y) / sd_y
        qp_Amat_list[[length(qp_Amat_list) + 1]] <- range_Amat
        qp_bvec_list[[length(qp_bvec_list) + 1]] <- range_bvec
        qp_meq_list[[length(qp_meq_list) + 1]] <- 0
        qp_type_list[[length(qp_type_list) + 1]] <- "qp_range_lower"
      }

      if (!any(is.null(qp_range_upper)) && nrow(X_block_up) > 0) {
        range_Amat <- -t(X_block_up)
        if (length(qp_range_upper) == 1) {
          range_Amat <- t(unique(t(range_Amat)))
          range_bvec <- rep(qp_range_upper, ncol(range_Amat))
        } else {
          range_bvec <- qp_range_upper
        }
        range_bvec <- -(range_bvec - mean_y) / sd_y
        qp_Amat_list[[length(qp_Amat_list) + 1]] <- range_Amat
        qp_bvec_list[[length(qp_bvec_list) + 1]] <- range_bvec
        qp_meq_list[[length(qp_meq_list) + 1]] <- 0
        qp_type_list[[length(qp_type_list) + 1]] <- "qp_range_upper"
      }

    } else if (!any(is.null(qp_range_upper)) &&
               !any(is.null(qp_range_lower)) &&
               nrow(X_block) > 0) {
      range_Amat <- cbind(t(X_block), -t(X_block))
      if (length(qp_range_lower) == 1) {
        if (length(qp_range_upper) == 1) {
          range_Amat <- t(unique(t(range_Amat)))
        }
        range_bvec_lower <- rep((qp_range_lower - mean_y) / sd_y,
                                ncol(range_Amat) / 2)
      } else {
        range_bvec_lower <- (qp_range_lower - mean_y) / sd_y
      }
      if (length(qp_range_upper) == 1) {
        range_bvec_upper <- rep(-(qp_range_upper - mean_y) / sd_y,
                                ncol(range_Amat) / 2)
      } else {
        range_bvec_upper <- -(qp_range_upper - mean_y) / sd_y
      }
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- range_Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- c(range_bvec_lower,
                                                     range_bvec_upper)
      qp_meq_list[[length(qp_meq_list) + 1]] <- 0
      qp_type_list[[length(qp_type_list) + 1]] <- "qp_range_both"

    } else if (!any(is.null(qp_range_upper)) && nrow(X_block) > 0) {
      range_Amat <- -t(X_block)
      if (length(qp_range_upper) == 1) {
        range_Amat <- t(unique(t(range_Amat)))
        range_bvec <- rep(qp_range_upper, ncol(range_Amat))
      } else {
        range_bvec <- qp_range_upper
      }
      range_bvec <- -(range_bvec - mean_y) / sd_y
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- range_Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- range_bvec
      qp_meq_list[[length(qp_meq_list) + 1]] <- 0
      qp_type_list[[length(qp_type_list) + 1]] <- "qp_range_upper"

    } else if (!any(is.null(qp_range_lower)) && nrow(X_block) > 0) {
      range_Amat <- t(X_block)
      if (length(qp_range_lower) == 1) {
        range_Amat <- t(unique(t(range_Amat)))
        range_bvec <- rep(qp_range_lower, ncol(range_Amat))
      } else {
        range_bvec <- qp_range_lower
      }
      range_bvec <- (range_bvec - mean_y) / sd_y
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- range_Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- range_bvec
      qp_meq_list[[length(qp_meq_list) + 1]] <- 0
      qp_type_list[[length(qp_type_list) + 1]] <- "qp_range_lower"
    }
  }

  ## Bundle the common derivative-matrix inputs once so each sign / order
  #  case only needs to specify what is changing.
  deriv_shared <- list(
    X_block = X_block,
    p_expansions = p_expansions,
    K = K,
    colnm_expansions = colnm_expansions,
    power1_cols = power1_cols,
    power2_cols = power2_cols,
    nonspline_cols = nonspline_cols,
    interaction_single_cols = interaction_single_cols,
    interaction_quad_cols = interaction_quad_cols,
    triplet_cols = triplet_cols,
    include_2way_interactions = include_2way_interactions,
    include_3way_interactions = include_3way_interactions,
    include_quadratic_interactions = include_quadratic_interactions,
    expansion_scales = expansion_scales,
    og_cols = og_cols
  )

  ## Positive first derivative
  target_pos1 <- .resolve_target(qp_positive_derivative)
  if (!identical(target_pos1, NA)) {
    args <- c(deriv_shared, list(
      sign_mult = +1,
      just_first = TRUE,
      target_vars = target_pos1,
      obs_map = obs_map,
      constraint_type = "qp_positive_derivative",
      obs_order = obs_order
    ))
    deriv_qp <- do.call(.build_deriv_qp, args)
    if (ncol(deriv_qp$Amat) > 0) {
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- deriv_qp$Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- deriv_qp$bvec
      qp_meq_list[[length(qp_meq_list) + 1]] <- deriv_qp$meq
      qp_type_list[[length(qp_type_list) + 1]] <- "qp_positive_derivative"
    }
  }

  ## Negative first derivative
  target_neg1 <- .resolve_target(qp_negative_derivative)
  if (!identical(target_neg1, NA)) {
    args <- c(deriv_shared, list(
      sign_mult = -1,
      just_first = TRUE,
      target_vars = target_neg1,
      obs_map = obs_map,
      constraint_type = "qp_negative_derivative",
      obs_order = obs_order
    ))
    deriv_qp <- do.call(.build_deriv_qp, args)
    if (ncol(deriv_qp$Amat) > 0) {
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- deriv_qp$Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- deriv_qp$bvec
      qp_meq_list[[length(qp_meq_list) + 1]] <- deriv_qp$meq
      qp_type_list[[length(qp_type_list) + 1]] <- "qp_negative_derivative"
    }
  }

  ## Positive second derivative
  target_pos2 <- .resolve_target(qp_positive_2ndderivative)
  if (!identical(target_pos2, NA)) {
    args <- c(deriv_shared, list(
      sign_mult = +1,
      just_first = FALSE,
      target_vars = target_pos2,
      obs_map = obs_map,
      constraint_type = "qp_positive_2ndderivative",
      obs_order = obs_order
    ))
    deriv_qp <- do.call(.build_deriv_qp, args)
    if (ncol(deriv_qp$Amat) > 0) {
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- deriv_qp$Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- deriv_qp$bvec
      qp_meq_list[[length(qp_meq_list) + 1]] <- deriv_qp$meq
      qp_type_list[[length(qp_type_list) + 1]] <- "qp_positive_2ndderivative"
    }
  }

  ## Negative second derivative
  target_neg2 <- .resolve_target(qp_negative_2ndderivative)
  if (!identical(target_neg2, NA)) {
    args <- c(deriv_shared, list(
      sign_mult = -1,
      just_first = FALSE,
      target_vars = target_neg2,
      obs_map = obs_map,
      constraint_type = "qp_negative_2ndderivative",
      obs_order = obs_order
    ))
    deriv_qp <- do.call(.build_deriv_qp, args)
    if (ncol(deriv_qp$Amat) > 0) {
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- deriv_qp$Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- deriv_qp$bvec
      qp_meq_list[[length(qp_meq_list) + 1]] <- deriv_qp$meq
      qp_type_list[[length(qp_type_list) + 1]] <- "qp_negative_2ndderivative"
    }
  }

  ## Monotonic increase (observation order, not per-variable)
  if (qp_monotonic_increase || qp_monotonic_decrease) {
    sign_mono <- if (qp_monotonic_increase) +1 else -1
    mono_type <- if (qp_monotonic_increase) {
      "qp_monotonic_increase"
    } else {
      "qp_monotonic_decrease"
    }

    X_block_mono <- X_block
    obs_order_mono <- obs_order
    if (obs_map$mode == "keyed") {
      obs_mono <- .lookup_qp_obs_nonvar(obs_map, mono_type)
      X_block_mono <- .subset_X_block(X_block_full, obs_order_full, obs_mono)
      if (!is.null(obs_mono)) {
        obs_order_mono <- obs_order_full[obs_order_full %in% obs_mono]
      }
    }

    inv_order <- order(obs_order_mono)
    X_block_orig <- X_block_mono[inv_order, , drop = FALSE]

    if (nrow(X_block_orig) >= 2) {
      value_constraints <- t(Reduce(
        'rbind',
        lapply(2:nrow(X_block_orig), function(i) {
          matrix(sign_mono * c(X_block_orig[i, ] - X_block_orig[i - 1, ]),
                 nrow = 1)
        })
      ))
      mono_Amat <- cbind(value_constraints)
      mono_Amat <- t(unique(t(mono_Amat)))
      mono_bvec <- rep(0, ncol(mono_Amat))
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- mono_Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- mono_bvec
      qp_meq_list[[length(qp_meq_list) + 1]] <- 0
      qp_type_list[[length(qp_type_list) + 1]] <- mono_type
    }
  }

  ## User-supplied constraint builders are evaluated last and then appended
  #  to the built-in derivative / range / monotonicity constraints.
  if (!is.null(qp_Amat_fxn) &
      !is.null(qp_bvec_fxn) &
      !is.null(qp_meq_fxn)) {

    custom_Amat <- qp_Amat_fxn(
      N = N_obs, p = p_expansions, K = K,
      X = X_block, colnm = colnm_expansions,
      scales = expansion_scales,
      deriv_fxn = all_derivatives_fxn, ...
    )
    custom_bvec <- qp_bvec_fxn(
      qp_Amat = custom_Amat, N = N_obs,
      p = p_expansions, K = K,
      X = X_block, colnm = colnm_expansions,
      scales = expansion_scales,
      deriv_fxn = all_derivatives_fxn, ...
    )
    custom_meq <- qp_meq_fxn(
      qp_Amat = custom_Amat, N = N_obs,
      p = p_expansions, K = K,
      X = X_block, colnm = colnm_expansions,
      scales = expansion_scales,
      deriv_fxn = all_derivatives_fxn, ...
    )

    ## Validate dimensions
    P_coef <- p_expansions * (K + 1)
    if (!is.null(custom_Amat) && nrow(custom_Amat) != P_coef) {
      stop(
        '\n \t qp_Amat_fxn returned a matrix with ', nrow(custom_Amat),
        ' rows, but the coefficient vector has length p*(K+1) = ',
        P_coef, '. The constraint matrix must have nrow = p*(K+1).',
        '\n \t Your function receives named arguments: ',
        'p (expansions per partition), K (number of interior knots), ',
        'N (observations), X (block-diagonal design), ',
        'scales (expansion scales), colnm (column names), ',
        'deriv_fxn (derivative function).\n')
    }
    if (!is.null(custom_bvec) && length(custom_bvec) != ncol(custom_Amat)) {
      stop(
        '\n \t qp_bvec_fxn returned a vector of length ', length(custom_bvec),
        ' but qp_Amat has ', ncol(custom_Amat), ' columns (constraints).',
        ' These must match.\n')
    }

    qp_Amat_list[[length(qp_Amat_list) + 1]] <- custom_Amat
    qp_bvec_list[[length(qp_bvec_list) + 1]] <- custom_bvec
    qp_meq_list[[length(qp_meq_list) + 1]] <- custom_meq
    qp_type_list[[length(qp_type_list) + 1]] <- "custom"
  }

  ## Preserve any user-supplied solve.QP objects by prepending them to
  # the built-in constraint lists.  This keeps the leading qp_meq
  # equalities in the user-specified order.
  if (!is.null(qp_Amat) & !is.null(qp_bvec) & !is.null(qp_meq)) {
    qp_Amat_list <- c(list(qp_Amat), qp_Amat_list)
    qp_bvec_list <- c(list(qp_bvec), qp_bvec_list)
    qp_meq_list <- c(list(qp_meq), qp_meq_list)
    qp_type_list <- c(list("user"), qp_type_list)
  }

  ## Reduce partition-local inequality columns before stacking the final
  #  solve.QP matrix. Columns spanning multiple partitions, such as
  #  monotonicity constraints, are left as supplied.
  if (qr_pivot_inequality_constraints) {
    reduced_qp <- lapply(seq_along(qp_Amat_list), function(i) {
      if (qp_type_list[[i]] %in%
          c("qp_monotonic_increase", "qp_monotonic_decrease",
            "qp_range_lower", "qp_range_upper", "qp_range_both")) {
        return(list(
          Amat = qp_Amat_list[[i]],
          bvec = qp_bvec_list[[i]],
          meq = qp_meq_list[[i]]
        ))
      }
      .reduce_qp_constraint_block(
        qp_Amat = qp_Amat_list[[i]],
        qp_bvec = qp_bvec_list[[i]],
        qp_meq = qp_meq_list[[i]],
        p_expansions = p_expansions,
        K = K,
        parallel_qr_qp = parallel_qr_qp,
        cl = cl
      )
    })
    qp_Amat_list <- lapply(reduced_qp, `[[`, "Amat")
    qp_bvec_list <- lapply(reduced_qp, `[[`, "bvec")
    qp_meq_list <- lapply(reduced_qp, `[[`, "meq")
  }

  ## If QP flags were supplied but no conformable columns were generated,
  # return the non-QP state rather than carrying a NULL constraint matrix
  # through the solver stack.
  if (length(qp_Amat_list) == 0) {
    if (include_warnings) {
      warning(
        "\n\t Quadratic programming constraints were requested, but no ",
        "conformable constraint columns were generated. Proceeding ",
        "without QP refinement.\n"
      )
    }
    return(list(
      qp_Amat = NULL,
      qp_bvec = NULL,
      qp_meq = 0L,
      quadprog = FALSE
    ))
  }

  ## Stack all active constraint sources into the single solve.QP interface.
  qp_Amat_out <- do.call(cbind, qp_Amat_list)
  qp_bvec_out <- do.call(c, qp_bvec_list)
  qp_meq_out <- sum(unlist(qp_meq_list))

  list(
    qp_Amat = qp_Amat_out,
    qp_bvec = qp_bvec_out,
    qp_meq = qp_meq_out,
    quadprog = TRUE
  )
}
