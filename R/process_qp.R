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
                            og_cols = NULL) {

  N_sub <- nrow(X_block)

  ## Recover the p-column expansion matrix from the block-diagonal design.
  # Each row of X_block belongs to exactly one partition k, and its
  # p_expansions nonzero entries sit in columns [(k-1)*p + 1] : [k*p].
  C_qp <- matrix(0, N_sub, p_expansions)
  colnames(C_qp) <- colnm_expansions
  partition_assignment <- integer(N_sub)
  for (i in 1:N_sub) {
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
  for (v in seq_along(deriv_list)) {
    deriv_mat <- deriv_list[[v]]
    for (i in 1:nrow(deriv_mat)) {
      k <- partition_assignment[i]
      col_P <- rep(0, p_expansions * (K + 1))
      col_P[((k - 1) * p_expansions + 1):(k * p_expansions)] <-
        sign_mult * deriv_mat[i, ]
      constraint_cols[[length(constraint_cols) + 1]] <- col_P
    }
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
#' @param qp_observations Optional integer vector of observation indices.
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
#' @param all_derivatives_fxn Function to compute derivatives from expansion
#'   matrices (the \code{all_derivatives} closure from \code{lgspline.fit}).
#' @param og_cols Optional character vector of original predictor column names.
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
#' set.seed(42)
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
                       all_derivatives_fxn = NULL,
                       og_cols = NULL,
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

  ## Build the block-diagonal design matrix for QP observations.
  # This is N_sub x P where P = p_expansions * (K+1).
  X_block <- Reduce("rbind", lapply(1:(K + 1), function(k) {
    if (nrow(X[[k]]) == 0) return(X[[k]])
    Reduce("cbind", lapply(1:(K + 1), function(j) {
      if (j == k) X[[k]] else 0 * X[[k]]
    }))
  }))

  if (any(!is.null(qp_observations))) {
    qp_observations <- try(c(qp_observations), silent = TRUE)
    if (any(inherits(qp_observations, 'try-error'))) {
      stop('\n \t qp_observations must be coercible to a numeric vector \n')
    }
    X_block <- X_block[c(unlist(order_list)) %in% qp_observations, , drop = FALSE]
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

    if (!any(is.null(qp_range_upper)) & !any(is.null(qp_range_lower))) {
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

    } else if (!any(is.null(qp_range_upper))) {
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

    } else if (!any(is.null(qp_range_lower))) {
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
      target_vars = target_pos1
    ))
    deriv_qp <- do.call(.build_deriv_qp, args)
    if (ncol(deriv_qp$Amat) > 0) {
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- deriv_qp$Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- deriv_qp$bvec
      qp_meq_list[[length(qp_meq_list) + 1]] <- deriv_qp$meq
    }
  }

  ## Negative first derivative
  target_neg1 <- .resolve_target(qp_negative_derivative)
  if (!identical(target_neg1, NA)) {
    args <- c(deriv_shared, list(
      sign_mult = -1,
      just_first = TRUE,
      target_vars = target_neg1
    ))
    deriv_qp <- do.call(.build_deriv_qp, args)
    if (ncol(deriv_qp$Amat) > 0) {
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- deriv_qp$Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- deriv_qp$bvec
      qp_meq_list[[length(qp_meq_list) + 1]] <- deriv_qp$meq
    }
  }

  ## Positive second derivative
  target_pos2 <- .resolve_target(qp_positive_2ndderivative)
  if (!identical(target_pos2, NA)) {
    args <- c(deriv_shared, list(
      sign_mult = +1,
      just_first = FALSE,
      target_vars = target_pos2
    ))
    deriv_qp <- do.call(.build_deriv_qp, args)
    if (ncol(deriv_qp$Amat) > 0) {
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- deriv_qp$Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- deriv_qp$bvec
      qp_meq_list[[length(qp_meq_list) + 1]] <- deriv_qp$meq
    }
  }

  ## Negative second derivative
  target_neg2 <- .resolve_target(qp_negative_2ndderivative)
  if (!identical(target_neg2, NA)) {
    args <- c(deriv_shared, list(
      sign_mult = -1,
      just_first = FALSE,
      target_vars = target_neg2
    ))
    deriv_qp <- do.call(.build_deriv_qp, args)
    if (ncol(deriv_qp$Amat) > 0) {
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- deriv_qp$Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- deriv_qp$bvec
      qp_meq_list[[length(qp_meq_list) + 1]] <- deriv_qp$meq
    }
  }

  ## Monotonic increase (observation order, not per-variable)
  if (qp_monotonic_increase || qp_monotonic_decrease) {
    sign_mono <- if(qp_monotonic_increase) +1 else -1

    ## Recover original observation order from partition-stacked X_block.
    #  order_list[[k]] gives the original row indices for partition k.
    #  Stacking them in partition order gives the permutation applied to
    #  build X_block; inverting it restores original observation order.
    obs_order    <- unlist(order_list)
    inv_order    <- order(obs_order)
    X_block_orig <- X_block[inv_order, , drop = FALSE]

    ## If qp_observations was applied, X_block may be a subset.
    #  In that case inv_order may be longer than nrow(X_block); guard.
    if(nrow(X_block_orig) != nrow(X_block)){
      ## Fallback: use X_block as-is (already filtered by qp_observations)
      X_block_orig <- X_block
    }

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
  }

  ## Preserve any user-supplied solve.QP objects by prepending them to
  # the built-in constraint lists.  This keeps the leading qp_meq
  # equalities in the user-specified order.
  if (!is.null(qp_Amat) & !is.null(qp_bvec) & !is.null(qp_meq)) {
    qp_Amat_list <- c(list(qp_Amat), qp_Amat_list)
    qp_bvec_list <- c(list(qp_bvec), qp_bvec_list)
    qp_meq_list <- c(list(qp_meq), qp_meq_list)
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
