## tune_Lambda.R
#  Tuning for the smoothing and ridge penalties via exact leave-one-out
#  or modified generalized cross-validation.
#
## Main pieces:
#    tune_Lambda()
#    .compute_gcvu() and .compute_gcvu_gradient()
#    .damped_bfgs() and .tune_grid_search()
#    small helpers for predictions, residuals, and meta-penalties


## Small helpers for outer tuning parallelism

.tuning_worker_count <- function(cl) {
  if (is.null(cl) || !inherits(cl, "cluster")) {
    return(1L)
  }
  as.integer(length(cl))
}


.tuning_eval_env <- function(env) {
  env_eval <- env
  env_eval$parallel <- FALSE
  env_eval$parallel_eigen <- FALSE
  env_eval$parallel_trace <- FALSE
  env_eval$parallel_aga <- FALSE
  env_eval$parallel_matmult <- FALSE
  env_eval$parallel_qr <- FALSE
  env_eval$parallel_unconstrained <- FALSE
  env_eval$cl <- NULL
  env_eval$chunk_size <- 1L
  env_eval$num_chunks <- 0L
  env_eval$rem_chunks <- 0L
  env_eval
}


.preserve_rng_state <- function(expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  force(expr)
}


.tuning_seed <- function(...) {
  vals <- abs(c(...))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) {
    return(1L)
  }
  seed <- sum(vals * seq_along(vals) * 1e6)
  as.integer(seed %% (.Machine$integer.max - 1L)) + 1L
}


.expand_tuning_grid <- function(initial_grid,
                                log_initial_wiggle,
                                log_initial_flat,
                                n_workers,
                                parallel_grideval) {
  if (!parallel_grideval || n_workers <= 6L) {
    return(initial_grid)
  }

  n_extra <- min(n_workers - 6L, 8L)
  if (n_extra <= 0L) {
    return(initial_grid)
  }

  wiggle_bounds <- range(exp(log_initial_wiggle)) * c(1e-1, 10)
  flat_bounds <- range(exp(log_initial_flat)) * c(1e-1, 10)
  seed <- .tuning_seed(log_initial_wiggle, log_initial_flat, n_workers, n_extra)

  extra_grid <- .preserve_rng_state({
    set.seed(seed)
    data.frame(
      wiggle = log(runif(n_extra, wiggle_bounds[1], wiggle_bounds[2])),
      flat = log(runif(n_extra, flat_bounds[1], flat_bounds[2]))
    )
  })

  rbind(initial_grid, extra_grid)
}


.tuning_adaptive_criterion_tol <- function(tol, criterion_value) {
  crit_scale <- max(abs(c(criterion_value)), sqrt(.Machine$double.eps),
                    na.rm = TRUE)
  max(tol, min(1e-5, 1e-5 * crit_scale))
}


.safe_tuning_outlist <- function(par,
                                 criterion_fxn,
                                 log_penalty_vec,
                                 env,
                                 include_warnings = FALSE,
                                 ...) {
  tryCatch({
    out <- criterion_fxn(c(unlist(par), log_penalty_vec),
                         log_penalty_vec,
                         env,
                         ...)
    crit <- out$criterion_value
    if (length(crit) != 1L || is.na(crit) || is.nan(crit) || !is.finite(crit)) {
      out$criterion_value <- Inf
    }
    out
  }, error = function(e) {
    if (include_warnings) {
      warning(conditionMessage(e))
    }
    list(criterion_value = Inf)
  })
}


.safe_tuning_value <- function(par,
                               criterion_fxn,
                               log_penalty_vec,
                               env,
                               include_warnings = FALSE,
                               ...) {
  .safe_tuning_outlist(par,
                       criterion_fxn,
                       log_penalty_vec,
                       env,
                       include_warnings,
                       ...)$criterion_value
}


.bfgs_damp_set <- function(damp, n_workers) {
  n_damps <- max(1L, min(n_workers, 8L))
  damp * 2^(-(seq_len(n_damps) - 1L))
}


## Residuals used in tuning criteria

#' Compute Residuals During Penalty Tuning
#'
#' @description
#' Computes the residuals used by the tuning criteria. Handles identity link,
#' general GLM link functions with pseudocount delta, custom deviance residual
#' functions, and observation weights.
#'
#' @param y List; response vectors by partition.
#' @param preds List; prediction vectors by partition.
#' @param delta Numeric; pseudocount for link function stabilization.
#' @param family GLM family object.
#' @param observation_weights List; observation weights by partition.
#' @param K Integer; number of interior knots (partitions - 1).
#' @param order_list List; observation indices per partition.
#' @param ... Additional arguments passed to \code{family$custom_dev.resids}.
#'
#' @return List of residual vectors, one per partition.
#'
#' @details
#' Three computation paths:
#' \enumerate{
#'   \item \strong{Identity link or no custom deviance}: Standard
#'     \eqn{g(y) - \hat{\eta}} residuals on the link scale, optionally
#'     weighted by observation weights for non-Gaussian families.
#'     For Gaussian identity link with heterogeneous weights, the weights
#'     have already been absorbed into X and y prior to this call.
#'   \item \strong{Custom deviance residuals}: Delegates to
#'     \code{family$custom_dev.resids(y, mu, order_indices, family,
#'     observation_weights, ...)}.
#' }
#'
#' @keywords internal
.compute_tuning_residuals <- function(y,
                                      preds,
                                      delta,
                                      family,
                                      observation_weights,
                                      K,
                                      order_list,
                                      ...) {
  if (paste0(family)[2] == "identity" ||
      is.null(family$custom_dev.resids)) {

    ## If not canonical Gaussian and weights are present, use them.
    ## Recall: for Gaussian identity, X and y were pre-weighted upstream.
    if (any(!is.null(observation_weights[[1]])) &&
        (paste0(family)[2] != "identity" ||
         paste0(family)[1] != "gaussian")) {
      residuals <- lapply(1:(K + 1), function(k) {
        (family$linkfun((y[[k]] + delta) / (1 + 2 * delta)) -
           (preds[[k]] + delta) / (1 + 2 * delta)) *
          c(observation_weights[[k]])
      })
    } else {
      residuals <- lapply(1:(K + 1), function(k) {
        family$linkfun((y[[k]] + delta) / (1 + 2 * delta)) -
          (preds[[k]] + delta) / (1 + 2 * delta)
      })
    }

  } else {
    ## Otherwise use the family's custom deviance residuals.
    residuals <- lapply(1:(K + 1), function(k) {
      family$custom_dev.resids(y[[k]],
                               family$linkinv(c(preds[[k]])),
                               order_list[[k]],
                               family,
                               observation_weights[[k]],
                               ...)
    })
  }

  return(residuals)
}


## Predictions used in tuning criteria

#' Compute Predictions During Penalty Tuning
#'
#' @description
#' Wrapper around \code{matmult_block_diagonal} for computing partition-wise
#' predictions \eqn{\mathbf{X}_{k} \boldsymbol{\beta}_{k}} during penalty
#' tuning.
#'
#' @param X List; design matrices by partition.
#' @param B List; coefficient vectors by partition.
#' @param K Integer; number of interior knots.
#' @param parallel Logical; use parallel computation.
#' @param cl Parallel cluster object.
#' @param chunk_size,num_chunks,rem_chunks Integer; parallel chunking parameters.
#'
#' @return List of prediction vectors, one per partition.
#'
#' @keywords internal
.compute_tuning_predictions <- function(X, B, K,
                                        parallel, cl,
                                        chunk_size, num_chunks,
                                        rem_chunks) {
  matmult_block_diagonal(X, B, K,
                         parallel, cl,
                         chunk_size, num_chunks, rem_chunks)
}


## Meta-penalty on the tuning parameters

#' Compute Regularization (Meta) Penalty on Penalty Parameters
#'
#' @description
#' Computes the regularization term that pulls predictor- and partition-specific
#' penalty parameters toward 1 on the raw (positive) scale. This acts as a
#' "meta-penalty" on the penalty magnitudes themselves.
#'
#' @param wiggle_penalty Numeric; current wiggle penalty on raw scale.
#' @param penalty_vec Numeric vector; current predictor/partition penalties
#'   on raw scale. May be empty (\code{c()}).
#' @param meta_penalty_coef Numeric; coefficient for the meta-penalty.
#' @param unique_penalty_per_predictor Logical; whether predictor-specific
#'   penalties are active.
#' @param unique_penalty_per_partition Logical; whether partition-specific
#'   penalties are active.
#'
#' @return Numeric scalar; the regularization penalty value.
#'
#' @details
#' The penalty takes the form:
#' \deqn{0.5 \times c_{\mathrm{meta}} \times \sum_{j} (\lambda_{j} - 1)^{2}
#'       + 0.5 \times 10^{-32} \times (\lambda_{w} - 1)^{2}}
#' where \eqn{\lambda_{j}} are predictor/partition penalties and
#' \eqn{\lambda_{w}} is the wiggle penalty.
#'
#' @keywords internal
.compute_meta_penalty <- function(wiggle_penalty,
                                  penalty_vec,
                                  meta_penalty_coef,
                                  unique_penalty_per_predictor,
                                  unique_penalty_per_partition) {
  if (unique_penalty_per_partition || unique_penalty_per_predictor) {
    0.5 * meta_penalty_coef * sum((penalty_vec - 1)^2) +
      0.5 * 1e-32 * ((wiggle_penalty - 1))^2
  } else {
    0.5 * 1e-32 * ((wiggle_penalty - 1))^2
  }
}


## Gradient of the meta-penalty

#' Compute Gradient of Regularization (Meta) Penalty
#'
#' @description
#' Computes the gradient of the meta-penalty with respect to the log-scale
#' penalty parameters, incorporating the exp parameterization chain rule.
#'
#' @param wiggle_penalty Numeric; current wiggle penalty on raw scale.
#' @param penalty_vec Numeric vector; current predictor/partition penalties
#'   on raw scale. May be empty (\code{c()}).
#' @param meta_penalty_coef Numeric; coefficient for the meta-penalty.
#' @param unique_penalty_per_predictor Logical; whether predictor-specific
#'   penalties are active.
#' @param unique_penalty_per_partition Logical; whether partition-specific
#'   penalties are active.
#'
#' @return Numeric vector; gradient of the meta-penalty on the log scale.
#'   Length equals 2 + length(penalty_vec).
#'
#' @details
#' Under exp parameterization \eqn{\lambda = \exp(\theta)}:
#' \deqn{\frac{\partial}{\partial \theta}
#'   \left[ 0.5 c (\exp(\theta) - 1)^{2} \right]
#'   = c (\lambda - 1) \lambda}
#'
#' @keywords internal
.compute_meta_penalty_gradient <- function(wiggle_penalty,
                                           penalty_vec,
                                           meta_penalty_coef,
                                           unique_penalty_per_predictor,
                                           unique_penalty_per_partition) {
  if (unique_penalty_per_partition || unique_penalty_per_predictor) {
    c(1e-32 * (wiggle_penalty - 1) * wiggle_penalty,
      0,
      meta_penalty_coef * (penalty_vec - 1) * penalty_vec)
  } else {
    c(1e-32 * (wiggle_penalty - 1) * wiggle_penalty,
      0)
  }
}


## Coefficient fit inside tuning

#' Fit Coefficients During Penalty Tuning: blockfit_solve or get_B
#'
#' @description
#' Dispatches to \code{blockfit_solve} when the blockfit conditions are met
#' (i.e. \code{env$use_blockfit} is TRUE), otherwise calls \code{get_B}.
#' On \code{blockfit_solve} failure, falls back to \code{get_B} automatically.
#'
#' @details
#' The blockfit condition mirrors \code{lgspline.fit}:
#' \code{blockfit && length(flat_cols) > 0 && K > 0}, pre-computed in
#' \code{tune_Lambda} and stored in \code{env$use_blockfit}.
#'
#' \code{return_G_getB} is set to TRUE by the callers so that
#' \code{B_list$G_list} contains the updated G matrices (after any GLM
#' weight iteration inside \code{get_B} or \code{blockfit_solve}).
#' These are needed immediately after this call for \code{AGAmult_wrapper}
#' and the trace computation.
#'
#' When correlation structure inputs are supplied, this wrapper does not
#' introduce a separate tuning-specific notation or solver path. Instead, the
#' same correlated coefficient estimator used in the final model fit is called
#' here inside each tuning-objective evaluation. In particular, any
#' structured-correlation Woodbury factorization
#' \eqn{\mathbf{G}^{-1} =
#' \mathbf{G}_{\mathrm{on}}^{-1} + \mathbf{G}_{\mathrm{off}}^{-1}}
#' is handled inside \code{get_B} and documented in \code{lgspline-details}.
#'
#' @param G_list List; eigendecomposition results from \code{compute_G_eigen},
#'   containing \code{Ghalf} and \code{GhalfInv}.
#' @param Lambda Matrix; current combined penalty matrix.
#' @param L_partition_list List; partition-specific penalty matrices.
#' @param env List; tuning environment from \code{.build_tuning_env}.
#' @param return_G_getB Logical; whether to return G inside the fit. Set to
#'   TRUE within \code{.compute_gcvu} and \code{.compute_gcvu_gradient} so
#'   that \code{B_list$G_list} carries the (possibly GLM-iterated) G matrices
#'   required for the subsequent \code{AGAmult_wrapper} and trace computations.
#' @param ... Additional arguments forwarded to the fitting routine.
#'
#' @return List; output of \code{blockfit_solve} or \code{get_B}, containing
#'   at minimum \code{$B} (coefficient list) and \code{$G_list}.
#'
#' @keywords internal
.fit_coefficients <- function(G_list,
                              Lambda,
                              L_partition_list,
                              env,
                              return_G_getB,
                              ...) {
  initial_active_ineq <- .tuning_seed_active_set(env)

  ## Attempt blockfit_solve when conditions are met
  if (env$use_blockfit) {
    B_list <- try({
      blockfit_solve(
        X                              = env$X,
        y                              = env$y,
        flat_cols                      = env$flat_cols,
        K                              = env$K,
        p_expansions                   = env$p_expansions,
        Lambda                         = Lambda,
        L_partition_list               = L_partition_list,
        unique_penalty_per_partition   = env$unique_penalty_per_partition,
        A                              = env$A,
        R_constraints                  = env$R_constraints,
        constraint_values              = env$constraint_value_vectors,
        X_gram                         = env$X_gram,
        Ghalf_full                     = G_list$Ghalf,
        GhalfInv_full                  = G_list$GhalfInv,
        family                         = env$family,
        order_list                     = env$order_list,
        glm_weight_function            = env$glm_weight_function,
        schur_correction_function      = env$schur_correction_function,
        need_dispersion_for_estimation = env$need_dispersion_for_estimation,
        dispersion_function            = env$dispersion_function,
        observation_weights            = env$observation_weights,
        homogenous_weights             = env$homogenous_weights,
        iterate                        = env$iterate,
        tol                            = env$tol,
        parallel_eigen                 = env$parallel & env$parallel_eigen,
        parallel_qr                    = env$parallel & env$parallel_qr,
        qr_pivot_smoothing_constraints =
          env$qr_pivot_smoothing_constraints,
        cl                             = env$cl,
        chunk_size                     = env$chunk_size,
        num_chunks                     = env$num_chunks,
        rem_chunks                     = env$rem_chunks,
        return_G_getB                  = return_G_getB,
        quadprog                       = env$quadprog,
        qp_Amat                        = env$qp_Amat,
        qp_bvec                        = env$qp_bvec,
        qp_meq                         = env$qp_meq,
        qp_score_function              = env$qp_score_function,
        keep_weighted_Lambda           = env$keep_weighted_Lambda,
        max_backfit_iter               = 100,
        Vhalf                          = env$Vhalf,
        VhalfInv                       = env$VhalfInv,
        initial_active_ineq            = initial_active_ineq,
        include_warnings               = env$include_warnings,
        verbose                        = env$verbose,
        ...
      )
    }, silent = TRUE)

    ## Fall back to get_B on failure
    if (inherits(B_list, "try-error")) {
      if (env$include_warnings) {
        warning("\n \t blockfit_solve failed during GCV tuning, ",
                "falling back to get_B \n")
      }
      B_list <- .fit_get_B(G_list, Lambda, L_partition_list,
                           env, return_G_getB, ...)
    }

    return(B_list)
  }

  ## Otherwise use the usual get_B path.
  .fit_get_B(G_list, Lambda, L_partition_list, env, return_G_getB, ...)
}


## Small get_B wrapper for tuning

#' Call get_B During Penalty Tuning
#'
#' @description
#' Internal wrapper that calls \code{get_B} with all arguments drawn from
#' the tuning environment \code{env}. Separated from \code{.fit_coefficients}
#' so the fallback path in \code{.fit_coefficients} is clean and does not
#' repeat the full argument list.
#'
#' This helper is the point at which penalty tuning inherits the full
#' constrained solver described in \code{lgspline-details}, including GLM
#' reweighting, dense correlation whitening, and the structured-correlation
#' Woodbury acceleration when available. No additional tuning-specific
#' approximation is introduced here beyond whatever \code{get_B} itself uses.
#'
#' @inheritParams .fit_coefficients
#' @keywords internal
.fit_get_B <- function(G_list,
                       Lambda,
                       L_partition_list,
                       env,
                       return_G_getB,
                       ...) {
  get_B(
    env$X,
    env$X_gram,
    Lambda,
    env$keep_weighted_Lambda,
    env$unique_penalty_per_partition,
    L_partition_list,
    env$A,
    env$Xy,
    env$y,
    env$K,
    env$p_expansions,
    env$R_constraints,
    G_list$Ghalf,
    G_list$GhalfInv,
    env$parallel & env$parallel_eigen,
    env$parallel & env$parallel_aga,
    env$parallel & env$parallel_matmult,
    env$parallel & env$parallel_qr,
    env$parallel & env$parallel_unconstrained,
    env$cl,
    env$chunk_size,
    env$num_chunks,
    env$rem_chunks,
    env$family,
    env$unconstrained_fit_fxn,
    env$iterate,
    env$qp_score_function,
    env$quadprog,
    env$qp_Amat,
    env$qp_bvec,
    env$qp_meq,
    prevB                          = NULL,
    prevUnconB                     = NULL,
    iter_count                     = 0,
    prev_diff                      = Inf,
    env$tol,
    env$constraint_value_vectors,
    env$order_list,
    env$glm_weight_function,
    env$schur_correction_function,
    env$need_dispersion_for_estimation,
    env$dispersion_function,
    env$observation_weights,
    env$homogenous_weights,
    return_G_getB,
    env$blockfit,
    env$just_linear_without_interactions,
    env$Vhalf,
    env$VhalfInv,
    gee_precomp = env$gee_precomp,
    initial_active_ineq = .tuning_seed_active_set(env),
    ...
  )
}


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# Subfunction: .compute_gcvu
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

#' Evaluate Modified GCV_u Criterion at a Given Penalty Configuration
#'
#' @description
#' Computes the modified GCV_u criterion for a given set of penalty
#' parameters. This is the objective function minimized during penalty tuning;
#' \code{gcv_gamma = 1} recovers ordinary GCV.
#'
#' @param par Numeric vector; log-scale penalty parameters. First two elements
#'   are log(wiggle_penalty) and log(flat_ridge_penalty). Remaining elements
#'   (if any) are log-scale predictor/partition penalties.
#' @param log_penalty_vec Numeric vector; log-scale predictor/partition
#'   penalties (passed separately for compatibility with the grid search).
#' @param env List; pre-computed objects and tuning configuration. Contains:
#'   \describe{
#'     \item{X, y, X_gram, Xy}{Design matrices, response, Gram matrices,
#'       cross-products.}
#'     \item{A, R_constraints, K, p_expansions, N_obs}{Constraint matrix and
#'       dimensions.}
#'     \item{smoothing_spline_penalty}{Integrated second-derivative penalty
#'       matrix.}
#'     \item{custom_penalty_mat}{Optional custom penalty matrix.}
#'     \item{colnm_expansions}{Character vector of expansion column names.}
#'     \item{unique_penalty_per_predictor, unique_penalty_per_partition}{Logicals.}
#'     \item{family}{GLM family object.}
#'     \item{delta}{Pseudocount for link function stabilization.}
#'     \item{gcv_gamma}{Modified-GCV inflation factor applied to the
#'       effective-degrees-of-freedom term in the denominator.}
#'     \item{meta_penalty}{Meta-penalty coefficient.}
#'     \item{order_list}{List of observation indices per partition.}
#'     \item{observation_weights, homogenous_weights}{Observation weighting.}
#'     \item{parallel flags, cl, chunk_size, num_chunks, rem_chunks}{Parallel
#'       config.}
#'     \item{unconstrained_fit_fxn, keep_weighted_Lambda, iterate}{Fitting
#'       config.}
#'     \item{qp_score_function, quadprog, qp_Amat, qp_bvec, qp_meq}{QP config.}
#'     \item{tol, sd_y}{Convergence tolerance and response scale.}
#'     \item{constraint_value_vectors}{Constraint values.}
#'     \item{glm_weight_function, schur_correction_function}{GLM functions.}
#'     \item{need_dispersion_for_estimation, dispersion_function}{Dispersion.}
#'     \item{blockfit, just_linear_without_interactions}{Blockfit config.}
#'     \item{use_blockfit}{Logical; pre-computed dispatch flag for blockfit.}
#'     \item{flat_cols}{Integer vector; pre-computed flat column indices.}
#'     \item{Vhalf, VhalfInv}{Correlation structure matrices.}
#'     \item{verbose, include_warnings}{Output control.}
#'   }
#' @param ... Additional arguments passed to fitting functions.
#'
#' @return List containing:
#' \describe{
#'   \item{GCV_u}{Numeric; GCV_u criterion value including meta-penalty.}
#'   \item{B}{List; fitted coefficient vectors by partition.}
#'   \item{GXX}{List; \eqn{\mathbf{G}_{k} \mathbf{X}_{k}^{\top}\mathbf{X}_{k}}
#'     matrices.}
#'   \item{G_list}{List; eigendecomposition results from
#'     \code{compute_G_eigen}.}
#'   \item{mean_W}{Numeric; \eqn{\mathrm{tr}(\mathbf{H})/N}, the average
#'     leverage entering the denominator of \eqn{\mathrm{GCV}_u}.}
#'   \item{sum_W}{Numeric; \eqn{\mathrm{tr}(\mathbf{H})}, the effective
#'     degrees of freedom.}
#'   \item{Lambda}{Matrix; combined penalty matrix
#'     \eqn{\boldsymbol{\Lambda}}.}
#'   \item{L1}{Matrix; baseline smoothness penalty component.}
#'   \item{L2}{Matrix; baseline ridge penalty component.}
#'   \item{L_predictor_list}{List; predictor-specific penalty matrices.}
#'   \item{L_partition_list}{List; partition-specific penalty matrices.}
#'   \item{numerator}{Numeric; sum of squared residuals.}
#'   \item{denominator}{Numeric; modified-GCV denominator
#'     \eqn{N(1 - \gamma\bar{W})^{2}}.}
#'   \item{residuals}{List; residual vectors by partition.}
#'   \item{denom_sq}{Numeric; squared denominator.}
#'   \item{AGAInv}{Matrix; \eqn{(\mathbf{A}^{\top}\mathbf{G}\mathbf{A})^{-1}}.}
#' }
#'
#' @keywords internal
.compute_gcvu <- function(par,
                          log_penalty_vec,
                          env,
                          ...) {
  verbose <- env$verbose
  gamma <- env$gcv_gamma

  if (verbose) cat("        gcvu_fxn start\n")

  ## Unpack penalties from log scale
  wiggle_penalty     <- exp(par[1])
  flat_ridge_penalty <- exp(par[2])
  if (env$unique_penalty_per_predictor || env$unique_penalty_per_partition) {
    penalty_vec <- exp(c(par[-c(1:2)]))
  } else {
    penalty_vec <- c()
  }

  ## Compute penalty matrix Lambda and components
  if (verbose) cat("        compute_Lambda\n")
  Lambda_list <- compute_Lambda(env$custom_penalty_mat,
                                env$smoothing_spline_penalty,
                                wiggle_penalty,
                                flat_ridge_penalty,
                                env$K,
                                env$p_expansions,
                                env$unique_penalty_per_predictor,
                                env$unique_penalty_per_partition,
                                penalty_vec,
                                env$colnm_expansions,
                                just_Lambda = FALSE)
  Lambda <- Lambda_list[[1]]
  L1     <- Lambda_list[[2]]
  L2     <- Lambda_list[[3]]

  ## Compute G matrices via eigendecomposition.
  #  keep_G = TRUE always: G is needed for the trace computation below
  #  regardless of whether get_B or blockfit_solve is used for fitting.
  if (verbose) cat("        compute_G_eigen\n")
  schur_corrections <- lapply(1:(env$K + 1), function(k) 0)
  G_list <- compute_G_eigen(env$X_gram,
                            Lambda,
                            env$K,
                            env$parallel & env$parallel_eigen,
                            env$cl,
                            env$chunk_size,
                            env$num_chunks,
                            env$rem_chunks,
                            env$family,
                            env$unique_penalty_per_partition,
                            Lambda_list$L_partition_list,
                            keep_G = TRUE,
                            schur_corrections)

  ## Fit coefficients: dispatches to blockfit_solve or get_B via
  #  .fit_coefficients, with automatic fallback on failure.
  #  return_G_getB = TRUE: B_list$G_list carries the (possibly GLM-iterated)
  #  G matrices needed for AGAmult_wrapper and trace computations below.
  if (verbose) cat("        gcvu_fxn fit coefficients\n")
  return_G_getB <- TRUE
  B_list <- .fit_coefficients(G_list, Lambda,
                              Lambda_list$L_partition_list,
                              env, return_G_getB, ...)
  .tuning_store_active_set(env, B_list$qp_info)
  G_list <- B_list$G_list
  B      <- B_list$B

  ## Compute (A^T G A)^{-1}
  if (verbose) cat("        gcvu_fxn AGAmult_wrapper\n")
  AGAInv <- invert(AGAmult_wrapper(G_list$G,
                                   env$A,
                                   env$K,
                                   env$p_expansions,
                                   env$R_constraints,
                                   env$parallel & env$parallel_aga,
                                   env$cl,
                                   env$chunk_size,
                                   env$num_chunks,
                                   env$rem_chunks) +
                     1e-16 * diag(ncol(env$A)))

  ## Compute G * X^T X
  if (verbose) cat("        gcvu_fxn matmult_block_diagonal for GXX\n")
  GXX <- matmult_block_diagonal(G_list$G,
                                env$X_gram,
                                env$K,
                                env$parallel & env$parallel_matmult,
                                env$cl,
                                env$chunk_size,
                                env$num_chunks,
                                env$rem_chunks)

  ## Compute trace of hat matrix
  if (verbose) cat("        gcvu_fxn compute_trace_UGXX_wrapper\n")
  sum_W <- compute_trace_H(G_list$G,
                           Lambda,
                           env$A,
                           AGAInv,
                           env$p_expansions,
                           env$R_constraints,
                           env$K,
                           env$parallel & env$parallel_trace,
                           env$cl,
                           env$chunk_size,
                           env$num_chunks,
                           env$rem_chunks,
                           env$unique_penalty_per_partition,
                           Lambda_list$L_partition_list)

  ## Predictions
  if (verbose) cat("        gcvu_fxn get predictions\n")
  preds <- .compute_tuning_predictions(env$X, B, env$K,
                                       FALSE,
                                       env$cl,
                                       env$chunk_size,
                                       env$num_chunks,
                                       env$rem_chunks)

  ## Residuals
  if (verbose) cat("        gcvu_fxn custom or default residuals\n")
  residuals <- .compute_tuning_residuals(env$y, preds, env$delta,
                                         env$family,
                                         env$observation_weights,
                                         env$K, env$order_list,
                                         ...)

  ## Compute GCV_u components
  if (verbose) cat("        gcvu_fxn GCVu operations\n")
  numerator   <- sum(unlist(residuals)^2)
  mean_W      <- sum_W / env$N_obs
  denominator <- env$N_obs * (1 - gamma * mean_W)^2
  denom_sq    <- denominator^2
  GCV_u       <- numerator / denominator

  ## Meta-penalty regularization
  if (verbose) cat("        gcvu_fxn penalization operations\n")
  mp <- .compute_meta_penalty(wiggle_penalty, penalty_vec,
                              env$meta_penalty,
                              env$unique_penalty_per_predictor,
                              env$unique_penalty_per_partition)

  if (verbose) cat("        done GCVu,", GCV_u, "\n")

  return(list(criterion_value = GCV_u + mp,
              GCV_u       = GCV_u + mp,
              B           = B,
              GXX         = GXX,
              G_list      = G_list,
              mean_W      = mean_W,
              sum_W       = sum_W,
              Lambda      = Lambda,
              L1          = L1,
              L2          = L2,
              L_predictor_list = Lambda_list$L_predictor_list,
              L_partition_list = Lambda_list$L_partition_list,
              numerator   = numerator,
              denominator = denominator,
              residuals   = residuals,
              denom_sq    = denom_sq,
              AGAInv      = AGAInv,
              qp_info     = B_list$qp_info))
}


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# Subfunction: .compute_gcvu_gradient
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

#' Compute Closed-Form Gradient of GCV_u Criterion
#'
#' @description
#' Computes the gradient of the modified GCV_u criterion with respect to the
#' log-scale penalty parameters using analytical derivatives of the hat matrix
#' trace and residual sum of squares.
#'
#' @param par Numeric vector; log-scale penalty parameters.
#' @param log_penalty_vec Numeric vector; log-scale predictor/partition
#'   penalties.
#' @param outlist List or NULL; pre-computed GCV_u components from
#'   \code{.compute_gcvu}. If NULL, they are computed internally.
#' @param env List; pre-computed objects and tuning configuration (same
#'   structure as in \code{.compute_gcvu}).
#' @param ... Additional arguments passed to fitting functions.
#'
#' @return List containing:
#' \describe{
#'   \item{GCV_u}{Numeric; GCV_u criterion value including meta-penalty.}
#'   \item{gradient}{Numeric vector; gradient on the log penalty scale.}
#'   \item{outlist}{List; GCV_u components (for reuse to avoid recomputation).}
#' }
#'
#' @details
#' The gradient is computed via:
#' \deqn{\frac{\partial \mathrm{GCV}_u}{\partial \theta}
#'   = \frac{1}{D^{2}} \left(
#'     \frac{\partial N}{\partial \theta} D
#'     - N \frac{\partial D}{\partial \theta}
#'   \right)}
#' where \eqn{N = \sum r_{i}^{2}} (numerator),
#' \eqn{D = n(1 - \gamma\bar{W})^{2}} (denominator), \eqn{\theta} is the
#' log-scale penalty parameter, and the chain rule
#' \eqn{d\lambda / d\theta = \lambda} (exp parameterization) is applied.
#'
#' After differentiating through \eqn{\mathbf{G}}, \eqn{\mathbf{U}}, and the
#' trace term once, the gradient is stored as
#' \eqn{\mathbf{M}_k = \partial \mathrm{GCV}_u /
#' \partial \boldsymbol{\Lambda}_k}. Each log-penalty derivative is then
#' \eqn{\sum_k \mathrm{tr}(\mathbf{M}_k\mathbf{L}_{j,k})}.
#'
#' @keywords internal
.compute_gcvu_gradient <- function(par,
                                   log_penalty_vec,
                                   outlist = NULL,
                                   env,
                                   ...) {
  verbose <- env$verbose
  gamma <- env$gcv_gamma

  if (verbose) cat("        gr_fxn start\n")

  ## Unpack penalties from log scale
  wiggle_penalty     <- exp(par[1])
  flat_ridge_penalty <- exp(par[2])
  if (env$unique_penalty_per_predictor || env$unique_penalty_per_partition) {
    penalty_vec <- exp(c(par[-c(1:2)]))
  } else {
    penalty_vec <- c()
  }

  ## Reparameterize
  lambda_1 <- wiggle_penalty
  lambda_2 <- flat_ridge_penalty

  if (verbose) {
    cat("        lambda_1, lambda_2: ", lambda_1, ", ", lambda_2, "\n")
  }

  ## Recompute the fitted GCV_u pieces only when they were not already
  #  passed in from the same parameter value.
  if (any(is.null(outlist))) {

    if (verbose) cat("        Lambda list\n")
    Lambda_list <- compute_Lambda(env$custom_penalty_mat,
                                  env$smoothing_spline_penalty,
                                  wiggle_penalty,
                                  flat_ridge_penalty,
                                  env$K,
                                  env$p_expansions,
                                  env$unique_penalty_per_predictor,
                                  env$unique_penalty_per_partition,
                                  penalty_vec,
                                  env$colnm_expansions,
                                  just_Lambda = FALSE)
    Lambda <- Lambda_list[[1]]
    L1     <- Lambda_list[[2]]
    L2     <- Lambda_list[[3]]

    if (verbose) cat("        G list\n")
    schur_corrections <- lapply(1:(env$K + 1), function(k) 0)
    G_list <- compute_G_eigen(env$X_gram,
                              Lambda,
                              env$K,
                              env$parallel & env$parallel_eigen,
                              env$cl,
                              env$chunk_size,
                              env$num_chunks,
                              env$rem_chunks,
                              env$family,
                              env$unique_penalty_per_partition,
                              Lambda_list$L_partition_list,
                              keep_G = TRUE,
                              schur_corrections)

    ## Fit coefficients: dispatches to blockfit_solve or get_B
    #  return_G_getB = TRUE: B_list$G_list carries the (possibly GLM-iterated)
    #  G matrices needed for AGAmult_wrapper and trace computations below.
    if (verbose) cat("        gr fxn fit coefficients\n")
    return_G_getB <- TRUE
    B_list <- .fit_coefficients(G_list, Lambda,
                                Lambda_list$L_partition_list,
                                env, return_G_getB, ...)
    .tuning_store_active_set(env, B_list$qp_info)
    G_list <- B_list$G_list
    B      <- B_list$B

    if (verbose) cat("        AGAmult_wrapper\n")
    AGAInv <- invert(AGAmult_wrapper(G_list$G,
                                     env$A,
                                     env$K,
                                     env$p_expansions,
                                     env$R_constraints,
                                     env$parallel & env$parallel_aga,
                                     env$cl,
                                     env$chunk_size,
                                     env$num_chunks,
                                     env$rem_chunks) +
                       1e-16 * diag(ncol(env$A)))

    if (verbose) cat("        GXX matmult_block_diagonal\n")
    GXX <- matmult_block_diagonal(G_list$G,
                                  env$X_gram,
                                  env$K,
                                  env$parallel & env$parallel_matmult,
                                  env$cl,
                                  env$chunk_size,
                                  env$num_chunks,
                                  env$rem_chunks)

    if (verbose) cat("        sum_W compute_trace_UGXX_wrapper\n")
    sum_W <- compute_trace_H(G_list$G,
                             Lambda,
                             env$A,
                             AGAInv,
                             env$p_expansions,
                             env$R_constraints,
                             env$K,
                             env$parallel & env$parallel_trace,
                             env$cl,
                             env$chunk_size,
                             env$num_chunks,
                             env$rem_chunks,
                             env$unique_penalty_per_partition,
                             Lambda_list$L_partition_list)

    if (verbose) cat("        gr fxn preds\n")
    preds <- .compute_tuning_predictions(env$X, B, env$K,
                                         FALSE,
                                         env$cl,
                                         env$chunk_size,
                                         env$num_chunks,
                                         env$rem_chunks)

    if (verbose) cat("        gr fxn residuals\n")
    residuals <- .compute_tuning_residuals(env$y, preds, env$delta,
                                           env$family,
                                           env$observation_weights,
                                           env$K, env$order_list,
                                           ...)

    if (verbose) cat("        gr fxn compute GCV_u\n")
    numerator   <- sum(unlist(residuals)^2)
    mean_W      <- sum_W / env$N_obs
    denominator <- env$N_obs * (1 - gamma * mean_W)^2
    denom_sq    <- denominator^2
    GCV_u       <- numerator / denominator

    mp <- .compute_meta_penalty(wiggle_penalty, penalty_vec,
                                env$meta_penalty,
                                env$unique_penalty_per_predictor,
                                env$unique_penalty_per_partition)

    if (verbose) cat("        gr fxn outlist\n")
    outlist <- list(criterion_value = GCV_u + mp,
                    GCV_u       = GCV_u + mp,
                    B           = B,
                    GXX         = GXX,
                    G_list      = G_list,
                    mean_W      = mean_W,
                    sum_W       = sum_W,
                    Lambda      = Lambda,
                    L1          = L1,
                    L2          = L2,
                    L_predictor_list = Lambda_list$L_predictor_list,
                    L_partition_list = Lambda_list$L_partition_list,
                    numerator   = numerator,
                    denominator = denominator,
                    residuals   = residuals,
                    denom_sq    = denom_sq,
                    AGAInv      = AGAInv,
                    qp_info     = B_list$qp_info)
  }

  ## Build the residual side of the criterion derivative.

  ## -2 * residuals^T * X
  if (verbose) cat("        neg2tresidX\n")
  neg2tresidX <- Reduce("cbind",
                        matmult_block_diagonal(
                          lapply(outlist$residuals,
                                 function(r) -2 * t(r)),
                          env$X,
                          env$K,
                          FALSE,
                          env$cl,
                          env$chunk_size,
                          env$num_chunks,
                          env$rem_chunks))

  ## Reuse the wiggle derivative for all penalty directions.
  #  L_*_list entries are current log-scale contributions to Lambda.
  cN <- 1 / outlist$denominator
  cW <- outlist$numerator * 2 * gamma * (1 - gamma * outlist$mean_W) /
    outlist$denom_sq
  xy_list <- env$Xy
  neg2_list <- .split_penalty_vector(c(neg2tresidX),
                                     env$K,
                                     env$p_expansions)
  block_F <- lapply(env$X_gram, function(XtX_k) cW * XtX_k)
  low_rank_F <- list(
    list(coef = 0.5 * cN, u = xy_list,   v = neg2_list),
    list(coef = 0.5 * cN, u = neg2_list, v = xy_list)
  )
  penalty_adjoint <- .compute_penalty_adjoint(
    outlist$G_list$G,
    env$A,
    outlist$AGAInv,
    env$K,
    env$p_expansions,
    block_F,
    low_rank_F,
    outlist$Lambda,
    outlist$L_partition_list,
    env$unique_penalty_per_partition
  )
  log_wiggle_gradient <- penalty_adjoint$wiggle_log_gradient
  log_flat_gradient <- sum(penalty_adjoint$M_sum * (lambda_1 * outlist$L2))

  if (verbose) cat("        Gradient start \n")
  gradient <- cbind(c(log_wiggle_gradient,
                      log_flat_gradient))

  if (env$unique_penalty_per_predictor) {
    predictor_penalty_gradient <- .compute_log_penalty_gradient(
      penalty_adjoint$M,
      outlist$L_predictor_list,
      M_sum = penalty_adjoint$M_sum
    )
    gradient <- cbind(c(c(gradient), predictor_penalty_gradient))
  }

  ## Partition-specific penalty gradients
  if (env$unique_penalty_per_partition) {
    partition_penalties <- penalty_vec[grep("partition", names(penalty_vec))]
    partition_penalty_gradient <- .compute_log_penalty_gradient(
      penalty_adjoint$M,
      outlist$L_partition_list,
      seq_along(partition_penalties)
    )
    gradient <- cbind(c(gradient, partition_penalty_gradient))
  }

  ## Add meta-penalty gradient
  regularizer <- .compute_meta_penalty_gradient(
    wiggle_penalty, penalty_vec, env$meta_penalty,
    env$unique_penalty_per_predictor,
    env$unique_penalty_per_partition)

  mp <- .compute_meta_penalty(wiggle_penalty, penalty_vec,
                              env$meta_penalty,
                              env$unique_penalty_per_predictor,
                              env$unique_penalty_per_partition)

  if (verbose) cat("        Gradient end \n")

  return(list(criterion_value = outlist$criterion_value,
              GCV_u    = outlist$GCV_u,
              gradient = gradient + regularizer,
              outlist  = outlist))
}


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# Subfunction: .damped_bfgs
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

#' Damped BFGS Optimizer for Penalty Tuning
#'
#' @description
#' Custom implementation of damped BFGS quasi-Newton optimization for
#' minimizing the selected tuning criterion. Uses step-size damping with
#' backtracking and Sherman-Morrison-Woodbury inverse Hessian updates.
#'
#' @param par Numeric vector; initial log-scale penalty parameters (first two
#'   elements are log(wiggle) and log(flat_ridge)).
#' @param log_penalty_vec Numeric vector; log-scale predictor/partition
#'   penalties appended to the optimization vector.
#' @param criterion_fxn Function; tuning-objective evaluation function with signature
#'   \code{function(par, log_penalty_vec, env, ...)}.
#' @param gr_fxn Function; gradient function with signature
#'   \code{function(par, log_penalty_vec, outlist, env, ...)}.
#' @param env List; tuning environment (passed through to
#'   \code{criterion_fxn} and \code{gr_fxn}).
#' @param tol Numeric; convergence tolerance for parameter change and the
#'   strict absolute criterion-change check. A small scale-aware plateau check
#'   is also used after the first few iterations.
#' @param max_iter Integer; maximum number of BFGS iterations (default 100).
#' @param ... Additional arguments passed to fitting functions.
#'
#' @return List containing the best parameter vector found, the corresponding
#'   criterion value, and the number of iterations performed.
#'
#' @details
#' The optimizer uses the following strategy:
#' \enumerate{
#'   \item Iterations 1-2: steepest descent with damping.
#'   \item Iteration 3+: BFGS quasi-Newton with inverse Hessian approximation
#'     updated via the standard secant condition. Falls back to identity
#'     matrix when the update is numerically unstable.
#'   \item Step acceptance: Armijo-like criterion (accept if
#'     \eqn{\mathrm{GCV}_{u}^{(\mathrm{new})} \leq
#'     \mathrm{GCV}_{u}^{(\mathrm{old})}}).
#'   \item Backtracking: damping factor halved on rejection; terminates
#'     when damp < \eqn{2^{-10}} (early iterations) or \eqn{2^{-12}}
#'     (later iterations).
#'   \item Convergence: after iteration 9, stop when either the strict
#'     \code{tol} checks are met or the accepted criterion changes remain
#'     below a small criterion-scale tolerance for several iterations.
#' }
#'
#' @keywords internal
.damped_bfgs <- function(par,
                         log_penalty_vec,
                         criterion_fxn,
                         gr_fxn,
                         env,
                         tol,
                         parallel_bfgs = TRUE,
                         max_iter = 100,
                         ...) {

  ## Initialize optimization parameters and storage
  lambda       <- cbind(par)
  old_lambda   <- lambda
  new_lambda   <- lambda
  n_params     <- length(lambda)
  Id           <- diag(n_params)
  initial_damp <- 0.5
  damp         <- initial_damp
  prev_gradient <- lambda * 0
  gradient      <- prev_gradient
  best_gradient <- prev_gradient
  old_gradient  <- prev_gradient
  outlist       <- NULL
  prev_outlist  <- NULL
  criterion_value <- Inf
  ridge         <- NULL
  dont_skip_gr  <- TRUE
  rho           <- NULL
  Inv           <- diag(n_params)
  best_lambda   <- lambda
  best_criterion_value <- criterion_value
  prev_lambda   <- old_lambda
  restart       <- TRUE
  n_workers     <- .tuning_worker_count(env$cl)
  criterion_env <- if (parallel_bfgs && !is.null(env$cl) && n_workers > 1L) {
    .tuning_eval_env(env)
  } else {
    env
  }
  small_change_count <- 0L
  small_change_patience <- 3L

  ## Evaluate and retain the grid-search start point before taking any
  # optimizer step. Otherwise a bad first gradient step can displace the
  # best-known solution even when the initialization was already better.
  outlist <- criterion_fxn(c(lambda[1], lambda[2], log_penalty_vec),
                           log_penalty_vec, criterion_env, ...)
  .tuning_store_active_set(env, outlist$qp_info)
  criterion_value <- outlist$criterion_value
  best_lambda <- lambda
  best_criterion_value <- criterion_value
  prev_outlist <- outlist

  ## Main optimization loop
  for (iter in 1:max_iter) {

    ## Compute gradient if needed
    if (dont_skip_gr) {
      result   <- gr_fxn(c(lambda[1], lambda[2], log_penalty_vec),
                         log_penalty_vec, outlist, criterion_env, ...)
      gradient <- cbind(result$gradient)
      outlist  <- result$outlist
      .tuning_store_active_set(env, outlist$qp_info)
    }

    ## First two iterations: steepest descent
    if (iter <= 2) {
      step_direction <- gradient

    } else {
      ## Add small ridge for stability if needed
      if (any(is.null(ridge))) ridge <- 1e-8 * diag(length(lambda))

      ## Update BFGS approximation
      if (dont_skip_gr) {

        ## Initial or restart BFGS approximation
        if (iter == 3 || restart) {
          diff_grad <- gradient - prev_gradient
          diff_lam  <- cbind(lambda - prev_lambda)
          denom     <- as.numeric(t(diff_grad) %**% diff_lam)
          Inv <- Inv +
            (denom + as.numeric(t(diff_grad) %**% Inv %**% diff_grad)) *
            (diff_lam %**% t(diff_lam)) / (denom^2) -
            (Inv %**% cbind(diff_grad) %**% t(diff_lam) +
               diff_lam %**% t(diff_grad) %**% Inv) / denom
          restart <- FALSE
        }

        ## Standard BFGS update
        if (iter > 3) {
          diff_grad <- gradient - prev_gradient
          diff_lam  <- cbind(lambda - prev_lambda)
          denom     <- as.numeric(t(diff_grad) %**% diff_lam)
          if (!is.na(denom) && abs(denom) > 1e-64) {
            rho   <- 1 / denom
            term1 <- Id - rho * (diff_lam %**% t(diff_grad))
            term2 <- Id - rho * (cbind(diff_grad) %**% t(diff_lam))
            Inv   <- term1 %**% Inv %**% term2 +
              rho * (diff_lam %**% t(diff_lam))
          } else {
            ## Reset if numerically unstable
            Inv     <- diag(length(gradient))
            restart <- TRUE
          }
        }
      }

      ## Compute BFGS step direction
      step_direction <- Inv %**% cbind(gradient)
    }

    ## Evaluate GCV at new point
    if (parallel_bfgs && !is.null(env$cl) && n_workers > 1L) {
      damp_set <- .bfgs_damp_set(damp, n_workers)
      candidate_list <- lapply(damp_set, function(damp_j) {
        list(damp = damp_j,
             lambda = lambda - damp_j * step_direction)
      })

      valid <- vapply(candidate_list, function(cand) {
        all(is.finite(exp(cand$lambda))) &&
          !any(is.nan(exp(cand$lambda))) &&
          !any(is.na(exp(cand$lambda)))
      }, logical(1))

      if (!any(valid)) {
        lambda       <- old_lambda
        gradient     <- old_gradient
        dont_skip_gr <- FALSE
        damp         <- min(damp_set) / 2
        if (damp < 2^-12 && iter > 9) {
          return(list(par        = c(best_lambda),
                      criterion_value = best_criterion_value,
                      gcv_u      = best_criterion_value,
                      iterations = iter))
        }
        next
      }

      env_eval <- .tuning_eval_env(env)
      prev_outlist <- outlist
      eval_results <- parallel::parLapply(
        env$cl,
        candidate_list[valid],
      function(cand,
               criterion_fxn,
               log_penalty_vec,
               env_eval,
               include_warnings,
               ...) {
          out <- tryCatch({
            res <- criterion_fxn(c(cand$lambda[1], cand$lambda[2],
                                   log_penalty_vec),
                                 log_penalty_vec,
                                 env_eval,
                                 ...)
            crit <- res$criterion_value
            if (length(crit) != 1L || is.na(crit) ||
                is.nan(crit) || !is.finite(crit)) {
              res$criterion_value <- Inf
            }
            res
          }, error = function(e) {
            if (include_warnings) {
              warning(conditionMessage(e))
            }
            list(criterion_value = Inf)
          })
          list(damp = cand$damp,
               lambda = cand$lambda,
               criterion_value = out$criterion_value,
               outlist = out)
        },
        criterion_fxn = criterion_fxn,
        log_penalty_vec = log_penalty_vec,
        env_eval = env_eval,
        include_warnings = env$include_warnings,
        ...
      )

      eval_values <- vapply(eval_results, function(res) {
        crit <- res$criterion_value
        if (length(crit) != 1L || is.na(crit) || !is.finite(crit)) {
          Inf
        } else {
          crit
        }
      }, numeric(1))

      if (iter <= 2) {
        improving <- which(is.finite(eval_values))
      } else {
        improving <- which(is.finite(eval_values) &
                             eval_values <= criterion_value)
      }

      if (length(improving) > 0L) {
        best_idx <- improving[which.min(eval_values[improving])[1]]
        best_eval <- eval_results[[best_idx]]
        new_lambda <- best_eval$lambda
        new_criterion_value <- best_eval$criterion_value
        outlist <- best_eval$outlist
      } else {
        dont_skip_gr <- FALSE
        outlist      <- prev_outlist
        damp         <- min(damp_set) / 2
        if (damp < 2^(-10) && iter > 0) {
          return(list(par        = c(best_lambda),
                      criterion_value = best_criterion_value,
                      gcv_u      = best_criterion_value,
                      iterations = iter))
        }
        next
      }

    } else {
      new_lambda <- lambda - damp * step_direction

      if (any(is.na(new_lambda))) {
        lambda       <- old_lambda
        gradient     <- old_gradient
        dont_skip_gr <- FALSE
        damp         <- damp / 2
        if (damp < 2^-12 && iter > 9) {
          return(list(par        = c(best_lambda),
                      criterion_value = best_criterion_value,
                      gcv_u      = best_criterion_value,
                      iterations = iter))
        }
        next
      } else {
        if (any(!is.finite(exp(new_lambda))) ||
            any(is.nan(exp(new_lambda))) ||
            any(is.na(exp(new_lambda)))) {
          new_lambda <- best_lambda
        }
        prev_outlist <- outlist
        outlist      <- criterion_fxn(c(new_lambda[1], new_lambda[2],
                                        log_penalty_vec),
                                      log_penalty_vec, criterion_env, ...)
        new_criterion_value <- outlist$criterion_value
        if (any(is.na(new_criterion_value))) {
          new_criterion_value <- criterion_value
        }
      }
    }

    ## Accept step if improvement or early iterations
    if (new_criterion_value <= criterion_value || iter <= 2) {
      ## Update solution history
      old_gradient  <- prev_gradient
      prev_gradient <- gradient
      dont_skip_gr  <- TRUE
      prev_outlist  <- outlist
      .tuning_store_active_set(env, outlist$qp_info)
      old_lambda    <- prev_lambda
      prev_lambda   <- lambda
      lambda        <- new_lambda
      damp          <- 1

      ## Track best solution globally, not just relative to the most recent
      # accepted iterate. Early damped steps can be exploratory and may be
      # accepted even when they are worse than the grid-search start.
      if (new_criterion_value <= best_criterion_value) {
        best_criterion_value <- new_criterion_value
        best_lambda <- lambda
      }

      ## Check convergence
      criterion_change <- abs(new_criterion_value - criterion_value)
      parameter_change <- max(abs(lambda - prev_lambda))
      adaptive_criterion_tol <- .tuning_adaptive_criterion_tol(
        tol, best_criterion_value
      )
      if (criterion_change < adaptive_criterion_tol) {
        small_change_count <- small_change_count + 1L
      } else {
        small_change_count <- 0L
      }

      if (((criterion_change < tol) ||
           (parameter_change < tol) ||
           (small_change_count >= small_change_patience)) &&
          (iter > 9)) {
        return(list(par        = c(best_lambda),
                    criterion_value = best_criterion_value,
                    gcv_u      = best_criterion_value,
                    iterations = iter))
      }
      criterion_value <- new_criterion_value
    } else {
      ## Reject step and backtrack
      dont_skip_gr <- FALSE
      outlist      <- prev_outlist
      damp         <- damp / 2
      if (damp < 2^(-10) && iter > 0) {
        return(list(par        = c(best_lambda),
                    criterion_value = best_criterion_value,
                    gcv_u      = best_criterion_value,
                    iterations = iter))
      }
    }
  }

  return(list(par        = c(best_lambda),
              criterion_value = best_criterion_value,
              gcv_u      = best_criterion_value,
              iterations = iter))
}


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# Subfunction: .tune_grid_search
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

#' Grid Search Initialization for Penalty Tuning
#'
#' @description
#' Evaluates the selected tuning criterion over a grid of initial wiggle and ridge
#' penalty values to find a good starting point for BFGS optimization.
#'
#' @param log_initial_wiggle Numeric vector; log-scale candidate values for
#'   the wiggle penalty.
#' @param log_initial_flat Numeric vector; log-scale candidate values for
#'   the flat ridge penalty.
#' @param log_penalty_vec Numeric vector; log-scale predictor/partition
#'   penalties.
#' @param criterion_fxn Function; tuning-objective evaluation function.
#' @param env List; tuning environment.
#' @param include_warnings Logical; whether to print warnings on failure.
#' @param ... Additional arguments passed to \code{criterion_fxn}.
#'
#' @return Numeric vector of length 2; the best (log_wiggle, log_flat) found.
#'
#' @keywords internal
.tune_grid_search <- function(log_initial_wiggle,
                              log_initial_flat,
                              log_penalty_vec,
                              criterion_fxn,
                              env,
                              include_warnings,
                              parallel_grideval = TRUE,
                              cl = NULL,
                              ...) {

  ## Create all combinations of the grid values
  initial_grid <- expand.grid(wiggle = log_initial_wiggle,
                              flat   = log_initial_flat)
  n_workers <- .tuning_worker_count(cl)
  initial_grid <- .expand_tuning_grid(initial_grid,
                                      log_initial_wiggle,
                                      log_initial_flat,
                                      n_workers,
                                      parallel_grideval)

  ## Function to safely evaluate the selected criterion
  safe_gcvu <- function(par) {
    .safe_tuning_value(par,
                       criterion_fxn,
                       log_penalty_vec,
                       env,
                       include_warnings,
                       ...)
  }

  ## Evaluate the criterion for each grid point
  if (parallel_grideval && !is.null(cl) && n_workers > 1L) {
    env_eval <- .tuning_eval_env(env)
    gcv_values <- unlist(parallel::parLapply(
      cl,
      seq_len(nrow(initial_grid)),
      function(i,
               initial_grid,
               criterion_fxn,
               log_penalty_vec,
               env_eval,
               include_warnings,
               ...) {
        tryCatch({
          out <- criterion_fxn(c(unlist(initial_grid[i, ]),
                                 log_penalty_vec),
                               log_penalty_vec,
                               env_eval,
                               ...)
          crit <- out$criterion_value
          if (length(crit) != 1L || is.na(crit) ||
              is.nan(crit) || !is.finite(crit)) {
            return(Inf)
          }
          crit
        }, error = function(e) {
          if (include_warnings) {
            warning(conditionMessage(e))
          }
          Inf
        })
      },
      initial_grid = initial_grid,
      criterion_fxn = criterion_fxn,
      log_penalty_vec = log_penalty_vec,
      env_eval = env_eval,
      include_warnings = include_warnings,
      ...
    ))
  } else {
    gcv_values <- apply(initial_grid, 1, safe_gcvu)
  }
  bads <- which(is.na(gcv_values) |
                  is.nan(gcv_values) |
                  !is.finite(gcv_values))

  if (length(bads) == length(gcv_values)) {
    stop("All tuning criteria for the initial grid were computed as NA,",
         " NaN, or non-finite: check your data for corrupt or missing values,",
         " try changing initial tuning grid, or try manual tuning instead.",
         " If you are setting no_intercept = TRUE, try experimenting with",
         " standardize_response = FALSE vs. TRUE.")
  } else if (length(bads) > 0) {
    gcv_values <- gcv_values[-c(bads)]
    initial_grid <- initial_grid[-c(bads), , drop = FALSE]
  }

  ## Find the best starting point
  best_index <- which.min(gcv_values)[1]
  best_start <- as.numeric(initial_grid[best_index, ])

  return(best_start)
}


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# Subfunction: .compute_tuning_delta
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

#' Compute Pseudocount Delta for Link Function Stabilization
#'
#' @description
#' Determines the pseudocount \eqn{\delta} used to stabilize link function
#' transformations during penalty tuning. For identity link or when
#' the response is naturally in the domain of the link function, returns 0.
#' Otherwise, finds the \eqn{\delta} that makes the transformed response
#' distribution most closely approximate a t-distribution.
#'
#' @param family GLM family object.
#' @param unl_y Numeric vector; unlisted response values (concatenated across
#'   partitions).
#' @param N_obs Integer; total sample size.
#' @param observation_weights List or NULL; observation weights by partition.
#' @param opt Logical; whether optimization is being performed.
#'
#' @return Numeric scalar; the pseudocount \eqn{\delta \geq 0}.
#'
#' @keywords internal
.compute_tuning_delta <- function(family, unl_y, N_obs,
                                  observation_weights, opt) {

  ## Check if delta is needed at all
  if (!opt) return(0)

  ## Identity link: no pseudocount needed
  if (paste0(family)[2] == "identity") return(0)

  ## Log link with all positive y: no pseudocount needed
  if (paste0(family)[2] == "log" && (min(unl_y) > 0)) return(0)

  ## Inverse/1/mu^2 link with no zeros: no pseudocount needed
  if (any(paste0(family)[2] %in% c("inverse", "1/mu^2")) &&
      (!any(unl_y == 0))) return(0)

  ## Logit link with no boundary values: no pseudocount needed
  if (paste0(family)[2] == "logit" && (!any(unl_y %in% c(0, 1)))) return(0)

  ## Otherwise, find optimal pseudocount via Brent optimization
  t_quants <- qt((seq(0, 1, len = N_obs + 2))[-c(1, N_obs + 2)], df = N_obs - 1)
  delta <- stats::optim(
    1 / 16,
    fn = function(par) {
      y_delta <- std(sort(family$linkfun((unl_y + par) / (1 + 2 * par))))
      if (length(unlist(observation_weights)) == length(unl_y)) {
        mean(abs(y_delta - t_quants) * unlist(observation_weights))
      } else {
        mean(abs(y_delta - t_quants))
      }
    },
    method = "Brent",
    lower = 1e-64,
    upper = 1
  )$par

  return(delta)
}


## Active-set warm start

.tuning_seed_active_set <- function(env) {
  if (is.null(env$active_set_cache) || is.null(env$qp_Amat)) {
    return(integer(0))
  }

  active_ineq <- env$active_set_cache$active_ineq
  if (is.null(active_ineq)) return(integer(0))

  active_ineq <- as.integer(active_ineq)
  active_ineq <- active_ineq[is.finite(active_ineq)]
  unique(active_ineq[active_ineq >= 1L & active_ineq <= ncol(env$qp_Amat)])
}


.tuning_store_active_set <- function(env, qp_info) {
  if (is.null(env$active_set_cache)) return(invisible(NULL))

  active_ineq <- if (!is.null(qp_info$active_ineq)) {
    qp_info$active_ineq
  } else {
    integer(0)
  }
  active_ineq <- as.integer(active_ineq)
  active_ineq <- active_ineq[is.finite(active_ineq)]
  if (!is.null(env$qp_Amat)) {
    active_ineq <- active_ineq[active_ineq >= 1L &
                                 active_ineq <= ncol(env$qp_Amat)]
  }

  env$active_set_cache$active_ineq <- unique(active_ineq)
  invisible(NULL)
}


.tuning_clear_active_set <- function(env) {
  if (!is.null(env$active_set_cache)) {
    env$active_set_cache$active_ineq <- integer(0)
  }
  invisible(NULL)
}


## Correlated tuning cache

.make_gee_tuning_precomp <- function(X, y, order_list, Vhalf, VhalfInv,
                                     observation_weights) {
  if (is.null(Vhalf) | is.null(VhalfInv)) {
    return(NULL)
  }

  perm <- unlist(order_list)
  VhalfInv_perm <- VhalfInv[perm, perm]
  X_block <- collapse_block_diagonal(X)
  y_block <- cbind(unlist(y))
  X_tilde <- VhalfInv_perm %**% X_block
  y_tilde <- VhalfInv_perm %**% y_block

  obs_wt <- c(unlist(observation_weights))
  if (length(obs_wt) == 0L) {
    obs_wt <- rep(1, nrow(X_tilde))
  } else if (length(obs_wt) == 1L) {
    obs_wt <- rep(obs_wt, nrow(X_tilde))
  }

  list(
    perm           = perm,
    VhalfInv_perm = VhalfInv_perm,
    X_block       = X_block,
    y_block       = y_block,
    X_tilde       = X_tilde,
    y_tilde       = y_tilde,
    obs_wt        = obs_wt,
    Gram_full     = crossprod(X_tilde, obs_wt * X_tilde),
    Xy_tilde      = crossprod(X_tilde, obs_wt * y_tilde)
  )
}


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# Subfunction: .build_tuning_env
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

#' Build Tuning Environment
#'
#' @description
#' Assembles a named list containing all pre-computed objects and configuration
#' needed by the tuning-objective evaluation and gradient functions during
#' penalty tuning.
#' This avoids deep nesting of closures and makes the dependencies explicit.
#'
#' @details
#' In addition to the standard tuning arguments, the environment stores two
#' pre-computed blockfit dispatch items:
#' \describe{
#'   \item{use_blockfit}{Logical; TRUE when \code{blockfit} is enabled,
#'     \code{flat_cols} is non-empty, and \code{K > 0}. Mirrors the dispatch
#'     logic in \code{lgspline.fit} so that the same fitting path is used
#'     during tuning as during the final fit.}
#'   \item{flat_cols}{Integer vector; column indices of non-interactive linear
#'     terms derived from \code{just_linear_without_interactions} and
#'     \code{colnm_expansions}. Pre-computed once here rather than re-derived
#'     at every tuning-objective evaluation.}
#' }
#'
#' @param flat_cols Integer vector; pre-computed flat column indices (passed
#'   in from \code{tune_Lambda} to avoid recomputation).
#' @param use_blockfit Logical; pre-computed dispatch flag (passed in from
#'   \code{tune_Lambda}).
#' @param gcv_gamma Numeric scalar, at least 1; modified-GCV multiplier for the
#'   effective-degrees-of-freedom term.
#'
#' @return Named list (the "tuning environment").
#'
#' @keywords internal
.build_tuning_env <- function(y, X, X_gram, Xy,
                              smoothing_spline_penalty,
                              A, R_constraints, K, p_expansions, N_obs,
                              custom_penalty_mat,
                              colnm_expansions,
                              unique_penalty_per_predictor,
                              unique_penalty_per_partition,
                              meta_penalty,
                              family, delta,
                              order_list,
                              observation_weights,
                              homogenous_weights,
                              parallel, parallel_eigen,
                              parallel_trace, parallel_aga,
                              parallel_matmult,
                              parallel_qr = FALSE,
                              qr_pivot_smoothing_constraints = TRUE,
                              parallel_unconstrained,
                              cl, chunk_size,
                              num_chunks, rem_chunks,
                              unconstrained_fit_fxn,
                              keep_weighted_Lambda,
                              iterate,
                              qp_score_function, quadprog,
                              qp_Amat, qp_bvec, qp_meq,
                              tol, sd_y, tuning_criterion, gcv_gamma,
                              constraint_value_vectors,
                              glm_weight_function,
                              schur_correction_function,
                              need_dispersion_for_estimation,
                              dispersion_function,
                              blockfit,
                              just_linear_without_interactions,
                              Vhalf, VhalfInv,
                              verbose, include_warnings,
                              flat_cols, use_blockfit) {

  gee_precomp <- .make_gee_tuning_precomp(X, y, order_list, Vhalf, VhalfInv,
                                          observation_weights)

  list(
    y                             = y,
    X                             = X,
    X_gram                        = X_gram,
    Xy                            = Xy,
    smoothing_spline_penalty      = smoothing_spline_penalty,
    A                             = A,
    R_constraints                 = R_constraints,
    K                             = K,
    p_expansions                  = p_expansions,
    N_obs                         = N_obs,
    custom_penalty_mat            = custom_penalty_mat,
    colnm_expansions              = colnm_expansions,
    unique_penalty_per_predictor  = unique_penalty_per_predictor,
    unique_penalty_per_partition  = unique_penalty_per_partition,
    meta_penalty                  = meta_penalty,
    family                        = family,
    delta                         = delta,
    order_list                    = order_list,
    observation_weights           = observation_weights,
    homogenous_weights            = homogenous_weights,
    parallel                      = parallel,
    parallel_eigen                = parallel_eigen,
    parallel_trace                = parallel_trace,
    parallel_aga                  = parallel_aga,
    parallel_matmult              = parallel_matmult,
    parallel_qr                   = parallel_qr,
    qr_pivot_smoothing_constraints = qr_pivot_smoothing_constraints,
    parallel_unconstrained        = parallel_unconstrained,
    cl                            = cl,
    chunk_size                    = chunk_size,
    num_chunks                    = num_chunks,
    rem_chunks                    = rem_chunks,
    unconstrained_fit_fxn         = unconstrained_fit_fxn,
    keep_weighted_Lambda          = keep_weighted_Lambda,
    iterate                       = iterate,
    qp_score_function             = qp_score_function,
    quadprog                      = quadprog,
    qp_Amat                       = qp_Amat,
    qp_bvec                       = qp_bvec,
    qp_meq                        = qp_meq,
    tol                           = tol,
    sd_y                          = sd_y,
    tuning_criterion              = tuning_criterion,
    gcv_gamma                     = gcv_gamma,
    constraint_value_vectors      = constraint_value_vectors,
    glm_weight_function           = glm_weight_function,
    schur_correction_function     = schur_correction_function,
    need_dispersion_for_estimation = need_dispersion_for_estimation,
    dispersion_function           = dispersion_function,
    blockfit                      = blockfit,
    just_linear_without_interactions = just_linear_without_interactions,
    Vhalf                         = Vhalf,
    VhalfInv                      = VhalfInv,
    verbose                       = verbose,
    include_warnings              = include_warnings,
    flat_cols                     = flat_cols,
    use_blockfit                  = use_blockfit,
    gee_precomp                   = gee_precomp,
    active_set_cache              = new.env(parent = emptyenv())
  )
}


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# Main function: tune_Lambda (exported)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

#' Tune Smoothing and Ridge Penalties via Leave-One-Out or Generalized Cross-Validation
#'
#' @description
#' Tunes the penalty matrix by minimizing either exact leave-one-out on the
#' transformed problem or modified GCV. This is the top-level tuning entry
#' point: it assembles \eqn{\boldsymbol{\Lambda}}, initializes from a small
#' grid search, and then refines the penalty parameters by quasi-Newton
#' optimization with either closed-form or finite-difference gradients.
#'
#' @param y List; response vectors by partition.
#' @param X List; design matrices by partition.
#' @param X_gram List; Gram matrices by partition.
#' @param smoothing_spline_penalty Matrix; integrated squared second derivative
#'   penalty.
#' @param A Matrix; smoothness constraints at knots.
#' @param K Integer; number of interior knots in 1-D, number of partitions - 1
#'   in higher dimensions.
#' @param p_expansions Integer; columns per partition.
#' @param N_obs Integer; total sample size.
#' @param opt Logical; TRUE to optimize penalties, FALSE to use initial values.
#' @param use_custom_bfgs Logical; TRUE for analytic gradient BFGS as natively
#'   implemented, FALSE for finite differences as implemented by
#'   \code{stats::optim()}. The analytic path is typically faster; however, for
#'   exact LOO tuning the leverage derivative can be numerically delicate in
#'   some problems, in which case \code{use_custom_bfgs = FALSE} remains a
#'   conservative fallback.
#' @param C Matrix; polynomial expansion matrix (used for initialization).
#'   This is the monomial expansion design used to derive starting values and
#'   is not the inequality-constraint matrix sometimes denoted by
#'   \eqn{\mathbf{C}} elsewhere in the package documentation.
#' @param colnm_expansions Character vector; column names of the expansion
#'   matrix.
#' @param wiggle_penalty,flat_ridge_penalty Fixed penalty values if provided.
#' @param initial_wiggle,initial_flat Numeric vectors; candidate values for grid
#'   search initialization on the raw (non-negative) scale. Converted to log
#'   scale internally for optimization.
#' @param unique_penalty_per_predictor,unique_penalty_per_partition Logical;
#'   allow predictor/partition-specific penalties.
#' @param penalty_vec Numeric vector; initial values for predictor/partition
#'   penalties on the raw (non-negative) scale. Converted to log scale
#'   internally for optimization. Use \code{c()} when no per-predictor or
#'   per-partition penalties are needed.
#' @param meta_penalty The "meta" ridge penalty, a regularization for
#'   predictor/partition penalties to pull them towards 1 on the raw scale.
#' @param family GLM family with optional custom tuning loss.
#' @param unconstrained_fit_fxn Function for unconstrained fitting.
#' @param keep_weighted_Lambda,iterate Logical controlling GLM fitting.
#' @param qp_score_function,quadprog,qp_Amat,qp_bvec,qp_meq Quadratic
#'   programming parameters (see arguments of \code{\link[lgspline]{lgspline}}).
#' @param tol Numeric; convergence tolerance. Used directly for strict
#'   absolute criterion and log-parameter checks, and as the floor for the
#'   adaptive plateau check in the custom BFGS path.
#' @param sd_y,delta Response standardization parameters.
#' @param tuning_criterion Character scalar selecting the tuning objective:
#'   \code{"loo"} for exact leave-one-out on the transformed tuning problem or
#'   \code{"gcv"} for generalized cross-validation. The LOO criterion itself is
#'   still evaluated exactly on that transformed problem; the caveat is that
#'   its fully analytic leverage derivative can be numerically sensitive in some
#'   datasets. For very large samples, \code{"gcv"} is often the more practical
#'   choice; as a rough guideline, it is recommended once the sample size is
#'   above about 250,000.
#' @param gcv_gamma Numeric scalar, at least 1, used only when
#'   \code{tuning_criterion = "gcv"}. Multiplies the effective degrees of
#'   freedom in the GCV denominator.
#' @param constraint_value_vectors List; constraint values.
#' @param parallel Logical; enable parallel computation.
#' @param parallel_eigen,parallel_trace,parallel_aga Logical; specific parallel
#'   flags.
#' @param parallel_matmult,parallel_qr,parallel_bfgs,parallel_grideval,
#'   qr_pivot_smoothing_constraints,parallel_unconstrained Logical; specific
#'   parallel flags. When \code{parallel_grideval = TRUE} or
#'   \code{parallel_bfgs = TRUE}, inner parallel flags such as
#'   \code{parallel_eigen} are ignored for those tuning fits.
#' @param cl Parallel cluster object.
#' @param chunk_size,num_chunks,rem_chunks Integer; parallel computation
#'   parameters.
#' @param shared_env Environment; shared variables exported to cluster workers.
#' @param custom_penalty_mat Optional custom penalty matrix.
#' @param order_list List; observation ordering by partition.
#' @param glm_weight_function,schur_correction_function Functions for GLM
#'   weights and corrections.
#' @param need_dispersion_for_estimation,dispersion_function Control dispersion
#'   estimation.
#' @param observation_weights Optional observation weights.
#' @param homogenous_weights Logical; TRUE if all weights equal.
#' @param blockfit Logical; when TRUE, the backfitting block decomposition
#'   (\code{blockfit_solve}) is used in place of \code{get_B} for the
#'   coefficient estimation step at each tuning-objective evaluation, provided
#'   \code{just_linear_without_interactions} is non-empty and \code{K > 0}.
#'   The dispatch condition is identical to that in \code{lgspline.fit} so
#'   that penalties are tuned under the same estimator that will be used for
#'   the final fit. Falls back to \code{get_B} automatically on failure.
#' @param just_linear_without_interactions Numeric; vector of columns for
#'   non-spline effects without interactions.
#' @param Vhalf,VhalfInv Square root and inverse square root correlation
#'   structure matrices. These are passed through to the coefficient-estimation
#'   step used inside each tuning-objective evaluation, so any dense or
#'   Woodbury-accelerated
#'   correlated solve is the same one used for the final fitted model.
#' @param verbose Logical; print progress.
#' @param include_warnings Logical; print warnings/try-errors.
#' @param ... Additional arguments passed to fitting functions.
#'
#' @return List containing:
#' \describe{
#'   \item{Lambda}{Baseline per-partition penalty matrix corresponding to
#'     \eqn{\lambda_w(\boldsymbol{\Lambda}_s +
#'     \lambda_r\boldsymbol{\Lambda}_r +
#'     \sum_{l=1}^{L}\lambda_{l,k}\boldsymbol{\Lambda}_{l,k})}.}
#'   \item{flat_ridge_penalty}{Optimized ridge penalty.}
#'   \item{wiggle_penalty}{Optimized smoothing penalty.}
#'   \item{other_penalties}{Optimized additional penalties
#'     \eqn{\lambda_{l,k}}.}
#'   \item{L_predictor_list}{Implementation lists storing additional
#'     penalty matrices \eqn{\boldsymbol{\Lambda}_{l,k}} for
#'     predictor-specific directions.}
#'   \item{L_partition_list}{Implementation lists storing additional
#'     penalty matrices \eqn{\boldsymbol{\Lambda}_{l,k}} for
#'     partition-specific directions.}
#' }
#'
#' @details
#' The tuning procedure consists of the following steps:
#' \enumerate{
#'   \item \strong{Preprocessing}: Convert raw-scale penalties to log scale,
#'     compute cross-products, determine pseudocount delta, ensure constraint
#'     matrix compatibility.
#'   \item \strong{Blockfit dispatch}: Pre-compute \code{flat_cols} and the
#'     \code{use_blockfit} flag so that every tuning-objective evaluation uses the same
#'     coefficient estimator as the final fit. \code{flat_cols} are identified
#'     by matching column names against
#'     \code{paste0("_", just_linear_without_interactions, "_")} in
#'     \code{colnm_expansions}.
#'   \item \strong{Grid search}: Evaluate the selected tuning criterion over a grid of
#'     (wiggle, ridge) penalty candidates to find a good starting point
#'     (see \code{.tune_grid_search}). When \code{parallel_grideval = TRUE}
#'     and a cluster is supplied, these candidate fits are spread across the
#'     workers, and if more than six workers are available the grid is
#'     enlarged with additional random raw-scale candidates drawn from
#'     \eqn{[\min/10,\max\times 10]} in each penalty direction.
#'   \item \strong{BFGS optimization}: Minimize the selected tuning criterion via either the custom
#'     damped BFGS with closed-form gradients (see \code{.damped_bfgs},
#'     \code{.compute_gcvu_gradient}, \code{.compute_loocv_gradient}) or base R's \code{stats::optim} with
#'     finite-difference gradients. The analytic paths differentiate the fitted
#'     criterion once to obtain
#'     \eqn{\mathbf{M}_k=\partial\ell/\partial\boldsymbol{\Lambda}_k};
#'     all log-penalty gradients then reduce to blockwise trace products
#'     with the corresponding current penalty direction.
#'     For LOO, the leverage derivative may be numerically delicate in some
#'     problems and can motivate the finite-difference fallback. After the
#'     first 10 iterations, the custom BFGS path also stops on a small
#'     scale-aware plateau in accepted criterion values. When
#'     \code{parallel_bfgs = TRUE} and a cluster is supplied, each damped step
#'     evaluates a batch of halved step sizes in parallel and takes the best
#'     improving candidate from that batch. If inequality constraints are
#'     present, the active set from the last accepted coefficient solve is used
#'     as the initial working set for the next tuning evaluation.
#'   \item \strong{Criterion choice}: \code{"loo"} uses exact per-observation
#'     leverages on the transformed tuning problem, while \code{"gcv"} uses the
#'     usual trace-based denominator with optional \code{gcv_gamma}
#'     inflation.
#'   \item \strong{Sample-size adjustment}: After optimization, divide the tuned
#'     penalties by \eqn{(N+1)/(N-1)} (equivalently multiply by
#'     \eqn{(N-1)/(N+1)}) for both GCV- and LOO-based tuning so that the
#'     final penalties are decreased at small sample sizes.
#'   \item \strong{Final Lambda}: Compute the final penalty matrix from
#'     optimized parameters via \code{compute_Lambda}.
#' }
#'
#' Parameterization: initial penalty values are accepted on the raw
#' (non-negative) scale and converted to natural log-scale internally,
#' i.e. raw_penalty = exp(theta), so that raw penalties are always positive.
#' The chain rule factor d(exp(theta))/d(theta) = exp(theta) = raw_penalty.
#'
#' If \eqn{\mathbf{L}_{j,k}} is the current contribution of
#' \eqn{\theta_j = \log\lambda_j} to partition \eqn{k}, then
#' \deqn{
#'   \frac{\partial \ell}{\partial \theta_j}
#'   =
#'   \sum_{k=0}^{K}
#'   \mathrm{tr}\{\mathbf{M}_k\mathbf{L}_{j,k}\}.
#' }
#' This avoids repeating the full derivative through
#' \eqn{\mathbf{G}}, \eqn{\mathbf{U}}, and \eqn{\mathbf{H}} for each penalty.
#'
#' The resulting penalty follows the same notation used in the paper:
#' \deqn{
#'   \boldsymbol{\Lambda}_k
#'   =
#'   \lambda_w\Bigl(\boldsymbol{\Lambda}_s +
#'   \lambda_r\boldsymbol{\Lambda}_r +
#'   \sum_{l=1}^{L}\lambda_{l,k}\boldsymbol{\Lambda}_{l,k}\Bigr),
#' }
#' with full penalty
#' \deqn{
#'   \boldsymbol{\Lambda}
#'   =
#'   \mathrm{blockdiag}(\boldsymbol{\Lambda}_0, \ldots,
#'   \boldsymbol{\Lambda}_K).
#' }
#' The current package stores these pieces using the implementation labels
#' \code{L1}, \code{L2}, \code{L_predictor_list}, and
#' \code{L_partition_list}, corresponding to
#' \eqn{\boldsymbol{\Lambda}_s}, \eqn{\boldsymbol{\Lambda}_r}, and the
#' additional \eqn{\boldsymbol{\Lambda}_{l,k}} components.
#'
#' If correlation structure inputs are supplied, each tuning-objective
#' evaluation calls the same constrained correlated solver used by the final
#' model fit. In the structured-correlation case, the Woodbury reduction
#' described in \code{lgspline-details} is inherited automatically through
#' \code{get_B}; there is no separate tuning-specific notation or solver.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[stats]{optim}} for Hessian-free optimization
#'   \item \code{\link{compute_Lambda}} for penalty matrix construction
#'   \item \code{\link{compute_G_eigen}} for eigendecomposition of penalized
#'     Gram matrices
#'   \item \code{\link{get_B}} for constrained coefficient estimation
#'   \item \code{\link{blockfit_solve}} for the backfitting block-decomposition
#'     estimator used when \code{blockfit = TRUE}
#' }
#'
#' @examples
#' \dontrun{
#' ## ## Example 1: Basic usage within lgspline ## ## ## ## ## ## ## ## ## ## ##
#' ## tune_Lambda is called internally by lgspline; direct calls are
#' #  for advanced users. Here we verify that the refactored version
#' #  produces identical output to the original.
#'
#' set.seed(1234)
#' t <- runif(200, -5, 5)
#' y <- sin(t) + rnorm(200, 0, 0.5)
#'
#' ## Fit with automatic penalty tuning (calls tune_Lambda internally)
#' fit1 <- lgspline(t, y, K = 3)
#' cat("Wiggle penalty:", fit1$penalties$wiggle_penalty, "\n")
#' cat("Ridge penalty:", fit1$penalties$flat_ridge_penalty, "\n")
#' cat("Trace (edf):", fit1$trace_XUGX, "\n")
#'
#' ## ## Example 2: Fixed penalties (no tuning) ## ## ## ## ## ## ## ## ## ## ##
#' fit2 <- lgspline(t, y, K = 3, opt = FALSE,
#'                  wiggle_penalty = 1e-4,
#'                  flat_ridge_penalty = 0.1)
#' cat("Fixed wiggle:", fit2$penalties$wiggle_penalty, "\n")
#'
#' ## ## Example 3: blockfit path #### ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#' ## When blockfit = TRUE and just_linear_without_interactions is non-empty,
#' #  tune_Lambda dispatches to blockfit_solve at each GCV evaluation,
#' #  ensuring penalties are tuned under the same estimator used in the final
#' #  fit. Verify that tuned penalties are consistent across both paths.
#'
#' set.seed(7)
#' n  <- 300
#' t1 <- runif(n, 0, 5)
#' t2 <- rnorm(n)
#' y2 <- sin(t1) + 0.5 * t2 + rnorm(n, 0, 0.3)
#' df <- data.frame(t1 = t1, t2 = t2)
#'
#' ## blockfit = TRUE uses blockfit_solve during tuning
#' fit_bf  <- lgspline(df, y2, K = 2, blockfit = TRUE,
#'                     just_linear_without_interactions = 2)
#' ## blockfit = FALSE uses get_B during tuning (original path)
#' fit_std <- lgspline(df, y2, K = 2, blockfit = FALSE,
#'                     just_linear_without_interactions = 2)
#'
#' cat("blockfit wiggle  :", fit_bf$penalties$wiggle_penalty, "\n")
#' cat("standard wiggle  :", fit_std$penalties$wiggle_penalty, "\n")
#' ## Penalties may differ slightly; predictions should be close.
#' cat("Max pred diff:", max(abs(fit_bf$ytilde - fit_std$ytilde)), "\n")
#'
#' ## ## Example 4: Verify refactored subfunctions ## ## ## ## ## ## ## ## ## ##
#' ## The internal .compute_meta_penalty should match hand calculation
#' mp <- lgspline:::.compute_meta_penalty(
#'   wiggle_penalty = 0.5,
#'   penalty_vec = c(predictor1 = 1.2, partition1 = 0.8),
#'   meta_penalty_coef = 1e-8,
#'   unique_penalty_per_predictor = TRUE,
#'   unique_penalty_per_partition = TRUE
#' )
#' expected <- 0.5 * 1e-8 * ((1.2 - 1)^2 + (0.8 - 1)^2) +
#'             0.5 * 1e-32 * (0.5 - 1)^2
#' stopifnot(abs(mp - expected) < 1e-20)
#' cat("Meta-penalty check passed.\n")
#'
#' ## ## Example 5: Verify gradient of meta-penalty ## ## ## ## ## ## ## ## ## #
#' gr <- lgspline:::.compute_meta_penalty_gradient(
#'   wiggle_penalty = 2.0,
#'   penalty_vec = c(predictor1 = 1.5),
#'   meta_penalty_coef = 1e-8,
#'   unique_penalty_per_predictor = TRUE,
#'   unique_penalty_per_partition = FALSE
#' )
#' # gr[1] should be 1e-32 * (2 - 1) * 2 = 2e-32
#' # gr[2] should be 0
#' # gr[3] should be 1e-8 * (1.5 - 1) * 1.5 = 7.5e-9
#' stopifnot(abs(gr[1] - 2e-32) < 1e-40)
#' stopifnot(gr[2] == 0)
#' stopifnot(abs(gr[3] - 7.5e-9) < 1e-17)
#' cat("Meta-penalty gradient check passed.\n")
#'
#' ## ## Example 6: Residual computation paths ## ##
#' ## Identity link ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#' r1 <- lgspline:::.compute_tuning_residuals(
#'   y = list(c(1, 2, 3)),
#'   preds = list(c(1.1, 1.9, 3.2)),
#'   delta = 0,
#'   family = gaussian(),
#'   observation_weights = list(NULL),
#'   K = 0,
#'   order_list = list(1:3)
#' )
#' stopifnot(max(abs(r1[[1]] - c(-0.1, 0.1, -0.2))) < 1e-10)
#' cat("Residual check passed.\n")
#' }
#'
#' @keywords internal
#' @export
tune_Lambda <- function(
    y,
    X,
    X_gram,
    smoothing_spline_penalty,
    A,
    K,
    p_expansions,
    N_obs,
    opt,
    use_custom_bfgs,
    C,
    colnm_expansions,
    wiggle_penalty,
    flat_ridge_penalty,
    initial_wiggle,
    initial_flat,
    unique_penalty_per_predictor,
    unique_penalty_per_partition,
    penalty_vec,
    meta_penalty,
    family,
    unconstrained_fit_fxn,
    keep_weighted_Lambda,
    iterate,
    qp_score_function, quadprog, qp_Amat, qp_bvec, qp_meq,
    tol,
    sd_y,
    delta,
    tuning_criterion,
    gcv_gamma,
    constraint_value_vectors,
    parallel,
    parallel_eigen,
    parallel_trace,
    parallel_aga,
    parallel_matmult,
    parallel_qr,
    parallel_bfgs,
    parallel_grideval,
    qr_pivot_smoothing_constraints,
    parallel_unconstrained,
    cl,
    chunk_size,
    num_chunks,
    rem_chunks,
    shared_env,
    custom_penalty_mat,
    order_list,
    glm_weight_function,
    schur_correction_function,
    need_dispersion_for_estimation,
    dispersion_function,
    observation_weights,
    homogenous_weights,
    blockfit,
    just_linear_without_interactions,
    Vhalf,
    VhalfInv,
    verbose,
    include_warnings,
    ...
) {

  if (verbose) cat("    Starting tuning\n")

  tuning_criterion <- match.arg(tuning_criterion, c("loo", "gcv"))

  if (tuning_criterion == "gcv" &&
      (!is.numeric(gcv_gamma) || length(gcv_gamma) != 1 ||
       !is.finite(gcv_gamma) || gcv_gamma < 1)) {
    stop("gcv_gamma must be a finite numeric scalar >= 1. ",
         "Set gcv_gamma = 1 to recover ordinary GCV.", call. = FALSE)
  }

  ## Step 1: Convert raw-scale inputs to log scale ## ##
  log_initial_wiggle <- log(initial_wiggle)
  log_initial_flat   <- log(initial_flat)
  if (length(penalty_vec) > 0) {
    log_penalty_vec <- log(penalty_vec)
  } else {
    log_penalty_vec <- c()
  }

  if (verbose) cat("    Xy\n")

  ## Step 2: Precompute cross-products and constants
  sN_obs <- sqrt(N_obs)
  Xy  <- vectorproduct_block_diagonal(X, y, K)
  Xyr <- Reduce("rbind", Xy)
  R_constraints <- ncol(A)
  Xt  <- lapply(X, t)
  unl_y <- unlist(y)

  ## Step 3: Compute pseudocount delta
  if (verbose) cat("    Getting pseudocount delta\n")

  if (is.null(delta) &&
      any(is.null(family$custom_dev.resids)) &&
      opt) {
    delta <- .compute_tuning_delta(family, unl_y, N_obs,
                                   observation_weights, opt)
  }
  if (is.null(delta)) delta <- 0

  if (verbose) cat("    GCV, gradient, and BFGS Function prep\n")

  ## Step 4: ensure the constraint matrix is always conformable
  if (any(is.null(A))) {
    A   <- cbind(rep(0, (K + 1) * p_expansions))
    A   <- cbind(A, A)
    R_constraints <- 2
  }

  ## Step 5: pre-compute the blockfit dispatch used during tuning
  #  by matching the flat linear columns the same way lgspline.fit does.
  flat_cols <- c()
  if (length(just_linear_without_interactions) > 0) {
    flat_cols <- which(colnm_expansions %in%
                         paste0("_", just_linear_without_interactions, "_"))
  }
  use_blockfit <- blockfit && length(flat_cols) > 0 && K > 0

  if (verbose && use_blockfit) {
    cat("    blockfit dispatch active during tuning (",
        length(flat_cols), "flat cols )\n")
  }

  ## Step 6: build the tuning environment
  env <- .build_tuning_env(
    y = y, X = X, X_gram = X_gram, Xy = Xy,
    smoothing_spline_penalty = smoothing_spline_penalty,
    A = A, R_constraints = R_constraints, K = K,
    p_expansions = p_expansions, N_obs = N_obs,
    custom_penalty_mat = custom_penalty_mat,
    colnm_expansions = colnm_expansions,
    unique_penalty_per_predictor = unique_penalty_per_predictor,
    unique_penalty_per_partition = unique_penalty_per_partition,
    meta_penalty = meta_penalty,
    family = family, delta = delta,
    order_list = order_list,
    observation_weights = observation_weights,
    homogenous_weights = homogenous_weights,
    parallel = parallel, parallel_eigen = parallel_eigen,
    parallel_trace = parallel_trace, parallel_aga = parallel_aga,
    parallel_matmult = parallel_matmult,
    parallel_qr = parallel_qr,
    qr_pivot_smoothing_constraints = qr_pivot_smoothing_constraints,
    parallel_unconstrained = parallel_unconstrained,
    cl = cl, chunk_size = chunk_size,
    num_chunks = num_chunks, rem_chunks = rem_chunks,
    unconstrained_fit_fxn = unconstrained_fit_fxn,
    keep_weighted_Lambda = keep_weighted_Lambda,
    iterate = iterate,
    qp_score_function = qp_score_function, quadprog = quadprog,
    qp_Amat = qp_Amat, qp_bvec = qp_bvec, qp_meq = qp_meq,
    tol = tol, sd_y = sd_y,
    tuning_criterion = tuning_criterion,
    gcv_gamma = gcv_gamma,
    constraint_value_vectors = constraint_value_vectors,
    glm_weight_function = glm_weight_function,
    schur_correction_function = schur_correction_function,
    need_dispersion_for_estimation = need_dispersion_for_estimation,
    dispersion_function = dispersion_function,
    blockfit = blockfit,
    just_linear_without_interactions = just_linear_without_interactions,
    Vhalf = Vhalf, VhalfInv = VhalfInv,
    verbose = verbose, include_warnings = include_warnings,
    flat_cols = flat_cols, use_blockfit = use_blockfit
  )

  ## Step 7: Optimization (if requested)
  if (opt) {
    criterion_fxn <- if (tuning_criterion == "loo") {
      .compute_loocv
    } else {
      .compute_gcvu
    }
    gr_fxn <- if (tuning_criterion == "loo") {
      .compute_loocv_gradient
    } else {
      .compute_gcvu_gradient
    }

    if (verbose) cat("    Starting grid search for initialization\n")

    ## Grid search for good starting values
    best_start <- .tune_grid_search(log_initial_wiggle,
                                    log_initial_flat,
                                    log_penalty_vec,
                                    criterion_fxn,
                                    env,
                                    include_warnings,
                                    parallel_grideval = parallel_grideval,
                                    cl = cl,
                                    ...)

    if (verbose) {
      cat("    Finished grid evaluations\n")
      cat("    Best from grid search: ", cbind(c(best_start)), "\n")
    }
    .tuning_clear_active_set(env)

    ## Run optimization
    if (use_custom_bfgs) {
      ## Custom damped BFGS with closed-form gradients
      res <- withCallingHandlers(
        try(.damped_bfgs(c(best_start, log_penalty_vec),
                         log_penalty_vec,
                         criterion_fxn,
                         gr_fxn,
                         env,
                         tol,
                         parallel_bfgs = parallel_bfgs,
                         max_iter = 100,
                         ...),
            silent = TRUE),
        warning = function(w) {
          if (include_warnings) warning(w) else invokeRestart("muffleWarning")
        },
        message = function(m) {
          if (include_warnings) message(m) else invokeRestart("muffleMessage")
        }
      )
      if (any(inherits(res, "try-error"))) {
        if (include_warnings) print(res)
        if (include_warnings) {
          warning("Custom BFGS implementation failed. Try use_custom_bfgs =",
                  " FALSE, or manual tuning. Resorting to best as selected",
                  " from grid search.")
        }
        par <- c(best_start, log_penalty_vec)
      } else {
        par <- res$par
      }

    } else {
      ## Base R BFGS with finite-difference gradients
      optim_env <- if (!is.null(cl) && .tuning_worker_count(cl) > 1L) {
        .tuning_eval_env(env)
      } else {
        env
      }
      res <- withCallingHandlers(
        try({
          safe_optim_fn <- function(par) {
            crit <- .safe_tuning_value(par,
                                       criterion_fxn,
                                       log_penalty_vec,
                                       optim_env,
                                       include_warnings = FALSE,
                                       ...)
            if (!is.finite(crit) || is.na(crit) || is.nan(crit)) {
              return(.Machine$double.xmax^0.25)
            }
            crit
          }
          optim(c(best_start, log_penalty_vec),
                fn = safe_optim_fn,
                method = "BFGS")
        }, silent = TRUE),
        warning = function(w) {
          if (include_warnings) warning(w) else invokeRestart("muffleWarning")
        },
        message = function(m) {
          if (include_warnings) message(m) else invokeRestart("muffleMessage")
        }
      )
      if (any(inherits(res, "try-error"))) {
        if (include_warnings) print(res)
        if (include_warnings) {
          warning("Base R BFGS failed. Try use_custom_bfgs = TRUE, or manual",
                  " tuning. Resorting to best as selected from grid search.")
        }
        par <- c(best_start, log_penalty_vec)
      } else {
        par <- res$par
      }
    }

    if (verbose) cat("    Finished tuning penalties\n")

    ## Decrease tuned penalties at small sample sizes by dividing the
    #  optimizer solution by (N+1)/(N-1) for both tuning criteria.
    sample_size_adjustment <- (N_obs - 1) / (N_obs + 1)
    wiggle_penalty     <- exp(par[1]) * sample_size_adjustment
    flat_ridge_penalty <- exp(par[2]) * sample_size_adjustment

    ## Update penalty vec for predictor- and partition-specific penalties
    if (length(log_penalty_vec) > 0) {
      penalty_vec[1:length(penalty_vec)] <-
        exp(c(par[-c(1, 2)])) * sample_size_adjustment
    }

  } else if (length(penalty_vec) == 0) {
    ## No per-predictor or per-partition penalties; nothing to do
  }
  ## When !opt and length(penalty_vec) > 0, penalty_vec is already on raw
  #  scale from input and used as-is

  ## Step 8: Compute final penalty matrix
  if (verbose) cat("    Final update\n")

  Lambda_list <- compute_Lambda(custom_penalty_mat,
                                smoothing_spline_penalty,
                                wiggle_penalty,
                                flat_ridge_penalty,
                                K,
                                p_expansions,
                                unique_penalty_per_predictor,
                                unique_penalty_per_partition,
                                penalty_vec,
                                colnm_expansions,
                                just_Lambda = FALSE)

  return(list("Lambda"            = Lambda_list$Lambda,
              "flat_ridge_penalty" = flat_ridge_penalty,
              "wiggle_penalty"     = wiggle_penalty,
              "other_penalties"    = penalty_vec,
              "L_predictor_list"   = Lambda_list$L_predictor_list,
              "L_partition_list"   = Lambda_list$L_partition_list))
}
