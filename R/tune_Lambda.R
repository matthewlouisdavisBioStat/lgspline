## tune_Lambda.R
## GCV-based tuning for the smoothing and ridge penalties.
##
## Main pieces:
##   tune_Lambda()
#    .compute_gcvu() and .compute_gcvu_gradient()
#    .damped_bfgs() and .tune_grid_search()
#    small helpers for predictions, residuals, and meta-penalties


## Residuals used in GCV_u

#' Compute Residuals for GCV Criterion During Penalty Tuning
#'
#' @description
#' Computes residuals used in the numerator of the GCV criterion. Handles
#' identity link, general GLM link functions with pseudocount delta, custom
#' deviance residual functions, and observation weights.
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


## Predictions used in GCV_u

#' Compute Predictions During Penalty Tuning
#'
#' @description
#' Wrapper around \code{matmult_block_diagonal} for computing partition-wise
#' predictions \eqn{\mathbf{X}_{k} \boldsymbol{\beta}_{k}} during GCV
#' penalty tuning.
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

#' Fit Coefficients During GCV Tuning: blockfit_solve or get_B
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
#' These are needed immediately after this call for \code{AGAmult_wrapper},
#' \code{GXX}, and the trace computation.
#'
#' When correlation structure inputs are supplied, this wrapper does not
#' introduce a separate tuning-specific notation or solver path. Instead, the
#' same correlated coefficient estimator used in the final model fit is called
#' here inside each GCV evaluation. In particular, any structured-correlation
#' Woodbury correction is handled inside \code{get_B} and documented in
#' \code{lgspline-details}.
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

#' Call get_B During GCV Tuning
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
  ## keep_G = TRUE always: G is needed for the trace computation below
  ## regardless of whether get_B or blockfit_solve is used for fitting.
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
  ## .fit_coefficients, with automatic fallback on failure.
  ## return_G_getB = TRUE: B_list$G_list carries the (possibly GLM-iterated)
  ## G matrices needed for AGAmult_wrapper and trace computations below.
  if (verbose) cat("        gcvu_fxn fit coefficients\n")
  return_G_getB <- TRUE
  B_list <- .fit_coefficients(G_list, Lambda,
                              Lambda_list$L_partition_list,
                              env, return_G_getB, ...)
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
                           env$parallel,
                           env$cl,
                           env$chunk_size,
                           env$num_chunks,
                           env$rem_chunks,
                           env$unique_penalty_per_partition,
                           Lambda_list$L_partition_list)

  ## Predictions
  if (verbose) cat("        gcvu_fxn get predictions\n")
  preds <- .compute_tuning_predictions(env$X, B, env$K,
                                       env$parallel & env$parallel_matmult,
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

  return(list(GCV_u       = GCV_u + mp,
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
              AGAInv      = AGAInv))
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
#' For predictor- and partition-specific penalties, a trace-ratio heuristic
#' is used:
#' \deqn{\frac{\partial \mathrm{GCV}_u}{\partial \lambda_{j}}
#'   \approx \frac{\mathrm{tr}(\mathbf{L}_{j})}{\mathrm{tr}(\boldsymbol{\Lambda})}
#'   \cdot \frac{\partial \mathrm{GCV}_u}{\partial \lambda_{w}}}
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
    ## return_G_getB = TRUE: B_list$G_list carries the (possibly GLM-iterated)
    ## G matrices needed for AGAmult_wrapper and trace computations below.
    if (verbose) cat("        gr fxn fit coefficients\n")
    return_G_getB <- TRUE
    B_list <- .fit_coefficients(G_list, Lambda,
                                Lambda_list$L_partition_list,
                                env, return_G_getB, ...)
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
                             env$parallel,
                             env$cl,
                             env$chunk_size,
                             env$num_chunks,
                             env$rem_chunks,
                             env$unique_penalty_per_partition,
                             Lambda_list$L_partition_list)

    if (verbose) cat("        gr fxn preds\n")
    preds <- .compute_tuning_predictions(env$X, B, env$K,
                                         env$parallel & env$parallel_matmult,
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
    outlist <- list(GCV_u       = GCV_u + mp,
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
                    AGAInv      = AGAInv)
  }

  ## Differentiate the same fitted objects above so we do not have to
  #  refit again just to get the gradient.

  ## Key intermediate for the derivative calculations.
  if (verbose) cat("        GhalfXy_temp_list \n")
  GhalfXy_temp_list <- compute_GhalfXy_temp_wrapper(
    outlist$G_list$G,
    outlist$G_list$Ghalf,
    env$A,
    outlist$AGAInv,
    env$Xy,
    env$p_expansions,
    env$K,
    env$parallel & env$parallel_aga,
    env$cl,
    env$chunk_size,
    env$num_chunks,
    env$rem_chunks)
  GhalfXy_temp <- GhalfXy_temp_list[[1]]
  AGAInvAGXy   <- GhalfXy_temp_list[[2]]

  ## dG/dlambda
  if (verbose) cat("        compute_dG_dlambda \n")
  dG_dlambda <- compute_dG_dlambda(outlist$G_list$G,
                                   outlist$Lambda,
                                   env$K,
                                   lambda_1,
                                   env$unique_penalty_per_partition,
                                   outlist$L_partition_list,
                                   env$parallel & env$parallel_matmult,
                                   env$cl,
                                   env$chunk_size,
                                   env$num_chunks,
                                   env$rem_chunks)

  ## dG^{1/2}/dlambda
  if (verbose) cat("        Compute dGhalf \n")
  dGhalf <- compute_dGhalf(dG_dlambda,
                           env$p_expansions,
                           env$K,
                           env$parallel & env$parallel_eigen,
                           env$cl,
                           env$chunk_size,
                           env$num_chunks,
                           env$rem_chunks)

  ## d(U G)/dlambda * X^T y
  if (verbose) cat("        compute_dG_u_dlambda_xy \n")
  dG_u_dlambda1_Xyr <- compute_dG_u_dlambda_xy(
    AGAInvAGXy,
    outlist$AGAInv,
    outlist$G_list$G,
    env$A,
    dG_dlambda,
    env$p_expansions,
    env$R_constraints,
    env$K,
    env$Xy,
    outlist$G_list$Ghalf,
    dGhalf,
    GhalfXy_temp,
    env$parallel & env$parallel_matmult,
    env$cl,
    env$chunk_size,
    env$num_chunks,
    env$rem_chunks)

  ## dW/dlambda (trace derivative)
  if (verbose) cat("        compute_dW_dlambda_wrapper \n")
  dW_dlambda <- compute_dW_dlambda_wrapper(
    outlist$G_list$G,
    env$A,
    outlist$GXX,
    outlist$G_list$Ghalf,
    dG_dlambda,
    dGhalf,
    outlist$AGAInv,
    env$p_expansions,
    env$K,
    env$parallel & env$parallel_matmult,
    env$cl,
    env$chunk_size,
    env$num_chunks,
    env$rem_chunks)

  ## -2 * residuals^T * X
  if (verbose) cat("        neg2tresidX\n")
  neg2tresidX <- Reduce("cbind",
                        matmult_block_diagonal(
                          lapply(outlist$residuals,
                                 function(r) -2 * t(r)),
                          env$X,
                          env$K,
                          env$parallel & env$parallel_matmult,
                          env$cl,
                          env$chunk_size,
                          env$num_chunks,
                          env$rem_chunks))

  ## Quotient rule for GCV_u derivative w.r.t. wiggle penalty
  dnumerator_dlambda1   <- c(neg2tresidX %**% dG_u_dlambda1_Xyr)
  ddenominator_dlambda1 <- 2 * (1 - gamma * outlist$mean_W) *
    (-gamma * dW_dlambda)
  dGCV_u_dlambda1 <- (dnumerator_dlambda1 * outlist$denominator -
                        outlist$numerator * ddenominator_dlambda1) /
    outlist$denom_sq

  ## Ridge gradient from the same trace-ratio approximation.
  dGCV_u_dlambda2 <- mean(diag(outlist$L2)) /
    mean(diag(outlist$Lambda)) *
    dGCV_u_dlambda1

  ## Apply the exp-parameterization chain rule.
  if (verbose) cat("        Gradient start \n")
  gradient <- cbind(c(dGCV_u_dlambda1 * lambda_1,
                      dGCV_u_dlambda2 * lambda_2))

  ## The remaining penalty gradients follow the wiggle derivative
  #  through the same trace-ratio approximation.
  if (env$unique_penalty_per_predictor) {
    predictor_penalties <- penalty_vec[grep("predictor", names(penalty_vec))]
    predictor_penalty_gradient <- sapply(
      seq_along(predictor_penalties), function(j) {
        mean(diag(outlist$L_predictor_list[[j]])) /
          mean(diag(outlist$Lambda)) *
          dGCV_u_dlambda1 *
          predictor_penalties[j]
      })
    gradient <- cbind(c(c(gradient), predictor_penalty_gradient))
  }

  ## Partition-specific penalty gradients
  if (env$unique_penalty_per_partition) {
    partition_penalties <- penalty_vec[grep("partition", names(penalty_vec))]
    partition_penalty_gradient <- sapply(
      seq_along(partition_penalties), function(j) {
        mean(diag(outlist$L_partition_list[[j]])) /
          mean(diag(outlist$Lambda +
                      outlist$L_partition_list[[j]])) *
          dGCV_u_dlambda1 *
          partition_penalties[j]
      })
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

  return(list(GCV_u    = outlist$GCV_u + mp,
              gradient = env$N_obs * gradient + regularizer,
              outlist  = outlist))
}


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# Subfunction: .damped_bfgs
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

#' Damped BFGS Optimizer for Modified-GCV Penalty Tuning
#'
#' @description
#' Custom implementation of damped BFGS quasi-Newton optimization for
#' minimizing the modified GCV_u criterion. Uses step-size damping with
#' backtracking and Sherman-Morrison-Woodbury inverse Hessian updates.
#'
#' @param par Numeric vector; initial log-scale penalty parameters (first two
#'   elements are log(wiggle) and log(flat_ridge)).
#' @param log_penalty_vec Numeric vector; log-scale predictor/partition
#'   penalties appended to the optimization vector.
#' @param gcvu_fxn Function; GCV_u evaluation function with signature
#'   \code{function(par, log_penalty_vec, env, ...)}.
#' @param gr_fxn Function; gradient function with signature
#'   \code{function(par, log_penalty_vec, outlist, env, ...)}.
#' @param env List; tuning environment (passed through to gcvu_fxn and
#'   gr_fxn).
#' @param tol Numeric; convergence tolerance for both GCV_u change and
#'   parameter change.
#' @param max_iter Integer; maximum number of BFGS iterations (default 100).
#' @param ... Additional arguments passed to fitting functions.
#'
#' @return List containing:
#' \describe{
#'   \item{par}{Numeric vector; best log-scale penalty parameters found.}
#'   \item{gcv_u}{Numeric; best GCV_u value achieved.}
#'   \item{iterations}{Integer; number of iterations performed.}
#' }
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
#' }
#'
#' @keywords internal
.damped_bfgs <- function(par,
                         log_penalty_vec,
                         gcvu_fxn,
                         gr_fxn,
                         env,
                         tol,
                         max_iter = 100,
                         ...) {

  ## Initialize optimization parameters and storage
  lambda       <- par
  old_lambda   <- lambda
  new_lambda   <- lambda
  n_params     <- length(lambda)
  Id           <- diag(n_params)
  initial_damp <- 0.5
  damp         <- initial_damp
  prev_gradient <- lambda * 0
  gradient      <- prev_gradient
  best_gradient <- Inf
  old_gradient  <- prev_gradient
  outlist       <- NULL
  prev_outlist  <- NULL
  gcv_u         <- Inf
  ridge         <- NULL
  dont_skip_gr  <- TRUE
  rho           <- NULL
  Inv           <- diag(n_params)
  best_lambda   <- lambda
  best_gcv_u    <- gcv_u
  prev_lambda   <- old_lambda
  restart       <- TRUE

  ## Main optimization loop
  for (iter in 1:max_iter) {

    ## Compute gradient if needed
    if (dont_skip_gr) {
      result   <- gr_fxn(c(lambda[1], lambda[2], log_penalty_vec),
                         log_penalty_vec, outlist, env, ...)
      gradient <- result$gradient
      outlist  <- result$outlist
    }

    ## First two iterations: steepest descent
    if (iter <= 2) {
      new_lambda <- lambda - damp * gradient

      ## Reset to best if numerical issues
      if (any(!is.finite(exp(new_lambda))) ||
          any(is.nan(exp(new_lambda))) ||
          any(is.na(exp(new_lambda)))) {
        new_lambda <- best_lambda
      }

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
      new_lambda <- lambda - damp * Inv %**% cbind(gradient)

      ## Reset to best if numerical issues
      if (any(!is.finite(exp(new_lambda))) ||
          any(is.nan(exp(new_lambda))) ||
          any(is.na(exp(new_lambda)))) {
        new_lambda <- best_lambda
      }
    }

    ## Evaluate GCV at new point
    if (any(is.na(new_lambda))) {
      ## Backtrack if invalid step
      lambda       <- old_lambda
      gradient     <- old_gradient
      dont_skip_gr <- FALSE
      damp         <- damp / 2
      if (damp < 2^-12 && iter > 9) {
        return(list(par        = best_lambda,
                    gcv_u      = best_gcv_u,
                    iterations = iter))
      }
      next
    } else {
      prev_outlist <- outlist
      outlist      <- gcvu_fxn(c(new_lambda[1], new_lambda[2],
                                 log_penalty_vec),
                               log_penalty_vec, env, ...)
      new_gcv_u    <- outlist$GCV_u
      if (any(is.na(new_gcv_u))) {
        new_gcv_u <- gcv_u
      }
    }

    ## Accept step if improvement or early iterations
    if (new_gcv_u <= gcv_u || iter <= 2) {
      ## Update solution history
      old_gradient  <- prev_gradient
      prev_gradient <- gradient
      dont_skip_gr  <- TRUE
      prev_outlist  <- outlist
      old_lambda    <- prev_lambda
      prev_lambda   <- lambda
      lambda        <- new_lambda
      damp          <- 1

      ## Track best solution
      if (new_gcv_u <= gcv_u || iter == 1) {
        best_gcv_u  <- new_gcv_u
        best_lambda <- lambda
      }

      ## Check convergence
      if (((abs(new_gcv_u - gcv_u) < tol) ||
           (max(abs(lambda - prev_lambda)) < tol)) &&
          (iter > 9)) {
        return(list(par        = best_lambda,
                    gcv_u      = best_gcv_u,
                    iterations = iter))
      }
      gcv_u <- new_gcv_u
    } else {
      ## Reject step and backtrack
      dont_skip_gr <- FALSE
      outlist      <- prev_outlist
      damp         <- damp / 2
      if (damp < 2^(-10) && iter > 0) {
        return(list(par        = best_lambda,
                    gcv_u      = best_gcv_u,
                    iterations = iter))
      }
    }
  }

  return(list(par        = best_lambda,
              gcv_u      = best_gcv_u,
              iterations = iter))
}


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# Subfunction: .tune_grid_search
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

#' Grid Search Initialization for Penalty Tuning
#'
#' @description
#' Evaluates the GCV_u criterion over a grid of initial wiggle and ridge
#' penalty values to find a good starting point for BFGS optimization.
#'
#' @param log_initial_wiggle Numeric vector; log-scale candidate values for
#'   the wiggle penalty.
#' @param log_initial_flat Numeric vector; log-scale candidate values for
#'   the flat ridge penalty.
#' @param log_penalty_vec Numeric vector; log-scale predictor/partition
#'   penalties.
#' @param gcvu_fxn Function; GCV_u evaluation function.
#' @param env List; tuning environment.
#' @param include_warnings Logical; whether to print warnings on failure.
#' @param ... Additional arguments passed to gcvu_fxn.
#'
#' @return Numeric vector of length 2; the best (log_wiggle, log_flat) found.
#'
#' @keywords internal
.tune_grid_search <- function(log_initial_wiggle,
                              log_initial_flat,
                              log_penalty_vec,
                              gcvu_fxn,
                              env,
                              include_warnings,
                              ...) {

  ## Create all combinations of the grid values
  initial_grid <- expand.grid(wiggle = log_initial_wiggle,
                              flat   = log_initial_flat)

  ## Function to safely evaluate gcv_u
  safe_gcvu <- function(par) {
    tryCatch({
      result <- gcvu_fxn(c(unlist(par), log_penalty_vec),
                         log_penalty_vec, env, ...)$GCV_u
      if (is.na(result) || is.nan(result)) {
        return(Inf)
      }
      return(result)
    }, error = function(e) {
      if (include_warnings) {
        return(Inf)
      }
      return(Inf)
    })
  }

  ## Evaluate GCV_u for each grid point
  gcv_values <- apply(initial_grid, 1, safe_gcvu)
  bads <- which(is.na(gcv_values) |
                  is.nan(gcv_values) |
                  !is.finite(gcv_values))

  if (length(bads) == length(gcv_values)) {
    stop("All GCV criteria for the initial tuning grid were computed as NA,",
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
#' transformations during GCV penalty tuning. For identity link or when
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


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# Subfunction: .build_tuning_env
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

#' Build Tuning Environment
#'
#' @description
#' Assembles a named list containing all pre-computed objects and configuration
#' needed by the GCV_u evaluation and gradient functions during penalty tuning.
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
#'     at every GCV evaluation.}
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
                              parallel_unconstrained,
                              cl, chunk_size,
                              num_chunks, rem_chunks,
                              unconstrained_fit_fxn,
                              keep_weighted_Lambda,
                              iterate,
                               qp_score_function, quadprog,
                               qp_Amat, qp_bvec, qp_meq,
                               tol, sd_y, gcv_gamma,
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
    use_blockfit                  = use_blockfit
  )
}


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# Main function: tune_Lambda (exported)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

#' Tune Smoothing and Ridge Penalties via Generalized Cross Validation
#'
#' @description
#' Optimizes smoothing spline and ridge regression penalties by minimizing the
#' GCV criterion. Uses BFGS optimization with analytical gradients or finite
#' differences. This is the top-level entry point that orchestrates grid search
#' initialization and quasi-Newton optimization via internal subfunctions.
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
#'   \code{stats::optim()}.
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
#' @param tol Numeric; convergence tolerance.
#' @param sd_y,delta Response standardization parameters.
#' @param gcv_gamma Numeric scalar, at least 1. Multiplies the effective
#'   degrees of freedom in the GCV denominator so that tuning minimizes the
#'   modified criterion
#'   \eqn{\sum r_i^2 / \{N(1 - \gamma \bar{W})^2\}}.
#'   \code{gcv_gamma = 1} recovers ordinary GCV; values above 1 penalize
#'   complexity more strongly and reduce the occasional severe undersmoothing
#'   of unmodified GCV. The default \code{1.4} follows Kim and Gu (2004,
#'   Section 4, equation 4.1), who recommend values in the range 1.2 to 1.4.
#' @param constraint_value_vectors List; constraint values.
#' @param parallel Logical; enable parallel computation.
#' @param parallel_eigen,parallel_trace,parallel_aga Logical; specific parallel
#'   flags.
#' @param parallel_matmult,parallel_unconstrained Logical; specific parallel
#'   flags.
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
#'   coefficient estimation step at each GCV evaluation, provided
#'   \code{just_linear_without_interactions} is non-empty and \code{K > 0}.
#'   The dispatch condition is identical to that in \code{lgspline.fit} so
#'   that penalties are tuned under the same estimator that will be used for
#'   the final fit. Falls back to \code{get_B} automatically on failure.
#' @param just_linear_without_interactions Numeric; vector of columns for
#'   non-spline effects without interactions.
#' @param Vhalf,VhalfInv Square root and inverse square root correlation
#'   structure matrices. These are passed through to the coefficient-estimation
#'   step used inside each GCV evaluation, so any dense or Woodbury-accelerated
#'   correlated solve is the same one used for the final fitted model.
#' @param verbose Logical; print progress.
#' @param include_warnings Logical; print warnings/try-errors.
#' @param ... Additional arguments passed to fitting functions.
#'
#' @return List containing:
#' \describe{
#'   \item{Lambda}{Final combined penalty matrix.}
#'   \item{flat_ridge_penalty}{Optimized ridge penalty.}
#'   \item{wiggle_penalty}{Optimized smoothing penalty.}
#'   \item{other_penalties}{Optimized predictor/partition penalties.}
#'   \item{L_predictor_list}{Predictor-specific penalty matrices.}
#'   \item{L_partition_list}{Partition-specific penalty matrices.}
#' }
#'
#' @details
#' The tuning procedure consists of the following steps:
#' \enumerate{
#'   \item \strong{Preprocessing}: Convert raw-scale penalties to log scale,
#'     compute cross-products, determine pseudocount delta, ensure constraint
#'     matrix compatibility.
#'   \item \strong{Blockfit dispatch}: Pre-compute \code{flat_cols} and the
#'     \code{use_blockfit} flag so that every GCV evaluation uses the same
#'     coefficient estimator as the final fit. \code{flat_cols} are identified
#'     by matching column names against
#'     \code{paste0("_", just_linear_without_interactions, "_")} in
#'     \code{colnm_expansions}.
#'   \item \strong{Grid search}: Evaluate GCV_u over a grid of
#'     (wiggle, ridge) penalty candidates to find a good starting point
#'     (see \code{.tune_grid_search}).
#'   \item \strong{BFGS optimization}: Minimize GCV_u via either the custom
#'     damped BFGS with closed-form gradients (see \code{.damped_bfgs},
#'     \code{.compute_gcvu_gradient}) or base R's \code{stats::optim} with
#'     finite-difference gradients.
#'   \item \strong{Modified GCV}: Use
#'     \eqn{N(1 - \gamma\bar{W})^2} in the denominator, where
#'     \eqn{\gamma = \code{gcv_gamma}}. Setting \eqn{\gamma > 1}
#'     follows Kim and Gu's modified GCV recommendation and penalizes
#'     effective degrees of freedom more strongly.
#'   \item \strong{Inflation}: Apply small inflation factor
#'     \eqn{((N+2)/(N-2))^{2}} to counteract in-sample bias toward
#'     underpenalization.
#'   \item \strong{Final Lambda}: Compute the final penalty matrix from
#'     optimized parameters via \code{compute_Lambda}.
#' }
#'
#' Parameterization: initial penalty values are accepted on the raw
#' (non-negative) scale and converted to natural log-scale internally,
#' i.e. raw_penalty = exp(theta), so that raw penalties are always positive.
#' The chain rule factor d(exp(theta))/d(theta) = exp(theta) = raw_penalty.
#'
#' The resulting penalty matrix follows the same decomposition used throughout
#' the paper and package documentation:
#' \deqn{
#'   \boldsymbol{\Lambda}
#'   =
#'   \lambda_{w}\mathbf{L}_{1}
#'   +
#'   \lambda_{r}\mathbf{L}_{2}
#'   +
#'   \sum_{j}\nu_{j}\mathbf{L}^{(\mathrm{pred})}_{j}
#'   +
#'   \sum_{k}\tau_{k}\mathbf{L}^{(\mathrm{part})}_{k},
#' }
#' where the predictor- and partition-specific sums are included only when the
#' corresponding options are active. Internally these components are returned as
#' \code{L1}, \code{L2}, \code{L_predictor_list}, and \code{L_partition_list}.
#'
#' If correlation structure inputs are supplied, each GCV evaluation calls the
#' same constrained correlated solver used by the final model fit. In the
#' structured-correlation case, the low-rank Woodbury correction described in
#' \code{lgspline-details} is therefore inherited automatically through
#' \code{get_B}; no separate tuning-specific notation is introduced here.
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
#' ## ## Example 1: Basic usage within lgspline ## ##
#' ## tune_Lambda is called internally by lgspline; direct calls are
#' ## for advanced users. Here we verify that the refactored version
#' ## produces identical output to the original.
#'
#' set.seed(42)
#' t <- runif(200, -5, 5)
#' y <- sin(t) + rnorm(200, 0, 0.5)
#'
#' ## Fit with automatic penalty tuning (calls tune_Lambda internally)
#' fit1 <- lgspline(t, y, K = 3)
#' cat("Wiggle penalty:", fit1$penalties$wiggle_penalty, "\n")
#' cat("Ridge penalty:", fit1$penalties$flat_ridge_penalty, "\n")
#' cat("Trace (edf):", fit1$trace_XUGX, "\n")
#'
#' ## ## Example 2: Fixed penalties (no tuning) ## ##
#' fit2 <- lgspline(t, y, K = 3, opt = FALSE,
#'                  wiggle_penalty = 1e-4,
#'                  flat_ridge_penalty = 0.1)
#' cat("Fixed wiggle:", fit2$penalties$wiggle_penalty, "\n")
#'
#' ## ## Example 3: blockfit path ## ##
#' ## When blockfit = TRUE and just_linear_without_interactions is non-empty,
#' ## tune_Lambda dispatches to blockfit_solve at each GCV evaluation,
#' ## ensuring penalties are tuned under the same estimator used in the final
#' ## fit. Verify that tuned penalties are consistent across both paths.
#'
#' set.seed(7)
#' n  <- 300
#' x1 <- runif(n, 0, 5)
#' x2 <- rnorm(n)
#' y2 <- sin(x1) + 0.5 * x2 + rnorm(n, 0, 0.3)
#' df <- data.frame(x1 = x1, x2 = x2)
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
#' ## ## Example 4: Verify refactored subfunctions ## ##
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
#' ## ## Example 5: Verify gradient of meta-penalty ## ##
#' gr <- lgspline:::.compute_meta_penalty_gradient(
#'   wiggle_penalty = 2.0,
#'   penalty_vec = c(predictor1 = 1.5),
#'   meta_penalty_coef = 1e-8,
#'   unique_penalty_per_predictor = TRUE,
#'   unique_penalty_per_partition = FALSE
#' )
#' ## gr[1] should be 1e-32 * (2 - 1) * 2 = 2e-32
#' ## gr[2] should be 0
#' ## gr[3] should be 1e-8 * (1.5 - 1) * 1.5 = 7.5e-9
#' stopifnot(abs(gr[1] - 2e-32) < 1e-40)
#' stopifnot(gr[2] == 0)
#' stopifnot(abs(gr[3] - 7.5e-9) < 1e-17)
#' cat("Meta-penalty gradient check passed.\n")
#'
#' ## ## Example 6: Residual computation paths ## ##
#' ## Identity link
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
    gcv_gamma,
    constraint_value_vectors,
    parallel,
    parallel_eigen,
    parallel_trace,
    parallel_aga,
    parallel_matmult,
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

  if (!is.numeric(gcv_gamma) || length(gcv_gamma) != 1 ||
      !is.finite(gcv_gamma) || gcv_gamma < 1) {
    stop("gcv_gamma must be a finite numeric scalar >= 1. ",
         "Set gcv_gamma = 1 to recover ordinary GCV.", call. = FALSE)
  }

  ## ## Step 1: Convert raw-scale inputs to log scale ## ##
  log_initial_wiggle <- log(initial_wiggle)
  log_initial_flat   <- log(initial_flat)
  if (length(penalty_vec) > 0) {
    log_penalty_vec <- log(penalty_vec)
  } else {
    log_penalty_vec <- c()
  }

  if (verbose) cat("    Xy\n")

  ## ## Step 2: Precompute cross-products and constants ## ##
  sN_obs <- sqrt(N_obs)
  Xy  <- vectorproduct_block_diagonal(X, y, K)
  Xyr <- Reduce("rbind", Xy)
  R_constraints <- ncol(A)
  Xt  <- lapply(X, t)
  unl_y <- unlist(y)

  ## ## Step 3: Compute pseudocount delta ## ##
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
    parallel_unconstrained = parallel_unconstrained,
    cl = cl, chunk_size = chunk_size,
    num_chunks = num_chunks, rem_chunks = rem_chunks,
    unconstrained_fit_fxn = unconstrained_fit_fxn,
    keep_weighted_Lambda = keep_weighted_Lambda,
    iterate = iterate,
    qp_score_function = qp_score_function, quadprog = quadprog,
    qp_Amat = qp_Amat, qp_bvec = qp_bvec, qp_meq = qp_meq,
    tol = tol, sd_y = sd_y,
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

  ## ## Step 7: Optimization (if requested) ## ##
  if (opt) {

    if (verbose) cat("    Starting grid search for initialization\n")

    ## Grid search for good starting values
    best_start <- .tune_grid_search(log_initial_wiggle,
                                    log_initial_flat,
                                    log_penalty_vec,
                                    .compute_gcvu,
                                    env,
                                    include_warnings,
                                    ...)

    if (verbose) {
      cat("    Finished grid evaluations\n")
      cat("    Best from grid search: ", cbind(c(best_start)), "\n")
    }

    ## Run optimization
    if (use_custom_bfgs) {
      ## Custom damped BFGS with closed-form gradients
      res <- withCallingHandlers(
        try(.damped_bfgs(c(best_start, log_penalty_vec),
                         log_penalty_vec,
                         .compute_gcvu,
                         .compute_gcvu_gradient,
                         env,
                         tol,
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
      res <- withCallingHandlers(
        try({
          optim(c(best_start, log_penalty_vec),
                fn = function(par) {
                  .compute_gcvu(par, log_penalty_vec, env, ...)$GCV_u
                },
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

    ## Inflate penalties to counteract in-sample bias toward underpenalization
    infl               <- ((N_obs + 2) / ((N_obs - 2)))^2
    wiggle_penalty     <- exp(par[1]) * infl
    flat_ridge_penalty <- exp(par[2]) * infl

    ## Update penalty vec for predictor- and partition-specific penalties
    if (length(log_penalty_vec) > 0) {
      penalty_vec[1:length(penalty_vec)] <- exp(c(par[-c(1, 2)])) * infl
    }

  } else if (length(penalty_vec) == 0) {
    ## No per-predictor or per-partition penalties; nothing to do
  }
  ## When !opt and length(penalty_vec) > 0, penalty_vec is already on raw
  ## scale from input and used as-is

  ## ## Step 8: Compute final penalty matrix ## ##
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
