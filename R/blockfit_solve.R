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
#' Update flat coefficients via pooled penalized regression on residuals:
#'
#' \deqn{
#' \mathbf{v}
#' =
#' \left(
#' \sum_{k=0}^{K}
#' \mathbf{X}_{\mathrm{flat}}^{(k)\top}
#' \mathbf{X}_{\mathrm{flat}}^{(k)}
#' +
#' \mathbf{\Lambda}_{\mathrm{flat}}
#' \right)^{-1}
#' \sum_{k=0}^{K}
#' \mathbf{X}_{\mathrm{flat}}^{(k)\top}
#' \left(
#' \mathbf{y}_k
#' -
#' \mathbf{Z}_k
#' \boldsymbol{\beta}_{\mathrm{spline}}^{(k)}
#' \right).
#' }
#' }
#'
#' For Gaussian identity link, backfitting is block coordinate descent on
#' a convex quadratic objective and converges rapidly.
#'
#' For GLMs without correlation structures, an IRLS outer loop wraps the
#' backfitting, with each IRLS step solving a weighted least squares
#' problem.
#'
#' When a working correlation structure is supplied via \code{Vhalf} and
#' \code{VhalfInv} (GEE estimation), the design matrices and response are
#' provided unwhitened; whitening is applied internally.
#'
#' For Gaussian identity link with GEE, a full-system closed-form solve
#' replaces naive backfitting because the whitened design is not
#' block-diagonal.
#'
#' For GLM GEE (non-identity link or non-Gaussian family), the Gaussian
#' backfitting solution provides a warm start. A damped sequential
#' quadratic programming (SQP) refinement loop handles nonlinear IRLS
#' iterations using the whitened design
#'
#' \deqn{
#' \mathbf{X}_{\mathrm{tilde}}
#' =
#' \mathbf{V}^{-1/2}
#' \mathbf{X}_{\mathrm{block}}.
#' }
#'
#' When \code{quadprog = TRUE}, a damped SQP refinement loop is applied
#' after backfitting convergence. The loop operates in full block-diagonal
#' form and enforces:
#'
#' \itemize{
#'   \item equality constraints via \eqn{\mathbf{A}}
#'   \item user-supplied inequality constraints via
#'     \code{qp_Amat}, \code{qp_bvec}, and \code{qp_meq}
#' }
#'
#' using \code{quadprog::solve.QP}.
#'
#' For GLM GEE with \code{quadprog = TRUE}, inequality constraints are
#' incorporated into the GEE SQP loop so that nonlinear link refinement
#' and inequality enforcement are handled simultaneously.
#'
#' After convergence, coefficients are reassembled into the standard
#' per-partition format (flat coefficients replicated across partitions)
#' for compatibility with downstream inference.
#'
#' The full constraint matrix \eqn{\mathbf{A}} from
#' \code{lgspline.fit} still contains flat-equality constraints; these are
#' trivially satisfied and the projection matrix \eqn{\mathbf{U}} handles
#' them correctly.
#'
#' @param X List of \eqn{K+1} design matrices
#'   (\eqn{N_k \times nc}). Unwhitened even when GEE.
#' @param y List of \eqn{K+1} response vectors. Unwhitened even when GEE.
#' @param flat_cols Integer vector indicating flat columns of
#'   \eqn{\mathbf{X}^{(k)}}.
#' @param K Integer; number of interior knots.
#' @param nc Integer; number of coefficients per partition.
#' @param Lambda \eqn{nc \times nc} penalty matrix.
#' @param L_partition_list List of partition-specific penalty matrices.
#' @param unique_penalty_per_partition Logical.
#' @param A Full \eqn{P \times nca} constraint matrix.
#' @param nca Integer; number of columns of \eqn{\mathbf{A}}.
#' @param constraint_values List of constraint right-hand sides.
#' @param X_gram List of Gram matrices
#'   \eqn{\mathbf{X}_k^{\top} \mathbf{X}_k}.
#' @param Ghalf_full,GhalfInv_full Lists of
#'   \eqn{\mathbf{G}^{1/2}} and \eqn{\mathbf{G}^{-1/2}} matrices.
#' @param family GLM family object.
#' @param order_list List of observation index vectors by partition.
#' @param glm_weight_function GLM weight function.
#' @param shur_correction_function Schur complement correction function.
#' @param need_dispersion_for_estimation Logical.
#' @param dispersion_function Dispersion estimation function.
#' @param observation_weights List of observation weights.
#' @param homogenous_weights Logical.
#' @param iterate Logical; if FALSE, single pass (no IRLS).
#' @param tol Convergence tolerance.
#' @param parallel_eigen,cl,chunk_size,num_chunks,rem_chunks Parallel arguments.
#' @param return_G_getB Logical.
#' @param quadprog Logical; apply SQP refinement if TRUE.
#' @param qp_Amat Inequality constraint matrix for
#'   \code{quadprog::solve.QP}.
#' @param qp_bvec Inequality constraint right-hand side.
#' @param qp_meq Number of leading equality constraints.
#' @param qp_score_function Score function for QP subproblem.
#' @param keep_weighted_Lambda Logical.
#' @param max_backfit_iter Integer.
#' @param Vhalf Square root of working correlation matrix.
#' @param VhalfInv Inverse square root of working correlation matrix.
#' @param include_warnings Logical.
#' @param verbose Logical.
#' @param ... Additional arguments passed to weight and dispersion functions.
#'
#' @return A list with elements:
#' \describe{
#'   \item{B}{
#'     List of \eqn{K+1} coefficient vectors
#'     (\eqn{nc \times 1}), flat coefficients replicated across partitions.
#'   }
#'   \item{G_list}{
#'     List containing \eqn{\mathbf{G}},
#'     \eqn{\mathbf{G}^{1/2}}, and
#'     \eqn{\mathbf{G}^{-1/2}},
#'     each a list of \eqn{K+1}
#'     \eqn{nc \times nc} matrices,
#'     or NULL if \code{return_G_getB} is FALSE.
#'
#'     \eqn{\mathbf{G}} satisfies
#'     \eqn{\mathbf{G}
#'     =
#'     \mathbf{G}^{1/2}
#'     (\mathbf{G}^{1/2})^{\top}} exactly.
#'   }
#'   \item{qp_info}{
#'     NULL when \code{quadprog = FALSE}. Otherwise a list with elements
#'     \code{solution} (final coefficient vector in block form),
#'     \code{lagrangian} (Lagrange multipliers from the last
#'     \code{solve.QP} call; one per column of the combined constraint
#'     matrix, ordered equalities first then inequalities),
#'     \code{active_constraints} (integer indices into the combined
#'     constraint matrix for constraints active at convergence, i.e.
#'     those with non-negligible multipliers or binding inequalities),
#'     \code{iact} (raw \code{iact} vector returned by
#'     \code{solve.QP}; 1-based indices of active inequality constraints),
#'     \code{info_matrix} (the penalized information matrix \code{Dmat}
#'     used in the final QP subproblem, unscaled),
#'     \code{Amat_combined} (combined equality/inequality constraint
#'     matrix passed to \code{solve.QP}),
#'     \code{bvec_combined} (corresponding right-hand side vector),
#'     \code{meq_combined} (number of leading equality constraints),
#'     \code{converged} (logical; TRUE if the loop exited via the
#'     tolerance criterion rather than the iteration/damp limit),
#'     and \code{final_deviance} (deviance value at convergence).
#'     Intended for use in constructing the variance-covariance matrix
#'     under active constraints.
#'   }
#' }
#'
#' @keywords internal
#' @export
blockfit_solve <- function(
    X,
    y,
    flat_cols,
    K,
    nc,
    Lambda,
    L_partition_list,
    unique_penalty_per_partition,
    A,
    nca,
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
  spline_cols <- setdiff(1:nc, flat_cols)
  nc_spline <- length(spline_cols)
  nr <- sum(sapply(y, length))
  is_gauss_id <- (paste0(family)[1] == 'gaussian' &
                    paste0(family)[2] == 'identity')

  ## Detect GEE: correlation structure present
  has_corr <- !is.null(Vhalf) & !is.null(VhalfInv)

  ## GLM + GEE requires SQP refinement even without inequality constraints,
  #  because the nonlinear link function interacts with the whitening
  #  transform. For Gaussian identity + GEE, the closed-form full-system
  #  solve is exact, so SQP is only needed when quadprog = TRUE.
  needs_gee_glm_sqp <- has_corr & !is_gauss_id

  ## Split design matrices into spline and flat blocks
  X_spline <- lapply(X, function(Xk) Xk[, spline_cols, drop = FALSE])
  X_flat   <- lapply(X, function(Xk) Xk[, flat_cols,   drop = FALSE])

  ## Split penalty matrices
  Lambda_spline <- Lambda[spline_cols, spline_cols]
  Lambda_flat   <- Lambda[flat_cols,   flat_cols]
  L_part_spline <- lapply(L_partition_list, function(Lk){
    if(is.numeric(Lk) && length(Lk) == 1 && Lk == 0) return(0)
    Lk[spline_cols, spline_cols]
  })

  ## Extract spline-only constraint matrix from A
  flat_rows_all   <- unlist(lapply(0:K, function(k) k * nc + flat_cols))
  spline_rows_all <- setdiff(1:(nc * (K + 1)), flat_rows_all)
  A_spline <- A[spline_rows_all, , drop = FALSE]
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

  ## Spline-only constraint values (e.g. from no_intercept)
  if(length(constraint_values) > 0){
    constraint_values_spline <- lapply(constraint_values, function(cv){
      cv[spline_cols, , drop = FALSE]
    })
  } else {
    constraint_values_spline <- constraint_values
  }

  ## Compute spline-only G^{1/2} for the Lagrangian projection
  X_gram_spline <- lapply(1:(K + 1), function(k){
    crossprod(X_spline[[k]])
  })
  schur_zero <- lapply(1:(K + 1), function(k) 0)

  G_sp <- compute_G_eigen(X_gram_spline,
                          Lambda_spline,
                          K,
                          parallel_eigen,
                          cl,
                          chunk_size, num_chunks, rem_chunks,
                          gaussian(),
                          unique_penalty_per_partition,
                          L_part_spline,
                          keep_G = TRUE,
                          schur_zero)
  Ghalf_sp <- G_sp$Ghalf

  ## G^{1/2} A for Lagrangian projection
  GhalfA_sp <- Reduce("rbind", lapply(1:(K + 1), function(k){
    rows <- ((k - 1) * nc_spline + 1):(k * nc_spline)
    Ghalf_sp[[k]] %**% A_spline[rows, , drop = FALSE]
  }))

  ## Flat Gram matrix (pooled) and penalized inverse
  XfXf <- Reduce("+", lapply(1:(K + 1), function(k){
    crossprod(X_flat[[k]])
  }))
  XfXf_pen_inv <- invert(XfXf + Lambda_flat)

  ## Lagrangian projection helper for the spline-only subproblem
  .lagrangian_project <- function(Xy_adj, Ghalf_cur, GhalfA_cur,
                                  cv_spline){
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

  ## Initialize coefficients
  beta_flat <- rep(0, nc_flat)
  beta_spline <- lapply(1:(K + 1), function(k) cbind(rep(0, nc_spline)))

  ## Initialize qp_info to NULL; populated only when quadprog = TRUE
  qp_info <- NULL

  ## For correlated Gaussian response
  if(is_gauss_id & has_corr){
    ## Gaussian identity + GEE.
    #  Full-system closed-form solve in whitened space. Cannot use
    #  naive partition-wise backfitting because the whitened design
    #  V^{-1/2} X_block is not block-diagonal.
    perm_bf <- unlist(order_list)
    VhalfInv_bf <- VhalfInv[perm_bf, perm_bf]

    X_block_bf <- collapse_block_diagonal(X)
    y_block_bf <- cbind(unlist(y))

    ## Full whitened design and response
    X_tilde_bf <- VhalfInv_bf %**% X_block_bf
    y_tilde_bf <- VhalfInv_bf %**% y_block_bf

    ## Full Gram + penalty
    Gram_bf <- crossprod(X_tilde_bf)

    if(unique_penalty_per_partition){
      Lambda_block_bf <- Reduce("rbind", lapply(1:(K+1), function(k){
        Reduce("cbind", lapply(1:(K+1), function(j){
          if(j == k) Lambda + L_partition_list[[k]] else 0 * Lambda
        }))
      }))
    } else {
      Lambda_block_bf <- Reduce("rbind", lapply(1:(K+1), function(k){
        Reduce("cbind", lapply(1:(K+1), function(j){
          if(j == k) Lambda else 0 * Lambda
        }))
      }))
    }

    G_bf_inv <- Gram_bf + Lambda_block_bf
    G_bf <- invert(G_bf_inv)

    ## G^{1/2} via eigendecomposition of full P x P G
    eig_bf <- eigen(G_bf, symmetric = TRUE)
    vals_bf <- eig_bf$values
    vals_bf[vals_bf <= 0] <- 0
    G_bf_half <- eig_bf$vectors %**%
      (t(eig_bf$vectors) * sqrt(vals_bf))

    Xy_bf <- crossprod(X_tilde_bf, y_tilde_bf)

    ## Lagrangian projection in full P-space using the full A
    #  Uses the G^{1/2}r* residual trick (see details)
    y_star_bf <- G_bf_half %**% Xy_bf
    X_star_bf <- G_bf_half %**% A

    comp_sc <- 1 / sqrt(K + 1)
    resids_bf <- .lm.fit(
      X_star_bf * comp_sc,
      y_star_bf * comp_sc
    )$residuals / comp_sc

    if(length(constraint_values) > 0){
      if(any(unlist(constraint_values) != 0)){
        cv_full <- Reduce("rbind", constraint_values)
        preds_bf <- X_star_bf %**%
          (invert(crossprod(X_star_bf) * comp_sc) %**%
             crossprod(A, cv_full * comp_sc))
        resids_bf <- resids_bf + c(preds_bf)
      }
    }

    beta_block_bf <- G_bf_half %**% cbind(resids_bf)

    ## Unpack into per-partition form
    result <- lapply(1:(K+1), function(k){
      cbind(beta_block_bf[(k-1)*nc + 1:nc])
    })

    ## Extract spline and flat coefficients for consistency
    beta_spline <- lapply(result, function(b) b[spline_cols, , drop=FALSE])
    beta_flat   <- result[[1]][flat_cols]

    VhalfInv_bf <- NULL
    X_tilde_bf  <- NULL

  } else if(is_gauss_id & !has_corr){
    ## Gaussian identity without correlation. Standard backfitting
    #  on the block-diagonal system

    ## Backfitting loop for alternating updates using Gaussian idnetity
    # - Iteratively update spline (`beta_spline`) and flat (`beta_flat`) coefficients.
    # - First, adjust `y` for current flat effects and update spline coefficients.
    # - Then, compute residuals to update flat coefficients.
    # - Track maximum changes and stop if all updates are below `tol`.
    for(bf_iter in 1:max_backfit_iter){

      ## Spline step: adjust y for flat contribution
      Xy_adj <- lapply(1:(K + 1), function(k){
        y_adj_k <- y[[k]] - X_flat[[k]] %**% cbind(beta_flat)
        crossprod(X_spline[[k]], y_adj_k)
      })
      beta_spline_new <- .lagrangian_project(Xy_adj,
                                             Ghalf_sp,
                                             GhalfA_sp,
                                             constraint_values_spline)

      ## Flat step: pooled penalized regression on residuals
      Xfr <- Reduce("+", lapply(1:(K + 1), function(k){
        r_k <- y[[k]] - X_spline[[k]] %**% beta_spline_new[[k]]
        crossprod(X_flat[[k]], r_k)
      }))
      beta_flat_new <- c(XfXf_pen_inv %**% Xfr)

      ## Convergence
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

  } else if(has_corr){
    ## GLM + GEE warm start.
    #  Gaussian backfitting on unwhitened data as approximate warm
    #  start. The correct GLM estimates are obtained by the subsequent
    #  GEE SQP refinement loop (needs_gee_glm_sqp = TRUE).

    for(bf_iter in 1:max_backfit_iter){

      Xy_adj <- lapply(1:(K + 1), function(k){
        y_adj_k <- y[[k]] - X_flat[[k]] %**% cbind(beta_flat)
        crossprod(X_spline[[k]], y_adj_k)
      })
      beta_spline_new <- .lagrangian_project(Xy_adj,
                                             Ghalf_sp,
                                             GhalfA_sp,
                                             constraint_values_spline)

      Xfr <- Reduce("+", lapply(1:(K + 1), function(k){
        r_k <- y[[k]] - X_spline[[k]] %**% beta_spline_new[[k]]
        crossprod(X_flat[[k]], r_k)
      }))
      beta_flat_new <- c(XfXf_pen_inv %**% Xfr)

      spline_diff <- max(abs(unlist(beta_spline_new) -
                               unlist(beta_spline)))
      flat_diff   <- max(abs(beta_flat_new - beta_flat))
      beta_spline <- beta_spline_new
      beta_flat   <- beta_flat_new

      if(max(spline_diff, flat_diff) < tol) break
    }

  } else {
    ## GLM without GEE IRWLS outer loop with backfitting
    #  d mu / d eta helper
    .get_mu_eta <- function(eta, fam){
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

    IRWLS_err  <- Inf
    damp_cnt  <- 0
    disp_temp <- 1

    for(IRWLS_iter in 1:100){

      ## Full linear predictor
      eta <- unlist(lapply(1:(K + 1), function(k){
        c(X_spline[[k]] %**% beta_spline[[k]] +
            X_flat[[k]] %**% cbind(beta_flat))
      }))
      mu  <- family$linkinv(eta)

      ## d mu / d eta
      mu_eta <- .get_mu_eta(eta, family)
      mu_eta <- ifelse(abs(mu_eta) < sqrt(.Machine$double.eps),
                       sqrt(.Machine$double.eps), mu_eta)

      ## Working response
      y_vec <- unlist(y)
      z <- eta + (y_vec - mu) / mu_eta

      ## Dispersion
      if(need_dispersion_for_estimation){
        disp_temp <- dispersion_function(mu = mu,
                                         y = y_vec,
                                         order_indices = unlist(order_list),
                                         family = family,
                                         observation_weights =
                                           unlist(observation_weights),
                                         Vhalf_inv = VhalfInv,
                                         ...)
      }

      ## IRWLS weights
      W <- c(glm_weight_function(mu, y_vec,
                                 unlist(order_list),
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

      ## Weighted spline Gram matrices and G
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

      ## Inner backfitting on the weighted working response
      for(bf_iter in 1:max_backfit_iter){

        ## Spline step (weighted)
        Xy_adj <- lapply(1:(K + 1), function(k){
          z_adj_k <- z_list[[k]] - X_flat[[k]] %**% cbind(beta_flat)
          crossprod(X_spline[[k]] * W_list[[k]], cbind(z_adj_k))
        })
        beta_spline_new <- .lagrangian_project(Xy_adj,
                                               Ghalf_sp_w,
                                               GhalfA_sp_w,
                                               constraint_values_spline)

        ## Flat step (weighted)
        XfWXf <- Reduce("+", lapply(1:(K + 1), function(k){
          crossprod(X_flat[[k]] * sqrt(W_list[[k]]))
        }))
        XfWXf_pen_inv_w <- invert(XfWXf + Lambda_flat)
        Xfr <- Reduce("+", lapply(1:(K + 1), function(k){
          r_k <- z_list[[k]] - X_spline[[k]] %**% beta_spline_new[[k]]
          crossprod(X_flat[[k]] * W_list[[k]], cbind(r_k))
        }))
        beta_flat_new <- c(XfWXf_pen_inv_w %**% Xfr)

        flat_diff <- max(abs(beta_flat_new - beta_flat))
        beta_spline <- beta_spline_new
        beta_flat   <- beta_flat_new

        if(flat_diff < tol) break
      }

      ## IRWLS convergence via deviance
      eta_new <- unlist(lapply(1:(K + 1), function(k){
        c(X_spline[[k]] %**% beta_spline[[k]] +
            X_flat[[k]] %**% cbind(beta_flat))
      }))
      mu_new <- family$linkinv(eta_new)

      if(!is.null(family$custom_dev.resids)){
        err_new <- mean(family$custom_dev.resids(
          y_vec, mu_new, unlist(order_list),
          family, unlist(observation_weights), ...))
      } else if(!is.null(family$dev.resids)){
        err_new <- mean(family$dev.resids(y_vec, mu_new, wt = 1))
      } else {
        err_new <- mean((y_vec - mu_new)^2)
      }

      if(verbose){
        cat('  \n Backfitting IRWLS iter', IRWLS_iter,
            '| deviance:', format(err_new, digits = 6), '\n')
      }

      ## This block handles updates to the IRWLS iteration error (`err_new`) and
      #   controls damping to prevent divergence or numerical instability.
      # - If `err_new` is NA or not finite, increment `damp_cnt` and break if too many
      #   consecutive failures occur (>= 10 iterations).
      # - If the new error improves on the previous (`err_new <= IRWLS_err`), update the
      #   IRWLS error, reset the damp counter, and check convergence:
      #   Stop if change in error is below `tol` after at least 5 iterations.
      # - If the new error worsens, increment `damp_cnt` and break if repeated >= 10
      if(is.na(err_new) | !is.finite(err_new)){
        damp_cnt <- damp_cnt + 1
        if(damp_cnt >= 10) break
      } else if(err_new <= IRWLS_err){
        prev_IRWLS_err <- IRWLS_err
        IRWLS_err <- err_new
        damp_cnt <- 0
        if(abs(prev_IRWLS_err - IRWLS_err) < tol & IRWLS_iter > 5) break
      } else {
        damp_cnt <- damp_cnt + 1
        if(damp_cnt >= 10) break
      }

      if(!iterate & IRWLS_iter >= 2) break
    }
  }

  ## Assemble full per-partition coefficients from backfitting solution
  #  (skip if is_gauss_id & has_corr since result was set directly above)
  if(!(is_gauss_id & has_corr)){
    result <- lapply(1:(K + 1), function(k){
      b <- rep(0, nc)
      b[spline_cols] <- beta_spline[[k]]
      b[flat_cols]   <- beta_flat
      cbind(b)
    })
  }

  ## GEE GLM SQP refinement
  #  X and y are unwhitened. Form X_tilde = V^{-1/2} X_block
  #  and y_tilde = V^{-1/2} y_block internally for correct info matrix.
  #  Original-scale linear predictor: XB_gee = X_block %*% beta_block.
  if(needs_gee_glm_sqp){

    if(verbose) cat("  GEE GLM SQP refinement\n")

    ## Permute correlation matrices to partition ordering
    perm <- unlist(order_list)
    Vhalf_perm    <- Vhalf[perm, perm]
    VhalfInv_perm <- VhalfInv[perm, perm]

    ## Construct full unwhitened block-diagonal design
    X_block <- Reduce("rbind", lapply(1:(K + 1), function(k){
      Reduce("cbind", lapply(1:(K + 1), function(j){
        if(nrow(X[[k]]) == 0){
          return(X[[k]])
        } else if(j == k) X[[k]] else 0 * X[[k]]
      }))
    }))
    beta_block <- cbind(unlist(result))

    ## Form full whitened design internally
    X_tilde_gee <- VhalfInv_perm %**% X_block
    y_block     <- cbind(unlist(y))
    y_tilde_gee <- VhalfInv_perm %**% y_block

    ## Construct full block-diagonal penalty matrix
    if(unique_penalty_per_partition){
      Lambda_block <- Reduce("rbind", lapply(1:(K + 1), function(k){
        Reduce("cbind", lapply(1:(K + 1), function(j){
          if(j == k) Lambda + L_partition_list[[k]] else 0 * Lambda
        }))
      }))
    } else {
      Lambda_block <- Reduce("rbind", lapply(1:(K + 1), function(k){
        Reduce("cbind", lapply(1:(K + 1), function(j){
          if(j == k) Lambda else 0 * Lambda
        }))
      }))
    }

    ## Combined constraint matrix: smoothness equalities + user inequalities
    if(quadprog && !is.null(qp_Amat)){
      qp_Amat_combined <- cbind(A, qp_Amat)
      qp_bvec_combined <- c(rep(0, ncol(A)), qp_bvec)
      qp_meq_combined  <- ncol(A) + qp_meq
    } else {
      ## No user inequalities: equality constraints only
      qp_Amat_combined <- A
      qp_bvec_combined <- rep(0, ncol(A))
      qp_meq_combined  <- ncol(A)
    }

    ## Incorporate nonzero constraint values into the equality RHS
    if(length(constraint_values) > 0){
      cv_vec <- Reduce("rbind", constraint_values)
      constr_rhs <- crossprod(A, cv_vec)
      if(length(constr_rhs) == ncol(A)){
        qp_bvec_combined[1:ncol(A)] <- c(constr_rhs)
      }
    }

    ## SQP iteration control
    gee_damp_cnt   <- 0
    gee_master_cnt <- 0
    gee_err <- Inf
    gee_converged <- FALSE

    ## Store last successful solve.QP output for qp_info
    last_qp_sol   <- NULL
    last_info_mat <- NULL

    ## Original-scale linear predictor
    #  (X_block is unwhitened, no Vhalf multiplication needed)
    XB_gee <- X_block %**% beta_block

    while(gee_err > tol &
          gee_damp_cnt < 10 &
          gee_master_cnt < 100){

      gee_master_cnt <- gee_master_cnt + 1
      damp <- 2^(-(gee_damp_cnt))

      ## Dispersion at current iterate (original scale)
      if(need_dispersion_for_estimation){
        dispersion_temp <- dispersion_function(
          mu = family$linkinv(XB_gee),
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

      ## GLM working weights on original scale
      W_gee <- c(glm_weight_function(
        family$linkinv(XB_gee),
        y_block,
        unlist(order_list),
        family,
        dispersion_temp,
        unlist(observation_weights),
        ...
      ))

      ## Schur correction: pass X_block and y_block (unwhitened)
      result_temp <- lapply(1:(K + 1), function(k){
        cbind(beta_block[(k - 1) * nc + 1:nc])
      })
      schur_correction <- schur_correction_function(
        list(X_block),
        list(y_block),
        list(cbind(unlist(result_temp))),
        dispersion_temp,
        list(unlist(order_list)),
        0,
        family,
        unlist(observation_weights),
        ...
      )
      if(!(any(unlist(schur_correction) != 0))){
        schur_correction_coll <- 0
      } else {
        schur_correction_coll <- collapse_block_diagonal(schur_correction)
      }

      ## Info matrix using full whitened design X_tilde_gee
      info <- crossprod(X_tilde_gee, W_gee * X_tilde_gee) +
        Lambda_block +
        schur_correction_coll
      sc <- sqrt(mean(abs(info)))

      ## Score: pass X_tilde_gee, y_tilde_gee, V^{-1/2} mu
      qp_score <- qp_score_function(
        X_tilde_gee,
        y_tilde_gee,
        VhalfInv_perm %**% cbind(family$linkinv(c(XB_gee))),
        unlist(order_list),
        dispersion_temp,
        VhalfInv_perm,
        unlist(observation_weights),
        ...
      )

      ## solve.QP step with combined equality and inequality constraints
      qp_sol <- try({quadprog::solve.QP(
        Dmat = info / sc,
        dvec = (qp_score -
                  Lambda_block %**% beta_block +
                  info %**% beta_block) / sc,
        Amat = qp_Amat_combined,
        bvec = qp_bvec_combined,
        meq  = qp_meq_combined
      )}, silent = TRUE)

      ## Fall back to zero-target if solve.QP fails
      if(any(inherits(qp_sol, 'try-error'))){
        if(verbose){
          cat("    GEE QP iter", gee_master_cnt,
              "- solve.QP failed, using current solution\n")
        }
        beta_new <- 0 * beta_block
      } else {
        ## Cache the raw solve.QP result and unscaled info for qp_info
        last_qp_sol   <- qp_sol
        last_info_mat <- info
        beta_new      <- qp_sol$solution
      }

      ## Accept or reject via damped update with deviance monitoring
      if(!iterate & gee_master_cnt > 1){
        ## Non-iterative: accept after first QP step
        beta_block <- beta_new
        gee_converged  <- TRUE
        gee_damp_cnt   <- 11
        gee_master_cnt <- 101
        gee_err <- tol - 1
      } else {
        ## Damped update
        beta_new <- (1 - damp) * beta_block + damp * beta_new
        ## Original-scale LP update
        XB_gee <- X_block %**% cbind(beta_new)

        ## Recompute dispersion and weights at damped point
        if(need_dispersion_for_estimation){
          dispersion_temp <- dispersion_function(
            mu = family$linkinv(XB_gee),
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
        W_gee <- c(glm_weight_function(
          family$linkinv(XB_gee),
          y_block,
          unlist(order_list),
          family,
          dispersion_temp,
          unlist(observation_weights),
          ...
        ))

        ## Evaluate deviance for convergence (original scale)
        if(!is.null(family$custom_dev.resids)){
          raw <- family$custom_dev.resids(
            y_block,
            family$linkinv(c(XB_gee)),
            unlist(order_list),
            family,
            unlist(observation_weights),
            ...
          )
          ## Weight residuals by 1/sqrt(W) before whitening.
          W_gee_safe <- pmax(W_gee, sqrt(.Machine$double.eps))
          err_new <- mean((
            VhalfInv_perm %**%
              cbind(sign(raw) * sqrt(abs(raw)) / sqrt(c(W_gee_safe)))
          )^2)
        } else if(is.null(family$dev.resids)){
          ## Compare against y_block (unwhitened)
          err_new <- mean(
            (unlist(observation_weights) *
               (y_block - cbind(family$linkinv(XB_gee))))^2
          )
        } else {
          ## Compare against y_block (unwhitened)
          err_new <- mean(
            unlist(observation_weights) *
              family$dev.resids(
                y_block,
                cbind(family$linkinv(XB_gee)),
                wt = 1
              )
          )
        }

        ## Step acceptance logic
        # - Track whether the new estimate `beta_new` improves the error (`err_new`).
        # - Increment damping counter if the step fails or produces invalid error.
        # - Accept the step if `err_new` decreases; reset damping counter and update `beta_block`.
        # - Trigger convergence if both coefficient change and error reduction are below `tol`
        #   after sufficient iterations.
        if(is.null(err_new) | is.na(err_new) | !is.finite(err_new)){
          gee_damp_cnt <- gee_damp_cnt + 1
        } else if(err_new <= gee_err){

          prev_gee_err <- gee_err
          gee_err <- err_new
          abs_diff <- max(abs(beta_new - beta_block))
          beta_block <- beta_new
          gee_damp_cnt <- 0

          if((abs_diff < tol) &
             (prev_gee_err - gee_err < tol) &
             (gee_master_cnt > 10)){
            gee_converged  <- TRUE
            gee_damp_cnt   <- 11
            gee_master_cnt <- 101
            gee_err <- tol - 1
          }

        } else {
          gee_damp_cnt <- gee_damp_cnt + 1
        }
      }

      if(verbose & (gee_master_cnt <= 100)){
        cat("    GEE QP iter", gee_master_cnt,
            "| deviance:", format(gee_err, digits = 6),
            "| damp:", format(damp, digits = 4), "\n")
      }

      ## Unpack into per-partition form
      result <- lapply(1:(K + 1), function(k){
        cbind(beta_block[1:nc + (k - 1) * nc])
      })
    }

    ## Assemble qp_info from last successful solve.QP call (GEE branch)
    if(!is.null(last_qp_sol)){
      lag_mult <- last_qp_sol$Lagrangian   # multipliers, one per constraint col
      iact_raw <- last_qp_sol$iact         # 1-based active constraint indices
      iact_raw <- iact_raw[iact_raw > 0]   # solve.QP pads with zeros
      # active_constraints: equality indices (always active) plus active
      # inequality indices shifted to the combined index space
      n_eq <- qp_meq_combined
      n_con <- ncol(qp_Amat_combined)
      eq_idx  <- seq_len(n_eq)
      ineq_active <- if(length(iact_raw) > 0) n_eq + iact_raw else integer(0)
      active_constraints <- sort(unique(c(eq_idx, ineq_active)))
      qp_info <- list(
        lagrangian  = lag_mult,
        active_constraints = active_constraints,
        iact   = iact_raw,
        Amat_active = qp_Amat_combined[, active_constraints, drop = FALSE],
        bvec_active = qp_bvec_combined[active_constraints],
        meq_active = qp_meq_combined[active_constraints],
        converged  = gee_converged
      )
    }

    ## Clean up big matrices
    Vhalf_perm <- NULL
    VhalfInv_perm <- NULL
    X_tilde_gee <- NULL
    y_tilde_gee <- NULL
  }

  ## Standard SQP refinement (non-GEE cases)
  #  For models without a correlation structure, SQP is used only when
  #  quadprog = TRUE to enforce user-supplied inequality constraints.
  #  For Gaussian identity + GEE + quadprog, the closed-form solve
  #  already handled inequality constraints above via the optional QP
  #  block in case (a-gee), so this block is skipped for that case.
  if(quadprog & !needs_gee_glm_sqp & !(is_gauss_id & has_corr)){

    if(verbose) cat("  QP refinement (SQP loop)\n")

    ## Construct full block-diagonal design matrix
    X_block <- Reduce("rbind", lapply(1:(K + 1), function(k){
      Reduce("cbind", lapply(1:(K + 1), function(j){
        if(nrow(X[[k]]) == 0){
          return(X[[k]])
        } else if(j == k) X[[k]] else 0 * X[[k]]
      }))
    }))
    beta_block <- cbind(unlist(result))

    ## Construct full block-diagonal penalty matrix
    if(unique_penalty_per_partition){
      Lambda_block <- Reduce("rbind", lapply(1:(K + 1), function(k){
        Reduce("cbind", lapply(1:(K + 1), function(j){
          if(j == k) Lambda + L_partition_list[[k]] else 0 * Lambda
        }))
      }))
    } else {
      Lambda_block <- Reduce("rbind", lapply(1:(K + 1), function(k){
        Reduce("cbind", lapply(1:(K + 1), function(j){
          if(j == k) Lambda else 0 * Lambda
        }))
      }))
    }

    y_block <- cbind(unlist(y))

    ## Damped SQP iteration control
    qp_damp_cnt   <- 0
    qp_master_cnt <- 0
    qp_err <- Inf
    qp_converged <- FALSE

    ## Store last successful solve.QP output for qp_info
    last_qp_sol   <- NULL
    last_info_mat <- NULL

    ## Combined constraint matrix (equalities + inequalities)
    qp_Amat_combined_std <- cbind(A, qp_Amat)
    qp_bvec_combined_std <- c(rep(0, ncol(A)), qp_bvec)
    qp_meq_combined_std  <- ncol(A) + qp_meq

    XB <- X_block %**% beta_block

    ## Damped SQP loop
    while(qp_err > tol & qp_damp_cnt < 10 & qp_master_cnt < 100){
      qp_master_cnt <- qp_master_cnt + 1
      damp <- 2^(-(qp_damp_cnt))

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

      ## GLM working weights
      W_qp <- c(glm_weight_function(family$linkinv(XB),
                                    y_block,
                                    unlist(order_list),
                                    family,
                                    dispersion_temp,
                                    unlist(observation_weights),
                                    ...))

      ## Schur correction in partition form
      result_temp <- lapply(1:(K + 1), function(k){
        cbind(beta_block[(k - 1) * nc + 1:nc])
      })
      schur_correction <-
        schur_correction_function(
          X,
          y,
          result_temp,
          dispersion_temp,
          order_list,
          K,
          family,
          observation_weights,
          ...
        )
      if(!(any(unlist(schur_correction) != 0))){
        schur_correction_coll <- 0
      } else {
        schur_correction_coll <- collapse_block_diagonal(schur_correction)
      }

      ## Penalized information matrix
      info <- crossprod(X_block, W_qp * X_block) +
        Lambda_block +
        schur_correction_coll

      sc <- sqrt(mean(abs(info)))

      ## Score from user-supplied QP score function
      qp_score <- qp_score_function(
        X_block,
        y_block,
        cbind(family$linkinv(XB)),
        unlist(order_list),
        dispersion_temp,
        NULL,
        unlist(observation_weights),
        ...
      )

      ## solve.QP step with both equality and inequality constraints
      qp_sol <- try({quadprog::solve.QP(
        Dmat = info / sc,
        dvec = (qp_score -
                  Lambda_block %**% beta_block +
                  info %**% beta_block) / sc,
        Amat = qp_Amat_combined_std,
        bvec = qp_bvec_combined_std,
        meq  = qp_meq_combined_std
      )}, silent = TRUE)

      ## Fall back to current solution if solve.QP fails
      if(any(inherits(qp_sol, 'try-error'))){
        if(verbose){
          cat("    QP iter", qp_master_cnt, "- solve.QP failed,",
              "using current solution\n")
        }
        beta_new <- 0 * beta_block
      } else {
        ## Cache raw solve.QP result and unscaled info for qp_info
        last_qp_sol   <- qp_sol
        last_info_mat <- info
        beta_new      <- qp_sol$solution
      }

      ## Accept or reject via damped update with deviance monitoring
      if(!iterate & qp_master_cnt > 1){
        ## Non-iterative: accept undamped after first QP step
        beta_block <- beta_new
        qp_converged  <- TRUE
        qp_damp_cnt   <- 11
        qp_master_cnt <- 101
        qp_err <- tol - 1
      } else {
        ## Damped update
        beta_new <- (1 - damp) * beta_block + damp * beta_new
        XB <- X_block %**% beta_new

        ## Evaluate deviance for convergence
        if(!is.null(family$custom_dev.resids) &
           is.null(family$dev.resids)){
          err_new <- mean(
            family$custom_dev.resids(y_block,
                                     family$linkinv(c(XB)),
                                     unlist(order_list),
                                     family,
                                     unlist(observation_weights),
                                     ...))
        } else if(is.null(family$dev.resids)){
          err_new <- mean(unlist(observation_weights) *
                            (y_block - family$linkinv(XB))^2)
        } else {
          err_new <-
            mean(unlist(observation_weights) *
                   family$dev.resids(y_block,
                                     family$linkinv(XB),
                                     wt = 1))
        }

        ## Step acceptance logic
        if(is.null(err_new) | is.na(err_new) | !is.finite(err_new)){
          qp_damp_cnt <- qp_damp_cnt + 1
        } else if(err_new <= qp_err){

          prev_qp_err <- qp_err
          qp_err <- err_new
          abs_diff <- max(abs(beta_new - beta_block))
          beta_block <- beta_new
          qp_damp_cnt <- 0

          if((abs_diff < tol) &
             (prev_qp_err - qp_err < tol) &
             (qp_master_cnt > 10)){
            qp_converged  <- TRUE
            qp_damp_cnt   <- 11
            qp_master_cnt <- 101
            qp_err <- tol - 1
          }

        } else {
          qp_damp_cnt <- qp_damp_cnt + 1
        }
      }

      if(verbose & (qp_master_cnt <= 100)){
        cat("    QP iter", qp_master_cnt,
            "| deviance:", format(qp_err, digits = 6),
            "| damp:", format(damp, digits = 4), "\n")
      }

      ## Unpack into per-partition form
      result <- lapply(1:(K + 1), function(k){
        cbind(beta_block[1:nc + (k - 1) * nc])
      })
    }

    ## Assemble qp_info from last successful solve.QP call (standard branch)
    if(!is.null(last_qp_sol)){
      lag_mult <- last_qp_sol$Lagrangian
      iact_raw <- last_qp_sol$iact
      iact_raw <- iact_raw[iact_raw > 0]
      n_eq  <- qp_meq_combined_std
      eq_idx  <- seq_len(n_eq)
      ineq_active <- if(length(iact_raw) > 0) n_eq + iact_raw else integer(0)
      active_constraints <- sort(unique(c(eq_idx, ineq_active)))
      qp_info <- list(
        solution           = c(beta_block),
        lagrangian         = lag_mult,
        active_constraints = active_constraints,
        iact               = iact_raw,
        info_matrix        = last_info_mat,
        Amat_combined      = qp_Amat_combined_std,
        bvec_combined      = qp_bvec_combined_std,
        meq_combined       = qp_meq_combined_std,
        converged          = qp_converged,
        final_deviance     = qp_err
      )
    }
  }

  ## Return G_list for downstream inference
  #  G is always returned as tcrossprod(Ghalf) to
  #  guarantee G = Ghalf %*% t(Ghalf) exactly for varcov computation.

  if(!return_G_getB){
    return(list(B = result, G_list = NULL, qp_info = qp_info))
  }

  if(is_gauss_id & !quadprog & !needs_gee_glm_sqp){
    ## G does not depend on coefficients for Gaussian identity;
    #  reuse full-dimensional matrices from lgspline.fit.
    #  G recomputed as tcrossprod(Ghalf) exactly.
    return(list(
      B = result,
      G_list = list(
        G        = lapply(Ghalf_full, function(mat) tcrossprod(mat)),
        Ghalf    = Ghalf_full,
        GhalfInv = GhalfInv_full
      ),
      qp_info = qp_info
    ))
  }

  ## GLM, post-QP, or GEE GLM: recompute full-dimensional G at
  #  converged estimates.
  if(!exists("dispersion_temp")){
    if(!exists("disp_temp")){
      dispersion_temp <- 1
    } else {
      dispersion_temp <- disp_temp
    }
  }

  ## For GEE GLM, recompute G using the per-partition approximation
  #  (diagonal blocks of V^{-1/2} X). This is the same pre-existing
  #  approximation used throughout the pipeline: cross-partition
  #  contributions from off-diagonal blocks of V^{-1/2} are ignored.
  #  This is for compatibility - the ultimate posterior variance
  #  covariance matrix sigma^2 UG will use the full dense matrix anyways.
  Xw <- lapply(1:(K + 1), function(k){
    if(nrow(X[[k]]) == 0) return(X[[k]])
    var_k <- glm_weight_function(
      family$linkinv(X[[k]] %**% cbind(c(result[[k]]))),
      y[[k]], order_list[[k]], family, dispersion_temp,
      observation_weights[[k]], ...)
    X[[k]] * c(sqrt(pmax(c(var_k), 0)))
  })
  X_gram_w <- lapply(1:(K + 1), function(k) crossprod(Xw[[k]]))
  schur_corr <- schur_correction_function(
    X, y, result, dispersion_temp,
    order_list, K, family, observation_weights, ...)
  G_list_out <- compute_G_eigen(X_gram_w,
                                Lambda, K,
                                parallel_eigen, cl,
                                chunk_size, num_chunks, rem_chunks,
                                family,
                                unique_penalty_per_partition,
                                L_partition_list,
                                keep_G = TRUE,
                                schur_corr)
  ## Guarantee G = tcrossprod(Ghalf) exactly
  #  so that downstream varcov (sigma^2 * UG) is computed consistently
  #  regardless of what compute_G_eigen returns for G.
  G_list_out$G <- lapply(G_list_out$Ghalf, function(mat) tcrossprod(mat))
  return(list(B = result, G_list = G_list_out, qp_info = qp_info))
}
