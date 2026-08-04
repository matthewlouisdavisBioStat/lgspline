## additive_lgspline.R

.additive_normalize_groups <- function(spline_groups, q_predictors,
                                       pred_colnames = NULL){
  if(is.null(spline_groups) || length(spline_groups) == 0L) return(NULL)
  if(!is.list(spline_groups)) spline_groups <- list(spline_groups)

  out <- lapply(spline_groups, function(grp){
    if(any(is.character(grp))){
      if(is.null(pred_colnames)){
        stop('\n \t Character entries in spline_groups require named predictors.\n')
      }
      grp <- unlist(lapply(grp, function(v){
        if(is.character(v)){
          idx <- which(pred_colnames == v)
          if(length(idx) == 0L) idx <- grep(v, pred_colnames, fixed = TRUE)
          idx
        } else {
          v
        }
      }))
    }
    unique(as.integer(grp[!is.na(grp)]))
  })

  out <- out[sapply(out, length) > 0L]
  if(length(out) == 0L) return(NULL)
  if(any(unlist(out) < 1L) || any(unlist(out) > q_predictors)){
    stop('\n \t spline_groups contains indices outside predictors columns.\n')
  }
  out
}


.additive_use_groups <- function(spline_groups, q_predictors,
                                 pred_colnames = NULL){
  groups <- .additive_normalize_groups(spline_groups, q_predictors,
                                       pred_colnames)
  !is.null(groups) && length(groups) > 1L
}


.additive_arg_for_term <- function(x, term_index, n_terms){
  if(is.null(x)) return(NULL)
  if(is.list(x) && !is.data.frame(x)){
    if(!is.null(x$additive_terms)){
      return(x$additive_terms[[term_index]])
    }
    if(length(x) == n_terms){
      return(x[[term_index]])
    }
  }
  if(!is.list(x) && length(x) == n_terms){
    return(x[[term_index]])
  }
  x
}


.additive_B_for_term <- function(B_predict, term_index, n_terms){
  if(is.null(B_predict)) return(NULL)
  if(is.list(B_predict) && !is.data.frame(B_predict)){
    if(!is.null(B_predict$additive_terms)){
      return(B_predict$additive_terms[[term_index]])
    }
    if(length(B_predict) == n_terms){
      return(B_predict[[term_index]])
    }
  }
  B_predict
}


.additive_is_per_term_arg <- function(x, n_terms){
  if(is.null(x)) return(TRUE)
  is.list(x) && !is.data.frame(x) &&
    (!is.null(x$additive_terms) || length(x) == n_terms)
}


.additive_empty_constraint <- function(x){
  is.null(x) || length(x) == 0L
}


.additive_term_cols <- function(group, anchor_cols, term_index){
  unique(c(group, anchor_cols))
}


.additive_raw_expansion_cols <- function(raw_names){
  lapply(raw_names, function(nm){
    toks <- regmatches(nm, gregexpr("_[0-9]+_", nm))[[1]]
    if(length(toks) == 0L) return(integer(0))
    unique(as.integer(gsub("_", "", toks)))
  })
}


.additive_anchor_only_expansion_cols <- function(raw_names, anchor_cols){
  if(length(anchor_cols) == 0L) return(integer(0))
  cols <- .additive_raw_expansion_cols(raw_names)
  which(vapply(cols, function(x){
    length(x) > 0L && all(x %in% anchor_cols)
  }, logical(1)))
}


.additive_local_factor_groups <- function(factor_groups, local_cols){
  if(is.null(factor_groups) || length(factor_groups) == 0L) return(NULL)

  out <- lapply(factor_groups, function(cols){
    match(cols, local_cols)
  })
  out <- lapply(out, function(cols) cols[!is.na(cols)])
  out <- out[sapply(out, length) > 1L]
  if(length(out) == 0L) return(NULL)
  out
}


.additive_local_exclude_expansions <- function(exclude_these_expansions,
                                              local_cols){
  if(is.null(exclude_these_expansions) ||
     length(exclude_these_expansions) == 0L){
    return(NULL)
  }

  out <- character(0)
  for(pattern in exclude_these_expansions){
    matches <- gregexpr("_[0-9]+_", pattern)
    toks <- regmatches(pattern, matches)[[1]]
    if(length(toks) == 0L){
      out <- c(out, pattern)
      next
    }
    global_cols <- as.integer(gsub("_", "", toks))
    if(!all(global_cols %in% local_cols)) next
    local_pattern <- pattern
    for(ii in seq_along(local_cols)){
      local_pattern <- gsub(paste0("_", local_cols[ii], "_"),
                            paste0("_", ii, "_"),
                            local_pattern,
                            fixed = TRUE)
    }
    out <- c(out, local_pattern)
  }
  if(length(out) == 0L) return(NULL)
  unique(out)
}


.additive_block_diag_rect <- function(matrix_list){
  matrix_list <- matrix_list[!sapply(matrix_list, is.null)]
  if(length(matrix_list) == 0L) return(matrix(0, 0, 0))

  nrs <- sapply(matrix_list, nrow)
  ncs <- sapply(matrix_list, ncol)
  out <- matrix(0, sum(nrs), sum(ncs))
  r0 <- 0L
  c0 <- 0L
  for(i in seq_along(matrix_list)){
    if(nrs[i] > 0L && ncs[i] > 0L){
      out[r0 + seq_len(nrs[i]), c0 + seq_len(ncs[i])] <- matrix_list[[i]]
    }
    r0 <- r0 + nrs[i]
    c0 <- c0 + ncs[i]
  }
  out
}


.additive_safe_linkfun <- function(mu, family){
  mu <- .default_glm_safe_mu(mu, family)
  out <- tryCatch(c(family$linkfun(mu)),
                  error = function(e) rep(NA_real_, length(mu)))
  bad <- !is.finite(out)
  if(any(bad)){
    out[bad] <- 0
  }
  out
}


.additive_term_eta <- function(term_fit, term_predictors, family,
                               B_predict = NULL){
  if(is.null(B_predict)){
    mu <- c(term_fit$predict(new_predictors = term_predictors))
  } else {
    mu <- c(term_fit$predict(new_predictors = term_predictors,
                             B_predict = B_predict))
  }
  .additive_safe_linkfun(mu, family)
}


.additive_rename_fit <- function(term_fit, og_cols){
  if(is.null(og_cols) || is.null(term_fit$B)) return(term_fit)
  term_fit$og_cols <- og_cols
  if(is.null(term_fit$.fit_call_args)) term_fit$.fit_call_args <- list()
  term_fit$.fit_call_args$og_cols <- og_cols

  old <- paste0('_', seq_along(og_cols), '_')
  new <- og_cols
  rename_one <- function(nm){
    for(ii in seq_along(old)){
      nm <- gsub(old[ii], new[ii], nm, fixed = TRUE)
    }
    nm
  }

  for(k in seq_along(term_fit$B)){
    nms <- sapply(names(term_fit$B[[k]]), rename_one)
    rownames(term_fit$B[[k]]) <- nms
    names(term_fit$B[[k]]) <- nms
  }
  if(!is.null(term_fit$B_raw)){
    for(k in seq_along(term_fit$B_raw)){
      nms <- sapply(names(term_fit$B_raw[[k]]), rename_one)
      rownames(term_fit$B_raw[[k]]) <- nms
      names(term_fit$B_raw[[k]]) <- nms
    }
  }
  if(!is.null(term_fit$A) && length(term_fit$A) > 0L &&
     !is.null(term_fit$p) && !is.null(term_fit$K) &&
     nrow(cbind(term_fit$A)) == term_fit$p * (term_fit$K + 1L)){
    rn <- sapply(names(term_fit$B[[1]]), rename_one)
    rownames(term_fit$A) <- paste0(
      rep(paste0('partition', seq_len(term_fit$K + 1L)),
          each = term_fit$p), '_', rn)
  }
  term_fit
}


.additive_make_term_predictors <- function(predictors, group, anchor_cols,
                                           offset_eta, term_index,
                                           pred_colnames){
  cols <- .additive_term_cols(group, anchor_cols, term_index)
  term_predictors <- cbind(predictors[, cols, drop = FALSE],
                           .additive_offset = c(offset_eta))
  if(!is.null(pred_colnames)){
    colnames(term_predictors) <- c(pred_colnames[cols], ".additive_offset")
  }
  term_predictors
}


.additive_prepare_new_predictors <- function(new_predictors,
                                             pred_colnames = NULL,
                                             include_warnings = FALSE){
  if(is.null(new_predictors)) return(NULL)

  if(!is.null(pred_colnames) &&
     (inherits(new_predictors, "data.frame") ||
      (is.matrix(new_predictors) && !is.null(colnames(new_predictors))))){
    np_df <- as.data.frame(new_predictors)
    np_colnames <- colnames(np_df)
    relevant_cols <- np_colnames[sapply(np_colnames, function(cn){
      if(cn %in% pred_colnames) return(TRUE)
      onehot_pattern <- paste0("^", cn, "_")
      any(grepl(onehot_pattern, pred_colnames))
    })]
    if(length(relevant_cols) > 0L &&
       length(relevant_cols) < length(np_colnames)){
      np_df <- np_df[, relevant_cols, drop = FALSE]
    }

    col_names <- colnames(np_df)
    cols_needing_encoding <- character(0)
    for(cn in col_names){
      if(cn %in% pred_colnames) next
      onehot_pattern <- paste0("^", cn, "_")
      if(any(grepl(onehot_pattern, pred_colnames))){
        cols_needing_encoding <- c(cols_needing_encoding, cn)
      }
    }

    if(length(cols_needing_encoding) > 0L){
      for(col_name in cols_needing_encoding){
        onehot_pattern <- paste0("^", col_name, "_")
        candidate_cols <- pred_colnames[grepl(onehot_pattern, pred_colnames)]
        prefix <- paste0(col_name, "_")
        levels_fit <- substr(candidate_cols,
                             nchar(prefix) + 1L,
                             nchar(candidate_cols))
        user_vals <- as.character(np_df[[col_name]])
        onehot_mat <- matrix(0L,
                             nrow = nrow(np_df),
                             ncol = length(levels_fit))
        colnames(onehot_mat) <- candidate_cols
        for(lv_idx in seq_along(levels_fit)){
          onehot_mat[user_vals == levels_fit[lv_idx], lv_idx] <- 1L
        }
        unrecognized <- unique(user_vals[!user_vals %in% levels_fit])
        if(length(unrecognized) > 0L && include_warnings){
          warning("\n\t predict: value(s) [",
                  paste(unrecognized, collapse = ", "),
                  "] in column '", col_name,
                  "' were not seen during fitting and will be treated as ",
                  "all-zero (reference level).\n")
        }

        col_idx <- which(colnames(np_df) == col_name)
        np_df <- cbind(np_df[, -col_idx, drop = FALSE],
                       as.data.frame(onehot_mat))
      }
    }

    if(any(colnames(np_df) %in% pred_colnames)){
      missing_cols <- setdiff(pred_colnames, colnames(np_df))
      if(length(missing_cols) > 0L){
        np_df[missing_cols] <- 0L
      }
      if(all(pred_colnames %in% colnames(np_df))){
        np_df <- np_df[, pred_colnames, drop = FALSE]
      }
      return(methods::as(cbind(np_df), 'matrix'))
    }
  }

  out <- try(methods::as(cbind(new_predictors), 'matrix'), silent = TRUE)
  if(inherits(out, 'try-error')){
    stop('\n\t new_predictors must be coercible to a matrix.\n')
  }
  if(!is.null(pred_colnames) && ncol(out) == length(pred_colnames)){
    colnames(out) <- pred_colnames
  }
  if(!is.null(pred_colnames) && ncol(out) != length(pred_colnames)){
    stop('\n\t new_predictors has incompatible columns after factor ',
         'encoding. Expected ', length(pred_colnames), ' columns.\n')
  }
  out
}


.additive_zero_offset_predictors <- function(term_predictors){
  out <- term_predictors
  out[, ncol(out)] <- 0
  out
}


.additive_dev <- function(y, mu, family, observation_weights,
                          order_indices, ...){
  wt <- .default_glm_observation_weights(observation_weights, length(y))
  if(!is.null(family$custom_dev.resids)){
    mean(family$custom_dev.resids(y, mu, order_indices, family, wt, ...))
  } else if(!is.null(family$dev.resids)){
    mean(wt * family$dev.resids(y, mu, wt = 1))
  } else {
    mean(wt * (y - mu)^2)
  }
}


.additive_split_original_beta <- function(beta, template_B){
  lens <- sapply(template_B, length)
  ends <- cumsum(lens)
  starts <- if(length(ends) == 1L) 1L else c(1L, ends[-length(ends)] + 1L)
  out <- lapply(seq_along(template_B), function(k){
    vals <- beta[starts[k]:ends[k]]
    names(vals) <- names(template_B[[k]])
    vals
  })
  names(out) <- names(template_B)
  out
}


.additive_draw_term_coefficients <- function(term_fit, sigmasq,
                                             enforce_qp_constraints = TRUE,
                                             max_slice_iterations = 1000L,
                                             include_warnings = TRUE){
  p_expansions <- term_fit$p
  K <- term_fit$K
  P_total <- p_expansions * (K + 1L)

  use_raw_path <- !is.null(term_fit$Ghalf) && !is.null(term_fit$B_raw)
  if(use_raw_path){
    use_U <- if(!is.null(term_fit$U)) term_fit$U else diag(P_total)
    Ghalf_bd <- matrix(0, P_total, P_total)
    for(k in seq_len(K + 1L)){
      rows <- ((k - 1L) * p_expansions + 1L):(k * p_expansions)
      Ghalf_bd[rows, rows] <- term_fit$Ghalf[[k]]
    }
    L_post <- (1 / term_fit$sd_y) * (use_U %**% Ghalf_bd)
    beta_mode <- unlist(term_fit$B_raw)

    has_ineq <- FALSE
    if(enforce_qp_constraints &&
       !is.null(term_fit$quadprog_list) &&
       !identical(term_fit$quadprog_list, list(NA)) &&
       !is.null(term_fit$quadprog_list$qp_Amat)){
      qp_Amat_full <- term_fit$quadprog_list$qp_Amat
      qp_bvec_full <- term_fit$quadprog_list$qp_bvec
      qp_meq_val <- term_fit$quadprog_list$qp_meq
      n_total <- ncol(qp_Amat_full)
      if(qp_meq_val < n_total){
        ineq_cols <- (qp_meq_val + 1L):n_total
        ineq_Amat <- qp_Amat_full[, ineq_cols, drop = FALSE]
        ineq_bvec <- qp_bvec_full[ineq_cols]
        has_ineq <- TRUE
      }
    }

    .check_feasible <- function(b){
      if(!has_ineq) return(TRUE)
      all(c(crossprod(ineq_Amat, cbind(b))) >=
            ineq_bvec - sqrt(.Machine$double.eps))
    }

    if(has_ineq){
      beta_cur <- beta_mode
      nu <- sqrt(sigmasq) * c(L_post %**% cbind(rnorm(ncol(L_post))))
      theta <- runif(1, 0, 2 * pi)
      theta_min <- theta - 2 * pi
      theta_max <- theta
      for(iter in seq_len(max_slice_iterations)){
        beta_prop <- beta_cur * cos(theta) + nu * sin(theta)
        if(.check_feasible(beta_prop)){
          beta_cur <- beta_prop
          break
        }
        if(theta < 0) theta_min <- theta else theta_max <- theta
        theta <- runif(1, theta_min, theta_max)
      }
      beta_raw <- beta_cur
    } else {
      z <- rnorm(ncol(L_post))
      beta_raw <- beta_mode + sqrt(sigmasq) * c(L_post %**% cbind(z))
    }

    coefs <- lapply(seq_len(K + 1L), function(k){
      inds <- ((k - 1L) * p_expansions + 1L):(k * p_expansions)
      bk <- beta_raw[inds] * term_fit$sd_y
      bk[1L] <- bk[1L] + term_fit$mean_y
      term_fit$backtransform_coefficients(bk)
    })
    names(coefs) <- names(term_fit$B)
    return(coefs)
  }

  if(!is.null(term_fit$varcovmat) && !is.null(term_fit$B)){
    beta_mode <- unlist(term_fit$B)
    vcov <- term_fit$varcovmat
    if(!is.null(term_fit$sigmasq_tilde) &&
       is.finite(term_fit$sigmasq_tilde) &&
       term_fit$sigmasq_tilde > 0){
      vcov <- vcov * sigmasq / term_fit$sigmasq_tilde
    }
    eig <- eigen(vcov, symmetric = TRUE)
    eig$values <- pmax(eig$values, 0)
    beta <- beta_mode +
      c(eig$vectors %**% cbind(sqrt(eig$values) * rnorm(length(beta_mode))))
    return(.additive_split_original_beta(beta, term_fit$B))
  }

  if(include_warnings){
    warning('\n\t Posterior coefficient covariance is unavailable for an ',
            'additive term; returning fitted coefficients for that term.\n')
  }
  term_fit$B
}


.additive_term_coef_names <- function(term_fit){
  nms <- names(unlist(term_fit$B_raw))
  if(is.null(nms) || length(nms) == 0L ||
     all(is.na(nms) | nms == "")){
    nms <- names(unlist(term_fit$B))
  }
  nms
}


.additive_offset_coef_cols <- function(term_fit){
  raw_names <- term_fit$raw_expansion_names
  if(!is.null(raw_names)){
    hits <- which(grepl(".additive_offset", raw_names, fixed = TRUE))
    if(length(hits) > 0L) return(hits)
  }

  if(!is.null(term_fit$B) && length(term_fit$B) > 0L){
    nms <- names(term_fit$B[[1L]])
    if(is.null(nms)){
      nms <- rownames(cbind(term_fit$B[[1L]]))
    }
    if(!is.null(nms)){
      hits <- which(grepl(".additive_offset", nms, fixed = TRUE))
      if(length(hits) > 0L) return(hits)
    }
  }

  nms <- .additive_term_coef_names(term_fit)
  if(is.null(nms) || is.null(term_fit$p)) return(integer(0))
  hits <- which(grepl(".additive_offset", nms, fixed = TRUE))
  if(length(hits) == 0L) return(integer(0))
  unique(((hits - 1L) %% term_fit$p) + 1L)
}


.additive_offset_coef_rows <- function(term_fit){
  offset_cols <- .additive_offset_coef_cols(term_fit)
  if(length(offset_cols) == 0L) return(integer(0))
  unlist(lapply(seq_len(term_fit$K + 1L), function(k){
    (k - 1L) * term_fit$p + offset_cols
  }), use.names = FALSE)
}


.additive_drop_offset_constraint_cols <- function(A, term_fit){
  if(is.null(A) || ncol(cbind(A)) == 0L) return(A)
  A <- cbind(A)
  offset_rows <- .additive_offset_coef_rows(term_fit)
  if(length(offset_rows) == 0L) return(A)
  keep <- colSums(abs(A[offset_rows, , drop = FALSE])) <=
    sqrt(.Machine$double.eps)
  A[, keep, drop = FALSE]
}


.additive_term_design_original <- function(term_fit){
  X_std <- collapse_block_diagonal(lapply(term_fit$X, term_fit$std_X))
  X_og <- X_std[term_fit$og_order, , drop = FALSE]
  offset_cols <- .additive_offset_coef_rows(term_fit)
  if(length(offset_cols) > 0L){
    X_og[, offset_cols] <- 0
  }
  X_og
}


.additive_new_design_original <- function(term_fit, term_predictors){
  ex <- term_fit$predict(new_predictors = term_predictors,
                         expansions_only = TRUE)
  N_new <- length(ex$partition_codes)
  P_term <- term_fit$p * (term_fit$K + 1L)
  X_full <- matrix(0, nrow = N_new, ncol = P_term)
  order_list_new <- term_fit$knot_expand_function(
    ex$partition_codes,
    ex$partition_bounds,
    N_new,
    cbind(seq_len(N_new)),
    term_fit$K)
  for(k in seq_len(term_fit$K + 1L)){
    if(nrow(ex$expansions[[k]]) == 0L) next
    rows <- c(order_list_new[[k]])
    cols <- ((k - 1L) * term_fit$p + 1L):(k * term_fit$p)
    X_full[rows, cols] <- ex$expansions[[k]]
  }
  offset_cols <- .additive_offset_coef_rows(term_fit)
  if(length(offset_cols) > 0L){
    X_full[, offset_cols] <- 0
  }
  X_full
}


.additive_term_flat_cols <- function(term_fit){
  flat_cols <- unique(c(term_fit$blockfit_flat_cols,
                        term_fit$additive_force_flat_cols))
  if(is.null(flat_cols) || length(flat_cols) == 0L ||
     is.null(term_fit$K) || term_fit$K == 0L){
    return(integer(0))
  }

  raw_names <- term_fit$raw_expansion_names
  if(!is.null(raw_names)){
    offset_cols <- grep(".additive_offset", raw_names, fixed = TRUE)
    if(length(offset_cols) == 0L){
      offset_cols <- .additive_offset_coef_cols(term_fit)
    }
    flat_cols <- setdiff(flat_cols, offset_cols)
  } else {
    flat_cols <- setdiff(flat_cols, .additive_offset_coef_cols(term_fit))
  }
  flat_cols <- as.integer(flat_cols[flat_cols >= 1L & flat_cols <= term_fit$p])
  unique(flat_cols)
}


.additive_term_lambda_full <- function(term_fit){
  has_part_pen <- length(term_fit$penalties$L_partition_list) ==
    (term_fit$K + 1L)
  flat_cols <- .additive_term_flat_cols(term_fit)
  offset_cols <- .additive_offset_coef_cols(term_fit)
  collapse_block_diagonal(
    lapply(seq_len(term_fit$K + 1L), function(k){
      L <- if(has_part_pen){
        term_fit$penalties$Lambda +
          term_fit$penalties$L_partition_list[[k]]
      } else {
        term_fit$penalties$Lambda
      }
      if(length(flat_cols) > 0L){
        L[flat_cols, flat_cols] <- L[flat_cols, flat_cols, drop = FALSE] /
          (term_fit$K + 1L)
      }
      if(length(offset_cols) > 0L){
        L[offset_cols, ] <- 0
        L[, offset_cols] <- 0
      }
      L
    })
  )
}


.additive_term_active_A <- function(term_fit){
  P_term <- term_fit$p * (term_fit$K + 1L)
  if(!is.null(term_fit$qp_info) &&
     !is.null(term_fit$qp_info$Amat_active) &&
     ncol(term_fit$qp_info$Amat_active) > 0L){
    A <- cbind(term_fit$qp_info$Amat_active)
  } else if(!is.null(term_fit$A) && ncol(cbind(term_fit$A)) > 0L){
    A <- cbind(term_fit$A)
  } else {
    A <- matrix(0, P_term, 0)
  }
  A <- .additive_drop_offset_constraint_cols(A, term_fit)
  if(nrow(A) != P_term){
    stop('\n\t Additive active constraint matrix has incompatible rows.\n')
  }
  A
}


.additive_term_flat_A <- function(term_fit){
  flat_cols <- .additive_term_flat_cols(term_fit)
  if(length(flat_cols) == 0L){
    return(matrix(0, term_fit$p * (term_fit$K + 1L), 0))
  }

  out <- matrix(0,
                nrow = term_fit$p * (term_fit$K + 1L),
                ncol = term_fit$K * length(flat_cols))
  cc <- 0L
  for(col in flat_cols){
    for(k in 2:(term_fit$K + 1L)){
      cc <- cc + 1L
      out[col, cc] <- -1
      out[(k - 1L) * term_fit$p + col, cc] <- 1
    }
  }
  out
}


.additive_term_zero_main_A <- function(term_fit){
  zero_cols <- term_fit$additive_zero_main_cols
  if(is.null(zero_cols) || length(zero_cols) == 0L){
    return(matrix(0, term_fit$p * (term_fit$K + 1L), 0))
  }
  zero_cols <- as.integer(zero_cols[zero_cols >= 1L &
                                      zero_cols <= term_fit$p])
  if(length(zero_cols) == 0L){
    return(matrix(0, term_fit$p * (term_fit$K + 1L), 0))
  }

  out <- matrix(0,
                nrow = term_fit$p * (term_fit$K + 1L),
                ncol = (term_fit$K + 1L) * length(zero_cols))
  cc <- 0L
  for(col in zero_cols){
    for(k in seq_len(term_fit$K + 1L)){
      cc <- cc + 1L
      out[(k - 1L) * term_fit$p + col, cc] <- 1
    }
  }
  out
}


.additive_term_equality_A <- function(term_fit){
  P_term <- term_fit$p * (term_fit$K + 1L)
  A <- if(!is.null(term_fit$A)) cbind(term_fit$A) else matrix(0, P_term, 0)
  A <- .additive_drop_offset_constraint_cols(A, term_fit)
  A <- cbind(A,
             .additive_term_flat_A(term_fit),
             .additive_term_zero_main_A(term_fit))
  if(nrow(A) != P_term){
    stop('\n\t Additive equality constraint matrix has incompatible rows.\n')
  }
  A
}


.additive_term_equality_rhs <- function(term_fit, A){
  if(is.null(A) || ncol(A) == 0L) return(numeric(0))
  cv <- term_fit$constraint_values
  if(is.null(cv) || length(cv) == 0L) return(rep(0, ncol(A)))

  cv_vec <- try(Reduce("rbind", cv), silent = TRUE)
  if(any(inherits(cv_vec, "try-error")) ||
     length(cv_vec) == 0L ||
     length(cv_vec) != nrow(A)){
    return(rep(0, ncol(A)))
  }

  rhs <- c(crossprod(A, cbind(cv_vec)))
  if(length(rhs) != ncol(A)) rep(0, ncol(A)) else rhs
}


.additive_combined_equalities <- function(term_fits){
  P_total <- sum(sapply(term_fits, function(fit) fit$p * (fit$K + 1L)))
  mats <- list()
  bvec <- c()
  row_start <- 1L

  for(j in seq_along(term_fits)){
    fit <- term_fits[[j]]
    P_term <- fit$p * (fit$K + 1L)
    A_term <- .additive_term_equality_A(fit)
    if(ncol(A_term) > 0L){
      full <- matrix(0, P_total, ncol(A_term))
      rows <- row_start:(row_start + P_term - 1L)
      full[rows, ] <- A_term
      mats[[length(mats) + 1L]] <- full
      bvec <- c(bvec, .additive_term_equality_rhs(fit, A_term))
    }
    row_start <- row_start + P_term
  }

  if(length(mats) == 0L){
    return(list(Amat = matrix(0, P_total, 0), bvec = numeric(0)))
  }
  list(Amat = do.call(cbind, mats), bvec = bvec)
}


.additive_combined_qp_constraints <- function(term_fits){
  P_total <- sum(sapply(term_fits, function(fit) fit$p * (fit$K + 1L)))
  eq_mats <- list()
  eq_bvec <- c()
  ineq_mats <- list()
  ineq_bvec <- c()
  row_start <- 1L

  for(j in seq_along(term_fits)){
    fit <- term_fits[[j]]
    P_term <- fit$p * (fit$K + 1L)
    qp <- fit$quadprog_list
    if(!is.null(qp) &&
       !identical(qp, list(NA)) &&
       !is.null(qp$qp_Amat) &&
       ncol(cbind(qp$qp_Amat)) > 0L){
      qp_Amat <- cbind(qp$qp_Amat)
      qp_bvec <- c(qp$qp_bvec)
      qp_meq <- if(is.null(qp$qp_meq)) 0L else as.integer(qp$qp_meq)
      rows <- row_start:(row_start + P_term - 1L)

      if(qp_meq > 0L){
        eq_cols <- seq_len(qp_meq)
        full_eq <- matrix(0, P_total, length(eq_cols))
        full_eq[rows, ] <- qp_Amat[, eq_cols, drop = FALSE]
        eq_mats[[length(eq_mats) + 1L]] <- full_eq
        eq_bvec <- c(eq_bvec, qp_bvec[eq_cols])
      }

      if(qp_meq < ncol(qp_Amat)){
        ineq_cols <- (qp_meq + 1L):ncol(qp_Amat)
        full_ineq <- matrix(0, P_total, length(ineq_cols))
        full_ineq[rows, ] <- qp_Amat[, ineq_cols, drop = FALSE]
        ineq_mats[[length(ineq_mats) + 1L]] <- full_ineq
        ineq_bvec <- c(ineq_bvec, qp_bvec[ineq_cols])
      }
    }
    row_start <- row_start + P_term
  }

  mats <- c(eq_mats, ineq_mats)
  if(length(mats) == 0L){
    return(list(Amat = matrix(0, P_total, 0),
                bvec = numeric(0),
                meq = 0L))
  }
  list(Amat = do.call(cbind, mats),
       bvec = c(eq_bvec, ineq_bvec),
       meq = length(eq_bvec))
}


.additive_reduce_constraint_columns <- function(Amat, bvec = NULL){
  nr <- if(is.null(Amat)) 0L else nrow(cbind(Amat))
  if(is.null(Amat) || ncol(cbind(Amat)) == 0L){
    out <- list(Amat = matrix(0, nr, 0))
    if(!is.null(bvec)) out$bvec <- numeric(0)
    return(out)
  }
  Amat <- cbind(Amat)
  qr_A <- qr(Amat, tol = sqrt(.Machine$double.eps))
  if(qr_A$rank == 0L){
    keep <- integer(0)
  } else {
    keep <- sort(qr_A$pivot[seq_len(qr_A$rank)])
  }
  out <- list(Amat = Amat[, keep, drop = FALSE])
  if(!is.null(bvec)) out$bvec <- c(bvec)[keep]
  out
}


.additive_combined_inequalities <- function(term_fits){
  P_total <- sum(sapply(term_fits, function(fit) fit$p * (fit$K + 1L)))
  mats <- list()
  bvec <- c()
  row_start <- 1L
  for(j in seq_along(term_fits)){
    fit <- term_fits[[j]]
    P_term <- fit$p * (fit$K + 1L)
    qp <- fit$quadprog_list
    if(!is.null(qp) &&
       !identical(qp, list(NA)) &&
       !is.null(qp$qp_Amat) &&
       ncol(cbind(qp$qp_Amat)) > 0L){
      qp_Amat <- cbind(qp$qp_Amat)
      qp_bvec <- c(qp$qp_bvec)
      qp_meq <- if(is.null(qp$qp_meq)) 0L else as.integer(qp$qp_meq)
      if(qp_meq < ncol(qp_Amat)){
        keep <- (qp_meq + 1L):ncol(qp_Amat)
        full <- matrix(0, P_total, length(keep))
        rows <- row_start:(row_start + P_term - 1L)
        full[rows, ] <- qp_Amat[, keep, drop = FALSE]
        mats[[length(mats) + 1L]] <- full
        bvec <- c(bvec, qp_bvec[keep])
      }
    }
    row_start <- row_start + P_term
  }
  if(length(mats) == 0L){
    return(list(Amat = matrix(0, P_total, 0), bvec = numeric(0)))
  }
  list(Amat = do.call(cbind, mats), bvec = bvec)
}


.additive_copy_coef_names <- function(coef, template){
  if(!is.null(rownames(template)) && length(rownames(template)) == nrow(coef)){
    rownames(coef) <- rownames(template)
  }
  if(!is.null(names(template)) && length(names(template)) == length(coef)){
    names(coef) <- names(template)
  }
  coef
}


.additive_split_combined_raw_beta <- function(beta_raw, term_fits){
  out <- vector("list", length(term_fits))
  names(out) <- names(term_fits)
  start <- 1L

  for(j in seq_along(term_fits)){
    fit <- term_fits[[j]]
    P_term <- fit$p * (fit$K + 1L)
    beta_term <- c(beta_raw[start:(start + P_term - 1L)])
    term_coefs <- lapply(seq_len(fit$K + 1L), function(k){
      inds <- ((k - 1L) * fit$p + 1L):(k * fit$p)
      template <- fit$B_raw[[k]]
      vals <- matrix(beta_term[inds], ncol = 1L)
      vals <- .additive_copy_coef_names(vals, template)
      vals
    })
    names(term_coefs) <- names(fit$B_raw)
    out[[j]] <- term_coefs
    start <- start + P_term
  }

  out
}


.additive_joint_refine <- function(term_fits, B_raw, y, family,
                                   weights, VhalfInv, sigmasq_tilde,
                                   fit_args, ...){
  N <- length(y)
  P_total <- sum(sapply(term_fits, function(fit) fit$p * (fit$K + 1L)))
  X_full <- do.call(cbind, lapply(term_fits,
                                  .additive_term_design_original))
  Lambda_full <- create_block_diagonal(
    lapply(term_fits, .additive_term_lambda_full)
  )

  eq <- .additive_combined_equalities(term_fits)
  qp <- .additive_combined_qp_constraints(term_fits)
  qp_eq_cols <- if(qp$meq > 0L) seq_len(qp$meq) else integer(0)
  qp_ineq_cols <- if(qp$meq < ncol(qp$Amat)){
    (qp$meq + 1L):ncol(qp$Amat)
  } else {
    integer(0)
  }
  eq_full <- cbind(eq$Amat, qp$Amat[, qp_eq_cols, drop = FALSE])
  b_eq_full <- c(eq$bvec, qp$bvec[qp_eq_cols])
  eq_reduced <- .additive_reduce_constraint_columns(eq_full, b_eq_full)
  ineq_full <- qp$Amat[, qp_ineq_cols, drop = FALSE]
  b_ineq_full <- qp$bvec[qp_ineq_cols]
  qp_Amat_combined <- cbind(eq_reduced$Amat, ineq_full)
  qp_bvec_combined <- c(eq_reduced$bvec, b_ineq_full)
  qp_meq_combined <- ncol(eq_reduced$Amat)

  beta_init <- cbind(unlist(B_raw))
  if(length(beta_init) != P_total){
    stop('\n\t Additive joint coefficient vector has incompatible length.\n')
  }

  y_raw <- cbind(y)
  if(!is.null(VhalfInv)){
    VhalfInv <- methods::as(VhalfInv, "matrix")
    X_design <- VhalfInv %**% X_full
    y_design <- VhalfInv %**% y_raw
    is_gee <- TRUE
    VhalfInv_perm <- VhalfInv
    deviance_fun <- .bf_gee_deviance
  } else {
    X_design <- X_full
    y_design <- y_raw
    is_gee <- FALSE
    VhalfInv_perm <- NULL
    deviance_fun <- .bf_deviance
  }

  sqp <- .bf_sqp_loop(
    X_design = X_design,
    y_design = y_design,
    X_block_raw = X_full,
    beta_init = beta_init,
    Lambda_block = Lambda_full,
    qp_Amat_combined = qp_Amat_combined,
    qp_bvec_combined = qp_bvec_combined,
    qp_meq_combined = qp_meq_combined,
    K = 0L,
    p_expansions = P_total,
    family = family,
    order_list = list(seq_len(N)),
    glm_weight_function = fit_args$glm_weight_function,
    schur_correction_function = fit_args$schur_correction_function,
    qp_score_function = fit_args$qp_score_function,
    need_dispersion_for_estimation = fit_args$need_dispersion_for_estimation,
    dispersion_function = fit_args$dispersion_function,
    observation_weights = list(c(weights)),
    iterate = fit_args$iterate_final_fit,
    tol = fit_args$tol,
    VhalfInv = VhalfInv,
    VhalfInv_perm = VhalfInv_perm,
    is_gee = is_gee,
    deviance_fun = deviance_fun,
    X_partitions = list(X_full),
    y_partitions = list(y_raw),
    verbose = fit_args$verbose,
    ...)

  B_raw_new <- .additive_split_combined_raw_beta(sqp$beta_block, term_fits)
  base_fit <- list(additive_terms = term_fits)
  B_new <- .additive_backtransform_combined_beta(c(sqp$beta_block), base_fit)

  for(j in seq_along(term_fits)){
    term_fits[[j]]$B_raw <- B_raw_new[[j]]
    term_fits[[j]]$B <- B_new[[j]]
  }

  qp_info <- .bf_assemble_qp_info(
    sqp$last_qp_sol,
    sqp$beta_block,
    qp_Amat_combined,
    qp_bvec_combined,
    qp_meq_combined,
    sqp$converged,
    sqp$final_deviance)

  if(!is.null(qp_info)){
    A_active <- qp_info$Amat_active
  } else if(qp_meq_combined > 0L){
    A_active <- qp_Amat_combined[, seq_len(qp_meq_combined), drop = FALSE]
  } else {
    A_active <- matrix(0, P_total, 0)
  }

  list(
    term_fits = term_fits,
    B = B_new,
    B_raw = B_raw_new,
    beta_raw = c(sqp$beta_block),
    qp_info = qp_info,
    A_active = A_active,
    converged = sqp$converged,
    final_deviance = sqp$final_deviance
  )
}


.additive_combined_system <- function(model_fit,
                                      VhalfInv = model_fit$VhalfInv,
                                      sigmasq_tilde = model_fit$sigmasq_tilde,
                                      use_glm_weights = TRUE,
                                      ...){
  term_fits <- model_fit$additive_terms
  N <- model_fit$N
  X_full <- do.call(cbind, lapply(term_fits,
                                  .additive_term_design_original))
  Lambda_full <- create_block_diagonal(
    lapply(term_fits, .additive_term_lambda_full)
  )
  P_total <- ncol(X_full)

  if(!is.null(model_fit$additive_A_active) &&
     nrow(cbind(model_fit$additive_A_active)) == P_total){
    A <- cbind(model_fit$additive_A_active)
  } else {
    A <- .additive_block_diag_rect(
      lapply(term_fits, .additive_term_active_A)
    )
  }
  A <- .additive_reduce_constraint_columns(A)$Amat

  weights <- model_fit$weights
  if(is.null(weights)) weights <- rep(1, N)
  weights <- c(weights)
  if(length(weights) == 1L) weights <- rep(weights, N)

  if(use_glm_weights){
    glm_weight_function <- model_fit$.fit_call_args$glm_weight_function
    if(is.null(glm_weight_function)){
      glm_weight_function <- default_glm_weight_function
    }
    W_glm <- c(glm_weight_function(
      model_fit$ytilde,
      model_fit$y,
      seq_len(N),
      model_fit$family,
      sigmasq_tilde,
      rep(1, N),
      ...
    ))
    W_glm <- pmax(W_glm, .Machine$double.eps)
  } else {
    W_glm <- rep(1, N)
  }

  if(!is.null(VhalfInv)){
    VhalfInv <- methods::as(VhalfInv, 'matrix')
    X_weighted <- VhalfInv %**% (X_full * sqrt(weights * W_glm))
  } else {
    X_weighted <- X_full * sqrt(weights * W_glm)
  }

  gram <- crossprod(X_weighted) + Lambda_full
  G <- invert(gram)
  Ghalf <- matinvsqrt(gram)

  if(!is.null(A) && ncol(A) > 0L){
    GA <- G %**% A
    U <- diag(P_total) -
      GA %**% tcrossprod(invert(crossprod(A, GA)), A)
  } else {
    A <- matrix(0, P_total, 0)
    U <- diag(P_total)
  }

  UGhalf <- U %**% Ghalf
  trace_XUGX <- sum((X_weighted %**% UGhalf)^2)

  inv_sd_y <- unlist(lapply(term_fits, function(fit){
    rep(1 / fit$sd_y, fit$p * (fit$K + 1L))
  }))
  unstd_scale <- unlist(lapply(term_fits, function(fit){
    rep(c(1, 1 / fit$expansion_scales), times = fit$K + 1L)
  }))

  ineq <- .additive_combined_inequalities(term_fits)

  list(
    X = X_full,
    X_weighted = X_weighted,
    Lambda = Lambda_full,
    A = A,
    G = G,
    Ghalf = Ghalf,
    U = U,
    trace_XUGX = trace_XUGX,
    inv_sd_y = inv_sd_y,
    unstd_scale = unstd_scale,
    ineq_Amat = ineq$Amat,
    ineq_bvec = ineq$bvec,
    VhalfInv = VhalfInv
  )
}


.additive_system_varcov <- function(system, sigmasq_tilde){
  UGhalf <- system$U %**% system$Ghalf
  vcov <- tcrossprod(UGhalf) * sigmasq_tilde
  d <- system$unstd_scale
  t(t(vcov * d) * d)
}


.additive_backtransform_combined_beta <- function(beta_raw, model_fit){
  term_fits <- model_fit$additive_terms
  out <- vector("list", length(term_fits))
  names(out) <- names(term_fits)
  start <- 1L
  for(j in seq_along(term_fits)){
    fit <- term_fits[[j]]
    P_term <- fit$p * (fit$K + 1L)
    beta_term <- beta_raw[start:(start + P_term - 1L)]
    term_coefs <- lapply(seq_len(fit$K + 1L), function(k){
      inds <- ((k - 1L) * fit$p + 1L):(k * fit$p)
      bk <- beta_term[inds] * fit$sd_y
      bk[1L] <- bk[1L] + fit$mean_y
      .additive_copy_coef_names(
        fit$backtransform_coefficients(bk),
        fit$B[[k]]
      )
    })
    names(term_coefs) <- names(fit$B)
    out[[j]] <- term_coefs
    start <- start + P_term
  }
  out
}


.additive_draw_combined_coefficients <- function(model_fit,
                                                 system,
                                                 sigmasq,
                                                 enforce_qp_constraints = TRUE,
                                                 max_slice_iterations = 1000L,
                                                 include_warnings = TRUE){
  beta_mode <- unlist(model_fit$B_raw)
  L_post <- t(t(system$U %**% system$Ghalf) * system$inv_sd_y)

  has_ineq <- enforce_qp_constraints &&
    !is.null(system$ineq_Amat) &&
    ncol(system$ineq_Amat) > 0L

  .check_feasible <- function(b){
    if(!has_ineq) return(TRUE)
    all(c(crossprod(system$ineq_Amat, cbind(b))) >=
          system$ineq_bvec - sqrt(.Machine$double.eps))
  }

  if(has_ineq){
    beta_cur <- beta_mode
    nu <- sqrt(sigmasq) * c(L_post %**% cbind(rnorm(ncol(L_post))))
    theta <- runif(1, 0, 2 * pi)
    theta_min <- theta - 2 * pi
    theta_max <- theta
    accepted <- FALSE
    for(iter in seq_len(max_slice_iterations)){
      beta_prop <- beta_cur * cos(theta) + nu * sin(theta)
      if(.check_feasible(beta_prop)){
        beta_cur <- beta_prop
        accepted <- TRUE
        break
      }
      if(theta < 0) theta_min <- theta else theta_max <- theta
      theta <- runif(1, theta_min, theta_max)
    }
    if(!accepted && include_warnings){
      warning('\n\t Slice sampler: no feasible point in ',
              max_slice_iterations, ' iterations, returning current state.\n')
    }
    beta_raw <- beta_cur
  } else {
    z <- rnorm(ncol(L_post))
    beta_raw <- beta_mode + sqrt(sigmasq) * c(L_post %**% cbind(z))
  }

  .additive_backtransform_combined_beta(beta_raw, model_fit)
}


.additive_generate_posterior_from_system <- function(model_fit,
                                                     system,
                                                     new_sigmasq_tilde,
                                                     new_predictors,
                                                     theta_1,
                                                     theta_2,
                                                     posterior_predictive_draw,
                                                     draw_dispersion,
                                                     include_posterior_predictive,
                                                     num_draws,
                                                     enforce_qp_constraints,
                                                     max_slice_iterations,
                                                     include_warnings,
                                                     ...){
  if(is.null(new_predictors)){
    new_predictors <- model_fit$predictors
  } else {
    new_predictors <- .additive_prepare_new_predictors(
      new_predictors,
      model_fit$og_cols,
      include_warnings)
  }

  .draw_sigmasq <- function(){
    if(!draw_dispersion) return(new_sigmasq_tilde)
    half_df <- 0.5 * (model_fit$N - isTRUE(model_fit$unbias_dispersion) *
                        system$trace_XUGX)
    shape <- theta_1 + half_df
    rate <- theta_2 + half_df * new_sigmasq_tilde
    if(shape <= 0){
      stop('\n\t Posterior inverse-gamma shape <= 0, increase theta_1.\n')
    }
    if(rate <= 0){
      stop('\n\t Posterior inverse-gamma rate <= 0, increase theta_2.\n')
    }
    s2 <- 1 / rgamma(1, shape, rate)
    if(is.nan(s2) || !is.finite(s2)){
      if(include_warnings){
        warning('\n\t Infinite/NaN dispersion draw, using MLE.\n')
      }
      s2 <- new_sigmasq_tilde
    }
    s2
  }

  .one_draw <- function(){
    s2 <- .draw_sigmasq()
    coefs <- .additive_draw_combined_coefficients(
      model_fit,
      system = system,
      sigmasq = s2,
      enforce_qp_constraints = enforce_qp_constraints,
      max_slice_iterations = max_slice_iterations,
      include_warnings = include_warnings)
    if(include_posterior_predictive){
      pm <- model_fit$predict(new_predictors = new_predictors,
                              B_predict = coefs)
      pp <- posterior_predictive_draw(length(pm), pm, sqrt(s2), ...)
      list(post_pred_draw = pp,
           post_draw_coefficients = coefs,
           post_draw_sigmasq = s2)
    } else {
      list(post_draw_coefficients = coefs,
           post_draw_sigmasq = s2)
    }
  }

  res <- lapply(seq_len(num_draws), function(m) .one_draw())
  if(num_draws == 1L) return(res[[1L]])

  out <- list(
    post_draw_coefficients = lapply(res, `[[`, "post_draw_coefficients"),
    post_draw_sigmasq = lapply(res, `[[`, "post_draw_sigmasq"))
  if(include_posterior_predictive){
    out$post_pred_draw <- Reduce("cbind",
                                 lapply(res, `[[`, "post_pred_draw"))
  }
  out
}


.additive_generate_posterior_correlation <- function(
    object,
    new_sigmasq_tilde = object$sigmasq_tilde,
    new_predictors = NULL,
    theta_1 = 0,
    theta_2 = 0,
    posterior_predictive_draw =
      function(N, mean, sqrt_dispersion, ...) {
        rnorm(N, mean, sqrt_dispersion)
      },
    draw_dispersion = TRUE,
    include_posterior_predictive = FALSE,
    num_draws = 1,
    enforce_qp_constraints = TRUE,
    correlation_param_mean = NULL,
    correlation_param_vcov_sc = NULL,
    correlation_VhalfInv_fxn = NULL,
    correlation_Vhalf_fxn = NULL,
    include_warnings = TRUE,
    max_slice_iterations = 1000L,
    ...
){
  if(is.null(correlation_param_mean)){
    correlation_param_mean <- object$VhalfInv_params_estimates
    if(is.null(correlation_param_mean)){
      stop(
        "\n\t Cannot draw additive correlation parameters: no point ",
        "estimates found and no 'correlation_param_mean' was supplied.\n"
      )
    }
  }
  correlation_param_mean <- c(correlation_param_mean)
  n_corr_par <- length(correlation_param_mean)

  if(is.null(correlation_param_vcov_sc)){
    correlation_param_vcov_sc <- object$VhalfInv_params_vcov
    if(is.null(correlation_param_vcov_sc)){
      stop(
        "\n\t Cannot draw additive correlation parameters: no inverse ",
        "Hessian found and no 'correlation_param_vcov_sc' was supplied.\n"
      )
    }
  }
  correlation_param_vcov_sc <- as.matrix(correlation_param_vcov_sc)
  if(nrow(correlation_param_vcov_sc) != n_corr_par ||
     ncol(correlation_param_vcov_sc) != n_corr_par){
    stop("\n\t 'correlation_param_vcov_sc' has incompatible dimensions.\n")
  }

  if(is.null(correlation_VhalfInv_fxn)){
    correlation_VhalfInv_fxn <- object$VhalfInv_fxn
    if(is.null(correlation_VhalfInv_fxn)){
      stop(
        "\n\t No VhalfInv_fxn found and no ",
        "'correlation_VhalfInv_fxn' was supplied.\n"
      )
    }
  }
  if(is.null(correlation_Vhalf_fxn)){
    correlation_Vhalf_fxn <- object$Vhalf_fxn
  }

  vcov_chol <- tryCatch({
    chol(correlation_param_vcov_sc)
  }, error = function(e){
    eig <- eigen(correlation_param_vcov_sc, symmetric = TRUE)
    min_eval <- max(abs(eig$values)) * sqrt(.Machine$double.eps)
    eig$values <- pmax(eig$values, min_eval)
    vcov_pd <- crossprod(t(eig$vectors) * sqrt(eig$values))
    if(include_warnings){
      warning(
        "\n\t correlation_param_vcov_sc is not positive definite; ",
        "projecting to nearest PD matrix.\n"
      )
    }
    chol(vcov_pd)
  })

  .one_corr_draw <- function(){
    max_corr_reject <- 50L
    for(attempt in seq_len(max_corr_reject)){
      phi_draw <- correlation_param_mean +
        c(crossprod(vcov_chol, rnorm(n_corr_par)))
      tr <- try({
        VhalfInv_draw <- correlation_VhalfInv_fxn(phi_draw)
        stopifnot(
          is.matrix(VhalfInv_draw),
          all(is.finite(VhalfInv_draw)),
          nrow(VhalfInv_draw) == object$N,
          ncol(VhalfInv_draw) == object$N
        )
        system <- .additive_combined_system(
          object,
          VhalfInv = VhalfInv_draw,
          sigmasq_tilde = new_sigmasq_tilde,
          use_glm_weights = TRUE,
          ...
        )
      }, silent = TRUE)
      if(!inherits(tr, "try-error")){
        one_draw <- .additive_generate_posterior_from_system(
          object,
          system = system,
          new_sigmasq_tilde = new_sigmasq_tilde,
          new_predictors = new_predictors,
          theta_1 = theta_1,
          theta_2 = theta_2,
          posterior_predictive_draw = posterior_predictive_draw,
          draw_dispersion = draw_dispersion,
          include_posterior_predictive = include_posterior_predictive,
          num_draws = 1L,
          enforce_qp_constraints = enforce_qp_constraints,
          max_slice_iterations = max_slice_iterations,
          include_warnings = include_warnings,
          ...
        )
        one_draw$post_draw_correlation_params <- phi_draw
        return(one_draw)
      }
    }
    stop(
      "\n\t Could not generate a valid additive correlation draw after ",
      max_corr_reject, " attempts.\n"
    )
  }

  res <- lapply(seq_len(num_draws), function(m) .one_corr_draw())
  if(num_draws == 1L) return(res[[1L]])

  out <- list(
    post_draw_coefficients = lapply(res, `[[`, "post_draw_coefficients"),
    post_draw_sigmasq = lapply(res, `[[`, "post_draw_sigmasq"),
    post_draw_correlation_params =
      lapply(res, `[[`, "post_draw_correlation_params"))
  if(include_posterior_predictive){
    out$post_pred_draw <- Reduce("cbind",
                                 lapply(res, `[[`, "post_pred_draw"))
  }
  out
}


.additive_plot_function <- function(model_fit,
                                    custom_response_lab = "y",
                                    custom_predictor_lab = NULL,
                                    custom_title = "Fitted Function",
                                    new_predictors = NULL,
                                    add = FALSE,
                                    vars = c(),
                                    n_grid = 200,
                                    ...){
  q <- model_fit$q
  og_cols <- model_fit$og_cols

  if(is.null(new_predictors)){
    if(length(vars) == 0L){
      vars <- if(q == 1L) 1L else model_fit$spline_groups[[1L]][1L]
    }
    var_idx <- if(is.character(vars)) match(vars[1L], og_cols) else vars[1L]
    vals <- seq(min(model_fit$predictors[, var_idx]),
                max(model_fit$predictors[, var_idx]),
                length.out = n_grid)
    new_predictors <- matrix(0, n_grid, q)
    colnames(new_predictors) <- og_cols
    new_predictors[, var_idx] <- vals
  } else {
    new_predictors <- as.matrix(new_predictors)
    var_idx <- if(length(vars) > 0L){
      if(is.character(vars)) match(vars[1L], og_cols) else vars[1L]
    } else {
      1L
    }
  }

  fit <- model_fit$predict(new_predictors = new_predictors)
  x <- new_predictors[, var_idx]
  ord <- order(x)
  xlab <- if(!is.null(custom_predictor_lab)){
    custom_predictor_lab
  } else if(!is.null(og_cols)){
    og_cols[var_idx]
  } else {
    "x"
  }

  if(add){
    graphics::lines(x[ord], fit[ord], ...)
  } else {
    graphics::plot(x[ord], fit[ord], type = "l",
                   xlab = xlab, ylab = custom_response_lab,
                   main = custom_title, ...)
  }
  invisible(NULL)
}


.fit_lgspline_additive <- function(fit_args, ...){
  dots <- list(...)

  predictors <- methods::as(fit_args$predictors, 'matrix')
  y <- c(fit_args$y)
  q_predictors <- ncol(predictors)
  N_obs <- nrow(predictors)
  pred_colnames <- colnames(predictors)
  groups <- .additive_normalize_groups(fit_args$spline_groups,
                                       q_predictors,
                                       pred_colnames)
  if(is.null(groups) || length(groups) <= 1L){
    stop('\n \t Additive fitting requires at least two spline groups.\n')
  }
  n_terms <- length(groups)

  global_equalities <-
    (!.additive_empty_constraint(fit_args$constraint_vectors) &&
       !.additive_is_per_term_arg(fit_args$constraint_vectors, n_terms)) ||
    (!.additive_empty_constraint(fit_args$constraint_values) &&
       !.additive_is_per_term_arg(fit_args$constraint_values, n_terms))
  global_inequalities <-
    (!is.null(fit_args$qp_Amat) &&
       !.additive_is_per_term_arg(fit_args$qp_Amat, n_terms)) ||
    (!is.null(fit_args$qp_bvec) &&
       !.additive_is_per_term_arg(fit_args$qp_bvec, n_terms))

  if(global_equalities || global_inequalities){
    stop('\n \t Global pre-built equality or QP matrices are not supported ',
         'for additive spline_groups because there is no single joined ',
         'coefficient vector. Supply per-term lists, use qp_Amat_fxn/',
         'qp_bvec_fxn/qp_meq_fxn, use built-in QP constraints, or fit a ',
         'joined spline term.\n')
  }

  critical_value <- if(!is.null(fit_args$critical_value)){
    fit_args$critical_value
  } else {
    qnorm(1 - 0.05 / 2)
  }
  spline_cols <- sort(unique(unlist(groups)))
  supplied_linear <- sort(unique(c(fit_args$just_linear_with_interactions,
                                   fit_args$just_linear_without_interactions)))
  anchor_cols <- sort(unique(c(setdiff(seq_len(q_predictors), spline_cols),
                               setdiff(supplied_linear, spline_cols))))

  explicit_pairs <- fit_args$additive_spline_interaction_pairs
  if(is.null(explicit_pairs)) explicit_pairs <- list()
  if(!is.list(explicit_pairs)) explicit_pairs <- list(explicit_pairs)
  explicit_pairs <- lapply(explicit_pairs, function(pair){
    unique(as.integer(pair[!is.na(pair)]))
  })
  explicit_pairs <- explicit_pairs[
    sapply(explicit_pairs, length) == 2L
  ]
  interaction_anchor_by_term <- vector("list", n_terms)
  for(j in seq_len(n_terms)) interaction_anchor_by_term[[j]] <- integer(0)
  for(pair in explicit_pairs){
    if(!all(pair %in% spline_cols)) next
    owner <- which(vapply(groups, function(grp) pair[1L] %in% grp,
                          logical(1)))[1L]
    if(is.na(owner)) next
    other <- pair[2L]
    if(other %in% groups[[owner]]) next
    interaction_anchor_by_term[[owner]] <- unique(c(
      interaction_anchor_by_term[[owner]], other))
  }
  term_anchor_cols_list <- lapply(seq_len(n_terms), function(j){
    sort(unique(c(anchor_cols, interaction_anchor_by_term[[j]])))
  })

  linear_with <- intersect(fit_args$just_linear_with_interactions, anchor_cols)
  linear_without <- intersect(fit_args$just_linear_without_interactions,
                              anchor_cols)
  if(length(linear_with) + length(linear_without) == 0L &&
     length(anchor_cols) > 0L){
    linear_without <- anchor_cols
  } else if(length(anchor_cols) > 0L){
    linear_without <- unique(c(linear_without,
                               setdiff(anchor_cols, linear_with)))
  }

  additive_max_iter <- if(!is.null(dots$additive_max_iter)){
    dots$additive_max_iter
  } else {
    50L
  }
  additive_tol <- if(!is.null(dots$additive_tol)){
    dots$additive_tol
  } else {
    fit_args$tol
  }

  term_eta <- matrix(0, N_obs, n_terms)
  term_fits <- vector("list", n_terms)
  names(term_fits) <- paste0("smooth", seq_len(n_terms))
  dev_history <- c()

  base_fit_args <- fit_args
  base_fit_args$spline_groups <- NULL
  base_fit_args$constraint_vectors <- cbind()
  base_fit_args$constraint_values <- cbind()
  base_fit_args$qp_Amat <- NULL
  base_fit_args$qp_bvec <- NULL
  base_fit_args$qp_meq <- 0
  base_fit_args$qp_Amat_fxn <- NULL
  base_fit_args$qp_bvec_fxn <- NULL
  base_fit_args$qp_meq_fxn <- NULL
  base_fit_args$factor_groups <- NULL

  for(iter in seq_len(if(fit_args$dummy_fit) 1L else additive_max_iter)){
    eta_old <- rowSums(term_eta)

    for(j in seq_len(n_terms)){
      offset_eta <- rowSums(term_eta[, -j, drop = FALSE])
      term_anchor_cols <- term_anchor_cols_list[[j]]
      local_cols <- .additive_term_cols(groups[[j]], term_anchor_cols, j)
      term_predictors <- .additive_make_term_predictors(
        predictors, groups[[j]], term_anchor_cols, offset_eta, j,
        pred_colnames)

      local_offset <- ncol(term_predictors)
      local_group <- match(groups[[j]], local_cols)
      local_anchor <- match(term_anchor_cols, local_cols)
      local_anchor <- local_anchor[!is.na(local_anchor)]
      local_interaction_anchor <- match(interaction_anchor_by_term[[j]],
                                        local_cols)
      local_interaction_anchor <-
        local_interaction_anchor[!is.na(local_interaction_anchor)]

      local_linear_with <- match(unique(c(linear_with,
                                          interaction_anchor_by_term[[j]])),
                                 local_cols)
      local_linear_with <- local_linear_with[!is.na(local_linear_with)]
      local_linear_without <- match(linear_without, local_cols)
      local_linear_without <- local_linear_without[!is.na(local_linear_without)]

      term_args <- base_fit_args
      term_args$predictors <- term_predictors
      term_args$y <- y
      term_args$K <- .additive_arg_for_term(fit_args$K, j, n_terms)
      term_args$custom_knots <-
        .additive_arg_for_term(fit_args$custom_knots, j, n_terms)
      term_args$make_partition_list <-
        .additive_arg_for_term(fit_args$make_partition_list, j, n_terms)
      term_args$previously_tuned_penalties <-
        .additive_arg_for_term(fit_args$previously_tuned_penalties,
                               j, n_terms)
      term_args$predictor_penalties <-
        .additive_arg_for_term(fit_args$predictor_penalties, j, n_terms)
      term_args$partition_penalties <-
        .additive_arg_for_term(fit_args$partition_penalties, j, n_terms)
      term_args$constraint_vectors <-
        .additive_arg_for_term(fit_args$constraint_vectors, j, n_terms)
      if(.additive_empty_constraint(term_args$constraint_vectors)){
        term_args$constraint_vectors <- cbind()
      }
      term_args$constraint_values <-
        .additive_arg_for_term(fit_args$constraint_values, j, n_terms)
      if(.additive_empty_constraint(term_args$constraint_values)){
        term_args$constraint_values <- cbind()
      }
      term_args$qp_Amat <- .additive_arg_for_term(fit_args$qp_Amat,
                                                  j, n_terms)
      term_args$qp_bvec <- .additive_arg_for_term(fit_args$qp_bvec,
                                                  j, n_terms)
      term_args$qp_meq <- .additive_arg_for_term(fit_args$qp_meq,
                                                 j, n_terms)
      if(is.null(term_args$qp_meq)) term_args$qp_meq <- 0
      term_args$qp_Amat_fxn <-
        .additive_arg_for_term(fit_args$qp_Amat_fxn, j, n_terms)
      term_args$qp_bvec_fxn <-
        .additive_arg_for_term(fit_args$qp_bvec_fxn, j, n_terms)
      term_args$qp_meq_fxn <-
        .additive_arg_for_term(fit_args$qp_meq_fxn, j, n_terms)
      term_args$offset <- local_offset
      term_args$just_linear_with_interactions <- local_linear_with
      term_args$just_linear_without_interactions <-
        unique(c(local_linear_without, local_offset))
      term_args$do_not_cluster_on_these <-
        unique(c(local_anchor, local_offset))
      term_args$no_intercept <- isTRUE(fit_args$no_intercept) || j > 1L
      term_args$include_2way_interactions <-
        isTRUE(fit_args$include_2way_interactions) &&
        length(unique(c(local_group, local_linear_with))) > 1L
      term_args$include_3way_interactions <-
        isTRUE(fit_args$include_3way_interactions) &&
        length(unique(c(local_group, local_linear_with))) > 2L
      term_args$include_quadratic_interactions <-
        isTRUE(fit_args$include_quadratic_interactions) &&
        length(unique(c(local_group, local_linear_with))) > 1L
      term_args$exclude_interactions_for <- c()
      term_args$exclude_these_expansions <-
        .additive_local_exclude_expansions(
          fit_args$exclude_these_expansions, local_cols)
      term_args$og_cols <- colnames(term_predictors)
      term_args$factor_groups <-
        .additive_local_factor_groups(fit_args$factor_groups, local_cols)

      fit_j <- do.call(lgspline.fit, c(term_args, dots))
      anchor_only_cols <- .additive_anchor_only_expansion_cols(
        fit_j$raw_expansion_names, local_anchor)
      interaction_anchor_only_cols <- .additive_anchor_only_expansion_cols(
        fit_j$raw_expansion_names, local_interaction_anchor)
      fit_j$additive_force_flat_cols <- which(
        fit_j$raw_expansion_names %in% paste0("_", local_linear_with, "_")
      )
      fit_j$additive_force_flat_cols <- unique(c(
        fit_j$additive_force_flat_cols,
        anchor_only_cols
      ))
      fit_j$additive_zero_main_cols <- if(j > 1L){
        unique(c(anchor_only_cols, interaction_anchor_only_cols))
      } else {
        interaction_anchor_only_cols
      }
      fit_j <- .additive_rename_fit(fit_j, colnames(term_predictors))
      term_fits[[j]] <- fit_j

      zero_predictors <- .additive_zero_offset_predictors(term_predictors)
      term_eta[, j] <- .additive_term_eta(fit_j, zero_predictors,
                                          fit_args$family)
    }

    eta_new <- rowSums(term_eta)
    mu_new <- fit_args$family$linkinv(eta_new)
    dev_history <- c(dev_history,
                     .additive_dev(y, mu_new, fit_args$family,
                                   fit_args$observation_weights,
                                   seq_len(N_obs), ...))

    if(fit_args$dummy_fit ||
       max(abs(eta_new - eta_old), na.rm = TRUE) < additive_tol){
      break
    }
  }

  eta <- rowSums(term_eta)
  ytilde <- fit_args$family$linkinv(eta)
  weights <- if(is.null(fit_args$observation_weights)){
    rep(1, N_obs)
  } else {
    c(fit_args$observation_weights)
  }
  if(length(weights) == 1L) weights <- rep(weights, N_obs)

  estimate_additive_sigmasq <- function(yfit){
    if(!fit_args$estimate_dispersion) return(1)
    if(fit_args$family$family == "gaussian" &&
       fit_args$family$link == "identity"){
      if(!is.null(fit_args$VhalfInv)){
        resid <- c(fit_args$VhalfInv %**% cbind(y - yfit))
      } else {
        resid <- y - yfit
      }
      return(mean(weights * resid^2))
    } else {
      return(fit_args$dispersion_function(
        mu = yfit,
        y = y,
        order_indices = seq_len(N_obs),
        family = fit_args$family,
        observation_weights = weights,
        VhalfInv = fit_args$VhalfInv,
        ...))
    }
  }

  sigmasq_tilde <- estimate_additive_sigmasq(ytilde)

  B <- lapply(term_fits, function(fit) fit$B)
  B_raw <- lapply(term_fits, function(fit) fit$B_raw)
  names(B) <- names(B_raw) <- names(term_fits)

  additive_fit_call_args <- list(
    og_cols = pred_colnames,
    glm_weight_function = fit_args$glm_weight_function,
    schur_correction_function = fit_args$schur_correction_function,
    qp_score_function = fit_args$qp_score_function,
    need_dispersion_for_estimation = fit_args$need_dispersion_for_estimation,
    dispersion_function = fit_args$dispersion_function,
    iterate_final_fit = fit_args$iterate_final_fit,
    tol = fit_args$tol
  )

  joint_qp_info <- NULL
  additive_A_active <- NULL
  if(!fit_args$dummy_fit){
    joint <- .additive_joint_refine(
      term_fits = term_fits,
      B_raw = B_raw,
      y = y,
      family = fit_args$family,
      weights = weights,
      VhalfInv = fit_args$VhalfInv,
      sigmasq_tilde = sigmasq_tilde,
      fit_args = fit_args,
      ...)

    term_fits <- joint$term_fits
    B <- joint$B
    B_raw <- joint$B_raw
    joint_qp_info <- joint$qp_info
    additive_A_active <- joint$A_active

    for(j in seq_len(n_terms)){
      term_predictors <- .additive_make_term_predictors(
        predictors, groups[[j]], term_anchor_cols_list[[j]],
        rep(0, N_obs), j,
        pred_colnames)
      term_eta[, j] <- .additive_term_eta(
        term_fits[[j]], term_predictors, fit_args$family,
        B_predict = B[[j]])
    }
    eta <- rowSums(term_eta)
    ytilde <- fit_args$family$linkinv(eta)
    sigmasq_tilde <- estimate_additive_sigmasq(ytilde)
  }

  A <- .additive_block_diag_rect(lapply(term_fits, function(fit){
    if(is.null(fit$A)) matrix(0, fit$P, 0) else cbind(fit$A)
  }))

  combined_base <- list(
    y = y,
    ytilde = ytilde,
    additive_terms = term_fits,
    B = B,
    B_raw = B_raw,
    N = N_obs,
    family = fit_args$family,
    estimate_dispersion = fit_args$estimate_dispersion,
    unbias_dispersion = fit_args$unbias_dispersion,
    weights = weights,
    VhalfInv = fit_args$VhalfInv,
    sigmasq_tilde = sigmasq_tilde,
    predictors = predictors,
    og_cols = pred_colnames,
    additive_A_active = additive_A_active,
    .fit_call_args = additive_fit_call_args
  )

  combined_system <- .additive_combined_system(
    combined_base,
    VhalfInv = fit_args$VhalfInv,
    sigmasq_tilde = sigmasq_tilde,
    use_glm_weights = TRUE,
    ...
  )
  trace_XUGX <- combined_system$trace_XUGX

  if(isTRUE(fit_args$unbias_dispersion) &&
     fit_args$estimate_dispersion &&
     N_obs > trace_XUGX){
    sigmasq_tilde <- sigmasq_tilde * N_obs / (N_obs - trace_XUGX)
    combined_base$sigmasq_tilde <- sigmasq_tilde
    combined_system <- .additive_combined_system(
      combined_base,
      VhalfInv = fit_args$VhalfInv,
      sigmasq_tilde = sigmasq_tilde,
      use_glm_weights = TRUE,
      ...
    )
    trace_XUGX <- combined_system$trace_XUGX
  }

  varcovmat <- NULL
  if(fit_args$return_varcovmat){
    varcovmat <- .additive_system_varcov(combined_system, sigmasq_tilde)
  }

  predict_function <- function(new_predictors = predictors,
                               parallel = FALSE,
                               cl = NULL,
                               chunk_size = NULL,
                               num_chunks = NULL,
                               rem_chunks = NULL,
                               B_predict = B,
                               take_first_derivatives = FALSE,
                               take_second_derivatives = FALSE,
                               expansions_only = FALSE,
                               se.fit = FALSE,
                               cv = 1.96){
    if(any(!is.null(new_predictors))){
      new_predictors <- .additive_prepare_new_predictors(
        new_predictors,
        pred_colnames,
        fit_args$include_warnings)
    } else {
      new_predictors <- predictors
    }

    N_new <- nrow(new_predictors)
    term_eta_new <- matrix(0, N_new, n_terms)
    expansions <- vector("list", n_terms)

    for(j in seq_len(n_terms)){
      offset_zero <- rep(0, N_new)
      term_predictors <- .additive_make_term_predictors(
        new_predictors, groups[[j]], term_anchor_cols_list[[j]],
        offset_zero, j,
        pred_colnames)
      if(expansions_only){
        expansions[[j]] <- term_fits[[j]]$predict(
          new_predictors = term_predictors,
          expansions_only = TRUE)
      } else {
        term_B <- .additive_B_for_term(B_predict, j, n_terms)
        term_eta_new[, j] <-
          .additive_term_eta(term_fits[[j]], term_predictors,
                             fit_args$family,
                             B_predict = term_B)
      }
    }

    if(expansions_only){
      names(expansions) <- names(term_fits)
      return(list(expansions = expansions,
                  spline_groups = groups))
    }

    eta_new <- rowSums(term_eta_new)
    final_preds <- fit_args$family$linkinv(eta_new)

    if(take_first_derivatives || take_second_derivatives){
      hbase <- sqrt(.Machine$double.eps)
      first <- vector("list", q_predictors)
      second <- vector("list", q_predictors)
      for(v in seq_len(q_predictors)){
        h <- hbase * pmax(abs(new_predictors[, v]), 1)
        plus <- new_predictors
        minus <- new_predictors
        plus[, v] <- plus[, v] + h
        minus[, v] <- minus[, v] - h
        f_plus <- predict_function(plus)
        f_minus <- predict_function(minus)
        first[[v]] <- (f_plus - f_minus) / (2 * h)
        if(take_second_derivatives){
          second[[v]] <- (f_plus - 2 * final_preds + f_minus) / (h^2)
        }
      }
      if(!is.null(pred_colnames)){
        names(first) <- pred_colnames
        names(second) <- pred_colnames
      }
      out <- list(
        preds = final_preds,
        first_deriv = if(q_predictors == 1L) first[[1L]] else first,
        second_deriv = if(take_second_derivatives){
          if(q_predictors == 1L) second[[1L]] else second
        } else {
          NULL
        })
      return(out)
    }

    if(se.fit){
      if(is.null(varcovmat)){
        if(fit_args$include_warnings){
          warning('\n\t se.fit = TRUE requires return_varcovmat = TRUE; ',
                  'returning fitted values with NA standard errors.\n')
        }
        return(list(fit = final_preds,
                    se.fit = rep(NA_real_, length(final_preds)),
                    lower = rep(NA_real_, length(final_preds)),
                    upper = rep(NA_real_, length(final_preds)),
                    cv = cv))
      }
      X_full_new <- do.call(cbind, lapply(seq_len(n_terms), function(j){
        offset_zero <- rep(0, N_new)
        term_predictors <- .additive_make_term_predictors(
          new_predictors, groups[[j]], term_anchor_cols_list[[j]],
          offset_zero, j,
          pred_colnames)
        .additive_new_design_original(term_fits[[j]], term_predictors)
      }))
      XV <- X_full_new %**% varcovmat
      var_eta <- pmax(rowSums(XV * X_full_new), 0)
      se_link <- sqrt(var_eta)
      lower <- fit_args$family$linkinv(eta_new - cv * se_link)
      upper <- fit_args$family$linkinv(eta_new + cv * se_link)
      return(list(fit = final_preds,
                  se.fit = se_link,
                  lower = lower,
                  upper = upper,
                  cv = cv))
    }

    final_preds
  }

  wald_univariate <- function(scale_vcovmat_by = 1,
                              cv = critical_value){
    est <- unlist(B)
    if(!is.null(varcovmat) && length(est) == nrow(varcovmat)){
      se <- sqrt(scale_vcovmat_by * diag(varcovmat))
    } else {
      se <- rep(NA_real_, length(est))
    }
    stat <- est / se
    tab <- cbind(
      Estimate = est,
      Std.Error = se,
      Statistic = stat,
      Lower = est - cv * se,
      Upper = est + cv * se,
      p.value = 2 * (1 - pnorm(abs(stat)))
    )
    colnames(tab)[6] <- "Pr(>|z|)"
    out <- list(coefficients = tab,
                est = est,
                se = se,
                stat = stat,
                interval_lb = tab[, "Lower"],
                interval_ub = tab[, "Upper"],
                pval = tab[, "Pr(>|z|)"],
                df.residual = Inf,
                test_type = "z",
                critical_value = cv,
                scale_vcovmat_by = scale_vcovmat_by,
                K = sapply(term_fits, `[[`, "K"),
                p = sapply(term_fits, `[[`, "p"),
                N = N_obs,
                family = fit_args$family)
    class(out) <- "wald_lgspline"
    out
  }

  generate_posterior <- function(new_sigmasq_tilde = sigmasq_tilde,
                                 new_predictors = predictors,
                                 theta_1 = 0,
                                 theta_2 = 0,
                                 posterior_predictive_draw =
                                   function(N_obs, mean, sqrt_dispersion, ...)
                                     rnorm(N_obs, mean, sqrt_dispersion),
                                 draw_dispersion = TRUE,
                                 include_posterior_predictive = FALSE,
                                 num_draws = 1,
                                 enforce_constraints = FALSE,
                                 enforce_qp_constraints = NULL,
                                 max_slice_iterations = 1000,
                                 draw_correlation = FALSE,
                                 include_warnings = fit_args$include_warnings,
                                 ...){
    if(!is.null(enforce_qp_constraints)){
      enforce_constraints <- enforce_qp_constraints
    }
    if(draw_correlation){
      stop('\n\t Use generate_posterior(object, draw_correlation = TRUE) ',
           'so the additive correlation sampler receives the fitted ',
           'correlation parameter fields.\n')
    }

    system <- .additive_combined_system(
      return_list,
      VhalfInv = return_list$VhalfInv,
      sigmasq_tilde = new_sigmasq_tilde,
      use_glm_weights = TRUE,
      ...
    )

    .additive_generate_posterior_from_system(
      return_list,
      system = system,
      new_sigmasq_tilde = new_sigmasq_tilde,
      new_predictors = new_predictors,
      theta_1 = theta_1,
      theta_2 = theta_2,
      posterior_predictive_draw = posterior_predictive_draw,
      draw_dispersion = draw_dispersion,
      include_posterior_predictive = include_posterior_predictive,
      num_draws = num_draws,
      enforce_qp_constraints = enforce_constraints,
      max_slice_iterations = max_slice_iterations,
      include_warnings = include_warnings,
      ...
    )
  }


  find_extremum <- function(vars = NULL,
                            quick_heuristic = TRUE,
                            initial = NULL,
                            B_predict = NULL,
                            minimize = FALSE,
                            stochastic = FALSE,
                            stochastic_draw = function(mu, sigma, ...){
                              rnorm(length(mu), mu, sigma)
                            },
                            sigmasq_predict = sigmasq_tilde,
                            custom_objective_function = NULL,
                            custom_objective_derivative = NULL,
                            ...){
    if(is.null(vars)){
      vars_idx <- seq_len(q_predictors)
    } else if(is.character(vars)){
      if(is.null(pred_colnames)){
        stop('\n\t Character vars require named predictor columns.\n')
      }
      vars_idx <- match(vars, pred_colnames)
      if(any(is.na(vars_idx))){
        stop('\n\t vars not found in predictor names.\n')
      }
    } else {
      vars_idx <- as.integer(vars)
    }
    vars_idx <- unique(vars_idx)
    if(any(vars_idx < 1L | vars_idx > q_predictors)){
      stop('\n\t vars contains indices outside predictor columns.\n')
    }

    base_vec <- colMeans(predictors, na.rm = TRUE)
    if(!is.null(initial)){
      if(!is.null(names(initial)) && !is.null(pred_colnames)){
        matched <- match(names(initial), pred_colnames)
        matched <- matched[!is.na(matched)]
        base_vec[matched] <- initial[names(initial) %in% pred_colnames]
      } else if(length(initial) == q_predictors){
        base_vec <- as.numeric(initial)
      } else if(length(initial) == length(vars_idx)){
        base_vec[vars_idx] <- as.numeric(initial)
      }
    }

    lower <- apply(predictors[, vars_idx, drop = FALSE], 2, min, na.rm = TRUE)
    upper <- apply(predictors[, vars_idx, drop = FALSE], 2, max, na.rm = TRUE)
    same <- !is.finite(lower) | !is.finite(upper) | lower == upper
    if(any(same)){
      lower[same] <- base_vec[vars_idx][same] - 1
      upper[same] <- base_vec[vars_idx][same] + 1
    }
    start <- pmin(pmax(base_vec[vars_idx], lower), upper)

    y_best <- if(minimize) min(y, na.rm = TRUE) else max(y, na.rm = TRUE)
    .eval <- function(par){
      x0 <- base_vec
      x0[vars_idx] <- par
      pred <- c(predict_function(new_predictors = matrix(x0, nrow = 1),
                                 B_predict = B_predict))
      val <- if(!is.null(custom_objective_function)){
        custom_objective_function(pred, sqrt(sigmasq_predict), y_best, ...)
      } else {
        pred
      }
      if(stochastic){
        val <- stochastic_draw(val, sqrt(sigmasq_predict), ...)
      }
      c(val)[1L]
    }
    .optim_obj <- function(par){
      val <- .eval(par)
      if(minimize) val else -val
    }

    opt <- stats::optim(
      par = start,
      fn = .optim_obj,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper)

    t_out <- base_vec
    t_out[vars_idx] <- opt$par
    if(!is.null(pred_colnames)) names(t_out) <- pred_colnames
    list(t = t_out, y = .eval(opt$par), convergence = opt$convergence)
  }

  make_partition_list <- list(
    additive_terms = lapply(term_fits, function(fit) fit$make_partition_list)
  )
  knots <- list(additive_terms = lapply(term_fits, function(fit) fit$knots))
  penalties <- list(additive_terms = lapply(term_fits, function(fit) fit$penalties))

  return_list <- list(
    y = y,
    ytilde = ytilde,
    eta = eta,
    additive_terms = term_fits,
    additive_term_eta = term_eta,
    additive_deviance = dev_history,
    X = list(additive_terms = lapply(term_fits, function(fit) fit$X)),
    A = A,
    B = B,
    B_raw = B_raw,
    K = sapply(term_fits, `[[`, "K"),
    p = sapply(term_fits, `[[`, "p"),
    q = q_predictors,
    P = sum(sapply(term_fits, `[[`, "P")),
    N = N_obs,
    penalties = penalties,
    knots = knots,
    partition_codes = list(additive_terms =
                             lapply(term_fits, function(fit) fit$partition_codes)),
    knot_expand_function = knot_expand_list,
    predict = predict_function,
    plot = function(model_fit_in = NULL, ...){
      if(is.null(model_fit_in)) model_fit_in <- return_list
      .additive_plot_function(model_fit_in, ...)
    },
    assign_partition = function(x){
      lapply(term_fits, function(fit) fit$assign_partition(x))
    },
    family = fit_args$family,
    estimate_dispersion = fit_args$estimate_dispersion,
    unbias_dispersion = fit_args$unbias_dispersion,
    mean_y = 0,
    sd_y = 1,
    og_order = seq_len(N_obs),
    order_list = list(seq_len(N_obs)),
    constraint_values = list(additive_terms =
                               lapply(term_fits, function(fit) fit$constraint_values)),
    constraint_vectors = list(additive_terms =
                                lapply(term_fits, function(fit) fit$constraint_vectors)),
    make_partition_list = make_partition_list,
    expansion_scales = list(additive_terms =
                              lapply(term_fits, function(fit) fit$expansion_scales)),
    raw_expansion_names = list(additive_terms =
                                 lapply(term_fits, function(fit) fit$raw_expansion_names)),
    return_varcovmat = fit_args$return_varcovmat,
    parallel_cluster_supplied = !is.null(fit_args$cl),
    weights = weights,
    VhalfInv = fit_args$VhalfInv,
    Vhalf = fit_args$Vhalf,
    G = combined_system$G,
    Ghalf = combined_system$Ghalf,
    U = combined_system$U,
    additive_Lambda = combined_system$Lambda,
    additive_A_active = combined_system$A,
    qp_info = list(additive_joint = joint_qp_info,
                   additive_terms =
                     lapply(term_fits, function(fit) fit$qp_info)),
    quadprog_list = list(additive_terms =
                           lapply(term_fits, function(fit) fit$quadprog_list)),
    sigmasq_tilde = sigmasq_tilde,
    trace_XUGX = trace_XUGX,
    varcovmat = varcovmat,
    wald_univariate = wald_univariate,
    generate_posterior = generate_posterior,
    find_extremum = find_extremum,
    critical_value = critical_value,
    spline_groups = groups,
    predictors = predictors,
    og_cols = pred_colnames,
    additive_anchor_cols = anchor_cols,
    .fit_call_args = additive_fit_call_args
  )

  class(return_list) <- c("additive_lgspline", "lgspline")
  return_list
}
