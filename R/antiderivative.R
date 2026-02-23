#' Compute the Anti-Derivative or Definite Integral of a Fitted lgspline
#'
#' @description
#' Given a fitted \code{lgspline} object, this function computes either the
#' indefinite anti-derivative (evaluated pointwise over new or training predictors)
#' or the definite integral over a specified rectangular domain. Two strategies
#' are available: an analytical power-rule method that operates directly on the
#' piecewise polynomial structure, and a numerical Gauss-Legendre quadrature method
#' that evaluates the fitted surface on a tensor-product grid via \code{predict()}.
#'
#' For low-degree polynomials on moderate domains such as a cubic spline of
#' \code{sin(t)} on \eqn{[-\pi, \pi]}, both methods agree closely and the
#' analytical route is fast and does not require choice of tuning parameters.
#' For higher-degree fits on large predictor ranges, especially with quartics,
#' the analytical approach can become severely ill-conditioned. The power rule
#' amplifies high-order monomial contributions disproportionately (a quartic term
#' at \eqn{x = 80} is scaled by \eqn{80^4/5 \approx 8 \times 10^6}), so
#' coefficients that cancel in the fitted surface may add constructively in the
#' integrated form. The numerical method evaluates moderate function values directly
#' and is immune to this amplification, making it the default for definite integrals.
#' The analytical method remains the default for indefinite anti-derivatives, where
#' a pointwise output vector is needed rather than a scalar.
#'
#' @section Analytical power-rule:
#' The power rule is applied term-by-term to each column of the design matrix
#' partition-by-partition. For integration with respect to variable \eqn{j},
#' each monomial column transforms as follows: the intercept gains a factor of
#' \eqn{x_j}; a linear term in \eqn{x_j} becomes \eqn{x_j^2/2}; quadratic, cubic,
#' and quartic powers of \eqn{x_j} become \eqn{x_j^3/3}, \eqn{x_j^4/4}, and
#' \eqn{x_j^5/5}; interaction terms involving \eqn{x_j} follow accordingly (e.g.,
#' \eqn{x_j x_m \to x_j^2 x_m / 2}); and any column not involving \eqn{x_j}
#' is treated as a constant and multiplied through by \eqn{x_j}. When integrating
#' over multiple variables the transformations are applied sequentially, yielding
#' the iterated multiple integral.
#'
#' Continuity across partitions is enforced by propagating an additive offset at
#' each internal knot boundary, and the global constant is pinned via
#' \code{boundary_value} at the anchor point \code{boundary_predictors} (defaulting
#' to zero for all integration variables).
#'
#' @section Numerical Gauss-Legendre quadrature:
#' The integration domain is discretised into a tensor-product grid of
#' Gauss-Legendre quadrature nodes (default: 50 per dimension, tunable via
#' \code{n_quad}). Predicted values at each node come from the model's
#' \code{predict()} method, which correctly handles partition assignment and
#' piecewise polynomial evaluation. The integral is the weighted sum of those
#' predicted values multiplied by the Jacobian \eqn{(b_j - a_j)/2} for each
#' dimension. For smooth polynomials, 30-50 nodes per dimension is typically
#' sufficient; highly partitioned models (large \eqn{K}) may benefit from more.
#'
#' @section Integration scale:
#' By default (\code{link_scale = FALSE}), integration is performed on the
#' \emph{response} scale \eqn{\mu}, which requires the numerical method for
#' non-identity links (because the inverse link is nonlinear and the analytical
#' power rule does not apply). Setting \code{link_scale = TRUE} integrates on the
#' \emph{link-transformed} (linear predictor) scale \eqn{\eta = g(\mu) = \mathbf{X}\boldsymbol{\beta}}
#' instead, which allows the analytical method even for non-identity families.
#' For the identity link the two scales coincide and \code{link_scale} has no effect.
#' If \code{link_scale = FALSE}, \code{method = "analytical"}, and the link is
#' not the identity, an error is raised.
#'
#' @param model_fit A fitted \code{lgspline} object.
#' @param new_predictors Default: NULL. Numeric matrix or vector of predictor
#'   values at which to evaluate the anti-derivative. When NULL the training data
#'   stored in \code{model_fit} are used. Must be coercible to a matrix with the
#'   same number of columns as the predictor matrix used during fitting.
#' @param vars Default: NULL. Character or integer vector identifying the
#'   variable(s) to integrate with respect to. When NULL all numeric predictor
#'   variables are integrated simultaneously.
#' @param initial_values Default: NULL. Numeric vector of length equal to the
#'   total number of predictors, supplying fixed values for predictors not among
#'   the integration variables. When NULL, non-integration predictors default to
#'   zero at the boundary (analytical) or to the midpoint of their training range
#'   (numerical).
#' @param lower Default: NULL. Lower bound(s) for definite integration. Must be
#'   supplied together with \code{upper}.
#' @param upper Default: NULL. Upper bound(s) for definite integration.
#' @param boundary_value Default: 0. Desired value of the anti-derivative at the
#'   anchor point (analytical method only).
#' @param boundary_predictors Default: NULL. Numeric vector specifying the anchor
#'   point (analytical method only). Length must equal the number of predictors.
#' @param B_predict Default: NULL. Optional list of coefficient vectors, one per
#'   partition. When NULL the fitted coefficients from \code{model_fit} are used.
#' @param link_scale Default: FALSE. Logical; when FALSE integration is on the
#'   response scale \eqn{\mu}. When TRUE integration is on the link-transformed
#'   (linear predictor) scale \eqn{\eta}. For non-identity links,
#'   \code{link_scale = FALSE} with \code{method = "analytical"} raises an error
#'   since the power rule cannot account for the inverse link transformation.
#' @param method Default: NULL. One of \code{"analytical"}, \code{"numerical"},
#'   or NULL for automatic selection (numerical for definite integrals, analytical
#'   for indefinite anti-derivatives).
#' @param n_quad Default: 50. Number of Gauss-Legendre quadrature nodes per
#'   integration dimension (numerical method). Total evaluation points scale as
#'   \code{n_quad^length(vars)}.
#' @param cl Default: NULL. Optional parallel cluster with respect to partitions.
#'
#' @return A single numeric scalar for definite integrals; a named list with one
#'   vector per integration variable plus a \code{"total"} element when
#'   \code{vars = NULL} (integrate all); or a numeric vector when integrating
#'   a subset of variables.
#'
#' @examples
#' set.seed(1234)
#' t <- seq(-pi, pi, length.out = 1000)
#' y <- sin(t) + rnorm(length(t), 0, 0.01)
#' fit <- lgspline(t, y, K = 4, opt = FALSE)
#'
#' ## Indefinite anti-derivative at training points
#' F_hat <- antiderivative(fit, boundary_value = -1)$total
#' plot(t, F_hat, type = 'l', main = 'Anti-derivative of fitted sin(t)')
#' points(t, -cos(t), col = 'grey', cex = 0.5)
#'
#' ## Definite integral over [-pi, pi] (should be ~0 for sin)
#' integral_val <- antiderivative(fit, lower = -pi, upper = pi)
#'
#' ## Force analytical method (fine here, small domain, low degree)
#' integral_ana <- antiderivative(fit, lower = -pi, upper = pi,
#'                                method = "analytical")
#'
#' ## 2D volcano example where numerical method is essential
#' data(volcano)
#' vlong <- cbind(
#'   rep(seq_len(nrow(volcano)), ncol(volcano)),
#'   rep(seq_len(ncol(volcano)), each = nrow(volcano)),
#'   as.vector(volcano)
#' )
#' colnames(vlong) <- c("Length", "Width", "Height")
#' fit_v <- lgspline(vlong[, 1:2], vlong[, 3], K = 18,
#'                   include_quadratic_interactions = TRUE, opt = FALSE)
#'
#' ## Numerical quadrature (default for definite) near correct
#' antiderivative(fit_v, lower = c(1, 1), upper = c(87, 61))
#'
#' ## Analytical method, likely wildly wrong due to amplification
#' antiderivative(fit_v, lower = c(1, 1), upper = c(87, 61),
#'                method = "analytical")
#'
#' @export
antiderivative <- function(
    model_fit,
    new_predictors    = NULL,
    vars              = NULL,
    initial_values    = NULL,
    lower             = NULL,
    upper             = NULL,
    boundary_value    = 0,
    boundary_predictors = NULL,
    B_predict         = NULL,
    link_scale        = FALSE,
    method            = NULL,
    n_quad            = 50L,
    cl                = NULL
){

  if(is.null(B_predict)) B_predict <- model_fit$B

  ## Validate bounds: both or neither
  has_bounds <- !is.null(lower) && !is.null(upper)
  if(xor(is.null(lower), is.null(upper))){
    stop(
      '\n\t Both lower and upper must be supplied for a definite integral, ',
      'or both must be NULL for the indefinite anti-derivative.\n'
    )
  }

  ## Pull structural info from the fitted model
  K         <- model_fit$K
  nc        <- model_fit$p
  col_names <- model_fit$raw_expansion_names
  p1        <- model_fit$power1_cols
  p2        <- model_fit$power2_cols
  p3        <- model_fit$power3_cols
  p4        <- model_fit$power4_cols
  s1        <- model_fit$interaction_single_cols
  sq        <- model_fit$interaction_quad_cols
  tr        <- model_fit$triplet_cols
  ns_cols   <- model_fit$nonspline_cols
  es        <- model_fit$expansion_scales
  q_total   <- model_fit$q

  og_cols <- tryCatch(
    get('og_cols', envir = environment(model_fit$predict)),
    error = function(e) NULL
  )

  ## Resolve integration variable indices
  if(is.character(vars)){
    coef_names <- rownames(model_fit$B[[1]])
    if(is.null(coef_names)){
      stop('\n\tUnable to resolve character vars: coefficient names missing.\n')
    }
    base_names   <- coef_names[coef_names != 'intercept']
    is_main_effect <- !grepl('\\^', names(base_names)) &
      !grepl('x', names(base_names))
    predictor_names <- base_names[is_main_effect]
    if(length(predictor_names) < model_fit$q){
      stop('\n\tCould not reliably identify predictor names from coefficients.\n')
    }
    predictor_names <- predictor_names[seq_len(model_fit$q)]
    match_idx <- match(vars, predictor_names)
    if(any(is.na(match_idx))){
      stop(
        '\n\tCharacter vars not found among predictors.\n',
        '\tAvailable predictors: ', paste(predictor_names, collapse = ', '), '\n'
      )
    }
    vars <- as.integer(match_idx)
  }

  if(is.null(vars)){
    vars_idx      <- model_fit$numerics
    integrate_all <- TRUE
  } else {
    if(is.character(vars)){
      if(is.null(og_cols)){
        stop(
          '\n\tCharacter vars requires named predictor columns. Use integer ',
          'column indices instead.\n'
        )
      }
      vars_idx <- match(vars, og_cols)
      if(any(is.na(vars_idx)))
        stop('\n\tSome entries of vars not found in predictor column names.\n')
    } else {
      vars_idx <- as.integer(vars)
    }
    integrate_all <- FALSE
  }
  n_integ <- length(vars_idx)

  if(has_bounds){
    lower <- as.numeric(lower)
    upper <- as.numeric(upper)
    if(length(lower) == 1L) lower <- rep(lower, n_integ)
    if(length(upper) == 1L) upper <- rep(upper, n_integ)
    if(length(lower) != n_integ || length(upper) != n_integ){
      stop(
        '\n\t lower and upper must each have length equal to the number of ',
        'integration variables (', n_integ, ').\n'
      )
    }
  }

  ## Determine link name for scale checks
  link_name <- tryCatch(model_fit$family$link, error = function(e) 'identity')
  if(is.null(link_name)) link_name <- 'identity'

  ## Choose method with automatic default
  if(is.null(method)){
    method <- if(has_bounds) 'numerical' else 'analytical'
  }
  method <- match.arg(method, c('analytical', 'numerical'))

  if(method == 'numerical' && !has_bounds){
    stop(
      "\n\tmethod = 'numerical' is only available for definite integrals ",
      "(supply lower and upper). For indefinite anti-derivatives use ",
      "method = 'analytical'.\n"
    )
  }

  ## Enforce link_scale / method / link consistency
  if(!link_scale && method == 'analytical' && link_name != 'identity'){
    stop(
      "\n\tlink_scale = FALSE with method = 'analytical' is not supported for ",
      "non-identity link functions ('", link_name, "'). Either set ",
      "link_scale = TRUE to integrate on the linear predictor scale, or use ",
      "method = 'numerical' (the default for definite integrals) to integrate ",
      "on the response scale.\n"
    )
  }

  ## Resolve inverse-link function for response-scale numerical integration
  inv_link <- if(!link_scale && link_name != 'identity'){
    tryCatch(model_fit$family$linkinv, error = function(e) identity)
  } else {
    identity
  }

  ## ================================================================
  ##  NUMERICAL METHOD: Gauss-Legendre quadrature via predict()
  ## ================================================================

  if(method == 'numerical' && has_bounds){

    ## Golub-Welsch algorithm for Gauss-Legendre nodes and weights on [-1, 1]
    .gauss_legendre <- function(n){
      i   <- seq_len(n - 1L)
      b   <- i / sqrt(4 * i^2 - 1)
      J   <- matrix(0, n, n)
      for(k in seq_along(b)){
        J[k, k + 1L] <- b[k]
        J[k + 1L, k] <- b[k]
      }
      eig     <- eigen(J, symmetric = TRUE)
      nodes   <- eig$values
      weights <- 2 * eig$vectors[1, ]^2
      ord     <- order(nodes)
      list(nodes = nodes[ord], weights = weights[ord])
    }

    gl <- .gauss_legendre(as.integer(n_quad))

    ## Map nodes and weights from [-1, 1] to [a_j, b_j]
    mapped_nodes   <- vector('list', n_integ)
    mapped_weights <- vector('list', n_integ)
    for(vi in seq_len(n_integ)){
      hw <- (upper[vi] - lower[vi]) / 2
      mp <- (upper[vi] + lower[vi]) / 2
      mapped_nodes[[vi]]   <- hw * gl$nodes + mp
      mapped_weights[[vi]] <- gl$weights * hw
    }

    ## Tensor-product grid and product weights
    grid           <- expand.grid(mapped_nodes)
    weight_grid    <- expand.grid(mapped_weights)
    product_weights <- apply(weight_grid, 1, prod)
    n_grid         <- nrow(grid)

    ## Fixed values for non-integration predictors
    if(is.null(initial_values)){
      train_raw <- tryCatch(model_fit$raw_X, error = function(e) NULL)
      if(!is.null(train_raw)){
        base_vec <- colMeans(matrix(
          c(apply(train_raw, 2, min), apply(train_raw, 2, max)),
          nrow = 2, byrow = TRUE
        ))
      } else {
        base_vec <- rep(0, q_total)
      }
    } else {
      base_vec <- as.numeric(initial_values)
      if(length(base_vec) != q_total){
        stop(
          '\n\tinitial_values must have length equal to the number of ',
          'predictor variables (', q_total, ').\n'
        )
      }
    }

    pred_mat <- matrix(rep(base_vec, each = n_grid), nrow = n_grid)
    for(vi in seq_len(n_integ)){
      pred_mat[, vars_idx[vi]] <- grid[, vi]
    }

    ## Evaluate on the grid with correct scale
    predicted <- model_fit$predict(new_predictors = pred_mat)

    if(link_scale){
      ## Integrate on link (linear predictor) scale
      if(link_name != 'identity'){
        # convert response to link scale if non-identity
        predicted <- tryCatch(
          model_fit$family$linkfun(predicted),
          error = function(e) stop("\nCannot convert to link scale\n")
        )
      }
    } else {
      ## Integrate on response scale
      # predicted is already on response scale from predict(), do nothing
      predicted <- predicted
    }
    total_val  <- sum(predicted * product_weights)

    return(total_val)
  }

  ## ================================================================
  ##  ANALYTICAL METHOD: power-rule anti-derivative
  ## ================================================================

  ## Build design matrix partitions
  if(is.null(new_predictors)){
    X_parts      <- model_fit$X
    order_list   <- model_fit$order_list
    N            <- model_fit$N
    use_training <- TRUE
    only_1       <- FALSE
  } else {
    new_predictors <- try(as.matrix(new_predictors), silent = TRUE)
    if(inherits(new_predictors, 'try-error'))
      stop("\n\t'new_predictors' must be coercible to a matrix.\n")
    only_1 <- (nrow(new_predictors) == 1L)
    if(only_1) new_predictors <- rbind(new_predictors, new_predictors)
    prep       <- model_fit$predict(new_predictors = new_predictors,
                                    expansions_only = TRUE)
    X_parts    <- prep$expansions
    order_list <- model_fit$knot_expand_function(
      prep$partition_codes,
      prep$partition_bounds,
      nrow(new_predictors),
      as.matrix(seq_len(nrow(new_predictors))),
      K
    )
    N            <- nrow(new_predictors)
    use_training <- FALSE
  }

  keep <- which(sapply(X_parts, nrow) > 0L)

  ## Extract raw x_j from design matrix partition (unstandardized linear column)
  .get_xj <- function(X_k, j){
    col_nm_j <- paste0('_', j, '_')
    all_lin  <- c(p1, ns_cols)
    ci_list  <- all_lin[col_names[all_lin] == col_nm_j]
    if(length(ci_list) == 0L) return(rep(0, nrow(X_k)))
    X_k[, ci_list[1L]]
  }

  ## Power-rule transformation for a single integration variable j
  .antideriv_one_var <- function(F_k, j, x_j_col){

    involves_j <- function(ci){
      grepl(paste0('_', j, '_'), col_names[ci], fixed = TRUE)
    }

    F_k[, 1L] <- F_k[, 1L] * x_j_col

    for(ci in unique(c(p1, ns_cols))){
      if(involves_j(ci)){
        F_k[, ci] <- F_k[, ci] * x_j_col / 2
      } else {
        F_k[, ci] <- F_k[, ci] * x_j_col
      }
    }

    for(ci in p2){
      if(involves_j(ci)){
        F_k[, ci] <- F_k[, ci] * x_j_col / 3
      } else {
        F_k[, ci] <- F_k[, ci] * x_j_col
      }
    }

    for(ci in p3){
      if(involves_j(ci)){
        F_k[, ci] <- F_k[, ci] * x_j_col / 4
      } else {
        F_k[, ci] <- F_k[, ci] * x_j_col
      }
    }

    for(ci in p4){
      if(involves_j(ci)){
        F_k[, ci] <- F_k[, ci] * x_j_col / 5
      } else {
        F_k[, ci] <- F_k[, ci] * x_j_col
      }
    }

    for(ci in s1[!(s1 %in% c(tr, sq))]){
      if(involves_j(ci)){
        F_k[, ci] <- F_k[, ci] * x_j_col / 2
      } else {
        F_k[, ci] <- F_k[, ci] * x_j_col
      }
    }

    for(ci in sq){
      if(involves_j(ci)){
        if(grepl(paste0('_', j, '_\\^2$'), col_names[ci])){
          F_k[, ci] <- F_k[, ci] * x_j_col / 3
        } else {
          F_k[, ci] <- F_k[, ci] * x_j_col / 2
        }
      } else {
        F_k[, ci] <- F_k[, ci] * x_j_col
      }
    }

    for(ci in tr){
      if(involves_j(ci)){
        F_k[, ci] <- F_k[, ci] * x_j_col / 2
      } else {
        F_k[, ci] <- F_k[, ci] * x_j_col
      }
    }

    F_k
  }

  ## Iterated anti-derivative: apply power rule sequentially over j_vec
  .antideriv_seq <- function(X_k, j_vec){
    F_k <- X_k
    for(j in j_vec){
      F_k <- .antideriv_one_var(F_k, j, .get_xj(X_k, j))
    }
    F_k
  }

  ## Boundary condition and per-partition continuity offsets
  if(!is.null(boundary_predictors) && length(boundary_predictors) != q_total){
    stop(
      '\n\tboundary_predictors must have length equal to the number of ',
      'predictor variables (', q_total, ').\n'
    )
  }

  {
    bnd_pred <- if(!is.null(boundary_predictors)){
      as.numeric(boundary_predictors)
    } else {
      rep(0, q_total)
    }
    for(j_anc in vars_idx) bnd_pred[j_anc] <- 0

    bnd_mat  <- rbind(bnd_pred, bnd_pred)
    prep_bnd <- model_fit$predict(new_predictors = bnd_mat,
                                  expansions_only = TRUE)
    keep_bnd <- which(sapply(prep_bnd$expansions, nrow) > 0L)
    k_bnd    <- keep_bnd[1L]
    X_anc    <- prep_bnd$expansions[[k_bnd]][1L, , drop = FALSE]

    train_X    <- model_fit$X
    train_keep <- which(sapply(train_X, nrow) > 0L)

    .offset_delta_train <- function(ki_tr, cur_off_tr){
      k_tr      <- train_keep[ki_tr]
      k_next_tr <- train_keep[ki_tr + 1L]
      x_bnd_tr  <- train_X[[k_tr]][nrow(train_X[[k_tr]]), , drop = FALSE]
      Fk_tr <- as.numeric(
        .antideriv_seq(x_bnd_tr, vars_idx) %*% B_predict[[k_tr]]
      ) + cur_off_tr
      Fn_tr <- as.numeric(
        .antideriv_seq(x_bnd_tr, vars_idx) %*% B_predict[[k_next_tr]]
      )
      Fk_tr - Fn_tr
    }

    offsets_tmp <- numeric(max(train_keep))
    cur_off_tmp <- 0
    for(ki_tr in seq_along(train_keep)){
      offsets_tmp[train_keep[ki_tr]] <- cur_off_tmp
      if(ki_tr < length(train_keep))
        cur_off_tmp <- .offset_delta_train(ki_tr, cur_off_tmp)
    }

    raw_at_anchor <- as.numeric(
      .antideriv_seq(X_anc, vars_idx) %*% B_predict[[k_bnd]]
    ) + offsets_tmp[k_bnd]

    offset_0 <- boundary_value - raw_at_anchor
  }

  ## Definite integral via analytical inclusion-exclusion (2^n corners)
  if(has_bounds){
    base_vec <- if(!is.null(initial_values)) as.numeric(initial_values) else rep(0, q_total)

    .eval_F_at <- function(pred_vec){
      pmat <- rbind(pred_vec, pred_vec)
      prep <- model_fit$predict(new_predictors = pmat, expansions_only = TRUE)
      kp   <- which(sapply(prep$expansions, nrow) > 0L)[1L]
      Xp   <- prep$expansions[[kp]][1L, , drop = FALSE]
      raw  <- as.numeric(.antideriv_seq(Xp, vars_idx) %*% B_predict[[kp]])
      off  <- offset_0 + if(kp <= length(offsets_tmp)) offsets_tmp[kp] else 0
      raw + off
    }

    n_corners <- 2L^n_integ
    total_val <- 0
    for(corner_idx in seq_len(n_corners)){
      bits    <- as.integer(intToBits(corner_idx - 1L))[seq_len(n_integ)]
      n_lower <- sum(bits == 0L)
      sgn     <- (-1L)^n_lower
      pred_vec <- base_vec
      for(vi in seq_len(n_integ)){
        pred_vec[vars_idx[vi]] <- if(bits[vi] == 1L) upper[vi] else lower[vi]
      }
      total_val <- total_val + sgn * .eval_F_at(pred_vec)
    }
    return(total_val)
  }

  ## Main loop: indefinite anti-derivative over all partitions
  out_vec <- numeric(N)

  ## Export parallel components
  if(!is.null(cl) && inherits(cl, 'cluster')){
    parallel::clusterExport(
      cl,
      varlist = c('X_parts', 'keep', 'B_predict', 'vars_idx',
                  'col_names', 'p1', 'p2', 'p3', 'p4',
                  's1', 'sq', 'tr', 'ns_cols', 'es',
                  'offsets_tmp', 'offset_0',
                  '.get_xj', '.antideriv_one_var', '.antideriv_seq'),
      envir = environment()
    )
    raw_vals <- parallel::parLapply(cl, seq_along(keep), function(ki){
      k <- keep[ki]
      as.numeric(.antideriv_seq(X_parts[[k]], vars_idx) %*% B_predict[[k]])
    })
    for(ki in seq_along(keep)){
      k           <- keep[ki]
      part_offset <- offset_0 + if(k <= length(offsets_tmp)) offsets_tmp[k] else 0
      out_vec[unlist(order_list[[k]])] <- raw_vals[[ki]] + part_offset
    }
  } else {
    for(ki in seq_along(keep)){
      k           <- keep[ki]
      part_offset <- offset_0 + if(k <= length(offsets_tmp)) offsets_tmp[k] else 0
      out_vec[unlist(order_list[[k]])] <-
        as.numeric(.antideriv_seq(X_parts[[k]], vars_idx) %*% B_predict[[k]]) +
        part_offset
    }
  }

  if(!use_training && only_1) out_vec <- out_vec[1L]

  ## If initial_values supplied with training data, build synthetic predictor matrix
  if(!is.null(initial_values) && use_training){

    init_vec <- as.numeric(initial_values)
    if(length(init_vec) != q_total){
      stop(
        '\n\tinitial_values must have length equal to the number of ',
        'predictor variables (', q_total, ').\n'
      )
    }

    N_train      <- model_fit$N
    out_vec_init <- numeric(N_train)

    synth_preds <- matrix(rep(init_vec, each = N_train), nrow = N_train, byrow = FALSE)
    for(k in keep){
      rows_k <- unlist(order_list[[k]])
      X_k    <- X_parts[[k]]
      for(j in vars_idx){
        synth_preds[rows_k, j] <- .get_xj(X_k, j)
      }
    }

    prep        <- model_fit$predict(new_predictors = synth_preds, expansions_only = TRUE)
    X_parts_syn <- prep$expansions
    order_syn   <- model_fit$knot_expand_function(
      prep$partition_codes, prep$partition_bounds,
      N_train, as.matrix(seq_len(N_train)), K
    )
    keep_syn <- which(sapply(X_parts_syn, nrow) > 0L)

    for(ki in seq_along(keep_syn)){
      k   <- keep_syn[ki]
      F_k <- .antideriv_seq(X_parts_syn[[k]], vars_idx)
      part_offset <- offset_0 + if(k <= length(offsets_tmp)) offsets_tmp[k] else 0
      out_vec_init[unlist(order_syn[[k]])] <-
        as.numeric(F_k %*% B_predict[[k]]) + part_offset
    }

    return(out_vec_init)
  }

  ## integrate_all: return marginal single-variable anti-derivatives + total
  if(integrate_all){
    var_labels <- if(!is.null(og_cols)){
      og_cols[vars_idx]
    } else {
      paste0('predictor_', vars_idx)
    }

    result_list        <- vector('list', n_integ + 1L)
    names(result_list) <- c(var_labels, 'total')

    for(vi in seq_len(n_integ)){
      j_v    <- vars_idx[vi]
      sv_out <- numeric(N)

      sv_offsets_tmp <- numeric(max(train_keep))
      sv_cur <- 0
      for(ki_tr in seq_along(train_keep)){
        k_tr <- train_keep[ki_tr]
        sv_offsets_tmp[k_tr] <- sv_cur
        if(ki_tr < length(train_keep)){
          k_next_tr <- train_keep[ki_tr + 1L]
          x_bnd_tr  <- model_fit$X[[k_tr]][nrow(model_fit$X[[k_tr]]), , drop = FALSE]
          Fk_sv <- as.numeric(
            .antideriv_one_var(x_bnd_tr, j_v, .get_xj(x_bnd_tr, j_v)) %*%
              B_predict[[k_tr]]
          ) + sv_cur
          Fn_sv <- as.numeric(
            .antideriv_one_var(x_bnd_tr, j_v, .get_xj(x_bnd_tr, j_v)) %*%
              B_predict[[k_next_tr]]
          )
          sv_cur <- Fk_sv - Fn_sv
        }
      }

      raw_sv_anchor <- as.numeric(
        .antideriv_one_var(X_anc, j_v, .get_xj(X_anc, j_v)) %*%
          B_predict[[k_bnd]]
      ) + sv_offsets_tmp[k_bnd]
      sv_offset_0 <- boundary_value - raw_sv_anchor

      for(ki in seq_along(keep)){
        k    <- keep[ki]
        X_k  <- X_parts[[k]]
        F_1v <- .antideriv_one_var(X_k, j_v, .get_xj(X_k, j_v))
        part_sv_offset <- sv_offset_0 + if(k <= length(sv_offsets_tmp)) sv_offsets_tmp[k] else 0
        sv_out[unlist(order_list[[k]])] <-
          as.numeric(F_1v %*% B_predict[[k]]) + part_sv_offset
      }
      if(!use_training && only_1) sv_out <- sv_out[1L]
      result_list[[vi]] <- sv_out
    }
    result_list[['total']] <- out_vec
    return(result_list)
  }

  out_vec
}
