#' Lagrangian Multiplier Smoothing Splines
#'
#' @docType _PACKAGE
#' @name lgspline-package
#' @aliases lgspline
#'
#' @import Rcpp RcppArmadillo parallel stats
#' @importFrom graphics plot points legend
#' @importFrom FNN get.knnx
#' @importFrom quadprog solve.QP
#' @importFrom RColorBrewer brewer.pal
#' @importFrom plotly plot_ly layout
#'
#' @keywords smoothing regression parametric constrained lagrangian
NULL

#' Fit Lagrangian Multiplier Smoothing Splines
#'
#' @description
#' Fits a Lagrangian multiplier smoothing spline by combining piecewise polynomial functions with explicit
#' smoothness constraints enforced solely using the method of Lagrangian multipliers.
#' The fitted function is penalized by the squared integrated second derivative.
#' Supports GLM families, quadratic programming, parallel processing, and provides comprehensive
#' tools for frequentist and Bayesian inference, optimization, and visualization.
#'
#' @param predictors Numeric matrix or data frame of predictor variables. Supports direct matrix input or formula interface when used with `data` argument. Must contain numeric predictors, with categorical variables pre-converted to numeric indicators.
#' @param y Numeric response variable vector representing the target outcome to be modeled.
#' @param formula Optional statistical formula for model specification, serving as an alternative to direct matrix input. Supports standard R formula syntax with special `spl()` function for defining spline terms.
#' @param response Alternative name for response variable, providing compatibility with different naming conventions. Takes precedence only if `y` is not supplied.
#' @param standardize_response Logical indicator controlling whether the response variable should be centered and scaled before model fitting. When TRUE, improves numerical stability and comparability of coefficients across different scales.
#' @param standardize_predictors_for_knots Logical flag determining whether predictor variables should be standardized before knot placement. Ensures consistent knot selection across different predictor scales.
#' @param standardize_expansions_for_fitting Logical switch to standardize polynomial basis expansions during model fitting. Provides computational stability during penalty tuning without affecting statistical inference, as design matrices are systematically backtransformed.
#' @param family Generalized linear model (GLM) distribution family specifying the error distribution and link function for model fitting. Defaults to Gaussian distribution with identity link. Supports custom family specifications, including user-defined link functions and optional custom tuning loss criteria.
#' @param glm_weight_function User-specified function for calculating generalized linear model weights used in computing the G matrix. Allows flexible weight specification based on fitted values, response, observation ordering, and family characteristics.
#' @param shur_correction_function Advanced function for computing Schur complements to properly account for uncertainty in dispersion *or other) parameter estimation, preventing potential underestimation of parameter uncertainty.
#' @param need_dispersion_for_estimation Logical indicator specifying whether a dispersion parameter is required for coefficient estimation, enabling more sophisticated modeling approaches.
#' @param dispersion_function Custom function for estimating the dispersion parameter, allowing flexible approaches to variance component estimation across different model specifications.
#' @param K Integer specifying the number of knot locations for spline partitions. When NULL, knots are automatically selected based on data characteristics to optimize model flexibility.
#' @param custom_knots Optional matrix providing user-specified knot locations in 1-D, allowing precise control over spline partitioning that overrides automatic knot selection.
#' @param cluster_on_indicators Logical flag determining whether indicator variables should be used for clustering knot locations, providing alternative clustering strategies.
#' @param make_partition_list Optional list allowing direct specification of custom partition assignments for advanced modeling scenarios.
#' @param previously_tuned_penalties Optional list of pre-computed penalty components from previous model fits, enabling efficient model refinement and cross-validation.
#' @param smoothing_spline_penalty Optional custom smoothing spline penalty matrix for fine-tuned complexity control.
#' @param opt Logical switch controlling whether model penalties should be automatically optimized via cross-validation.
#' @param use_custom_bfgs Logical indicator selecting between a custom damped-BFGS optimization method with analytical gradients or base R's BFGS implementation with finite-difference gradient approximation.
#' @param delta Numeric pseudocount used for stabilizing optimization in non-identity link function scenarios.
#' @param tol Numeric convergence tolerance controlling the precision of optimization algorithms.
#' @param log_initial_wiggle Numeric vector of initial grid points for wiggle penalty optimization, specified on the natural logarithmic scale.
#' @param log_initial_flat Numeric vector of initial grid points for ridge penalty optimization, specified on the natural logarithmic scale.
#' @param wiggle_penalty Numeric penalty controlling the integrated squared second derivative, governing function smoothness.
#' @param flat_ridge_penalty Numeric flat ridge penalty for additional regularization.
#' @param unique_penalty_per_partition Logical flag allowing unique complexity penalties for each spline partition.
#' @param unique_penalty_per_predictor Logical flag permitting unique complexity penalties for individual predictors.
#' @param penalty_ridge Numeric "meta-penalty" applied to predictor and partition penalties during tuning, encouraging penalty values close to unity.
#' @param predictor_penalties Optional list of custom penalties specified per predictor.
#' @param partition_penalties Optional list of custom penalties specified per partition.
#' @param include_quadratic_terms Logical switch to include squared predictor terms in basis expansion.
#' @param include_cubic_terms Logical switch to include cubic predictor terms in basis expansion.
#' @param include_quartic_terms Logical switch to include quartic predictor terms in basis expansion.
#' @param include_2way_interactions Logical switch to include two-way interactions between predictors.
#' @param include_3way_interactions Logical switch to include three-way interactions between predictors.
#' @param include_quadratic_interactions Logical switch to include quadratic interaction terms.
#' @param just_linear_with_interactions Integer vector specifying columns to retain linear terms with interactions.
#' @param just_linear_without_interactions Integer vector specifying columns to retain only linear terms without interactions.
#' @param exclude_interactions_for Integer vector indicating columns to exclude from all interaction terms.
#' @param exclude_these_expansions Character vector specifying basis expansions to be excluded from the model.
#' @param custom_basis_fxn Optional user-defined function for generating custom basis expansions, following specific column naming conventions.
#' @param include_constrain_fitted Logical switch to constrain fitted values at knot points.
#' @param include_constrain_first_deriv Logical switch to constrain first derivatives at knot points.
#' @param include_constrain_second_deriv Logical switch to constrain second derivatives at knot points.
#' @param include_constrain_interactions Logical switch to constrain interaction terms at knot points.
#' @param cl Parallel processing cluster object for distributed computation.
#' @param chunk_size Integer specifying custom fixed chunk size for parallel processing.
#' @param parallel_eigen Logical flag to enable parallel processing for eigenvalue decomposition computations.
#' @param parallel_trace Logical flag to enable parallel processing for trace computation.
#' @param parallel_aga Logical flag to enable parallel processing for specific matrix operations.
#' @param parallel_matmult Logical flag to enable parallel processing for block-diagonal matrix multiplication.
#' @param parallel_unconstrained Logical flag to enable parallel processing for unconstrained maximum likelihood estimation.
#' @param parallel_find_neighbors Logical flag to enable parallel processing for neighbor identification.
#' @param parallel_penalty Logical flag to enable parallel processing for penalty matrix construction.
#' @param parallel_make_constraint Logical flag to enable parallel processing for constraint matrix generation.
#' @param unconstrained_fit_fxn Custom function for fitting unconstrained models per partition.
#' @param keep_weighted_Lambda Logical flag to retain generalized linear model weights in constraints using Tikhonov parameterization.
#' @param iterate_tune Logical switch to use iterative optimization during penalty tuning.
#' @param iterate_final_fit Logical switch to use iterative optimization for final model fitting.
#' @param qp_Amat Constraint matrix for quadratic programming formulation.
#' @param qp_bvec Constraint vector for quadratic programming formulation.
#' @param qp_meq Number of equality constraints in quadratic programming setup.
#' @param qp_monotonic_increase Logical flag to constrain the function to be monotonically increasing.
#' @param qp_monotonic_decrease Logical flag to constrain the function to be monotonically decreasing.
#' @param qp_range_upper Numeric upper bound for constrained fitted values.
#' @param qp_range_lower Numeric lower bound for constrained fitted values.
#' @param qp_Amat_fxn Custom function for generating constraint matrix in quadratic programming.
#' @param qp_bvec_fxn Custom function for generating constraint vector in quadratic programming.
#' @param qp_meq_fxn Custom function for determining equality constraints in quadratic programming.
#' @param constraint_value_vectors Matrix of constraint values for sum constraints.
#' @param constraint_vectors Matrix of vectors for sum constraints.
#' @param return_G Logical switch to return the unscaled variance-covariance matrix.
#' @param return_Ghalf Logical switch to return the square root of the variance-covariance matrix.
#' @param return_U Logical switch to return the constraint matrix.
#' @param estimate_dispersion Logical flag to estimate the dispersion parameter.
#' @param return_varcovmat Logical switch to return the scaled variance-covariance matrix.
#' @param custom_penalty_mat Optional custom penalty matrix for individual partitions.
#' @param cluster_args Named vector of arguments controlling clustering procedures.
#' @param dummy_dividor Small numeric constant to prevent division by zero in computational routines.
#' @param dummy_adder Small numeric constant to prevent division by zero in computational routines.
#' @param verbose Logical flag to print general progress messages during model fitting.
#' @param verbose_tune Logical flag to print detailed progress messages during penalty tuning.
#' @param expansions_only Logical switch to return only basis expansions without full model fitting, useful for advanced modeling tasks.
#' @param observation_weights Numeric vector of observation-specific weights for generalized least squares estimation.
#' @param do_not_cluster_on_these Vector specifying predictor columns to exclude from clustering procedures.
#' @param neighbor_tolerance Numeric tolerance for determining neighboring partitions in multivariate spline contexts.
#' @param null_constraint Alternative parameterization of constraint values.
#' @param critical_value Numeric value used for constructing confidence intervals.
#' @param data Optional data frame providing context for formula-based model specification.
#' @param weights Alternative name for observation weights, maintained for interface compatibility.
#' @param no_intercept Logical flag to remove intercept, constraining it to a fixed value.
#' @param VhalfInv Matrix representing custom response covariance structure for advanced modeling.
#' @param include_warnings Logical switch to control display of warning messages during model fitting.
#' @param ... Additional arguments passed to the unconstrained model fitting function.
#'
#' @return A list containing various model components including fitted values,
#' coefficients, penalties, prediction and plotting functions, and inference tools.
#'
#' @examples
#' \dontrun{
#' ## ## ## 1-D Example
#' ## 1D Data generating function
#' set.seed(1234)
#' x <- seq(-9, 9, length.out = 1000)
#' slinky <- function(x) {
#'   (50 * cos(x * 2) +-2 * x ^ 2 + (0.25 * x) ^ 4 + 80)
#'
#' }
#' coil <- function(x) {
#'   (100 * cos(x * 2) +-1.5 * x ^ 2 + (0.1 * x) ^ 4 + (0.05 * x ^ 3) + (-0.01 *
#'                                                                         x ^ 5) +
#'      (0.00002*x^6) -(0.000001*x^7) + 100)
#' }
#' exponential_log <- function(x) {
#'   unlist(c(sapply(x, function(xx) {
#'     if (xx <= 1) {
#'       100 * (exp(xx) - exp(1))
#'
#'     } else {
#'       100 * (log(xx))
#'
#'     }
#'
#'   })))
#'
#' }
#' scaled_abs_gamma <- function(x) {
#'   2*sqrt(gamma(abs(x)))
#' }
#'
#' ## Composite function
#' fxn <- function(x)(slinky(x) +
#'                      coil(x) +
#'                      exponential_log(x) +
#'                      scaled_abs_gamma(x))
#'
#' ## Bind together with random noise
#' dat <- cbind(x, fxn(x) + rnorm(length(x), 0, 50))
#' colnames(dat) <- c('x', 'y')
#' x <- dat[,'x']
#' y <- dat[,'y']
#'
#'
#' ## Fit Model
#' model_fit <- lgspline(x, y)
#' model_fit <- lgspline(y ~ ., as.data.frame(dat))
#'
#' ## Basic Functionality
#' predict(model_fit)
#' leave_one_out(model_fit)
#' plot(model_fit)
#' points(dat, cex = 0.35)
#' coef(model_fit)
#' summary(model_fit)
#' generate_posterior(model_fit)
#' find_extremum(model_fit)
#'
#' ## Incorporate range constraints, custom knots, keep penalization identical
#' # across partitions
#' model_fit <- lgspline(y ~ spl(x),
#'                       standardize_response = FALSE,
#'                       unique_penalty_per_partition = FALSE,
#'                       custom_knots = cbind(c(-2, -1, 0, 1, 2)),
#'                       data = data.frame(x = x, y = y)[order(x),],
#'                       qp_range_lower = -250,
#'                       qp_range_upper = 400)
#'
#' ## Plotting the constraints and knots
#' plot(model_fit,
#'      custom_title = 'Fitted Function Constrained to Lie Between (-250, 400)',
#'      cex.main = 0.75)
#' # knot locations
#' abline(v = model_fit$knots)
#' # lower bound from quadratic program
#' abline(h = -250, lty = 2)
#' # upper bound from quadratic program
#' abline(h = 400, lty = 2)
#' # observed data
#' points(x, y, cex = 0.24)
#'
#' ## ## ## ## Fit to volcano dataset
#' ## Prep
#' set.seed(1234)
#' data('volcano')
#' volcano_long <-
#'   Reduce('rbind', lapply(1:nrow(volcano), function(i){
#'     t(sapply(1:ncol(volcano), function(j){
#'       c(i, j, volcano[i,j])
#'     }))
#'   }))
#' colnames(volcano_long) <- c('Length', 'Width', 'Height')
#'
#' ## Fit, with 50 partitions
#' model_fit <- lgspline(volcano_long[,c(1, 2)],
#'                       volcano_long[,3],
#'                       include_quartic_terms = TRUE,
#'                       K = 49,
#'                       opt = FALSE,
#'                       return_U = FALSE,
#'                       return_varcov = FALSE,
#'                       estimate_variance = TRUE,
#'                       return_G = FALSE,
#'                       include_constrain_second_deriv = FALSE,
#'                       unique_penalty_per_predictor = FALSE,
#'                       unique_penalty_per_partition = FALSE,
#'                       wiggle_penalty = 1e-8,
#'                       flat_ridge_penalty = 1e-2,
#'                       return_Ghalf = FALSE)
#'
#' ## Predictions on new data with interactive visual + formulas
#' new_input <- expand.grid(seq(min(volcano_long[,1]),
#'                              max(volcano_long[,1]),
#'                              length.out = 250),
#'                          seq(min(volcano_long[,2]),
#'                              max(volcano_long[,2]),
#'                              length.out = 250))
#' model_fit$plot(new_predictors = new_input,
#'                show_formulas = TRUE,
#'                custom_response_lab = "Height",
#'                custom_title = 'Volcano 3-D Map',
#'                digits = 2)
#'
#' ## ## ## ## Trees example (Advanced techniques)
#' ## Custom l1-regularization in lgspline
#' data('trees')
#' set.seed(1234)
#'
#' ## L1-regularization constraint function on standardized coefficients
#' # Bound all coefficients to be less than a certain value (l1_bound) in absolute
#' # magnitude such that | B^{(j)}_k | < lambda for all j = 1....P
#' l1_constraint_matrix <- function(p, K) {
#'   ## Total number of coefficients
#'   P <- p * (K + 1)
#'
#'   ## Create diagonal matrices for L1 constraint
#'   # First matrix: lamdba > -bound
#'   # Second matrix: -lambda > -bound
#'   first_diag <- diag(P)
#'   second_diag <- -diag(P)
#'
#'   ## Combine matrices
#'   l1_Amat <- cbind(first_diag, second_diag)
#'
#'   return(l1_Amat)
#' }
#'
#' ## By default, bounds absolute value of coefficients to be < 10
#' l1_bound_vector <- function(qp_Amat,
#'                             scales,
#'                             l1_bound) {
#'
#'   ## Combine matrices
#'   l1_bvec <- rep(-l1_bound, ncol(qp_Amat)) * c(1, scales)
#'
#'   return(l1_bvec)
#' }
#'
#' ## Fit model, using predictor-response formulation, assuming
#' # Gamma-distributed response, and custom quadratic-programming constraints
#' model_fit <- lgspline(
#'   Volume ~ spl(Girth) + Height*Girth,
#'   data = with(trees, cbind(Girth, Height, Volume)),
#'   include_quartic_terms = TRUE,
#'   family = Gamma(link = 'log'),
#'   K = 1,
#'   qp_Amat_fxn = function(N, p, K, X_block, colnm, scales) {
#'     mat <- l1_constraint_matrix(p, K)
#'     mat
#'   },
#'   qp_bvec_fxn = function(qp_Amat, N, p, K, X_block, colnm, scales) {
#'     vec <- l1_bound_vector(qp_Amat, scales, 10)
#'     vec
#'   },
#'   qp_meq_fxn = function(qp_Amat, N, p, K, X_block, colnm, scales) 0
#' )
#'
#' ## Notice, interaction effect are constant across partitions as is
#' # effect of Height alone
#' print(summary(model_fit))
#'
#' ## Penalties
#' print(model_fit$penalties)
#'
#' ## Plot results
#' plot(model_fit, custom_predictor_lab1 = 'Girth',
#'      custom_predictor_lab2 = 'Height',
#'      custom_response_lab = 'Volume',
#'      custom_title = 'Girth and Height Predicting Volume of Trees',
#'      show_formulas = TRUE)
#'
#' ## Verify magnitude of unstandardized coefficients does not exceed l1 bound (10)
#' print(max(abs(unlist(model_fit$B))))
#'
#' ## Thompson-sampling step
#' # Draw from posterior
#' postdraw <- generate_posterior(model_fit,
#'                                draw_dispersion = FALSE)
#' # Stochastically optimize the fitted function for "t",
#' # the points to evaluate next
#' find_extremum(
#'   model_fit,
#'   B_predict = postdraw$post_draw_coefficients,
#'   stochastic = TRUE
#' )
#'
#' ## Find coordinates where prediction is closest to median
#' find_extremum(
#'   model_fit,
#'   minimize = TRUE,
#'   custom_objective_function = function(mu, sigma, ybest){
#'     0.5*(mu - median(trees$Volume))^2
#'   },
#'   custom_objective_gradient = function(mu, sigma, ybest, d_mu){
#'     (mu - median(trees$Volume)) * d_mu
#'   }
#' )
#' print(median(trees$Volume))
#'
#'
#' ## ## ## ## Compare inference to survreg for special case of K = 0,
#' # only linear predictors, no penalties
#'
#' ## The concept here is that when using models where dispersion is
#' # being estimated and is required for estimating beta coefficients,
#' # we use a shur complement correction function to adjust our
#' # variance-covariance matrix for inference
#' require(survival)
#' df <- data.frame(na.omit(pbc[,c('time',
#'                                 'trt','stage','hepato','bili','status')]))
#'
#' ## Adjusted shur-correction function to match survreg,
#' # with N/(N-P) upscaling component
#' weibull_shur_correction2 <- function(X,
#'                                      y,
#'                                      B,
#'                                      dispersion,
#'                                      order_list,
#'                                      K,
#'                                      family,
#'                                      observation_weights,
#'                                      status){
#'   lapply(1:(K+1), function(k){
#'     if(nrow(X[[k]]) < 1){
#'       return(0)
#'     } else {
#'       mu <- family$linkinv(c(X[[k]] %*% B[[k]]))
#'       s <- status[order_list[[k]]]
#'       obs <- y[[k]]
#'       z <- (log(obs) - log(mu))/sqrt(dispersion)
#'       exp_z <- exp(z)
#'       zexp_z <- z*exp_z
#'       weights <- c(observation_weights[[k]])
#'
#'       ## Correction via Shur complement
#'       # I = ( I_bb I_bs^{T} )
#'       #     ( I_bs I_ss     )
#'       # for b = beta, s = dispersion (scale)
#'       #I_bb <- t(X[[k]]) %*% cbind(weights * exp_z * X[[k]])
#'       I_bs <- t(X[[k]]) %*% cbind(weights * zexp_z * sqrt(dispersion))
#'       I_ss <- -sum(
#'         weights * (
#'           (s + 2*s*z + zexp_z + exp_z * z^2)
#'         )
#'       ) * (nrow(X[[k]])- ncol(X[[k]]))/
#'         nrow(X[[k]]) # = N/(N-P), included that isn't in native
#'       compl <- I_bs %*% matrix(-1/I_ss) %*% t(I_bs)
#'       # Shur complement correction to pass on to compute_G_eigen()
#'       return(compl)
#'     }
#'   })
#' }
#'
#' ## Weibull AFT using lgspline
#' model_fit <- lgspline(time ~ trt + stage + hepato + bili,
#'                       df,
#'                       family = weibull_family(),
#'                       need_dispersion_for_estimation = TRUE,
#'                       dispersion_function = weibull_dispersion_function,
#'                       glm_weight_function = weibull_glm_weight_function,
#'                       shur_correction_function = weibull_shur_correction2,
#'                       unconstrained_fit_fxn = unconstrained_fit_weibull,
#'                       opt = FALSE,
#'                       wiggle_penalty = 0,
#'                       flat_ridge_penalty = 0,
#'                       K = 0,
#'                       status = pbc$status!=0)
#' print(summary(model_fit))
#'
#'
#' ## Survreg results
#' survreg_fit <- survreg(Surv(time, status!=0) ~ trt + stage + hepato + bili,
#'                        df)
#' print(summary(survreg_fit))
#'
#' ## Dispersion = scale^2
#' print(c(sqrt(model_fit$sigmasq_tilde), survreg_fit$scale))
#'
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[quadprog]{solve.QP}} for quadratic programming optimization
#'   \item \code{\link[plotly]{plot_ly}} for interactive plotting
#' }
#'
#' @export
#' @rdname lgspline-package
#' @aliases lgspline
#' @concept smoothing spline
#' @docType function
lgspline <- function(
    predictors = NULL,
    y = NULL,
    formula = NULL,
    response = NULL,
    standardize_response = TRUE,
    standardize_predictors_for_knots = TRUE,
    standardize_expansions_for_fitting = TRUE,
    family = gaussian(),
    glm_weight_function = function(mu,
                                   y,
                                   order_indices,
                                   family,
                                   dispersion,
                                   observation_weights,
                                   ...){
      if(any(!is.null(observation_weights))){
        family$variance(mu) * observation_weights
      } else {
        family$variance(mu)
      }
    },
    shur_correction_function = function(X,
                                        y,
                                        B,
                                        dispersion,
                                        order_list,
                                        K,
                                        family,
                                        observation_weights,
                                        ...){
      lapply(1:(K+1), function(k)0)
    },
    need_dispersion_for_estimation = FALSE,
    dispersion_function = function(mu,
                                   y,
                                   order_indices,
                                   family,
                                   observation_weights,
                                   ...){ 1 },
    K = NULL,
    custom_knots = NULL,
    cluster_on_indicators = FALSE,
    make_partition_list = NULL,
    previously_tuned_penalties = NULL,
    smoothing_spline_penalty = NULL,
    opt = TRUE,
    use_custom_bfgs = TRUE,
    delta = NULL,
    tol = 10*sqrt(.Machine$double.eps),
    log_initial_wiggle = c(-25, 20, -15, -10, -5),
    log_initial_flat = c(-14, -7),
    wiggle_penalty = 2e-7,
    flat_ridge_penalty = 2e-8,
    unique_penalty_per_partition = TRUE,
    unique_penalty_per_predictor = TRUE,
    penalty_ridge = 1e-8,
    predictor_penalties = NULL,
    partition_penalties = NULL,
    include_quadratic_terms = TRUE,
    include_cubic_terms = TRUE,
    include_quartic_terms = FALSE,
    include_2way_interactions = TRUE,
    include_3way_interactions = TRUE,
    include_quadratic_interactions = TRUE,
    just_linear_with_interactions = NULL,
    just_linear_without_interactions = NULL,
    exclude_interactions_for = NULL,
    exclude_these_expansions = NULL,
    custom_basis_fxn = NULL,
    include_constrain_fitted = TRUE,
    include_constrain_first_deriv = TRUE,
    include_constrain_second_deriv = TRUE,
    include_constrain_interactions = TRUE,
    cl = NULL,
    chunk_size = NULL,
    parallel_eigen = TRUE,
    parallel_trace = FALSE,
    parallel_aga = FALSE,
    parallel_matmult = FALSE,
    parallel_unconstrained = FALSE,
    parallel_find_neighbors = FALSE,
    parallel_penalty = FALSE,
    parallel_make_constraint = FALSE,
    unconstrained_fit_fxn = unconstrained_fit_default,
    keep_weighted_Lambda = FALSE,
    iterate_tune = TRUE,
    iterate_final_fit = TRUE,
    qp_Amat = NULL,
    qp_bvec = NULL,
    qp_meq = 0,
    qp_monotonic_increase = FALSE,
    qp_monotonic_decrease = FALSE,
    qp_range_upper = NULL,
    qp_range_lower = NULL,
    qp_Amat_fxn = NULL,
    qp_bvec_fxn = NULL,
    qp_meq_fxn = NULL,
    constraint_value_vectors = cbind(),
    constraint_vectors = cbind(),
    return_G = TRUE,
    return_Ghalf = TRUE,
    return_U = TRUE,
    estimate_dispersion = TRUE,
    return_varcovmat = TRUE,
    custom_penalty_mat = NULL,
    cluster_args = c(custom_centers = NA, nstart = 10),
    dummy_dividor = 0.00000000000000000000012345672152894,
    dummy_adder = 0.000000000000000002234567210529,
    verbose = FALSE,
    verbose_tune = FALSE,
    expansions_only = FALSE,
    observation_weights = NULL,
    do_not_cluster_on_these = c(),
    neighbor_tolerance = 1 + 1e-16,
    null_constraint = NULL,
    critical_value = qnorm(1-0.05/2),
    data = NULL,
    weights = NULL,
    no_intercept = FALSE,
    VhalfInv = NULL,
    include_warnings = TRUE,
    ...
  ){

  if(verbose){
    cat('Pre-Processing\n')
  }

  ## If returning expansions only, set VhalfInv to NULL no matter what
  if(expansions_only){
    VhalfInv <- NULL
  }

  ## Standardize response is not compatible with VhalfInv
  if(!is.null(VhalfInv) & standardize_response & include_warnings){
    warning('VhalfInv provided but standardize response is TRUE: standardize_response = TRUE is not compatible with the argument for VhalfInv, so standardize_response = FALSE will be imposed. You can work around this by pre-normalizing the response before inputting into the function.')
    standardize_response <- FALSE
  }

  ## Update naming conventions, if first argument is a formula and second is a
  # data frame, assumed by R-like interfaces for user convenience.
  if(any(!is.null(predictors)) & any(!is.null(y))){
    if(any(class(y) == 'data.frame' & inherits(predictors, "formula"))){
      data <- y
    }
  }

  ## Update naming conventions, if response supplied in place of "y" for user
  # convenience.
  if(any(is.null(predictors)) & any(!is.null(formula))){
    predictors <- formula
  } else if(any(is.null(predictors))){
    stop('Predictors argument is NULL without formula supplied. Either supply a formula to predictors OR formula argument, or a data frame to predictors argument, or a matrix of numeric predictors to the predictor argument.')
  }

  ## Update naming conventions, if response supplied in place of "y" for user
  # convenience.
  if(any(is.null(y)) & any(!is.null(response))){
    y <- response
    rm(response)
  }

  ## Weights is just an R-friendly argument to be passed to observation_weights
  # actually used by the function.
  if(any(is.null(observation_weights)) &
     any(!is.null(weights))){
    observation_weights <- weights
    weights <- NULL
  }

  ## Check cluster args for compatibility
  if(any(!is.na(cluster_args[[1]]))){
    ncluster <- try({nrow(cluster_args[[1]])},
                     silent = TRUE)
    if(any(class(ncluster) == 'try-error')){
      stop('custom_centers should be a matrix, do not include any other arguments within the cluster_args function if you include custom_centers.')
    }
    if(!is.null(K)){
      if(ncluster != (K+1) & include_warnings){
        warning('K must be equal to number of custom_centers minus 1. Updating K for compatibility.')
        K <- ncluster - 1
      }
    } else {
      K <- ncluster - 1
    }
  }

  ## Check data and formula argument
  if(!is.null(data) & !inherits(predictors, "formula")) {
    stop("If submitting data argument, formula must be supplied and variables must match. Otherwise, use predictors and y (or response) arguments directly.\n",
         "Example: lgspline(y ~ spl(x1, x2) + x3 + x4*x5, data = my_data)\n",
         "Example: lgspline(y ~ ., data = my_data)\n")
  }

  ## Handle formula interface
  if(inherits(predictors, "formula")) {
    ## Check data argument
    if(is.null(data)) {
      stop("When using formula interface, data argument must be provided. ",
           "Example: lgspline(y ~ spl(x1, x2) + x3 + x4*x5, data = my_data)")
    }

    ## Try to coerce data to data.frame
    tryCatch({
      data <- as.data.frame(data)
    }, error = function(e) {
      stop("Could not coerce data argument to data.frame. ",
           "Please provide data in a format coercible to data.frame. ",
           "Examples: data.frame, tibble, or matrix.")
    })

    ## Check column names are present
    if(any(is.null(colnames(data)))){
      stop('Column names of data must be supplied if formula is supplied, and column names must match what is provided in formula.')
    }

    ## Stringed formula
    form_paste0 <- paste0(predictors)

    ## If formula is y ~ ., replace with y ~ spl(x1, x2, ...)
    if(form_paste0[1] == '~' &
       (gsub(' ', '', form_paste0[3]) %in% c('.', '0+.','.+0')) &
       length(form_paste0) == 3){
         if(form_paste0[2] == '.'){
           stop('. ~ . is not valid for this function, specify the response directly, like y ~ .')
         } else if(gsub(' ', '', form_paste0[3]) == '0+.' |
                   gsub(' ', '', form_paste0[3]) == '.+0'){
           ## No intercept
           predictors <- as.formula(
             paste0(form_paste0[2],
                    ' ~ 0 + spl(',
                    paste(colnames(data)[-which(colnames(data) ==
                                                  form_paste0[2])],
                          collapse = ', '),
                    ')')
           )
         } else {
           ## With intercept
           predictors <- as.formula(
             paste0(form_paste0[2],
                    ' ~ spl(',
                    paste(colnames(data)[-which(colnames(data) ==
                                                form_paste0[2])],
                          collapse = ', '),
                    ')')
           )
         }
    }

    ## Check numeric columns
    non_numeric <- c(1:ncol(data))[!sapply(data, is.numeric)]
    if(length(non_numeric) > 0) {
      for(v in non_numeric){
        one_hot <- create_onehot(data[[v]])
        data <- cbind(data, one_hot[,-1,drop=FALSE]) # dummy-intercept coding
        data[[v]] <- NULL
      }
    }

    ## Parse formula
    terms <- terms(predictors)
    term_labels <- attr(terms, "term.labels")

    ## Check for no intercept specification in formula
    if(inherits(predictors, "formula")) {
      formula_text <- gsub(" ", "", Reduce(paste, deparse(predictors)))
      if(grepl("\\+0|0\\+", formula_text)) {
        no_intercept <- TRUE
      }
    }

    ## Get the factors matrix which shows interaction structure
    factors <- attr(terms, "factors")

    ## Initialize term containers
    spline_terms <- character()
    linear_no_int <- character()
    linear_with_int <- character()

    ## First pass to identify spline terms and extract their variables
    for(term in term_labels) {
      if(grepl("^spl\\(.*\\)$", term)) {
        vars <- gsub("^spl\\((.*)\\)$", "\\1", term)
        spline_terms <- c(spline_terms, trimws(strsplit(vars, ",")[[1]]))
      }
    }

    ## Identify interaction structure using factors matrix
    var_names <- rownames(factors)[-1] # Remove responses

    ## Extract all types of variables in formula
    all_formula_vars <- unique(c(spline_terms, var_names))
    formula_cols <- which(colnames(data) %in% all_formula_vars)

    ## Non-spline variables that appear in interactions (from factors matrix)
    interaction_terms <- term_labels[attr(terms, "order") > 1]
    nonspline_interact_vars <- unique(unlist(lapply(interaction_terms,
                                                    function(term) {
      if(grepl(":", term)) {
        vars <- strsplit(term, ":")[[1]]
        vars[!vars %in% spline_terms]
      }
    })))

    ## Add explicit interactions for spline terms
    if(length(spline_terms) > 1) {
      # Generate all possible interactions between spline terms
      spline_interactions <- combn(spline_terms, 2, simplify=FALSE)
      spline_triplets <- if(length(spline_terms) >= 3) {
        combn(spline_terms, 3, simplify=FALSE)
      } else {
        list()
      }

      # Add 2-way interactions
      for(pair in spline_interactions) {
        interaction_terms <- c(interaction_terms, paste(pair, collapse=":"))
      }

      # Add 3-way interactions
      for(triplet in spline_triplets) {
        interaction_terms <- c(interaction_terms, paste(triplet, collapse=":"))
      }
    }

    ## Get indices of variables in raw expansions (after response removed)
    resp_ind <- which(colnames(data) == paste0(terms[[2]]))
    var_positions <- match(var_names, colnames(data[,-resp_ind]))

    ## Get allowed interaction pairs
    allowed_pairs <- lapply(interaction_terms[grepl(":", interaction_terms)],
                            function(term) {
      vars <- strsplit(term, ":")[[1]]
      match(vars, colnames(data[,-resp_ind]))
    })

    ## Generate exclusion patterns for non-spline vars
    exclude_patterns <- c()

    ## Get all possible 2-way interactions between ANY variables
    vars <- colnames(data[,-resp_ind])
    for(ii in seq_along(vars)) {
      for(jj in seq_along(vars)) {
        if(ii != jj) {
          pattern <- get_interaction_patterns(c(vars[ii], vars[jj]))

          ## Skip if not in formula
          if(any(!(c(vars[ii], vars[jj]) %in% all_formula_vars))){
            next
          }

          ## Skip if any variable in exclude_interactions_for
          if(!is.null(exclude_interactions_for)){
            if(any(c(ii, jj) %in% exclude_interactions_for)){
              next
            }
          }

          ## Only keep spline-spline interactions and explicit interactions
          if(!all(c(vars[ii], vars[jj]) %in% spline_terms) &
             length(spline_terms) > 0 &
             all(c(ii, jj) %in% formula_cols)){
            exclude_patterns <- c(exclude_patterns, pattern)
          }
        }
      }
    }

    ## Add all possible 3-way interactions to exclusions
    if(include_3way_interactions) {
      for(ii in seq_along(vars)) {
        for(jj in seq_along(vars)) {
          for(kk in seq_along(vars)) {
            if(ii != jj && jj != kk && ii != kk) {
              triplet_vars <- c(vars[ii], vars[jj], vars[kk])
              pattern <- get_interaction_patterns(triplet_vars)

              ## Skip if not in formula
              if(any(!(triplet_vars %in% all_formula_vars))){
                next
              }

              ## Skip if any variable in exclude_interactions_for
              if(!is.null(exclude_interactions_for)){
                if(any(c(ii, jj, kk) %in% exclude_interactions_for)){
                  next
                }
              }

              ## Skip if only spline terms
              if(all(triplet_vars %in% spline_terms)){
                next
              }

              ## If ANY var is a spline term but not ALL are spline terms,
              # we should exclude this interaction
              if(any(triplet_vars %in% spline_terms) &&
                 !all(triplet_vars %in% spline_terms)) {
                exclude_patterns <- c(exclude_patterns, pattern)
                next
              }

              ## For non-spline terms, allow explicitly specified interactions
              if(!any(triplet_vars %in% spline_terms)) {
                # Get interaction terms that could involve these variables
                relevant_terms <- interaction_terms[grepl(paste(triplet_vars,
                                                                collapse="|"),
                                                          interaction_terms)]

                # Check if this exact triplet exists in any order
                matches_interaction <- any(sapply(strsplit(relevant_terms, ":"),
                                                  function(term) {
                  length(term) == 3 && all(sort(triplet_vars) == sort(term))
                }))

                if(!matches_interaction) {
                  exclude_patterns <- c(exclude_patterns, pattern)
                }
              }
            }
          }
        }
      }
    }

    ## Remove explicitly allowed interactions from exclusions
    for(term in interaction_terms) {
      vars <- strsplit(term, ":")[[1]]
      allowed <- get_interaction_patterns(vars)
      exclude_patterns <- setdiff(exclude_patterns, allowed)
    }

    ## Convert to positional notation
    vars <- colnames(data[,-resp_ind])
    for(i in seq_along(vars)) {
      exclude_patterns <- gsub(vars[i], paste0("_", i, "_"), exclude_patterns)
    }

    ## Append to custom exclusions
    if(!is.null(exclude_these_expansions)){
      exclude_these_expansions <- c(exclude_these_expansions,
                                    exclude_patterns)
    } else if (length(exclude_patterns) > 0){
      exclude_these_expansions <- exclude_patterns
    }

    ## For each variable, determine if linear with or without interactions
    for(var in var_names) {
      if(length(spline_terms) > 0){
        if(var %in% spline_terms) next # Skip spline terms
      }

      ## Check if this variable appears in any interactions
      var_terms <- which(factors[var,] > 0)
      if(length(var_terms) > 0) {
        ## If any term containing this variable has order > 1,
        # it's in an interaction
        if(any(attr(terms, "order")[var_terms] > 1)) {
          linear_with_int <- c(linear_with_int, var)
        } else {
          linear_no_int <- c(linear_no_int, var)
        }
      } else {
        linear_no_int <- c(linear_no_int, var)
      }
    }

    ## Remove duplicates and ensure proper separation
    linear_with_int <- unique(linear_with_int)
    linear_no_int <- setdiff(unique(linear_no_int),
                             c(spline_terms, linear_with_int))
    linear_no_int <- linear_no_int[!(substr(linear_no_int,
                                           1,
                                           4) == 'spl(')]

    ## Create predictors matrix and response for compatibility
    predictors <- data[, formula_cols, drop=FALSE]
    y <- data[,paste0(terms[[2]])]

    ## Convert variable names to column indices
    new_just_linear_without_interactions <- match(linear_no_int,
                                              colnames(predictors))
    new_just_linear_with_interactions <- match(linear_with_int,
                                           colnames(predictors))
    new_just_linear_without_interactions <- new_just_linear_without_interactions[
      !(new_just_linear_without_interactions %in% spline_terms)
    ]
    new_just_linear_with_interactions <- new_just_linear_with_interactions[
      !(new_just_linear_with_interactions %in% spline_terms)
    ]
    if(is.null(just_linear_with_interactions)){
      just_linear_with_interactions <- new_just_linear_with_interactions
    } else{
      just_linear_with_interactions <- unique(c(
        just_linear_with_interactions,
        new_just_linear_with_interactions
      ))
    }
    if(is.null(just_linear_without_interactions)){
      just_linear_without_interactions <- new_just_linear_without_interactions
    } else {
      just_linear_without_interactions <- unique(c(
        just_linear_without_interactions,
        new_just_linear_without_interactions
      ))
    }
  }

  ## Not a formula - try to coerce to matrix
  tryCatch({
    predictors <- as.matrix(predictors)
  }, error = function(e) {
    stop("Could not coerce predictors to matrix. ",
         "predictors must be either a formula or an object coercible to matrix.",
         "Examples:\n",
         "  Formula: lgspline(y ~ spl(x1, x2) + x3, data = my_data)\n",
         "  Matrix:  lgspline(predictors = Tmat, y = y)")
  })

  ## Check numeric type
  if(any(!is.numeric(predictors))){
    stop("predictors matrix must be numeric. ",
         "Please convert categorical variables to numeric indicators.")
  }

  ## Check response for missings
  if(any(is.na(y) | is.nan(y) | !is.finite(y))){
    stop("NA, NaN, or infinite value detected in response.")
  }

  ## Original predictor names
  og_cols <- colnames(predictors)
  if(!any(is.null(og_cols))){
    replace_colnames <- TRUE
  } else {
    replace_colnames <- FALSE
  }

  ## Alternative parameterization of null constraint for ease of use
  if(any(!is.null(null_constraint)) &
     length(constraint_vectors) > 0 &
     length(constraint_value_vectors) == 0){
     constraint_value_vectors <-
       constraint_vectors %**%
       invert(gramMatrix(cbind(constraint_vectors))) %**%
       cbind(c(null_constraint))
  }

  ## Check nrow of input predictors and matrix coersion
  t <- try({if(nrow(as(predictors,'matrix')) < 3){
    stop('Need at least 3 observations to fit model')
  }}, silent = TRUE)
  if(class(t) == 'try-error'){
    stop('Cannot coerce predictors to a matrix')
  }

  ## Check if no spline terms - if so, set K = 0
  if(length(unique(c(just_linear_with_interactions,
              just_linear_without_interactions))) == ncol(predictors)){
    K <- 0
  }

  ## Check if custom knots is not missing, that it can be
  # coerced to a matrix
  if(any(!(is.null(custom_knots)))){
    custom_knots <- try(cbind(custom_knots),
                        silent = TRUE)
    if(any(class(custom_knots) == 'try-error') & include_warnings){
      warning('custom_knots must be a matrix, or should be coercible to it. custom_knots will be ignored.')
      custom_knots <- NULL
    }
  }

  ## Model fit procedure called
  model_fit <- try({lgspline.fit(predictors,
                                 y,
                                 standardize_response,
                                 standardize_predictors_for_knots,
                                 standardize_expansions_for_fitting,
                                 family,
                                 glm_weight_function,
                                 shur_correction_function,
                                 need_dispersion_for_estimation,
                                 dispersion_function,
                                 K,
                                 custom_knots,
                                 cluster_on_indicators,
                                 make_partition_list,
                                 previously_tuned_penalties,
                                 smoothing_spline_penalty,
                                 opt,
                                 use_custom_bfgs,
                                 delta,
                                 tol,
                                 log_initial_wiggle,
                                 log_initial_flat,
                                 wiggle_penalty,
                                 flat_ridge_penalty,
                                 unique_penalty_per_partition,
                                 unique_penalty_per_predictor,
                                 penalty_ridge,
                                 predictor_penalties,
                                 partition_penalties,
                                 include_quadratic_terms,
                                 include_cubic_terms,
                                 include_quartic_terms,
                                 include_2way_interactions,
                                 include_3way_interactions,
                                 include_quadratic_interactions,
                                 just_linear_with_interactions,
                                 just_linear_without_interactions,
                                 exclude_interactions_for,
                                 exclude_these_expansions,
                                 custom_basis_fxn,
                                 include_constrain_fitted,
                                 include_constrain_first_deriv,
                                 include_constrain_second_deriv,
                                 include_constrain_interactions,
                                 cl,
                                 chunk_size,
                                 parallel_eigen,
                                 parallel_trace,
                                 parallel_aga,
                                 parallel_matmult,
                                 parallel_unconstrained,
                                 parallel_find_neighbors,
                                 parallel_penalty,
                                 parallel_make_constraint,
                                 unconstrained_fit_fxn,
                                 keep_weighted_Lambda,
                                 iterate_tune,
                                 iterate_final_fit,
                                 qp_Amat,
                                 qp_bvec,
                                 qp_meq,
                                 qp_monotonic_increase,
                                 qp_monotonic_decrease,
                                 qp_range_upper,
                                 qp_range_lower,
                                 qp_Amat_fxn,
                                 qp_bvec_fxn,
                                 qp_meq_fxn,
                                 constraint_value_vectors,
                                 constraint_vectors,
                                 return_G,
                                 return_Ghalf,
                                 return_U,
                                 estimate_dispersion,
                                 return_varcovmat,
                                 custom_penalty_mat,
                                 cluster_args,
                                 dummy_dividor,
                                 dummy_adder,
                                 verbose,
                                 verbose_tune,
                                 expansions_only,
                                 observation_weights,
                                 do_not_cluster_on_these,
                                 neighbor_tolerance,
                                 no_intercept,
                                 VhalfInv,
                                 include_warnings,
                                 ...
  )}, silent = TRUE)

  ## Return try error if model fails to to be fit
  if(any(class(model_fit) == 'try-error') & include_warnings){
    warning("Model fitting error: try verbose = TRUE, checking for NAs, adjusting starting tuning grid, or K. If using parallel options, check your cluster, submit to cl argument if valid, and make sure base R parallel package is loaded.")
    return(model_fit)
  }

  ## Return expansions only and associated components
  if(expansions_only){
    return(c(model_fit, list(og_cols = colnames(predictors))))
  }

  ## Rename elements of B_raw and B according to actual column names
  if(replace_colnames){
    ## perform for B and rownames A
    og_colnames_match <- cbind(og_cols, paste0('_', 1:ncol(predictors), '_'))
    new_names <- sapply(names(model_fit$B[[1]]), function(nm){
      for(ii in 1:nrow(og_colnames_match)){
        nm <- gsub(og_colnames_match[ii,2], og_colnames_match[ii,1], nm)
      }
      nm
    })

    for(k in 1:(model_fit$K+1)){
      rownames(model_fit$B[[k]]) <- new_names
      names(model_fit$B[[k]]) <- new_names
    }
    rownames(model_fit$A) <- paste0(rep(paste0('partition',
                                               1:(model_fit$K+1)),
                                        each = model_fit$p),
                                    "_",
                                    new_names)
  }


  ## Inference using Wald:
  # score Test/LR Test can be obtained using these components as well
  if(return_varcovmat){

    ## Univariate inference
    wald_univariate <- function(scale_vcovmat_by = 1,
                                cv = critical_value){
      return_list <- list(
        est = unlist(model_fit$B),
        se = sqrt(scale_vcovmat_by * diag(model_fit$varcovmat))
      )
      return_list$stat <- return_list$est/return_list$se
      return_list$interval_lb <- return_list$se*
        (return_list$stat - cv)
      return_list$interval_ub <- return_list$se*
        (return_list$stat + cv)

      ## If normal errors, use exact t-test
      # Otherwise use Wald's N(0,1) approximation
      if(!(any(!(paste0(family)[1:4] == paste0(gaussian())[1:4])))){
        return_list$pval <- 2*(1-pt(abs(return_list$stat),
                                    df = model_fit$N - model_fit$trace_XUGX))
      } else {
        return_list$pval <- 2*(1-pnorm(abs(return_list$stat)))
      }

      return(return_list)
    }
  } else {
    wald_univariate <- function(scale_vcovmat_by = 1,
                                cv = critical_value){
      NULL
    }
  }

  ## Univariate inference
  model_fit$wald_univariate <- wald_univariate

  ## Function for generating draws from posterior/posterior predictive
  model_fit$generate_posterior <- function(new_sigmasq_tilde =
                                           model_fit$sigmasq_tilde,
                                           new_predictors = predictors,
                                           theta_1 = 0,
                                           theta_2 = 0,
                                           posterior_predictive_draw =
                                             function(N,
                                                      mean,
                                                      sqrt_dispersion,
                                                      ...)rnorm(
                                                        N, mean, sqrt_dispersion
                                                      ),
                                           draw_dispersion = TRUE,
                                           include_posterior_predictive = FALSE,
                                           num_draws = 1,
                                           ...){

    ## Check compatibility, that new_predictors should be a matrix
    if(any(!is.null(new_predictors))){
      new_predictors <- try(as(cbind(new_predictors), 'matrix'), silent = TRUE)
      if(any(class(new_predictors) == 'try-error')){
        stop('New predictors should be able to be coerced into matrix form.')
      }
    }

    ## Quick and dirty way around this in R.....
    only_1 <- FALSE
    if(nrow(new_predictors) == 1){
      only_1 <- TRUE
      new_predictors <- rbind(new_predictors, new_predictors)
    }

    ## Helpful components
    nc <- model_fit$p # number of cubic expansions (P when K = 0)
    K <- model_fit$K # number of partitions - 1
    nr <- model_fit$N # number of observations in-sample

    res <- lapply(1:num_draws,function(m){
      ## Draw a dispersion parameter, if applicable, from InvG distribution
      if(draw_dispersion){
        shape <- theta_1 + 0.5*nr
        rate <- theta_2 + 0.5*sum((model_fit$y - model_fit$ytilde)^2)
        if(shape <= 0){
          stop("Posterior inverse-gamma shape is <= 0, increase theta_1 argument to draw a dispersion parameter.")
        }
        if(rate <= 0){
          stop("Posterior inverse-gamma rate is <= 0, increase theta_2 argument to draw a dispersion parameter.")
        }
        post_draw_sigmasq <-
          1/rgamma(1,
                   shape,
                   rate)

        ## If degenerate or infinite, default to the point-estimate provided
        # and provide warning
        if((is.nan(post_draw_sigmasq) | !is.finite(post_draw_sigmasq)) & include_warnings){
          warning("Infinite or NaN posterior draw of dispersion detected.")
          post_draw_sigmasq <- new_sigmasq_tilde
        }
      } else {
        post_draw_sigmasq <- new_sigmasq_tilde
      }

      ## Draw posterior "errors" of beta coefficients under smoothing constraints
      # Purpose here is to generate draws on the standardized y-scale,
      # using standardized X in the design matrix
      # we generate posterior draws of beta on model-scale,
      # before backtransforming
      # Unscaled by dispersion
      # = UG^{1/2}z
      post_draw_coefficients_err <-
        (1/model_fit$sd_y) * # un-scale the dispersion that was drawn
        sqrt(post_draw_sigmasq) * # sqrt-dispersion that was drawn
        (model_fit$U %**%
           cbind(Reduce("c", lapply(1:(K+1),function(k){
             c(model_fit$Ghalf[[k]] %**%
                 cbind(rnorm(nc)))
           })))) # UG^{1/2}z

      ## Add to B MAPs we've already fit, then backtransform for raw scale
      post_draw_coefficients <- lapply(1:(K+1),function(k){
        raw_draw <- post_draw_coefficients_err[1:nc +(k-1)*nc] +
          model_fit$B_raw[[k]]

        ## Un-scale, based on centered-and-scaled y
        raw_draw <- raw_draw * model_fit$sd_y # multiply by sd of y

        ## Add mean of y to all intercepts
        raw_draw[1] <- raw_draw[1] + model_fit$mean_y

        ## Backtransform for un-standardized predictors
        return(model_fit$backtransform_coefficients(raw_draw))

      })

      ## Return posterior predictive draws
      if(include_posterior_predictive){
        ## Posterior-predictive mean
        post_pred_mean <- model_fit$predict(
          new_predictors,
          B_predict = post_draw_coefficients)

        ## Posterior-predictive realization
        post_pred_draw <- posterior_predictive_draw(length(post_pred_mean),
                                post_pred_mean,
                                sqrt(post_draw_sigmasq),
                                ...)

        return(list(post_pred_draw = post_pred_draw,
                    post_draw_coefficients = post_draw_coefficients,
                    post_draw_sigmasq = post_draw_sigmasq))

      ## Return posterior coefficient draws, and sigma sq
      } else {
        return(list(
          post_draw_coefficients = post_draw_coefficients,
          post_draw_sigmasq = post_draw_sigmasq
        ))
      }
    })
    if(num_draws == 1){
      if(only_1){
        res[[1]][[1]] <- res[[1]][[1]][1]
      }
      return(res[[1]])
    }

    ## Combine results
    post_draw_coefficients <- lapply(res, `[[`, "post_draw_coefficients")
    post_draw_sigmasq <- lapply(res, `[[`, "post_draw_sigmasq")
    if(include_posterior_predictive){
      post_pred_draw <-
        Reduce("cbind", lapply(res, `[[`, "post_pred_draw"))

      if(only_1){
        post_pred_draw <- post_pred_draw[1,,drop=FALSE]
      }

      return(list(post_pred_draw = post_pred_draw,
                  post_draw_coefficients = post_draw_coefficients,
                  post_draw_sigmasq = post_draw_sigmasq))
    }
    return(list(
      post_draw_coefficients = post_draw_coefficients,
      post_draw_sigmasq = post_draw_sigmasq
    ))
  }

  ## Find global maximum/minimum
  model_fit$find_extremum <- function(
    quick_heuristic = TRUE, # only search top-performing partition
    initial = NULL, # initial values, useful for fixing binary predictors which aren't optimized
    parallel = FALSE, # run in parallel
    cl = NULL, # cluster
    B_predict = NULL, # custom coefficients, if desired
    minimize = FALSE, # minimize vs. maximize
    stochastic = FALSE, # add noise to candidates proposed by L-BFGS-B
    stochastic_draw = function(mu,
                               sigma){N <- length(mu)
                               rnorm(
                                 N, mu, sigma
                               )},
    sigmasq_tilde = model_fit$sigmasq_tilde, # Variance for stochastic optimization
    custom_objective_function = NULL,# custom function for maximizing/minimizing with args mean (mu), std dev (sigma), and best-observed
    custom_objective_gradient = NULL # custom gradient of function for maximizing/minimizing with args mean (mu), std dev (sigma), best observed thus far, and  chain-rule derivative of fitted function (x')^{t}b to pass through
    # example for expected improvement:
    # custom_objective_function = function(mu, sigma, y_best) {
    #   d <- y_best - mu
    #   d * pnorm(d/sigma) + sigma * dnorm(d/sigma)
    # }
    # custom_objective_gradient = function(mu, sigma, y_best, d_mu) {
    #   d <- y_best - mu
    #   z <- d/sigma
    #   d_z <- -d_mu/sigma
    #   pnorm(z)*(-d_mu) + d*dnorm(z)*d_z - sigma*z*dnorm(z)*d_z
    # }
    ){
    sigma_tilde = sqrt(sigmasq_tilde)
    ## Switch for maximizing or minimizing function
    # Since optim() default minimizes functions, this is -1 for maximize
    # Needed for implementation details that simply using
    # optim() option is insufficient for
    min_or_max <- 2*(minimize-0.5)

    ## Re-assign predictions if B_predict offered
    if(any(is.null(B_predict))){
      B_predict <- model_fit$B
    } else {
      model_fit$ytilde <-
        model_fit$predict(
          predictors,
          B_predict = B_predict
        )
    }

    ## Use all partitions by default
    partitions <- 1:(model_fit$K+1)

    ## If any NaN, return randomly selected value and predicted performance
    if(any(is.nan(model_fit$ytilde))){
      dummy_draw <- c(sapply(1:ncol(predictors), function(j){
        runif(1, min(predictors[,j]), max(predictors[,j]))
      }))
      dummy_y <- model_fit$predict(rbind(dummy_draw))
      return(list(
        t = dummy_draw,
        y = dummy_y
      ))
    }

    ## Only use partition with best fitted, by default
    if(quick_heuristic){
      best_fitted <- which.max(model_fit$ytilde * (-min_or_max))
      partitions <- which(sapply(1:(model_fit$K+1),function(k){
        best_fitted %in% model_fit$order_list[[k]]
      }))
      parallel <- FALSE
    }


    ## Find top-performing predictors value
    if(parallel & any(!is.null(cl))){
      ## Create shared environment in global environment
      assign("shared_env", new.env(), envir = .GlobalEnv)

      ## Prepare shared variables
      shared_vars <- list(
        B_predict = B_predict,
        custom_objective_function = custom_objective_function,
        custom_objective_gradient = custom_objective_gradient,
        sigma_tilde = sigma_tilde,
        stochastic = stochastic,
        stochastic_draw = stochastic_draw,
        get_polynomial_expansions = get_polynomial_expansions,
        take_interaction_2ndderivative = take_interaction_2ndderivative,
        take_derivative = take_derivative,
        make_derivative_matrix = make_derivative_matrix,
        knot_expand_list = knot_expand_list,
        initial = initial,
        family = family,
        nc = nc,
        nr = nr
      )

      ## Assign shared variables to shared environment
      for(nm in names(shared_vars)) {
        assign(nm, shared_vars[[nm]], envir = shared_env)
      }

      ## Export shared environment
      clusterExport(cl, "shared_env")

      ## Setup each cluster node with necessary functions
      clusterEvalQ(cl, {
        library(lgspline)
        `%**%` <- efficient_matrix_mult
        ## Load shared environment objects
        list2env(as.list(shared_env), .GlobalEnv)
      })

      ## Run this in parallel if desired
      # remove empty partitions first
      partitions_keep <- c(c(), which(sapply(partitions, function(k){
        nrow(model_fit$X[[k]])
      }) > 0))
      best_per_partition <- parLapply(cl, partitions[partitions_keep],
                                      function(k, predictors, model_fit
        ){

        if(any(!is.null(initial))){
          predictors_vals <- initial
        } else {
          ## Extract best fitted value for initialization
          yk <- model_fit$X[[k]] %**% B_predict[[k]]
          best <- which.max(-yk*min_or_max)
          predictors_vals <- predictors[model_fit$order_list[[k]][best],
                                        , drop=FALSE]
        }

        ## Quasi-newton optimization
        opt <- stats::optim(
          predictors_vals,
          fn = function(par){
            if(!is.null(custom_objective_function)){
              ## Prediction
              pred <- model_fit$predict(new_predictors = rbind(c(par)),
                                        parallel = FALSE,
                                        cl = NULL,
                                        chunk_size = NULL,
                                        num_chunks = NULL,
                                        rem_chunks = NULL,
                                        B_predict = B_predict)
              if(stochastic){
                ## Add random noise if desired
                pred <- stochastic_draw(pred, sigma_tilde)
              }
              ## Throw into custom objective if desired
              min_or_max*custom_objective_function(pred,
                                                   sigma_tilde,
                                                   max(-y*min_or_max))
            } else {
              ## Otherwise, no custom objective
              pred <- model_fit$predict(new_predictors = rbind(c(par)),
                                        parallel = FALSE,
                                        cl = NULL,
                                        chunk_size = NULL,
                                        num_chunks = NULL,
                                        rem_chunks = NULL,
                                        B_predict = B_predict)
              if(stochastic){
                pred <- stochastic_draw(pred, sigma_tilde)
              }
              min_or_max*pred
            }
          },
          gr = function(par){
            if(!is.null(custom_objective_gradient)) {
              ## Repeat for gradient
              pred <- model_fit$predict(new_predictors = rbind(c(par)),
                                        parallel = FALSE,
                                        cl = NULL,
                                        chunk_size = NULL,
                                        num_chunks = NULL,
                                        rem_chunks = NULL,
                                        B_predict = B_predict)
              gr <- model_fit$predict(new_predictors = rbind(c(par)),
                                      parallel = FALSE,
                                      cl = NULL,
                                      chunk_size = NULL,
                                      num_chunks = NULL,
                                      rem_chunks = NULL,
                                      B_predict = B_predict,
                                      take_first_derivatives = TRUE)$first_deriv
              gr_par <- rep(0, length(par))
              gr_raw <- min_or_max*custom_objective_gradient(pred,
                                                             sigma_tilde,
                                                             max(-y*min_or_max),
                                                             gr)
              gr_par[model_fit$numerics] <- gr_raw
              gr_par
            } else {
              gr_par <- rep(0, length(par))
              gr_raw <- min_or_max*model_fit$predict(
                new_predictors = rbind(c(par)),
                parallel = FALSE,
                cl = NULL,
                chunk_size = NULL,
                num_chunks = NULL,
                rem_chunks = NULL,
                B_predict = B_predict,
                take_first_derivatives = TRUE)$first_deriv
              gr_par[model_fit$numerics] <- gr_raw
              gr_par
            }
          },
          method = 'L-BFGS-B',
          lower = apply(predictors, 2, min),
          upper = apply(predictors, 2, max)
        )

        return(t(cbind(c(opt$par))))

      },
      predictors,
      model_fit
     )
    } else {

      ## Go through each partition, optimize the cubic function within
      # remove empty partitions first
      partitions_keep <- c(c(), which(sapply(partitions, function(k){
        nrow(model_fit$X[[k]])
      }) > 0))
      best_per_partition <- lapply(partitions[partitions_keep], function(k){


        if(any(!is.null(initial))){
          predictors_vals <- initial
        } else {
          ## Extract best fitted value for initialization
          yk <- model_fit$X[[k]] %**% B_predict[[k]]
          best <- which.max(-yk*min_or_max)
          predictors_vals <- predictors[model_fit$order_list[[k]][best],
                                        , drop=FALSE]
        }

        ## Quasi-newton optimization
        opt <- stats::optim(
          predictors_vals,
          fn = function(par){
            if(!is.null(custom_objective_function)){
              ## Prediction
              pred <- model_fit$predict(new_predictors = rbind(c(par)),
                                        parallel = FALSE,
                                        cl = NULL,
                                        chunk_size = NULL,
                                        num_chunks = NULL,
                                        rem_chunks = NULL,
                                        B_predict = B_predict)
              if(stochastic){
                ## Add random noise if desired
                pred <- stochastic_draw(pred, sigma_tilde)
              }
              ## Throw into custom objective if desired
              min_or_max*custom_objective_function(pred,
                                             sigma_tilde,
                                             max(-y*min_or_max))
            } else {
              ## Otherwise, no custom objective
              pred <- model_fit$predict(new_predictors = rbind(c(par)),
                                           parallel = FALSE,
                                           cl = NULL,
                                           chunk_size = NULL,
                                           num_chunks = NULL,
                                           rem_chunks = NULL,
                                           B_predict = B_predict)
              if(stochastic){
                pred <- stochastic_draw(pred, sigma_tilde)
              }
              min_or_max*pred
            }
          },
          gr = function(par){
            if(!is.null(custom_objective_gradient)) {
              ## Repeat for gradient
              pred <- model_fit$predict(new_predictors = rbind(c(par)),
                                        parallel = FALSE,
                                        cl = NULL,
                                        chunk_size = NULL,
                                        num_chunks = NULL,
                                        rem_chunks = NULL,
                                        B_predict = B_predict)
              gr <- model_fit$predict(new_predictors = rbind(c(par)),
                                      parallel = FALSE,
                                      cl = NULL,
                                      chunk_size = NULL,
                                      num_chunks = NULL,
                                      rem_chunks = NULL,
                                      B_predict = B_predict,
                                      take_first_derivatives = TRUE)$first_deriv
              gr_par <- rep(0, length(par))
              gr_raw <- min_or_max*custom_objective_gradient(pred,
                                                       sigma_tilde,
                                                       max(-y*min_or_max),
                                                       gr)
              gr_par[model_fit$numerics] <- gr_raw
              gr_par
            } else {
              gr_par <- rep(0, length(par))
              gr_raw <- min_or_max*model_fit$predict(
                new_predictors = rbind(c(par)),
                parallel = FALSE,
                cl = NULL,
                chunk_size = NULL,
                num_chunks = NULL,
                rem_chunks = NULL,
                B_predict = B_predict,
                take_first_derivatives = TRUE)$first_deriv
              gr_par[model_fit$numerics] <- gr_raw
              gr_par
            }
          },
          method = 'L-BFGS-B',
          lower = apply(predictors, 2, min),
          upper = apply(predictors, 2, max)
        )

        return(t(cbind(c(opt$par))))

      })
    }

    ## Find the global optimum out of all optimal-per-partitions
    best_per_partition <- Reduce("rbind", best_per_partition)
    preds <- model_fit$predict(new_predictors = best_per_partition,
                               B_predict = B_predict)
    global_max <- which.max(-min_or_max*preds)

    ## L-BFGS-B struggle to find optimum on borders,
    # but these are present in data at least in 1-D
    # not performed for custom functions
    if((-min_or_max*preds[global_max] < max(-min_or_max*model_fit$ytilde)) &
         is.null(custom_objective_function) &
         is.null(custom_objective_gradient)){
      return(list(
        t = predictors[which.max(-min_or_max*model_fit$ytilde),,drop=FALSE],
        y = -min_or_max*max(-min_or_max*model_fit$ytilde)
      ))
    }

    ## Otherwise, return the optimized value
    extr <- best_per_partition[global_max, ,drop=FALSE]
    colnames(extr) <- colnames(predictors)
    return(list(
      t = extr,
      y = preds[global_max]
    ))
  }


  ## One-dimensional plotting function
  plot_lgspline_1d <- function(modfit,
                               show_formulas,
                               digits,
                               legend_pos,
                               custom_ylab,
                               custom_predictor_lab,
                               custom_formula_lab,
                               custom_title,
                               text_size_formula,
                               xlim1d,
                               ylim1d,
                               ...) {

    ## For preventing stack issues
    model_fit <- modfit
    drop(modfit)

    ## Linear term and name
    xvals <- lapply(model_fit$X, function(x) x[,2])

    ## For customizing xlab and legend predictor label
    v1 <- colnames(model_fit$X[[1]])[2]
    if(is.null(custom_predictor_lab)){
      if(replace_colnames){
        custom_predictor_lab <- og_cols[as.numeric(substr(v1, 2, nchar(v1)-1))]
      } else {
        custom_predictor_lab <- v1
      }
    }

    ## Fitted values re-organized into list format
    y_fitted <- lapply(1:(model_fit$K + 1), function(k){
      model_fit$ytilde[model_fit$order_list[[k]]]
    })

    ## Rainbow gradient
    cols = rainbow(model_fit$K+1)

    ## Xlab defaults to actual variable name if subitted as NULL
    if(is.null(custom_predictor_lab)){
      xlab <- names(model_fit$B[[1]])[2]
    } else {
      xlab <- custom_predictor_lab
    }

    ## Default xlim/ylim preventing stack issues
    if(is.null(ylim1d)){
      ylim <- c(min(unlist(y_fitted),
                   model_fit$y), max(unlist(y_fitted),
                                     model_fit$y))
    } else {
      ylim <- ylim1d
    }
    if(is.null(xlim1d)){
      xlim <- c(min(unlist(xvals)), max(unlist(xvals)))
    } else {
      xlim <- xlim1d
    }

    ## Basic plot
    plot(xvals[[1]],
         y_fitted[[1]],
         ylim = ylim,
         xlim = xlim,
         xlab = xlab,
         ylab = custom_ylab,
         col = cols[1],
         main = custom_title,
         ...)

    ## Add in other partitions
    if(model_fit$K >= 1){
      for(k in 2:(model_fit$K + 1)){
        points(xvals[[k]],
               y_fitted[[k]],
               xlab = xlab,
               ylab = custom_ylab,
               col = cols[k],
               main = custom_title,
               ...)
      }
    }

    ## Add formulas if requested - using existing names
    if(show_formulas) {
      formulas <- sapply(1:(model_fit$K+1), function(k) {
        coefs <- round(model_fit$B[[k]], digits)
        names(coefs) <- gsub(
          names(model_fit$B[[k]])[2],
          xlab,
          names(coefs)
        )
        names(coefs) <- gsub("\\^2", "²", names(coefs))
        names(coefs) <- gsub("\\^3", "³", names(coefs))
        names(coefs) <- gsub("\\^4", "⁴", names(coefs))
        names(coefs) <- gsub(v1, custom_predictor_lab, names(coefs))
        paste0(custom_formula_lab, " = ", paste(coefs, names(coefs),
                                                collapse = " + "))
      })
      formulas <- gsub('intercept', '', formulas)
      formulas <- gsub('  ', ' ', formulas)
      legend(legend_pos,
             legend = formulas,
             col = cols,
             lwd = 2,
             cex = text_size_formula)
    }
  }

  ## Two-dimensional plotting function
  plot_lgspline_2d <- function(modfit,
                               show_formulas,
                               digits,
                               custom_zlab,
                               custom_formula_lab,
                               custom_predictor_lab1,
                               custom_predictor_lab2,
                               custom_title,
                               text_size_formula,
                               ...) {
    model_fit <- modfit

    ## Modification such that when plotting a categorical + spline effect,
    # we do not plot spline effect vs. spline effect^2, based on how
    # the polynomial expansions are arranged
    if(length(model_fit$nonspline_cols) > 0){
      if(length(model_fit$nonspline_cols) == 2){
        xvals1 <-
          lapply(model_fit$X, function(x) x[,model_fit$nonspline_cols[1]])
        v1 <- colnames(model_fit$X[[1]])[model_fit$nonspline_cols[1]]
        xvals2 <-
          lapply(model_fit$X, function(x) x[,model_fit$nonspline_cols[2]])
        v2 <- colnames(model_fit$X[[1]])[model_fit$nonspline_cols[2]]
      } else {
        xvals1 <-
          lapply(model_fit$X, function(x) x[,model_fit$power1_cols[1]])
        v1 <- colnames(model_fit$X[[1]])[model_fit$power1_cols[1]]
        xvals2 <-
          lapply(model_fit$X, function(x) x[,model_fit$nonspline_cols[1]])
        v2 <- colnames(model_fit$X[[1]])[model_fit$nonspline_cols[1]]
      }
    } else {
      xvals1 <-
        lapply(model_fit$X, function(x) x[,model_fit$power1_cols[1]])
      v1 <- colnames(model_fit$X[[1]])[model_fit$power1_cols[1]]
      xvals2 <-
        lapply(model_fit$X, function(x) x[,model_fit$power1_cols[2]])
      v2 <- colnames(model_fit$X[[1]])[model_fit$power1_cols[2]]
    }

    ## For customizing formula and xlab names
    if(is.null(custom_predictor_lab1)){
      if(replace_colnames){
        custom_predictor_lab1 <- og_cols[as.numeric(substr(v1, 2, nchar(v1)-1))]
      } else {
        custom_predictor_lab1 <- v1
      }
    }
    if(is.null(custom_predictor_lab2)){
      if(replace_colnames){
        custom_predictor_lab2 <- og_cols[as.numeric(substr(v2, 2, nchar(v2)-1))]
      } else {
        custom_predictor_lab2 <- v2
      }
    }

    ## For swapping out custom labels from formulas
    if(any(is.null(og_cols))){
      og_cols <- c('_1_', '_2_')
    }

    ## Fitted values in block-diagonal order
    y_fitted <- lapply(1:(model_fit$K+1), function(k) {
      model_fit$ytilde[model_fit$order_list[[k]]]
    })

    ## Combine data for plotting
    plot_data <- data.frame(
      x = unlist(xvals1),
      y = unlist(xvals2),
      z = unlist(y_fitted),
      partition = factor(rep(1:(model_fit$K+1), sapply(xvals1, length)))
    )

    ## Create formulas for hover text if requested
    if(show_formulas) {
      formulas <- sapply(1:(model_fit$K+1), function(k) {
        coefs <- round(model_fit$B[[k]], digits)
        names(coefs) <- gsub("\\^2", "²", names(coefs))
        names(coefs) <- gsub("\\^3", "³", names(coefs))
        names(coefs) <- gsub("\\^4", "⁴", names(coefs))
        names(coefs) <- gsub(og_cols[1], custom_predictor_lab1, names(coefs))
        names(coefs) <- gsub(og_cols[2], custom_predictor_lab2, names(coefs))
        paste0(custom_formula_lab, " = ", paste(coefs, names(coefs),
                                                collapse = " + "))
      })
      formulas <- gsub('intercept', '', formulas)
      formulas <- gsub('  ', ' ', formulas)
      plot_data$formula <- rep(formulas, sapply(xvals1, length))
    }

    ## Show formulas or not
    if(show_formulas){
      text <- ~formula
    } else {
      text <- NULL
    }

    ## Create plotly plot
    p <- plotly::layout(
      plotly::plot_ly(plot_data,
                      x = ~x,
                      y = ~y,
                      z = ~z,
                      color = ~partition,
                      colors = colorRampPalette(
                        RColorBrewer::brewer.pal(8, "Spectral"))(model_fit$K+1),
                      type = "scatter3d",
                      mode = "markers",
                      text = text,
                      connectgaps = TRUE,
                      hoverinfo = if(show_formulas) "text" else "x+y+z+name",
                      hoverlabel = list(font = list(size = text_size_formula)),
                      ...
      ),
      scene = list(
        xaxis = list(title = custom_predictor_lab1),
        yaxis = list(title = custom_predictor_lab2),
        zaxis = list(title = custom_zlab)
      ),
      title = custom_title
    )

    return(p)
  }

  ## Wrapper
  model_fit$plot <- function(model_fit_in = model_fit,
                             show_formulas = FALSE,
                             digits = 4,
                             legend_pos = "topright",
                             custom_response_lab = "y",
                             custom_predictor_lab = NULL,
                             custom_predictor_lab1 = NULL,
                             custom_predictor_lab2 = NULL,
                             custom_formula_lab = NULL,
                             custom_title = "Fitted Function",
                             text_size_formula = NULL,
                             new_predictors = NULL,
                             xlim = NULL,
                             ylim = NULL,
                             ...){

    ## Check compatibility, that new_predictors should be a matrix
    if(any(!is.null(new_predictors))){
      new_predictors <- try(as(cbind(new_predictors), 'matrix'), silent = TRUE)
      if(any(class(new_predictors) == 'try-error')){
        stop('New predictors should be able to be coerced into matrix form.')
      }
    }

    ## Default text_size_formula depends on q
    if(is.null(text_size_formula)){
      text_size_formula <- ifelse(model_fit_in$q == 1,
                                0.8,
                                8)
    }

    ## Default custom_formula_lab = g(E[y]) for g, a link function
    if(is.null(custom_formula_lab)){
      if(paste0(model_fit_in$family)[2] == 'identity' &
         paste0(model_fit_in$family)[1] == 'gaussian'){
        custom_formula_lab <- custom_response_lab
      } else {
        custom_formula_lab <- paste0(model_fit_in$family$link,
                                     '(E[',
                                     custom_response_lab,
                                     '])')
      }
    }

    ## Reset model-fit components for new predictors
    if(any(!is.null(new_predictors))){
      ## Get basis and knot expansions
      prep <- model_fit_in$predict(new_predictors = new_predictors,
                                   just_expansions = TRUE)
      model_fit_in$X <- prep$expansions

      ## Get order of y by partition
      model_fit_in$order_list <- model_fit_in$knot_expand_function(
                                                  prep$partition_codes,
                                                  prep$partition_bounds,
                                                  nrow(new_predictors),
                                                  cbind(1:nrow(new_predictors)),
                                                  model_fit_in$K)

      ## Make new prediction
      model_fit_in$ytilde <-
        model_fit_in$predict(new_predictors = new_predictors)
    }

    ## 1-D plotting
    if(model_fit_in$q == 1){
      plot_lgspline_1d(model_fit_in,
                       show_formulas,
                       digits,
                       legend_pos,
                       custom_response_lab,
                       custom_predictor_lab,
                       custom_formula_lab,
                       custom_title,
                       text_size_formula,
                       xlim,
                       ylim,
                       ...)
    ## 2-D plotting
    } else if(model_fit_in$q == 2){
      plot_lgspline_2d(model_fit_in,
                       show_formulas,
                       digits,
                       custom_response_lab,
                       custom_formula_lab,
                       custom_predictor_lab1,
                       custom_predictor_lab2,
                       custom_title,
                       text_size_formula,
                       ...)
    } else if(include_warnings){
      warning("No default plotting functions implemented for q > 2")
    }
  }

  ## Important information
  model_fit$critical_value <- critical_value

  ## Set S3 class
  class(model_fit) <- "lgspline"
  return(model_fit)
}

#' Low-Level Fitting for Lagrangian Smoothing Splines
#'
#' @description
#' The core function for fitting Lagrangian smoothing splines with
#' less user-friendliness.
#'
#' @inheritParams lgspline
#'
#' @keywords internal
#' @export
lgspline.fit <- function(predictors,
                         y = NULL,
                         standardize_response = TRUE,
                         standardize_predictors_for_knots = TRUE,
                         standardize_expansions_for_fitting = TRUE,
                         family = gaussian(),
                         glm_weight_function = function(mu,
                                                        y,
                                                        order_indices,
                                                        family,
                                                        dispersion,
                                                        observation_weights,
                                                        ...){
                           if(any(!is.null(observation_weights))){
                             family$variance(mu) * observation_weights
                           } else {
                             family$variance(mu)
                           }
                         },
                         shur_correction_function = function(X,
                                                             y,
                                                             B,
                                                             dispersion,
                                                             order_list,
                                                             K,
                                                             family,
                                                             observation_weights,
                                                             ...){
                           lapply(1:(K+1), function(k)0)
                         },
                         need_dispersion_for_estimation = FALSE,
                         dispersion_function = function(mu,
                                                        y,
                                                        order_indices,
                                                        family,
                                                        observation_weights,
                                                        ...) { 1 },
                         K = NULL,
                         custom_knots = NULL,
                         cluster_on_indicators = FALSE,
                         make_partition_list = NULL,
                         previously_tuned_penalties = NULL,
                         smoothing_spline_penalty = NULL,
                         opt = TRUE,
                         use_custom_bfgs = TRUE,
                         delta = NULL,
                         tol = 10*sqrt(.Machine$double.eps),
                         log_initial_wiggle = c(-25, 20, -15, -10, -5),
                         log_initial_flat = c(-14, -7),
                         wiggle_penalty = 2e-7,
                         flat_ridge_penalty = 2e-8,
                         unique_penalty_per_partition = TRUE,
                         unique_penalty_per_predictor = TRUE,
                         penalty_ridge = 1e-8,
                         predictor_penalties = NULL,
                         partition_penalties = NULL,
                         include_quadratic_terms = TRUE,
                         include_cubic_terms = TRUE,
                         include_quartic_terms = FALSE,
                         include_2way_interactions = TRUE,
                         include_3way_interactions = TRUE,
                         include_quadratic_interactions = TRUE,
                         just_linear_with_interactions = NULL,
                         just_linear_without_interactions = NULL,
                         exclude_interactions_for = NULL,
                         exclude_these_expansions = NULL,
                         custom_basis_fxn = NULL,
                         include_constrain_fitted = TRUE,
                         include_constrain_first_deriv = TRUE,
                         include_constrain_second_deriv = TRUE,
                         include_constrain_interactions = TRUE,
                         cl = NULL,
                         chunk_size = NULL,
                         parallel_eigen = TRUE,
                         parallel_trace = FALSE,
                         parallel_aga = FALSE,
                         parallel_matmult = FALSE,
                         parallel_unconstrained = FALSE,
                         parallel_find_neighbors = FALSE,
                         parallel_penalty = FALSE,
                         parallel_make_constraint = FALSE,
                         unconstrained_fit_fxn = unconstrained_fit_default,
                         keep_weighted_Lambda = FALSE,
                         iterate_tune = TRUE,
                         iterate_final_fit = TRUE,
                         qp_Amat = NULL,
                         qp_bvec = NULL,
                         qp_meq = 0,
                         qp_monotonic_increase = FALSE,
                         qp_monotonic_decrease = FALSE,
                         qp_range_upper = NULL,
                         qp_range_lower = NULL,
                         qp_Amat_fxn = NULL,
                         qp_bvec_fxn = NULL,
                         qp_meq_fxn = NULL,
                         constraint_value_vectors = cbind(),
                         constraint_vectors = cbind(),
                         return_G = TRUE,
                         return_Ghalf = TRUE,
                         return_U = TRUE,
                         estimate_dispersion = TRUE,
                         return_varcovmat = TRUE,
                         custom_penalty_mat = NULL,
                         cluster_args = c(custom_centers = NA, nstart = 10),
                         dummy_dividor = 0.00000000000000000000012345672152894,
                         dummy_adder = 0.000000000000000002234567210529,
                         verbose = FALSE,
                         verbose_tune = FALSE,
                         expansions_only = FALSE,
                         observation_weights = NULL,
                         do_not_cluster_on_these = c(),
                         neighbor_tolerance = 1 + 1e-16,
                         no_intercept = FALSE,
                         VhalfInv = NULL,
                         include_warnings = TRUE,
                         ...){

  if(verbose){
    cat("Starting\n")
  }

  ## Do not cluster on these should include all linear terms
  if(!is.null(just_linear_with_interactions)){
    do_not_cluster_on_these <- unique(c(do_not_cluster_on_these,
                                 just_linear_with_interactions))
  }
  if(!is.null(just_linear_without_interactions)){
    do_not_cluster_on_these <- unique(c(do_not_cluster_on_these,
                                  just_linear_without_interactions))
  }

  ## Accept raw predictors (the T matrix) and get dimensions
  predictors <- as(predictors,'matrix')
  qcols <- ncol(predictors)
  nr <- nrow(predictors)

  ## Return error message if any terms are > q
  vecdummy <- c(1,
          just_linear_with_interactions,
          just_linear_without_interactions,
          exclude_interactions_for)
  if(any(
    c(1,
      just_linear_with_interactions,
      just_linear_without_interactions,
      exclude_interactions_for) > qcols
    )){
    print(c(1,
            just_linear_with_interactions,
            just_linear_without_interactions,
            exclude_interactions_for))
    stop('Elements in just_linear_with_interactions, just_linear_without_interactions, and/or exclude_interactions_for are greater than the number of columns of predictors matrix.')
  }

  ## Original y vector of response
  y_og <- y

  ## Initialize all variables as numeric by default
  numerics <- 1:qcols

  ## Separate some variables based on desired polynomial expansions
  if(any(is.null(just_linear_with_interactions))){
    just_linear_with_interactions <- c()
  }
  if(any(is.null(just_linear_without_interactions))){
    just_linear_without_interactions <- c()
  }
  if(any(is.null(exclude_interactions_for))){
    exclude_interactions_for <- c()
  }
  numerics <- numerics[!(numerics %in% c(just_linear_with_interactions,
                                         just_linear_without_interactions))]
  intercept <- 1

  ## No interaction terms, set the corresponding options to FALSE
  # only one interaction term available = no interactions
  if(length(exclude_interactions_for) >= (qcols - 1)){
    include_2way_interactions = FALSE
    include_3way_interactions = FALSE
    include_quadratic_interactions = FALSE
  }
  ## With only two terms available for interactions, exclude 3-ways
  if(length(exclude_interactions_for) >= (qcols - 2)){
    include_3way_interactions = FALSE
  }

  if(verbose){
    cat("Polynomial Expansions\n")
  }

  ## Get cubic expansions for design matrix predictors
  C <- get_polynomial_expansions(predictors,
                                 numerics,
                                 just_linear_with_interactions,
                                 just_linear_without_interactions,
                                 exclude_interactions_for,
                                 include_quadratic_terms,
                                 include_cubic_terms,
                                 include_quartic_terms,
                                 include_2way_interactions,
                                 include_3way_interactions,
                                 include_quadratic_interactions,
                                 exclude_these_expansions,
                                 custom_basis_fxn,
                                 ...)

  ## Number of cubic expansions per-partition (little p = nc)
  nc <- ncol(C)

  ## In 1-D, set K = number of constraints given by custom knots
  # if not null
  if(any(!(is.null(custom_knots))) & qcols == 1){
    if(is.null(K)){
      K <- nrow(custom_knots)
    }
  }

  ## Default K
  orig_null <- FALSE
  if(is.null(K)){
    orig_null <- TRUE
    K <- round(max(min(24/(1 +
                           1*(qcols > 1) +
                           1*(paste0(family)[1] != 'gaussian' |
                              paste0(family)[2] != 'identity')),
                        nr/nc),
                 0)/(1 +
                 1*(qcols > 1) +
                 1*(paste0(family)[1] != 'gaussian' |
                    paste0(family)[2] != 'identity')))
  }
  if(K == 0){
    unique_penalty_per_partition <- FALSE
  }

  ## Catch error where we need to cluster on some variables, but we have none
  # allowed
  if(length(do_not_cluster_on_these) == length(qcols) &
     length(qcols) > 1 &
     K > 0){
    stop("Must include at least 1 variable to cluster on if multiple variables are present.")
  }

  ## K can't be greater than max of number of observations or q
  # for kmeans clustering purposes
  if(K >= max(c(nr, qcols))) {
    if(include_warnings){
      warning('Max (N, q) too samll for K. K = max(N, q) - 2 will be used.')
    }
    K <- max(max(c(nr, qcols)) - 2, 0)
  }

  ## Detect if parallel, and K > 0
  if(any(!(is.null(cl))) & K > 0){
    if(any(class(cl) == 'cluster')){
      parallel <- TRUE
      ncores <- length(cl)

      ## if K was not inserted as an argument, multiply minimum by 50
      if(orig_null){
        K <- K*ncores
      }

      ## extract the chunk sizes, number of chunks, and odd-out remaining chunks
      if(is.null(chunk_size)){
        chunk_size <- max(1, ceiling((K + 1) / (4 * ncores)))
      }
      num_chunks <- (K+1) %/% chunk_size
      rem_chunks <- (K+1) %% chunk_size
    } else {
      parallel <- FALSE
    }
  } else{
    parallel <- FALSE
  }

  if(verbose){
    cat("Standardization\n")
  }

  ## Standardize outcome for identity link
  if(paste0(family)[2] == 'identity' &
     paste0(family)[1] == 'gaussian' &
     length(unique(y)) > 1 &
     standardize_response){

    mean_y <- mean(y)
    sd_y <- try(sd(y),silent = TRUE)
    if(any(class(sd_y) == 'try-error')){
      sd_y <- 1
    }
    y <- (y - mean_y)/sd_y
  } else {
    sd_y <- 1
    mean_y <- 0
  }

  ## For cardinal knot-placement,
  # Scale between (0, 1) one-dimension, or standardize N(0,1) higher dim.
  if(standardize_predictors_for_knots){
    if(length(numerics) == 1){
      minns <- apply(predictors, 2, min)
      maxxs <- apply(predictors, 2, max)
      for (j in 1:qcols) {
        predictors[, j] <-
          (predictors[, j] - minns[j] + dummy_adder) /
          (maxxs[j] - minns[j] + dummy_dividor)
      }
    } else {
      means <- apply(predictors, 2, mean)
      sds <- apply(predictors, 2, function(x)tryCatch(sd(x),
                                                      error = function(err)1))
      for (j in 1:qcols) {
        predictors[, j] <-
          (predictors[, j] - means[j] + dummy_adder) /
          (sds[j] + dummy_dividor)
      }
    }
  } else {
    if(length(numerics) == 1){
      maxxs <- rep(1, qcols)
      minns <- rep(0, qcols)
      dummy_adder <- 0
      dummy_divider <- 0
    } else {
      means <- rep(0, qcols)
      sds <- rep(1, qcols)
      dummy_adder <- 0
      dummy_divider <- 0
    }
  }


  ## Transform function for cardinal knot placement
  transf <- function(X) {
    if(length(numerics) == 1){
      for (j in 1:ncol(X)) {
        X[, j] <-
          (X[, j] - minns[j] + dummy_adder) /
          (maxxs[j] - minns[j] + dummy_dividor)
      }
    } else {
      for (j in 1:ncol(X)) {
        X[, j] <-
          (X[, j] - means[j] + dummy_adder) /
          (sds[j] + dummy_dividor)
      }
    }
    X
  }

  ## Inverse transform function for cardinal knot placement
  inv_transf <- function(Xsc) {
    if(length(numerics) == 1){
      for (j in 1:ncol(Xsc)) {
        Xsc[, j] <-
          (Xsc[, j] *
             (maxxs[j] - minns[j] + dummy_dividor) + minns[j] - dummy_adder)
      }
    } else {
      for (j in 1:ncol(Xsc)) {
        Xsc[, j] <-
          (Xsc[, j] *
             (sds[j] + dummy_dividor) + means[j] - dummy_adder)
      }
    }
    Xsc
  }

  if(verbose){
    cat("Get Knots\n")
  }

  if(length(numerics) == 1 & qcols == 1){
    partitions <- NULL

    ## Needed for determining knot locations,
    # compute l1-norms of rows of standardized columns
    partition_codes <- rowMeans(predictors)
    if(any(!(is.null(custom_knots)))){
      ## Scaled by std devs, shift between (0, 1)
      knot_values <- transf(custom_knots)
      partition_bounds <- sort(rowMeans(knot_values))
      kvb <- partition_bounds

    } else if(K > 0){
      ## Knot values at partition_bounds
      kvb <- seq(0,1,length.out = K + 2)[-c(1, K + 2)]
      knot_values <- cbind(kvb)[,rep(1, qcols),drop=FALSE]

      ## Mapped to a single value, since all quantiles are equal
      partition_bounds <- kvb
    }

    ## Compatibility with knot expand function when K = 0
    if(K == 0){
      partition_bounds <- c()
    }
  } else {

    ## If custom knots, replace the partition knots
    if(!any(is.null(make_partition_list))){
      partitions <- make_partition_list

    } else {
      ## Get partitions based on kmeans clustering
      partitions <- make_partitions(predictors,
                                    cluster_args,
                                    cluster_on_indicators,
                                    K,
                                    parallel & parallel_find_neighbors,
                                    cl,
                                    do_not_cluster_on_these,
                                    neighbor_tolerance)
    }
    knot_values <- rbind(partitions$knots)

    ## For compatibility when no knots
    if(K == 0){
      partition_codes <- c()
      partition_bounds <- c()

      ## Code the partitions as an arbitrary monotonic transform
    } else {
      if(verbose){
        cat("Assign partitions\n")
      }
      partition_codes <- partitions$assign_partition(predictors)
      partition_bounds <- 1:nrow(partitions$centers)
    }
  }
  if(verbose){
    cat("Expansion Standardize\n")
  }

  ## Back transform to raw scale,
  # now that knots have been established on standardized scale
  predictors <- inv_transf(predictors)

  ## Index columns of C by variable type for penalization purposes later
  intercept_col <- 1
  colnm_C <- colnames(C)
  power1_cols <- 2:(length(numerics) + 1)
  if(length(numerics) == 0){
    power2_cols <- c()
    power1_cols <- c()
    include_constrain_second_deriv <- FALSE
  } else {
    power2_cols <- which(substr(colnm_C, nchar(colnm_C)-1,
                                nchar(colnm_C)) == '^2')
  }
  power3_cols <- which(substr(colnm_C, nchar(colnm_C)-1,
                              nchar(colnm_C)) == '^3')
  power4_cols <- which(substr(colnm_C, nchar(colnm_C)-1,
                              nchar(colnm_C)) == '^4')
  interaction_cols <- grep("_x_", colnm_C)
  if(length(numerics) > 2 & length(interaction_cols) > 0){
    triplet_cols <- interaction_cols[
      which(sapply(colnm_C[interaction_cols], function(col){
        grepl('_x_',substr(col, regexpr('_x_',col)[[1]]+3,nchar(col)))
      }))]
  } else {
    triplet_cols <- c()
  }
  quad_cols <- which(substr(colnm_C, nchar(colnm_C)-1, nchar(colnm_C)) == "^2")
  interaction_quad_cols <- intersect(
    interaction_cols,quad_cols
  )
  interaction_single_cols <- interaction_cols[!(interaction_cols %in% c(
    triplet_cols, interaction_quad_cols
  ))]

  ## Append non-spline terms
  nonspline_cols <- c(
    which(colnm_C %in% c(paste0("_", just_linear_with_interactions, "_"),
                         paste0("_", just_linear_without_interactions, "_")))
  )
  nonspline_cols <- nonspline_cols[!(nonspline_cols %in%
                                           c(power1_cols,
                                             interaction_single_cols,
                                             interaction_quad_cols,
                                             triplet_cols))]

  ## Standardize columns of C using expansion/(q0.69 - q0.31)
  # This is a p-1 length vector, it excludes the intercept
  C_scales <- apply(C[,-intercept_col,drop=FALSE], 2, function(x){
    1
    if(length(unique(x)) >= 2){
      ## Near sigma for a normal distribution
      # (i.e. this is close to 1 for N(0,1))
      abs(quantile(x, 0.69) - quantile(x, 0.31))
    } else {
      1
    }
  })
  C_scales[C_scales == 0] <- 1
  names(C_scales) <- colnm_C[-intercept_col, drop=FALSE]
  # Set back to 1 if not desired
  if(!standardize_expansions_for_fitting){
    C_scales <- 0*C_scales + 1
  }

  ## Function to un-standardize columns of C
  std_X <- function(unstd_X_in){
    sweep(unstd_X_in, 2, c(1, C_scales), "/")
  }
  max_C <- apply(C, 2, max)
  min_C <- apply(C, 2, min)
  max_min_C <- rbind(c(max_C), c(min_C))
  C <- std_X(C)
  unstd_X <- function(std_X_in){
    sweep(std_X_in, 2, c(1, C_scales), "*")
  }

  ## If no intercept enforced, include constraint on A indicating this
  if(no_intercept & length(constraint_vectors) < 1){
    constr <- sapply(1:(K+1), function(k){
      vec <- rep(0, nc*(K+1))
      vec[nc*(k-1) + 1] <- 1
      vec
    })
    constraint_vectors <- cbind(constr)
    constraint_value_vectors <- 0*constraint_vectors
  } else if(no_intercept){
    constr <- sapply(1:(K+1), function(k){
      vec <- rep(0, nc*(K+1))
      vec[nc*(k-1) + 1] <- 1
      vec
    })
    constraint_vectors <- cbind(constraint_vectors,
                                constr)
    constraint_value_vectors <- cbind(constraint_value_vectors,
                                      0*constr
                                    )
  }

  ## Adjust coefficients after un-standardizing
  backtransform_coefficients <- function(coef) {
    # Extract intercept and slope coefficients
    intercept <- coef[intercept_col]
    slopes <- coef[-intercept_col]

    # Back-transform slope coefficients
    backtransformed_slopes <- slopes / C_scales

    # Combine intercept and back-transformed slopes
    cbind(c(intercept, backtransformed_slopes))
  }

  ## Adjust coefficients for future standardizing
  forwtransform_coefficients <- function(coef) {
    # Extract intercept and slope coefficients
    intercept <- coef[intercept_col]
    slopes <- coef[-intercept_col]

    # Back-transform slope coefficients
    backtransformed_slopes <- slopes / C_scales

    # Combine intercept and back-transformed slopes
    cbind(c(intercept, backtransformed_slopes))
  }

  if(verbose){
    cat("Knot Expand\n")
  }

  ## Get knot expansions
  X <- knot_expand_list(partition_codes,
                        partition_bounds,
                        nr,
                        C,
                        K)

  ## Assign y to their partitions
  y <- knot_expand_list(partition_codes,
                        partition_bounds,
                        nr,
                        cbind(y),
                        K)

  ## If custom variance-covariance structure specified
  if(!is.null(VhalfInv)){
    VhalfInv <- try(as(VhalfInv,'matrix'), silent = TRUE)
    if(any(class(VhalfInv) == 'try-error')){
      if(include_warnings){
        warning('VhalfInv cannot be converted to a N by N matrix, it will not be considered here.')
      }
      VhalfInv <- NULL
    } else if(any(unique(dim(VhalfInv)) != nr)){
      if(include_warnings){
        warning('VhalfInv should be a N by N matrix. It will not be considered here.')
      }
      VhalfInv <- NULL
    } else {
      y_expand_og <- y
      y <- knot_expand_list(partition_codes,
                            partition_bounds,
                            nr,
                            VhalfInv %**% cbind(y_og),
                            K)
    }
  }

  ## Get observation weight expansions
  if(any(!is.null(observation_weights))){
    ## Coerce to N x 1 vector if not already
    if(nrow(cbind(observation_weights)) != nr |
       ncol(cbind(observation_weights)) != 1){
      stop('Observation weights must be an N x 1 vector.')
    }
    observation_weights_og <- observation_weights
    homogenous_weights <- (length(unique(observation_weights_og)) == 1)
    observation_weights <-
      knot_expand_list(partition_codes,
                       partition_bounds,
                       nr,
                       cbind(observation_weights),
                       K)
  } else {
    observation_weights_og <- rep(1, nr)
    observation_weights <- lapply(1:(K+1), function(k)cbind(rep(1,
                                                        length(y[[k]]))))
    homogenous_weights <- TRUE
  }

  ## Save the original ordering to each partition
  order_list <- knot_expand_list(partition_codes,
                                 partition_bounds,
                                 nr,
                                 cbind(1:nr),
                                 K)
  og_order <- order(unlist(order_list))

  ## Return derivatives per-partition of an expanded matrix
  all_derivatives <- function(X){
    lapply(X, function(C){
      make_derivative_matrix(
        nc,
        C,
        power1_cols,
        interaction_single_cols,
        interaction_quad_cols,
        triplet_cols,
        K,
        include_2way_interactions,
        include_3way_interactions,
        include_quadratic_interactions,
        colnm_C,
        C_scales
      )
    })
  }

  if(verbose){
    cat("2nd Derivative Penalty\n")
  }


  ## Compute integrated squared second derivative of fitted function
  # evaluated over bounds of the support
  # can be replaced with arbitrary p by p matrix if desired
  if(!(!(any(is.null(smoothing_spline_penalty))))){
    ## Compute the gram matrix for the squared integrated second derivative
    # Standard standardization
    max_min_C <- std(max_min_C)

    smoothing_spline_penalty <-
      get_2ndDerivPenalty_wrapper(K,
                                  colnm_C,
                                  max_min_C,
                                  power1_cols,
                                  power2_cols,
                                  power3_cols,
                                  power4_cols,
                                  interaction_single_cols,
                                  interaction_quad_cols,
                                  triplet_cols,
                                  nonspline_cols,
                                  nc,
                                  parallel & parallel_penalty,
                                  cl)
    colnames(smoothing_spline_penalty) <- colnames(C)
  }

  if(verbose){
    cat("Constraint Matrix\n")
  }


  ## Making a constraint matrix
  A <- 0
  if((K > 0 & length(numerics) == 1 & length(nonspline_cols) == 0)){

    ## Knot basis expansions
    CKnots <- get_polynomial_expansions(inv_transf(cbind(knot_values)),
                                        numerics,
                                        just_linear_with_interactions,
                                        just_linear_without_interactions,
                                        exclude_interactions_for,
                                        include_quadratic_terms,
                                        include_cubic_terms,
                                        include_quartic_terms,
                                        include_2way_interactions,
                                        include_3way_interactions,
                                        include_quadratic_interactions,
                                        exclude_these_expansions,
                                        custom_basis_fxn,
                                        ...)
    if(K == 1){
      CKnots <- rbind(CKnots)
    }
    if(nrow(CKnots) < K){
      CKnots <- rbind(CKnots, matrix(0, K - nrow(CKnots), ncol = ncol(CKnots)))
    }

    ## Constraint matrix A
    A <- make_constraint_matrix(nc,
                                CKnots,
                                power1_cols,
                                nonspline_cols,
                                interaction_single_cols,
                                interaction_quad_cols,
                                triplet_cols,
                                K,
                                include_constrain_fitted,
                                include_constrain_first_deriv,
                                include_constrain_second_deriv,
                                include_constrain_interactions,
                                include_2way_interactions,
                                include_3way_interactions,
                                include_quadratic_interactions,
                                colnm_C,
                                C_scales)
    ## apply standardization to the rows of A,
    # once constraints un-standardized are derived

    if(length(constraint_vectors) > 0 & length(constraint_value_vectors > 0)){
      A <- cbind(A, constraint_vectors)
    }
    A <- sweep(A, 1, rep(c(1, C_scales), K+1), "/")
    if(any(!is.finite(A))) stop(paste0('A is not finite', C_scales))

    ## Otherwise, if we do have knots.....
  } else if(K > 0){


    ## how many chunks/individual matrices will A be composed of
    chunk <- nrow(knot_values) %/% K
    rem <- nrow(knot_values) %% K # don't forget straggling rows

    ## permute knot values
    knot_values_perm <- knot_values[1:nrow(knot_values),,drop=FALSE]

    if(parallel & parallel_make_constraint){
      A <- Reduce("cbind",
            parLapply(cl,
             1:chunk,
             function(i){
               ## Select the knot values in the chunk
               knot_values_chunk <- knot_values_perm[1:K + (i-1)*K,,drop=FALSE]

               ## Get polynomial expansions of knot quantile values
               CKnots_chunk <- rbind(get_polynomial_expansions(
                 inv_transf(knot_values_chunk),
                 numerics,
                 just_linear_with_interactions,
                 just_linear_without_interactions,
                 exclude_interactions_for,
                 include_quadratic_terms,
                 include_cubic_terms,
                 include_quartic_terms,
                 include_2way_interactions,
                 include_3way_interactions,
                 include_quadratic_interactions,
                 exclude_these_expansions,
                 custom_basis_fxn,
                 ...))
               rownames(CKnots_chunk) <- rownames(knot_values_chunk)

               ## Constraint matrix A
               make_constraint_matrix(nc,
                                      CKnots_chunk,
                                      power1_cols,
                                      nonspline_cols,
                                      interaction_single_cols,
                                      interaction_quad_cols,
                                      triplet_cols,
                                      K,
                                      include_constrain_fitted,
                                      include_constrain_first_deriv,
                                      include_constrain_second_deriv,
                                      include_constrain_interactions,
                                      include_2way_interactions,
                                      include_3way_interactions,
                                      include_quadratic_interactions,
                                      colnm_C,
                                      C_scales)
                                   }))
    } else {

      for(i in 1:chunk){

        ## Permute knot_quantile_value_combinations
        knot_values_chunk <- knot_values_perm[1:K + (i-1)*K,,drop=FALSE]

        ## Get polynomial expansions of knot quantile values
        CKnots_chunk <- rbind(
          get_polynomial_expansions(inv_transf(knot_values_chunk),
                                    numerics,
                                    just_linear_with_interactions,
                                    just_linear_without_interactions,
                                    exclude_interactions_for,
                                    include_quadratic_terms,
                                    include_cubic_terms,
                                    include_quartic_terms,
                                    include_2way_interactions,
                                    include_3way_interactions,
                                    include_quadratic_interactions,
                                    exclude_these_expansions,
                                    custom_basis_fxn,
                                    ...))
        rownames(CKnots_chunk) <- rownames(knot_values_chunk)

        ## Constraint matrix A
        A <- cbind(A, make_constraint_matrix(nc,
                                             CKnots_chunk,
                                             power1_cols,
                                             nonspline_cols,
                                             interaction_single_cols,
                                             interaction_quad_cols,
                                             triplet_cols,
                                             K,
                                             include_constrain_fitted,
                                             include_constrain_first_deriv,
                                             include_constrain_second_deriv,
                                             include_constrain_interactions,
                                             include_2way_interactions,
                                             include_3way_interactions,
                                             include_quadratic_interactions,
                                             colnm_C,
                                             C_scales))
        if(i == 1){
          ## Remove 0 column
          A <- A[,-1,drop=FALSE]
        }
      }
    }
    if(rem > 0){
      ## Permute knot_quantile_value_combinations
      knot_values_chunk <-
        knot_values_perm[rev(c(nrow(knot_values_perm):1)[1:rem]),,drop=FALSE]


      ## Get polynomial expansions of knot quantile values
      temp_dat <- inv_transf(knot_values_chunk)
      only_1 <- FALSE
      if(nrow(temp_dat) == 1){
        only_1 <- TRUE
        temp_dat <- rbind(temp_dat, temp_dat)
      }
      CKnots_chunk <- rbind(
        get_polynomial_expansions(temp_dat,
                                  numerics,
                                  just_linear_with_interactions,
                                  just_linear_without_interactions,
                                  exclude_interactions_for,
                                  include_quadratic_terms,
                                  include_cubic_terms,
                                  include_quartic_terms,
                                  include_2way_interactions,
                                  include_3way_interactions,
                                  include_quadratic_interactions,
                                  exclude_these_expansions,
                                  custom_basis_fxn,
                                  ...))
      if(only_1){
        CKnots_chunk <- CKnots_chunk[1,,drop=FALSE]
      }
      rownames(CKnots_chunk) <- rownames(knot_values_chunk)
      dummy <- matrix(0, nrow = K - rem, ncol = ncol(CKnots_chunk))
      rownames(dummy) <- paste0(sample(1:nrow(dummy)), '_', 2:(nrow(dummy)+1))
      CKnots_chunk <- rbind(CKnots_chunk, dummy)

      ## Constraint matrix A
      A <- cbind(A, make_constraint_matrix(nc,
                                           CKnots_chunk,
                                           power1_cols,
                                           nonspline_cols,
                                           interaction_single_cols,
                                           interaction_quad_cols,
                                           triplet_cols,
                                           K,
                                           include_constrain_fitted,
                                           include_constrain_first_deriv,
                                           include_constrain_second_deriv,
                                           include_constrain_interactions,
                                           include_2way_interactions,
                                           include_3way_interactions,
                                           include_quadratic_interactions,
                                           colnm_C,
                                           C_scales))
    }

    ## Remove all 0 columns
    A <- A[,which(apply(abs(A), 2, sum) > 1e-16),drop=FALSE]

    ## Bind other constraints, standardize
    if(length(constraint_vectors) > 0){
      A <- cbind(A, constraint_vectors)
    }
    A <- sweep(A, 1, rep(c(1, C_scales), K+1), "/")
    if(any(!is.finite(A))) stop(paste0('A is not finite', C_scales))
    if(any(is.na(A))) stop(paste0('A is na somewhere', C_scales))

  } else {
    ## If missing constraints, apply custom constraints if desired only
    ## or do not include A at all
    if(length(constraint_vectors) > 0){
      A <- cbind(constraint_vectors)
      A <- sweep(A, 1, rep(c(1, C_scales), K+1), "/")
    } else {
      A <- NULL
    }
  }
  if(!(any(is.null(A)))){
    nca <- ncol(A)
  }

  ## Convert non-0 null vectors to (K+1) list of corresponding partitions
  # Adjust for intercept being shifted by mean y
  if(length(constraint_value_vectors) > 0){

    constraint_value_vectors <- lapply(1:(K+1),function(k){
      vec <- cbind(constraint_value_vectors)[1:nc + (k-1)*nc,,drop=FALSE]
      vec[1,] <- (vec[1,] - mean_y)/sd_y
      vec * c(1, C_scales)
    })

  }

  ## With only one predictor, we really only need one penalty
  if(ncol(predictors) == 1){
    unique_penalty_per_predictor <- FALSE
  }

  if(verbose){
    cat("Penalty and SQP Setup\n")
  }


  ## Getting unique penalties for predictors/partitions, if not specified
  log_penalty_vec <- c()
  if(unique_penalty_per_predictor & any(is.null(predictor_penalties))){

    ## Initialize
    predictor_penalties <- sapply(colnm_C[c(power1_cols,
                                            nonspline_cols)],
                                  function(j)rnorm(1, 0, 0.00001))
    names(predictor_penalties) <- paste0('predictor',
                                         colnm_C[c(power1_cols,
                                                   nonspline_cols)])
    log_penalty_vec <- c(log_penalty_vec, predictor_penalties)

  } else if(unique_penalty_per_predictor){
    if(length(unique_penalty_per_predictor) !=
       length(c(power1_cols, nonspline_cols))){
      stop('Custom predictor_penalties are not the same length as number of predictors in model.')
    }
    names(predictor_penalties) <- paste0('predictor',
                                         colnm_C[c(power1_cols,
                                                   nonspline_cols)])
    log_penalty_vec <- c(log_penalty_vec, predictor_penalties)
  }
  if(unique_penalty_per_partition & any(is.null(partition_penalties))){
    ## Initialize
    partition_penalties <- sapply(1:(K+1),
                                  function(j)rnorm(1, 0, 0.00001))
    names(partition_penalties) <- paste0('partition', 1:(K+1))
    log_penalty_vec <- c(log_penalty_vec, partition_penalties)
  } else if(unique_penalty_per_partition){
    names(partition_penalties) <- paste0('partition', 1:(K+1))
    log_penalty_vec <- c(log_penalty_vec, partition_penalties)
  }

  ## Update quadprog variable, if the correct arguments are made
  if(qp_monotonic_decrease |
     qp_monotonic_increase |
     any(!(is.null(qp_range_upper))) |
     any(!(is.null(qp_range_lower))) |
     (any(!is.null(qp_Amat_fxn)) &
      any(!is.null(qp_bvec_fxn)) &
      any(!is.null(qp_meq_fxn))) |
     (any(!is.null(qp_Amat)) &
      any(!is.null(qp_bvec)) &
      any(!is.null(qp_meq)))){
    quadprog <- TRUE
  } else {
    quadprog <- FALSE
  }

  if(verbose){
    cat("Parallel Setup\n")
  }

  ## Export components for parallel processing
  if(parallel & !is.null(cl)) {
    library(parallel)

    ## Create shared environment in global environment
    assign("shared_env", new.env(), envir = .GlobalEnv)

    ## Assign key variables to shared environment
    shared_vars <- list(
      A = A,
      nca = ncol(A),
      K = K,
      nc = nc,
      snr = sqrt(nr),
      chunk_size = chunk_size,
      num_chunks = num_chunks,
      rem_chunks = rem_chunks,
      parallel = parallel,
      X = X,
      log_penalty_vec = log_penalty_vec,
      unique_penalty_per_partition = unique_penalty_per_partition,
      keep_weighted_Lambda = keep_weighted_Lambda,
      custom_penalty_mat = custom_penalty_mat,
      glm_weight_function = glm_weight_function,
      shur_correction_function = shur_correction_function,
      unconstrained_fit_fxn = unconstrained_fit_fxn,
      observation_weights = observation_weights
    )

    for(nm in names(shared_vars)) {
      assign(nm, shared_vars[[nm]], envir = shared_env)
    }

    ## Export shared environment
    clusterExport(cl, "shared_env")

    ## Setup each cluster node with necessary functions
    clusterEvalQ(cl, {
      library(lgspline)
      `%**%` <- efficient_matrix_mult
      ## Load shared environment objects
      list2env(as.list(shared_env), .GlobalEnv)
    })
  } else {
    shared_env <- NULL
  }

  ## X^{T}WX
  ## Account for weights
  if(!is.null(VhalfInv)){
    X_og <- X
    ord <- unlist(order_list)
    VhalfInvX <- VhalfInv[ord,ord] %**% collapse_block_diagonal(X)
    VhalfInvX <- VhalfInvX[og_order,,drop=FALSE]
    X <- lapply(1:(K+1), function(k){
      VhalfInvX[order_list[[k]], (k-1)*nc +1:nc]
    })
  }
  if(( (paste0(family)[1] == 'gaussian' &
        paste0(family)[2] == 'identity')) &
     !homogenous_weights){
    X <- lapply(1:(K+1), function(k){
      X[[k]] * c(sqrt(observation_weights[[k]]))
    })
  }
  X_gram <- compute_gram_block_diagonal(X,
                                        parallel,
                                        cl,
                                        chunk_size,
                                        num_chunks,
                                        rem_chunks)
  ## Switch back after
  if(( (paste0(family)[1] == 'gaussian' &
        paste0(family)[2] == 'identity')) &
     !homogenous_weights){
    X <- lapply(1:(K+1), function(k){
      X[[k]] / c(sqrt(observation_weights[[k]]))
    })
  }
  if(!is.null(VhalfInv)){
    X <- X_og
    rm(X_og)
  }

  ## Quadprog setup
  if(quadprog){

    ## Initialize empty constraint lists
    qp_Amat_list <- list()
    qp_bvec_list <- list()
    qp_meq_list <- list()

    ## Big-matrix components (not memory efficient anymore)
    X_block <- Reduce("rbind", lapply(1:(K+1), function(k){
      dummy <- 0*X[[k]]
      Reduce("cbind",lapply(1:(K+1),function(j){
        if(nrow(X[[k]]) == 0){
          return(X[[k]])
        } else if(j == k) X[[k]] else 0*X[[k]]
      }))
    }))

    ## Constraints on range of fitted values
    if(!(any(is.null(qp_range_upper))) | !any(is.null(qp_range_lower))){

      ## Both upper and lower
      if(!(any(is.null(qp_range_upper))) & !any(is.null(qp_range_lower))){
        qp_Amat <- cbind(t(X_block), -t(X_block))
        if(length(qp_range_lower) == 1){
          ## If only single-bounds are given,
          # then use unique values of qp Amatrix
          if(length(qp_range_upper == 1)){
            qp_Amat <- t(unique(t(qp_Amat)))
          }
          qp_bvec <- rep(qp_range_lower, ncol(qp_Amat)/2)
        } else {
          qp_bvec <- qp_range_lower
        }
        ## Don't forget y is standardized
        qp_bvec_lower <- (qp_bvec - mean_y) / sd_y

        if(length(qp_range_upper) == 1){
          qp_bvec <- rep(qp_range_upper, ncol(qp_Amat)/2)
        } else {
          qp_bvec <- qp_range_upper
        }
        ## Don't forget y is standardized
        qp_bvec_upper <- -(qp_bvec - mean_y) / sd_y

        ## Combine
        qp_bvec <- c(qp_bvec_lower, qp_bvec_upper)

        ## Append
        qp_Amat_list[[length(qp_Amat_list) + 1]] <- qp_Amat
        qp_bvec_list[[length(qp_bvec_list) + 1]] <- qp_bvec
        qp_meq_list[[length(qp_meq_list) + 1]] <- 0

        ## Just upper
      } else if(!(any(is.null(qp_range_upper)))){
        qp_Amat <- -t(X_block)
        if(length(qp_range_upper) == 1){
          qp_Amat <- t(unique(t(qp_Amat)))
          qp_bvec <- rep(qp_range_upper, ncol(qp_Amat))
        } else {
          qp_bvec <- qp_range_upper
        }
        ## Don't forget y is standardized
        qp_bvec <- -(qp_bvec - mean_y) / sd_y

        ## Append
        qp_Amat_list[[length(qp_Amat_list) + 1]] <- qp_Amat
        qp_bvec_list[[length(qp_bvec_list) + 1]] <- qp_bvec
        qp_meq_list[[length(qp_meq_list) + 1]] <- 0

        ## Just lower
      } else if(!any(is.null(qp_range_lower))){
        qp_Amat <- t(X_block)
        if(length(qp_range_lower) == 1){
          qp_Amat <- t(unique(t(qp_Amat)))
          qp_bvec <- rep(qp_range_lower, ncol(qp_Amat))
        } else {
          qp_bvec <- qp_range_lower
        }

        ## Don't forget y is standardized
        qp_bvec <- (qp_bvec - mean_y) / sd_y

        ## Append
        qp_Amat_list[[length(qp_Amat_list) + 1]] <- qp_Amat
        qp_bvec_list[[length(qp_bvec_list) + 1]] <- qp_bvec
        qp_meq_list[[length(qp_meq_list) + 1]] <- 0
      }

    }

    ## Monotonic increasing constraint
    if(qp_monotonic_increase){

      ## First, create constraints for fitted values
      value_constraints <- t(Reduce('rbind', lapply(2:nr, function(i) {
        matrix(c(X_block[i,] - X_block[i-1,]), nrow = 1)
      })))

      ## Compute first derivative matrix
      derivs <- make_derivative_matrix(
        nc,  # number of columns
        X_block,  # design matrix
        power1_cols,  # linear term columns
        interaction_single_cols,  # single interaction columns
        interaction_quad_cols,  # quadratic interaction columns
        triplet_cols,  # triplet interaction columns
        K,  # number of knots
        include_2way_interactions,
        include_3way_interactions,
        include_quadratic_interactions,
        colnm_C,  # column names
        C_scales,  # scaling
        just_first_derivatives = TRUE
      )

      ## Extract first derivatives for each variable
      first_derivative_constraints <- Reduce("rbind",
       lapply(derivs$first_derivative, function(deriv_matrix) {
         ## Ensure non-negative first derivatives for monotonic increasing
         t(Reduce('rbind', lapply(1:nrow(deriv_matrix), function(i) {
           matrix(c(deriv_matrix[i,]), nrow = 1)  # enforce non-negativity
         })))
       })
      )

      ## Combine value and derivative constraints
      qp_Amat <- cbind(value_constraints, first_derivative_constraints)
      qp_Amat <- t(unique(t(qp_Amat)))
      qp_bvec <- rep(0, ncol(qp_Amat))

      ## Append
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- qp_Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- qp_bvec
      qp_meq_list[[length(qp_meq_list) + 1]] <- 0

      ## Monotonic decreasing constraint
    } else if(qp_monotonic_decrease){

      ## First, create constraints for fitted values
      value_constraints <- -t(Reduce('rbind', lapply(2:nr, function(i) {
        matrix(c(X_block[i,] - X_block[i-1,]), nrow = 1)
      })))

      ## Compute first derivative matrix
      derivs <- make_derivative_matrix(
        nc,  # number of columns
        X_block,  # design matrix
        power1_cols,  # linear term columns
        interaction_single_cols,  # single interaction columns
        interaction_quad_cols,  # quadratic interaction columns
        triplet_cols,  # triplet interaction columns
        K,  # number of knots
        include_2way_interactions,
        include_3way_interactions,
        include_quadratic_interactions,
        colnm_C,  # column names
        C_scales,  # scaling
        just_first_derivatives = TRUE
      )

      ## Extract first derivatives for each variable
      first_derivative_constraints <- Reduce("rbind",
       lapply(derivs$first_derivative, function(deriv_matrix) {
         ## Ensure non-positive first derivatives for monotonic decreasing
         t(Reduce('rbind', lapply(1:nrow(deriv_matrix), function(i) {
           -matrix(c(deriv_matrix[i,]), nrow = 1)  # enforce non-positive
         })))
       })
      )

      ## Combine value and derivative constraints
      qp_Amat <- cbind(value_constraints, first_derivative_constraints)
      qp_Amat <- t(unique(t(qp_Amat)))
      qp_bvec <- rep(0, ncol(qp_Amat))

      ## Append
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- qp_Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- qp_bvec
      qp_meq_list[[length(qp_meq_list) + 1]] <- 0

    }
    ## Custom constraints
    if(!is.null(qp_Amat_fxn) &
       !is.null(qp_bvec_fxn) &
       !is.null(qp_meq_fxn)){
      qp_Amat <- qp_Amat_fxn(nr, nc, K, X_block, colnm_C, C_scales)
      qp_bvec <- qp_bvec_fxn(qp_Amat, nr, nc, K, X_block, colnm_C, C_scales)
      qp_meq <- qp_meq_fxn(qp_Amat, nr, nc, K, X_block, colnm_C, C_scales)

      ## Append
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- qp_Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- qp_bvec
      qp_meq_list[[length(qp_meq_list) + 1]] <- qp_meq
    }

    ## Combined and overwrite constraint matrices/vectors
    qp_Amat <- do.call(cbind, qp_Amat_list)
    qp_bvec <- do.call(c, qp_bvec_list)
    qp_meq <- sum(unlist(qp_meq_list))

    ## Reduce memory constraint
    X_block <- NULL
  }

  ## Return basis-expansions only and partitioned y without model fitting
  # and associated components
  if(expansions_only){
    return(list(
      X = X,
      y = y,
      A = A,
      penalties  =   compute_Lambda(custom_penalty_mat,
                                    smoothing_spline_penalty,
                                    wiggle_penalty,
                                    flat_ridge_penalty,
                                    K,
                                    nc,
                                    unique_penalty_per_predictor,
                                    unique_penalty_per_partition,
                                    exp(log_penalty_vec),
                                    colnm_C,
                                    just_Lambda = FALSE),
      order_list = order_list,
      og_order = og_order,
      C_scales = C_scales,
      colnm_C = colnm_C,
      K = K,
      knots = knots,
      partitions = partitions,
      partition_codes = partition_codes,
      partition_bounds = partition_bounds
    ))
  }

  if(verbose){
    cat("Incorporate Variance Structure\n")
  }

  ## Incorporate user-specified variance-covariance structure on response,
  # if desired
  if(!is.null(VhalfInv)){
    X_og <- X
    ord <- unlist(order_list)
    VhalfInvX <- VhalfInv[ord,ord] %**% collapse_block_diagonal(X)
    VhalfInvX <- VhalfInvX[og_order,,drop=FALSE]
    X <- lapply(1:(K+1), function(k){
      VhalfInvX[order_list[[k]], (k-1)*nc +1:nc]
    })
  }

  if(verbose){
    cat("Tune Smoothing Spline Penalty\n")
  }

  ## This is to incorporate weights efficiently for linear regression outcomes
  # Remember to back-transform later
  if(( (paste0(family)[1] == 'gaussian' &
        paste0(family)[2] == 'identity')) &
     !homogenous_weights){
    X <- lapply(1:(K+1), function(k){
       X[[k]] * c(sqrt(observation_weights[[k]]))
    })
  }
  if(((paste0(family)[1] == 'gaussian' &
       paste0(family)[2] == 'identity')) &
     !homogenous_weights){
    y <- lapply(1:(K+1), function(k){
      y[[k]] * c(sqrt(observation_weights[[k]]))
    })
  }

  ## Model components
  if(!(!any(is.null(previously_tuned_penalties)))){

    ## Prior precision Lambda
    tL <- tune_Lambda(
      y,
      X,
      X_gram,
      smoothing_spline_penalty,
      A,
      K,
      nc,
      nr,
      opt,
      use_custom_bfgs,
      C,
      colnm_C,
      wiggle_penalty,
      flat_ridge_penalty,
      log_initial_wiggle,
      log_initial_flat,
      unique_penalty_per_predictor,
      unique_penalty_per_partition,
      log_penalty_vec,
      penalty_ridge,
      family,
      unconstrained_fit_fxn,
      keep_weighted_Lambda,
      iterate_tune,
      quadprog,
      qp_Amat,
      qp_bvec,
      qp_meq,
      tol,
      sd_y,
      delta,
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
      shur_correction_function,
      need_dispersion_for_estimation,
      dispersion_function,
      observation_weights,
      homogenous_weights,
      verbose_tune,
      ...)
  } else {

    ## Use previously-submitted Lambda
    tL <- previously_tuned_penalties
    rm(previously_tuned_penalties)

  }
  flat_ridge_penalty <- tL$flat_ridge_penalty
  wiggle_penalty <- tL$wiggle_penalty

  if(verbose){
    cat("Prep for final fitting\n")
  }

  ## Final fit
  if(K == 0){
    ## ensuring compatibility with no A
    if(any(is.null(A))){
      ## for compatibility, albeit inefficient
      A <- cbind(rep(0, (K+1)*nc))
      A <- cbind(A, A)
      nca <- 2
    }
  }
  Xy <- vectorproduct_block_diagonal(X, y, K)
  shur_corrections <- lapply(1:(K+1), function(k)0)
  G_list <- compute_G_eigen(X_gram,
                            tL$Lambda,
                            K,
                            parallel & parallel_eigen,
                            cl,
                            chunk_size,
                            num_chunks,
                            rem_chunks,
                            family,
                            unique_penalty_per_partition,
                            tL$L_partition_list,
                            keep_G = (return_G |
                                      return_U |
                                      estimate_dispersion |
                                      return_varcovmat),
                            shur_corrections)

  if(verbose){
    cat('Last fit\n')
  }
  ## Get coefficient and correlation matrix estimates
  return_G_getB <- TRUE
  B_list <-  try({get_B(
              X,
              X_gram,
              tL$Lambda,
              keep_weighted_Lambda,
              unique_penalty_per_partition,
              tL$L_partition_list,
              A,
              Xy,
              y,
              K,
              nc,
              nca,
              G_list$Ghalf,
              G_list$GhalfInv,
              parallel & parallel_eigen,
              parallel & parallel_aga,
              parallel & parallel_matmult,
              parallel & parallel_unconstrained,
              cl,
              chunk_size,
              num_chunks,
              rem_chunks,
              family,
              unconstrained_fit_fxn,
              iterate_final_fit,
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
              shur_correction_function,
              need_dispersion_for_estimation,
              dispersion_function,
              observation_weights,
              homogenous_weights,
              return_G_getB,
              ...)}, silent = TRUE)
            if(any(class(B_list) == 'try-error')){
              print(B_list)
              stop('Failure in fitting final model')
            }
  B <- B_list$B
  G_list <- B_list$G_list

  ## This is backtransforming from earlier,
  # if we have Gaussian weighted response
  if(( (paste0(family)[1] == 'gaussian' &
        paste0(family)[2] == 'identity')) &
     !homogenous_weights){
    X <- lapply(1:(K+1), function(k){
      X[[k]] / c(sqrt(observation_weights[[k]]))
    })
  }
  if(((paste0(family)[1] == 'gaussian' &
       paste0(family)[2] == 'identity')) &
     !homogenous_weights){
    y <- lapply(1:(K+1), function(k){
      y[[k]] / c(sqrt(observation_weights[[k]]))
    })
  }

  ## Get original design matrix now, and original y expansions, after fitting
  if(!is.null(VhalfInv)){
    X <- X_og
    y <- y_expand_og
    rm(X_og)
  }

  if(verbose){
    cat("After Fitting Processing \n")
  }

  ## For assigning out of data clusters
  if(length(c(numerics, nonspline_cols)) > 1 & K > 0){
    assign_partition <- partitions$assign_partition
  } else if(length(numerics) == 1 & K > 0){
    assign_partition <- function(x)rowMeans(cbind(x))
  } else {
    assign_partition <- function(x)0.5
  }

  ## Raw coefficients, useful for incorporation into Bayesian techniques
  B_raw <- B

  ## Un-scale, based on centered-and-scaled y
  B <- lapply(B, function(b)b * sd_y) # multiply all by sd of y

  ## Then add mean of y to all intercepts
  B <- lapply(1:(K+1), function(k){
    b <- B[[k]]
    b[1] <-
      b[1] + mean_y
    b
  })

  ## Rename B coefficients for interpretability,
  # adjust for unstandardized predictors
  B <- lapply(1:(K+1),function(k){
    B[[k]] <- backtransform_coefficients(B[[k]])
    names(B[[k]]) <- colnm_C
    B[[k]]
  })
  names(B) <- paste0('partition',1:(K+1))

  ## Predict function for new data
  predict_function <- function(new_predictors = predictors,
                               parallel = FALSE,
                               cl = NULL,
                               chunk_size = NULL,
                               num_chunks = NULL,
                               rem_chunks = NULL,
                               B_predict = B,
                               take_first_derivatives = FALSE,
                               take_second_derivatives = FALSE,
                               just_expansions = FALSE){

      ## Check compatibility, that new_predictors should be a matrix
      if(any(!is.null(new_predictors))){
        new_predictors <- try(as(cbind(new_predictors), 'matrix'), silent = TRUE)
        if(any(class(new_predictors) == 'try-error')){
          stop('New predictors should be able to be coerced into matrix form.')
        }
      }

      ## Avoid R rbind issue with 1 row only for certain internal functions
      if(nrow(new_predictors) == 1){
        new_predictors <- rbind(new_predictors, new_predictors)
        only_1 <- TRUE
      } else {
        only_1 <- FALSE
      }

      ## Accept predictors as matrix
      new_predictors <- transf(as(new_predictors, 'matrix'))


      ## Needed for determining knot locations,
      # compute l2-norms of rows of standardized columns
      partition_codes_new <- assign_partition(new_predictors)


      ## Back transform, now that knots have been established
      new_predictors <- inv_transf(new_predictors)


      ## Cubic expansions
      C_new <- get_polynomial_expansions(new_predictors,
                                         numerics,
                                         just_linear_with_interactions,
                                         just_linear_without_interactions,
                                         exclude_interactions_for,
                                         include_quadratic_terms,
                                         include_cubic_terms,
                                         include_quartic_terms,
                                         include_2way_interactions,
                                         include_3way_interactions,
                                         include_quadratic_interactions,
                                         exclude_these_expansions,
                                         custom_basis_fxn,
                                         ...)

        ## Knot expansions
        X_new <- knot_expand_list(
          partition_codes_new,
          partition_bounds,
          length(partition_codes_new),
          C_new,
          K
        )

        ## If just the expansions are desired
        if(just_expansions){
          if(only_1){
            partition_codes_new <- partition_codes_new[1]
            C_new <- C_new[1, , drop=FALSE]
            X_new <- lapply(X_new,function(x){
              if(!any(is.null(x))){
                if(nrow(x) == 2){
                  return(x[1,,drop=FALSE])
                } else {
                  x
                }
              } else{
                x
              }
            })
          }
          return(list("expansions" = X_new,
                      "partition_codes" = partition_codes_new,
                      "partition_bounds" = partition_bounds))
        }

        ## Re-order predictions after
        order_list <- knot_expand_list(
          partition_codes_new,
          partition_bounds,
          length(partition_codes_new),
          cbind(1:nrow(C_new)),
          K)

        ## Only use relevant blocks
        keep_blocks <- which(sapply(1:(K+1),function(k){
          nrow(X_new[[k]]) > 0
        }))
        order_list <- order_list[keep_blocks]

        ## Predictions
        preds <-
          unlist(
            matmult_block_diagonal(
              X_new[keep_blocks],
              B_predict[keep_blocks],
              length(keep_blocks) - 1,
              parallel,
              cl,
              chunk_size,
              num_chunks,
              rem_chunks))[order(unlist(order_list))]
        if(only_1){
          preds <- preds[1]
        }
        final_preds <- family$linkinv(preds)

        ## If returning derivatives
        if(take_first_derivatives | take_second_derivatives){
          derivs <- make_derivative_matrix(
            nc,
            C_new,
            power1_cols,
            interaction_single_cols,
            interaction_quad_cols,
            triplet_cols,
            K,
            include_2way_interactions,
            include_3way_interactions,
            include_quadratic_interactions,
            colnm_C,
            C_scales,
            !take_second_derivatives)

          if(only_1){
            partition_codes_new <- partition_codes_new[1]
          }


          if(take_first_derivatives){

            Cprime_new <- Reduce("rbind",
                                 lapply(1:length(derivs$first_derivative),
                                        function(var){
                                          d <- derivs$first_derivative[[var]]
                                          if(only_1){
                                            return(d[1,,drop=FALSE])
                                          } else{
                                            return(d)
                                          }
                                        }))


            ## Knot expansions
            Xprime_new <- knot_expand_list(
              partition_codes_new,
              partition_bounds,
              length(partition_codes_new),
              Cprime_new,
              K
            )

            ## Derivative of predictions
            preds_prime <-
              unlist(
                matmult_block_diagonal(
                  Xprime_new[keep_blocks],
                  B_predict[keep_blocks],
                  length(keep_blocks) - 1,
                  parallel,
                  cl,
                  chunk_size,
                  num_chunks,
                  rem_chunks))

            final_preds_prime <- family$linkinv(preds_prime)
          } else {
            final_preds_prime <- NULL
          }

          ## If returning second derivatives
          if(take_second_derivatives){

            Cdprime_new <- Reduce("rbind",
                                  lapply(1:length(derivs$first_derivative),
                                         function(var){
                                           d <- derivs$second_derivative[[var]]
                                           if(only_1){
                                             return(d[1,,drop=FALSE])
                                           } else{
                                             return(d)
                                           }
                                         }))


            Xdprime_new <- knot_expand_list(
              partition_codes_new,
              partition_bounds,
              length(partition_codes_new),
              Cdprime_new,
              K
            )

            ## Second derivative of predictions
            preds_dprime <-
              unlist(
                matmult_block_diagonal(
                  Xdprime_new[keep_blocks],
                  B_predict[keep_blocks],
                  length(keep_blocks) - 1,
                  parallel,
                  cl,
                  chunk_size,
                  num_chunks,
                  rem_chunks))

            final_preds_dprime <- family$linkinv(preds_dprime)
          } else {
            final_preds_dprime <- NULL
          }

          return(list(
            preds = final_preds,
            first_deriv = final_preds_prime,
            second_deriv = final_preds_dprime
          ))

        } else {
          return(final_preds)
        }
  }

  ## Get fitted values
  ytilde <- predict_function()

  ## Clean knots, back transform to raw-scale
  if(K == 0){
    knots <- NULL
  } else {
    knots <- inv_transf(knot_values)
    if(length(numerics) == 1 & length(nonspline_cols) == 0){
      rownames(knots) <- paste0(1:K, '_', 2:(K+1))
    }
  }

  ## List of items to return
  return_list <- list("y" = y_og,
                      "ytilde" = ytilde,
                      "X" = X,
                      "A" = A,
                      "B" = B,
                      "B_raw" = B_raw,
                      "K" = K,
                      "p" = nc,
                      "q" = ncol(predictors),
                      "P" = (K+1)*nc,
                      "N" = nr,
                      "penalties" = tL,
                      "knot_scale_transf" = transf,
                      "knot_scale_inv_transf" = inv_transf,
                      "knots" = knots,
                      "partition_codes" = partition_codes,
                      "knot_expand_function" = knot_expand_list,
                      "predict" = predict_function,
                      "assign_partition" = assign_partition,
                      "family" = family,
                      "estimate_dispersion" = estimate_dispersion,
                      "backtransform_coefficients" = backtransform_coefficients,
                      "forwtransform_coefficients" = forwtransform_coefficients,
                      "mean_y" = mean_y,
                      "sd_y" = sd_y,
                      "og_order" = og_order,
                      "order_list" = order_list,
                      "constraint_value_vectors" = constraint_value_vectors,
                      "constraint_vectors" = constraint_vectors,
                      "make_partition_list" = partitions,
                      "C_scales" = C_scales,
                      "take_derivative" = take_derivative,
                      "take_interaction_2ndderivative" =
                        take_interaction_2ndderivative,
                      "get_all_derivatives_insample" = function(expansions){
                        all_derivatives(expansions)},
                      "numerics" = numerics,
                      "power1_cols" = power1_cols,
                      "power2_cols" = power2_cols,
                      "power3_cols" = power3_cols,
                      "power4_cols" = power4_cols,
                      "quad_cols" = quad_cols,
                      "interaction_single_cols" = interaction_single_cols,
                      "interaction_quad_cols" = interaction_quad_cols,
                      "triplet_cols" = triplet_cols,
                      "nonspline_cols" = nonspline_cols,
                      "return_varcovmat" = return_varcovmat,
                      "raw_expansion_names" = colnm_C,
                      "std_X" = std_X,
                      "unstd_X" = unstd_X,
                      "any_parallel" = parallel,
                      "weights" = observation_weights)

  if(verbose){
    cat("Optional Components\n")
  }

  ## We need U and sigma^2 to compute sigma^2*UG
  if(return_varcovmat){
    return_U <- TRUE
    estimate_dispersion <- TRUE
  }

  ## Option is offered to not return these matrices to save memory/time

  ## Return scaled variance-covariance matrix components of coefficients
  # Note: these are on centered-and-scaled y, standardized-X scale
  # Backtransforms are needed to get the varcov on raw scale
  # An option provided below
  if(return_G){
    return_list$G <- G_list$G
  }
  if(return_Ghalf){
    return_list$Ghalf <- G_list$Ghalf
  }
  if(return_U){
    if(verbose){
      cat("U\n")
    }
    if(K == 0 & length(constraint_value_vectors) == 0){
      return_list$U <- diag(nc*(K+1))
      ## ensuring compatibility with no A
      if(any(is.null(A))){
        ## for compatibility, albeit inefficient
        A <- cbind(rep(0, (K+1)*nc))
        A <- cbind(A, A)
        nca <- 2
      }
    } else {
      return_list$U <- get_U(
        G_list$G,
        A,
        K,
        nc,
        nca
      )
    }
  }
  ## Estimate sigma^2
  if(estimate_dispersion){
    if(verbose){
      cat("Variance Est \n")
    }
    ## Compute trace of XUGX^{T} = trace of UGX^{T}X
    if(K == 0){
      trace_XUGX <- sum(unlist(lapply(
        matmult_block_diagonal(
          G_list$G,
          X_gram,
          K,
          parallel =
            FALSE,
          cl = NULL,
          chunk_size,
          num_chunks,
          rem_chunks),
        diag)))
    } else {

      ## No support of parallelism in sub-functions here due to errors,
      # relative efficiency of existing code,
      # and the fact that this is an optional step post-fitting
      trace_XUGX <- compute_trace_UGXX_wrapper(
        G_list$G,
        A,
        # GX^{T}X
        matmult_block_diagonal(G_list$G,
                               X_gram,
                               K,
                               parallel = parallel & parallel_matmult,
                               cl = cl,
                               chunk_size,
                               num_chunks,
                               rem_chunks),
        G_list$Ghalf,
        # (A^{T}GA)^{-1}
        invert(AGAmult_wrapper(G_list$G,
                               A,
                               K,
                               nc,
                               nca,
                               parallel = parallel & parallel_aga,
                               cl = cl,
                               chunk_size,
                               num_chunks,
                               rem_chunks)),
        nc,
        K,
        parallel = FALSE,
        cl = cl,
        chunk_size,
        num_chunks,
        rem_chunks)
    }
    if(trace_XUGX < 0 & include_warnings){
      warning('Trace of XUGX^{T} is < 0, which most often indicates a failure of convergence when fitting (i.e. the constrained maximum likelihood estimate was not found). Try re-fitting, different knot locations, greater penalties, or a less complicated model. Alteratively, try to recompute the trace manually using XUGUX^{T} instead.')
    }

    ## Estimating exponential dispersion or variance
    if(paste0(family)[1] == 'gaussian'){
      ## Dispersion estimate (variance for Gaussian family)
      return_list$sigmasq_tilde <-
        sum(observation_weights_og * (y_og - ytilde)^2 / (nr - trace_XUGX))
    } else {
      ## Dispersion estimate (using custom function)
      return_list$sigmasq_tilde <- dispersion_function(
                          ytilde,
                          y_og,
                          1:length(y_og), # this is original order!
                          family,
                          observation_weights_og,
                          ...)
    }

    ## Effective degrees of freedom is the trace, when we have penalization
    return_list$trace_XUGX <- trace_XUGX

  } else {
    ## Otherwise, return 1 for dispersion
    return_list$sigmasq_tilde <- 1
  }

  if(return_varcovmat){
    if(verbose){
      cat("VarCov Mat \n")
    }

    ## Use UGU^{T} parameterization rather than just UG for numeric stability
    return_list$varcovmat <-
      matmult_U(return_list$U, G_list$G, nc, K) %**%
      t(return_list$U)

    ## Un-standardize
    d <- rep(c(1, 1/C_scales), each = K + 1)
    return_list$varcovmat <-
      return_list$sigmasq_tilde *
      t(t(return_list$varcovmat * d) * d)

    ## Replace < 0 diagonals with 0
    if(any(diag(return_list$varcovmat) < 0) & include_warnings){
      warning("Variance-covariance matrix has diagonal elements < 0, model most likely did not converge when fitting. Try re-fitting, a simpler model, changing knot locations, or increasing the penalties.")
      for(ij in 1:nrow(return_list$varcovmat)){
        return_list$varcovmat[ij,ij] <- max(0,
                                            return_list$varcovmat[ij, ij])
      }
    }
  }

  ## Afterwards, update X to be unstandardized
  return_list$X <- lapply(return_list$X,
                          unstd_X)

  return(return_list)
}

#' Compute leave-one-out cross-validated predictions for Gaussian response/identity link
#'
#' @description
#' Computes the leave-one-out cross-validated predictions from a model fit,
#' assuming Gaussian-distributed response with identity link.
#'
#' @param model_fit A fitted Lagrangian smoothing spline model
#'
#' @return A vector of leave-one-out cross-validated predictions
#'
#' @examples
#' \dontrun{
#'
#' ## Basic usage with Gaussian response
#' set.seed(1234)
#' x <- matrix(rnorm(50), ncol=1)
#' y <- sin(x) + rnorm(50, 0, .25)
#' fit <- lgspline(x, y)
#' loo <- leave_one_out(fit)
#' plot(x , loo, col = 'blue')
#' points(x, y, col = 'red')
#'
#' }
#'
#' @export
leave_one_out <- function(model_fit){

  ## Collapse X into a single block-diagonal matrix, unscaled
  X_block <- collapse_block_diagonal(
    lapply(model_fit$X, model_fit$std_X)
  )

  ## UG efficiently multiplied together
  UG <-  matmult_U(model_fit$U, model_fit$G, model_fit$p,  model_fit$K)

  ## The expensive operation
  diag_XUGX <- rowSums((X_block %**% UG) * X_block)

  ## Order it correctly
  diag_XUGX <- diag_XUGX[unlist(model_fit$og_order)]

  ## LOO predictions
  leave_one_out <-
    model_fit$y -
    1/(1 - diag_XUGX) *
    (model_fit$y - model_fit$ytilde)

  return(leave_one_out)
}
