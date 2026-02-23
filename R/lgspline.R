#' Lagrangian Multiplier Smoothing Splines
#'
#' @docType package
#' @keywords internal
#' @name lgspline-package
#' @rdname lgspline-package
#' @aliases lgspline-package
#'
#' @import Rcpp RcppArmadillo methods stats
#' @importFrom graphics plot points legend
#' @importFrom FNN get.knnx
#' @importFrom quadprog solve.QP
#' @importFrom RColorBrewer brewer.pal
#' @importFrom plotly plot_ly
#'
#' @keywords smoothing regression parametric constrained lagrangian multiplier
"_PACKAGE"

#' Fit Lagrangian Multiplier Smoothing Splines
#'
#' @description
#'
#' A comprehensive software package for fitting a variant of smoothing splines
#' as a constrained optimization problem, avoiding the need to algebraically
#' disentangle a spline basis after fitting, and allowing for interpretable
#' interactions and non-spline effects to be included.
#'
#' \code{lgspline} fits piecewise polynomial regression splines constrained to be smooth where
#' they meet, penalized by the squared, integrated, second-derivative of the
#' estimated function with respect to predictors.
#'
#' The method of Lagrangian multipliers is used to enforce the following:
#' \itemize{
#'   \item Equivalent fitted values at knots
#'   \item Equivalent first derivatives at knots, with respect to predictors
#'   \item Equivalent second derivatives at knots, with respect to predictors
#' }
#'
#' The coefficients are penalized by a closed-form of the traditional
#' cubic smoothing spline penalty, as well as tunable modifications that allow
#' for unique penalization of multiple predictors and partitions.
#'
#' This package supports model fitting for multiple spline and non-spline effects, GLM families,
#' Weibull accelerated failure time (AFT) models, correlation structures, quadratic
#' programming constraints, and extensive customization for user-defined models.
#'
#' In addition, parallel processing capabilities and comprehensive
#' tools for visualization, frequentist, and Bayesian inference are provided.
#'
#' @usage lgspline(predictors = NULL, y = NULL, formula = NULL, response = NULL,
#'                 standardize_response = TRUE, standardize_predictors_for_knots = TRUE,
#'                 standardize_expansions_for_fitting = TRUE, family = gaussian(),
#'                 glm_weight_function = function(mu, y, order_indices, family, dispersion,
#'                                                observation_weights, ...) {
#'                   if(any(!is.null(observation_weights))){
#'                     family$variance(mu) * observation_weights
#'                   } else {
#'                     family$variance(mu)
#'                   }
#'                 }, schur_correction_function = function(X, y, B, dispersion, order_list, K,
#'                                                        family, observation_weights, ...) {
#'                   lapply(1:(K+1), function(k) 0)
#'                 }, need_dispersion_for_estimation = FALSE,
#'                 dispersion_function = function(mu, y, order_indices, family,
#'                                                observation_weights, VhalfInv,
#'                                                ...) { 1 },
#'                 K = NULL, custom_knots = NULL, cluster_on_indicators = FALSE,
#'                 make_partition_list = NULL, previously_tuned_penalties = NULL,
#'                 smoothing_spline_penalty = NULL, opt = TRUE, use_custom_bfgs = TRUE,
#'                 delta = NULL, tol = 10*sqrt(.Machine$double.eps),
#'                 initial_wiggle = exp(c(-25, -14, -3, -0.2, 1)),
#'                 initial_flat = c(0.1, 10), wiggle_penalty = 2e-07,
#'                 flat_ridge_penalty = 0.5, unique_penalty_per_partition = TRUE,
#'                 unique_penalty_per_predictor = TRUE, meta_penalty = 1e-08,
#'                 predictor_penalties = NULL, partition_penalties = NULL,
#'                 include_quadratic_terms = TRUE, include_cubic_terms = TRUE,
#'                 include_quartic_terms = NULL, include_2way_interactions = TRUE,
#'                 include_3way_interactions = TRUE, include_quadratic_interactions = FALSE,
#'                 offset = c(), just_linear_with_interactions = NULL,
#'                 just_linear_without_interactions = NULL, exclude_interactions_for = NULL,
#'                 exclude_these_expansions = NULL, custom_basis_fxn = NULL,
#'                 include_constrain_fitted = TRUE, include_constrain_first_deriv = TRUE,
#'                 include_constrain_second_deriv = TRUE,
#'                 include_constrain_interactions = TRUE, cl = NULL, chunk_size = NULL,
#'                 parallel_eigen = TRUE, parallel_trace = FALSE, parallel_aga = FALSE,
#'                 parallel_matmult = FALSE, parallel_unconstrained = TRUE,
#'                 parallel_find_neighbors = FALSE, parallel_penalty = FALSE,
#'                 parallel_make_constraint = FALSE,
#'                 unconstrained_fit_fxn = unconstrained_fit_default,
#'                 keep_weighted_Lambda = FALSE, iterate_tune = TRUE,
#'                 iterate_final_fit = TRUE, blockfit = FALSE,
#'                 qp_score_function = function(X, y, mu, order_list, dispersion,
#'                                              VhalfInv, observation_weights, ...) {
#'                   if(!is.null(observation_weights)) {
#'                     crossprod(X, cbind((y - mu)*observation_weights))
#'                   } else {
#'                     crossprod(X, cbind(y - mu))
#'                   }
#'                 }, qp_observations = NULL, qp_Amat = NULL, qp_bvec = NULL, qp_meq = 0,
#'                 qp_positive_derivative = FALSE, qp_negative_derivative = FALSE,
#'                 qp_positive_2ndderivative = FALSE, qp_negative_2ndderivative = FALSE,
#'                 qp_monotonic_increase = FALSE, qp_monotonic_decrease = FALSE,
#'                 qp_range_upper = NULL, qp_range_lower = NULL, qp_Amat_fxn = NULL,
#'                 qp_bvec_fxn = NULL, qp_meq_fxn = NULL, constraint_values = cbind(),
#'                 constraint_vectors = cbind(), return_G = TRUE, return_Ghalf = TRUE,
#'                 return_U = TRUE, estimate_dispersion = TRUE, unbias_dispersion = NULL,
#'                 return_varcovmat = TRUE, return_lagrange_multipliers = FALSE,
#'                 custom_penalty_mat = NULL,
#'                 cluster_args = c(custom_centers = NA, nstart = 10),
#'                 dummy_dividor = 1.2345672152894e-22,
#'                 dummy_adder = 2.234567210529e-18,
#'                 verbose = FALSE, verbose_tune = FALSE,
#'                 dummy_fit = FALSE, auto_encode_factors = TRUE,
#'                 observation_weights = NULL, do_not_cluster_on_these = c(),
#'                 neighbor_tolerance = 1 + 1e-08, null_constraint = NULL,
#'                 critical_value = qnorm(1 - 0.05/2), data = NULL, weights = NULL,
#'                 no_intercept = FALSE, correlation_id = NULL, spacetime = NULL,
#'                 correlation_structure = NULL, VhalfInv = NULL, Vhalf = NULL,
#'                 VhalfInv_fxn = NULL, Vhalf_fxn = NULL, VhalfInv_par_init = c(),
#'                 REML_grad = NULL, custom_VhalfInv_loss = NULL, VhalfInv_logdet = NULL,
#'                 include_warnings = TRUE, ...)
#'
#' @details
#' A flexible and interpretable implementation of smoothing splines including:
#' \itemize{
#'   \item Multiple predictors and interaction terms
#'   \item Various GLM families and link functions
#'   \item Correlation structures for longitudinal/clustered data
#'   \item Shape constraints via quadratic programming
#'   \item Parallel computation for large datasets
#'   \item Comprehensive inference tools
#' }
#'
#' @param predictors Default: NULL. Numeric matrix or data frame of predictor variables. Supports direct matrix input or formula interface when used with \code{data} argument. Must contain numeric predictors, with categorical variables pre-converted to numeric indicators (or set \code{auto_encode_factors = TRUE} to encode them automatically when using the formula interface).
#' @param y Default: NULL. Numeric response variable vector representing the target/outcome/dependent variable etc. to be modeled.
#' @param formula Default: NULL. Optional statistical formula for model specification, serving as an alternative to direct matrix input. Supports standard R formula syntax with special \code{spl()} function for defining spline terms.
#' @param response Default: NULL. Alternative name for response variable, providing compatibility with different naming conventions. Takes precedence only if \code{y} is not supplied.
#' @param standardize_response Default: TRUE. Logical indicator controlling whether the response variable should be centered by mean and scaled by standard deviation before model fitting. When TRUE, tends to improve numerical stability. Only offered for identity link functions (when \code{family$link == 'identity'}).
#' @param standardize_predictors_for_knots Default: TRUE. Logical flag determining whether predictor variables should be standardized before knot placement. Ensures consistent knot selection across different predictor scales, and that no one predictor dominates in terms of influence on knot placement. For all expansions (x), standardization corresponds to dividing by the difference in 69th and 31st percentiles = x / (quantile(x, 0.69) - quantile(x, 0.31)).
#' @param standardize_expansions_for_fitting Default: TRUE. Logical switch to standardize polynomial basis expansions during model fitting. Provides computational stability during penalty tuning without affecting statistical inference, as design matrices, variance-covariance matrices, and coefficient estimates are systematically backtransformed after fitting to account for the standardization. If TRUE, \eqn{\mathbf{U}} and \eqn{\mathbf{G}} will remain on the transformed scale, and \code{B_raw} as returned will correspond to the coefficients fitted on the expansion-standardized scale.
#' @param family Default: \code{gaussian()}. Generalized linear model (GLM) distribution family specifying the error distribution and link function for model fitting. Defaults to Gaussian distribution with identity link. Supports custom family specifications, including user-defined link functions and optional custom tuning loss criteria. Minimally requires 1) family name (family) 2) link name (link) 3) linkfun (link function) 4) linkinv (link function inverse) 5) variance (mean variance relationship function, \eqn{\mathrm{Var}(Y|\mu)}).
#' @param glm_weight_function Default: function that returns \code{family$variance(mu) * observation_weights} if weights exist, \code{family$variance(mu)} otherwise. Codes the mean-variance relationship of a GLM or GLM-like model, the diagonal \eqn{\mathbf{W}} matrix of \eqn{\mathbf{X}^{\top}\mathbf{W}\mathbf{X}} that appears in the information. This can be replaced with a user-specified function. It is used for updating \eqn{\mathbf{G} = (\mathbf{X}^{\top}\mathbf{W}\mathbf{X} + \boldsymbol{\Lambda})^{-1}} after obtaining constrained estimates of coefficients. This is not used for fitting unconstrained models, but for iterating between updates of \eqn{\mathbf{U}}, \eqn{\mathbf{G}}, and beta coefficients afterwards.
#' @param schur_correction_function Default: function that returns list of zeros. Advanced function for computing Schur complements \eqn{\mathbf{S}} to add to \eqn{\mathbf{G}} to properly account for uncertainty in dispersion or other nuisance parameter estimation. The effective information becomes \eqn{\mathbf{G}^{*} = (\mathbf{G}^{-1} + \mathbf{S})^{-1}}. The incorrect spelling using "shur" is from a previous version (< 1.0), and will return a warning, and ignored.
#' @param need_dispersion_for_estimation Default: FALSE. Logical indicator specifying whether a dispersion parameter is required for coefficient estimation. This is not needed for canonical regular exponential family models, but is often needed otherwise (such as fitting Weibull AFT models).
#' @param dispersion_function Default: function that returns 1. Custom function for estimating the dispersion parameter. Unless \code{need_dispersion_for_estimation} is TRUE, this will not affect coefficient estimates.
#' @param K Default: NULL. Integer specifying the number of knot locations for spline partitions. This can intuitively be considered the total number of partitions minus 1.
#' @param custom_knots Default: NULL. Optional matrix providing user-specified knot locations in 1-D.
#' @param cluster_on_indicators Default: FALSE. Logical flag determining whether indicator variables should be used for clustering knot locations.
#' @param make_partition_list Default: NULL. Optional list allowing direct specification of custom partition assignments. The intent is that the \code{make_partition_list} returned by one model can be supplied here to keep the same knot locations for another.
#' @param previously_tuned_penalties Default: NULL. Optional list of pre-computed penalty components from previous model fits.
#' @param smoothing_spline_penalty Default: NULL. Optional custom smoothing spline penalty matrix for fine-tuned complexity control.
#' @param opt Default: TRUE. Logical switch controlling whether model penalties should be automatically optimized via generalized cross-validation. Turn this off if \code{previously_tuned_penalties} are supplied AND desired, otherwise, the \code{previously_tuned_penalties} will be ignored.
#' @param use_custom_bfgs Default: TRUE. Logical indicator selecting between a native implementation of damped-BFGS optimization method with closed-form gradients or base R's BFGS implementation with finite-difference approximation of gradients. The custom implementation is more efficient when gradients are available, since function output and gradients are returned simultaneously, instead of separately, effectively cutting the number of large matrix operations in half.
#' @param delta Default: NULL. Numeric pseudocount used for stabilizing optimization in non-identity link function scenarios.
#' @param tol Default: \code{10*sqrt(.Machine$double.eps)}. Numeric convergence tolerance controlling the precision of optimization algorithms.
#' @param initial_wiggle Default: \code{exp(c(-25, -14, -3, -0.2, 1))}. Numeric vector of initial grid points for wiggle penalty optimization, specified on the raw (non-negative) scale. Conversion to log scale for the internal optimizer is handled automatically.
#' @param initial_flat Default: \code{c(0.1, 10)}. Numeric vector of initial grid points for ridge penalty optimization, specified on the raw scale; multiply by wiggle penalty to get the total ridge penalty (this parameter is a ratio of ridge to wiggle).
#' @param wiggle_penalty Default: 2e-7. Numeric penalty controlling the integrated squared second derivative, governing function smoothness. Applied to both smoothing spline penalty (alone) and is multiplied by \code{flat_ridge_penalty} for penalizing linear terms.
#' @param flat_ridge_penalty Default: 0.5. Numeric flat ridge penalty for additional regularization on only intercepts and linear terms (will not affect interactions or quadratic/cubic/quartic terms by default). If \code{custom_penalty_mat} is supplied, the penalty will be for the custom penalty matrix instead. This penalty is multiplied with \code{wiggle_penalty} to obtain the total ridge penalty — hence, by default, the ridge penalization on linear terms is half of the magnitude of non-linear terms.
#' @param unique_penalty_per_partition Default: TRUE. Logical flag allowing the magnitude of the smoothing spline penalty to differ across partitions.
#' @param unique_penalty_per_predictor Default: TRUE. Logical flag allowing the magnitude of the smoothing spline penalty to differ between predictors.
#' @param meta_penalty Default: 1e-8. Numeric "meta-penalty" applied to predictor and partition penalties during tuning. The minimization of GCV is modified to be a penalized minimization problem, with penalty \eqn{0.5 \times \mathrm{meta\_penalty} \times (\sum \log(\mathrm{penalty}))^{2}}, such that penalties are pulled towards 1 on the absolute scale and thus, no multiplicative effect.
#' @param predictor_penalties Default: NULL. Optional vector of custom penalties specified per predictor, on the raw (positive) scale.
#' @param partition_penalties Default: NULL. Optional vector of custom penalties specified per partition, on the raw (positive) scale.
#' @param include_quadratic_terms Default: TRUE. Logical switch to include squared predictor terms in basis expansions.
#' @param include_cubic_terms Default: TRUE. Logical switch to include cubic predictor terms in basis expansions.
#' @param include_quartic_terms Default: NULL. Logical switch to include quartic predictor terms in basis expansions. This is highly recommended for fitting models with multiple predictors to avoid over-specified constraints. When NULL (by default), will internally set to FALSE if only one predictor is present, and TRUE otherwise.
#' @param include_2way_interactions Default: TRUE. Logical switch to include linear two-way interactions between predictors.
#' @param include_3way_interactions Default: TRUE. Logical switch to include three-way interactions between predictors.
#' @param include_quadratic_interactions Default: FALSE. Logical switch to include linear-quadratic interaction terms.
#' @param offset Default: Empty vector. When non-missing, this is a vector of column indices/names to include as offsets. \code{lgspline} will automatically introduce constraints such that the coefficient for offset terms are 1.
#' @param just_linear_with_interactions Default: NULL. Integer vector specifying columns to retain linear terms with interactions.
#' @param just_linear_without_interactions Default: NULL. Integer vector specifying columns to retain only linear terms without interactions. These columns are eligible for blockfitting when \code{blockfit = TRUE}; their coefficients are pooled across partitions via backfitting rather than the Lagrangian projection.
#' @param exclude_interactions_for Default: NULL. Integer vector indicating columns to exclude from all interaction terms.
#' @param exclude_these_expansions Default: NULL. Character vector specifying basis expansions to be excluded from the model. These must be named columns of the data, or in the form \code{"_1_"}, \code{"_2_"}, \code{"_1_x_2_"}, \code{"_2_^2"} etc. where \code{"1"} and \code{"2"} indicate column indices of predictor matrix input.
#' @param custom_basis_fxn Default: NULL. Optional user-defined function for generating custom basis expansions. See \code{\link{get_polynomial_expansions}}.
#' @param include_constrain_fitted Default: TRUE. Logical switch to constrain fitted values at knot points.
#' @param include_constrain_first_deriv Default: TRUE. Logical switch to constrain first derivatives at knot points.
#' @param include_constrain_second_deriv Default: TRUE. Logical switch to constrain second derivatives at knot points.
#' @param include_constrain_interactions Default: TRUE. Logical switch to constrain interaction terms at knot points.
#' @param cl Default: NULL. Parallel processing cluster object for distributed computation (use \code{parallel::makeCluster()}).
#' @param chunk_size Default: NULL. Integer specifying custom fixed chunk size for parallel processing.
#' @param parallel_eigen Default: TRUE. Logical flag to enable parallel processing for eigenvalue decomposition computations.
#' @param parallel_trace Default: FALSE. Logical flag to enable parallel processing for trace computation.
#' @param parallel_aga Default: FALSE. Logical flag to enable parallel processing for specific matrix operations involving \eqn{\mathbf{G}} and \eqn{\mathbf{A}}.
#' @param parallel_matmult Default: FALSE. Logical flag to enable parallel processing for block-diagonal matrix multiplication.
#' @param parallel_unconstrained Default: TRUE. Logical flag to enable parallel processing for unconstrained maximum likelihood estimation for models other than identity-link Gaussian. Ignored for identity link Gaussian, as the step of obtaining unconstrained estimates is not explicitly performed. Users can set \code{parallel_matmult = TRUE} to ensure parallel processing for fitting identity link Gaussian models.
#' @param parallel_find_neighbors Default: FALSE. Logical flag to enable parallel processing for neighbor identification (which partitions are neighbors).
#' @param parallel_penalty Default: FALSE. Logical flag to enable parallel processing for penalty matrix construction.
#' @param parallel_make_constraint Default: FALSE. Logical flag to enable parallel processing for constraint matrix generation.
#' @param unconstrained_fit_fxn Default: \code{\link{unconstrained_fit_default}}. Custom function for fitting unconstrained models per partition.
#' @param keep_weighted_Lambda Default: FALSE. Logical flag to retain generalized linear model weights in penalty constraints using Tikhonov parameterization. It is advised to turn this to TRUE when fitting non-canonical GLMs. The default \code{\link{unconstrained_fit_default}} by default assumes canonical GLMs for setting up estimating equations; this is not valid with non-canonical GLMs. With \code{keep_weighted_Lambda = TRUE}, the Tikhonov parameterization binds \eqn{\boldsymbol{\Lambda}^{1/2}}, the square-root penalty matrix, to the design matrix \eqn{\mathbf{X}_{k}} for each partition k, and \code{family$linkinv(0)} to the response vector \eqn{\mathbf{y}_{k}} for each partition before finding unconstrained estimates using base R's \code{glm.fit} function. The potential issue is that the weights of the information matrix will appear in the penalty, such that the effective penalty is \eqn{\boldsymbol{\Lambda}_{\mathrm{eff}} = \mathbf{L}^{1/2}\mathbf{W}\mathbf{L}^{1/2}} rather than just \eqn{\mathbf{L}^{1/2}\mathbf{L}^{1/2}}. If FALSE, this approach will only be used to supply initial values to a native implementation of damped Newton-Raphson for fitting GLM models (see \code{\link{damped_newton_r}} and \code{\link{unconstrained_fit_default}}). For Gamma with log-link, this is fortunately a non-issue since the mean-variance relationship is essentially stabilized, so \code{keep_weighted_Lambda = TRUE} is strongly advised.
#' @param iterate_tune Default: TRUE. Logical switch to use iterative optimization during penalty tuning. If FALSE, \eqn{\mathbf{G}} and \eqn{\mathbf{U}} are constructed from unconstrained \eqn{\boldsymbol{\beta}} estimates when tuning.
#' @param iterate_final_fit Default: TRUE. Logical switch to use iterative optimization for final model fitting. If FALSE, \eqn{\mathbf{G}} and \eqn{\mathbf{U}} are constructed from unconstrained \eqn{\boldsymbol{\beta}} estimates when fitting the final model after tuning.
#' @param blockfit Default: FALSE. Logical switch to use backfitting for models with mixed spline and non-interactive linear ("flat") terms. When \code{TRUE} and the model contains columns identified by \code{just_linear_without_interactions}, the fitting is routed through \code{blockfit_solve}, which splits the design matrix into a spline block (spline and interaction effects) and a flat block (non-interactive linear terms). Coefficients for the flat block are pooled (shared identically) across all \eqn{K+1} partitions via backfitting, while spline-block coefficients are estimated per partition via the Lagrangian projection. This is more efficient and numerically stable when there are many non-spline effects without interactions relative to the number of spline effects. Blockfit option requires all four conditions: (1) \code{blockfit = TRUE}, (2) flat columns are non-empty, (3) \code{K > 0}, and (4) no correlation structure is active (\code{Vhalf} and \code{VhalfInv} both \code{NULL}). QP inequality constraints are supported within \code{blockfit_solve} via a post-backfitting SQP refinement loop. If any condition fails or \code{blockfit_solve} throws an error, the standard \code{get_B} path executes as a fallback.
#' @param qp_score_function Default: \eqn{\mathbf{X}^{\top}(\mathbf{y} - \mathrm{E}[\mathbf{y}])}, where \eqn{\mathrm{E}[\mathbf{y}] = \boldsymbol{\mu}}. A function returning the score of the log-likelihood for optimization (excluding penalization/priors involving \eqn{\boldsymbol{\Lambda}}), which is needed for the formulation of quadratic programming problems, when \code{blockfit = TRUE}, and correlation-structure fitting for GLMs, all relying on \code{\link[quadprog]{solve.QP}}. Accepts arguments \code{"X, y, mu, order_list, dispersion, VhalfInv, observation_weights, ..."} in order.
#' @param qp_observations Default: NULL. Numeric vector of observations to apply constraints to for monotonic and range quadratic programming constraints. Useful for saving computational resources.
#' @param qp_Amat Default: NULL. Constraint matrix for quadratic programming formulation. The \code{Amat} argument of \code{\link[quadprog]{solve.QP}}.
#' @param qp_bvec Default: NULL. Constraint vector for quadratic programming formulation. The \code{bvec} argument of \code{\link[quadprog]{solve.QP}}.
#' @param qp_meq Default: 0. Number of equality constraints in quadratic programming setup. The \code{meq} argument of \code{\link[quadprog]{solve.QP}}.
#' @param qp_positive_derivative,qp_positive_2ndderivative,qp_monotonic_increase Default: FALSE. Logical flags to constrain the function to have positive first derivatives/second derivatives/be monotonically increasing using quadratic programming with respect to the order (ascending rows) of the input data set.
#' @param qp_negative_derivative,qp_negative_2ndderivative,qp_monotonic_decrease Default: FALSE. Logical flags to constrain the function to have negative first derivatives/second derivatives/be monotonically decreasing using quadratic programming with respect to the order (ascending rows) of the input data set.
#' @param qp_range_upper Default: NULL. Numeric upper bound for constrained fitted values using quadratic programming.
#' @param qp_range_lower Default: NULL. Numeric lower bound for constrained fitted values using quadratic programming.
#' @param qp_Amat_fxn Default: NULL. Custom function for generating Amat matrix in quadratic programming.
#' @param qp_bvec_fxn Default: NULL. Custom function for generating bvec vector in quadratic programming.
#' @param qp_meq_fxn Default: NULL. Custom function for determining meq equality constraints in quadratic programming.
#' @param constraint_values Default: \code{cbind()}. Matrix of constraint values for sum constraints. The constraint enforces \eqn{\mathbf{C}^{\top}(\boldsymbol{\beta} - \mathbf{c}) = \mathbf{0}} in addition to smoothing constraints, where \eqn{\mathbf{C}} = \code{constraint_vectors} and \eqn{\mathbf{c}} = \code{constraint_values}.
#' @param constraint_vectors Default: \code{cbind()}. Matrix of vectors for sum constraints. The constraint enforces \eqn{\mathbf{C}^{\top}(\boldsymbol{\beta} - \mathbf{c}) = \mathbf{0}} in addition to smoothing constraints, where \eqn{\mathbf{C}} = \code{constraint_vectors} and \eqn{\mathbf{c}} = \code{constraint_values}.
#' @param return_G Default: TRUE. Logical switch to return the unscaled variance-covariance matrix without smoothing constraints (\eqn{\mathbf{G}}), which is block diagonal of p x p matrices for each of the K+1 partitions. When a correlation structure is present, the returned \eqn{\mathbf{G}_{\mathrm{correct}}} uses the full whitened Gram \eqn{\mathbf{X}^{\top}\mathbf{V}^{-1}\mathbf{X}} rather than the block-diagonal \eqn{\mathbf{X}^{\top}\mathbf{X}}.
#' @param return_Ghalf Default: TRUE. Logical switch to return the matrix square root of the unscaled variance-covariance matrix without smoothing constraints (\eqn{\mathbf{G}^{1/2}}). When a correlation structure is active, \eqn{\mathbf{G}^{1/2}} is computed from the full penalized GLS information \eqn{(\mathbf{X}^{\top}\mathbf{V}^{-1}\mathbf{X} + \boldsymbol{\Lambda})^{-1}}.
#' @param return_U Default: TRUE. Logical switch to return the constraint projection matrix \eqn{\mathbf{U}}.
#' @param estimate_dispersion Default: TRUE. Logical flag to estimate the dispersion parameter after fitting.
#' @param unbias_dispersion Default: NULL. Logical switch to multiply final dispersion estimates by \eqn{N/(N - \mathrm{trace}(\mathbf{H}))} which, in the case of Gaussian-distributed errors with identity link function, provides unbiased estimates of variance. When NULL (by default), is set to TRUE for Gaussian + identity link and FALSE otherwise.
#' @param return_varcovmat Default: TRUE. Logical switch to return the variance-covariance matrix of the estimated coefficients. This is needed for performing Wald inference. When a correlation structure is active, the variance-covariance matrix is computed as \eqn{\sigma^{2}(\mathbf{U}\mathbf{G}_{\mathrm{correct}}^{1/2})(\mathbf{U}\mathbf{G}_{\mathrm{correct}}^{1/2})^{\top}} where \eqn{\mathbf{G}_{\mathrm{correct}} = (\mathbf{X}^{\top}\mathbf{V}^{-1}\mathbf{X} + \boldsymbol{\Lambda})^{-1}}. This is the correct version, since the block-diagonal version constructed ignores off-block diagonal components, a computational convenience not available for generic correlation structures.
#' @param return_lagrange_multipliers Default: FALSE. Logical switch to return the Lagrangian multiplier vector associated with the smoothness constraints \eqn{ = (\mathbf{A}^{\top}\mathbf{G}\mathbf{A})^{-1}\mathbf{A}^{\top}(\hat{\boldsymbol{\beta}} - \boldsymbol{\beta_0})}. Requires \code{return_G = TRUE} and \code{K > 0} or non-trivial constraints. For inequality constraints, the multipliers as returned by \code{\link{solve.QP}} are used.
#' @param custom_penalty_mat Default: NULL. Optional \eqn{p \times p} custom penalty matrix for individual partitions to replace the default ridge penalty applied to linear-and-intercept terms only. This can be interpreted as proportional to the prior correlation matrix of coefficients for non-spline effects, and will appear in the penalty matrix for all partitions. It is recommended to first run the function using \code{dummy_fit = TRUE} so you have an idea of where the expansions appear in each partition, what "p" is, and you can carefully customize your penalty matrix after.
#' @param cluster_args Default: \code{c(custom_centers = NA, nstart = 10)}. Named vector of arguments controlling clustering procedures. If the first argument is not NA, this will be treated as custom cluster centers and all other arguments ignored. Otherwise, default base R k-means clustering will be used with all other arguments supplied to \code{kmeans} (for example, by default, the "nstart" argument as provided). Custom centers must be a \eqn{K \times q} matrix with one column for each predictor in order of their appearance in input predictor/data, and one row for each center.
#' @param dummy_dividor Default: 0.00000000000000000000012345672152894. Small numeric constant to prevent division by zero in computational routines.
#' @param dummy_adder Default: 0.000000000000000002234567210529. Small numeric constant to prevent division by zero in computational routines.
#' @param verbose Default: FALSE. Logical flag to print general progress messages during model fitting (does not include during tuning).
#' @param verbose_tune Default: FALSE. Logical flag to print detailed progress messages during penalty tuning specifically.
#' @param dummy_fit Default: FALSE. Logical switch to run the full pipeline (including penalty construction) but set coefficients to zero, allowing users to inspect the design matrix structure, penalty matrices, and partitioning information before fitting. Replaces the deprecated \code{expansions_only} argument. An early-return path in \code{lgspline.fit} skips the expensive fitting steps (\code{compute_G_eigen}, \code{get_B}, trace computation, variance-covariance matrix) while retaining penalty and partitioning information for downstream use such as setting up custom constraints.
#' @param auto_encode_factors Default: TRUE. Logical switch to automatically one-hot encode factor or character variables when using the formula interface. Encoded variables are appended to \code{just_linear_without_interactions} unless the original column appears in an interaction term. Interactions will be excluded for auto-encoded variables by default; users can manually one-hot encode their variables and include in formula to achieve interactions instead, or use the design matrix format and indicate the column indices of interactions. Dummy-intercept is avoided, since one-hot encoding is suitable for interpretation of penalized fixed effects as "random" effects.Turning off is useful for saving memory when working with datasets with many character variables that don't end up being used in the model.
#' @param observation_weights Default: NULL. Numeric vector of observation-specific weights for generalized least squares estimation. The argument \code{weights} is an alias for this. This is the equivalent to the usual \code{weights} in \code{glm} and \code{lm}.
#' @param do_not_cluster_on_these Default: \code{c()}. Vector specifying predictor columns to exclude from clustering procedures, in addition to the non-spline effects by default. Accepts either numeric column indices or character column names; character entries are resolved to numeric indices via exact match (with \code{grep()} fallback) against predictor column names in both the formula and direct-matrix paths. If column names are unavailable in the direct-matrix path, a warning is issued and unresolvable character entries are dropped.
#' @param neighbor_tolerance Default: \code{1 + 1e-8}. Numeric tolerance for determining neighboring partitions using k-means clustering. Greater values means more partitions are likely to be considered neighbors. Intended for internal use only (modify at your own risk!).
#' @param null_constraint Default: NULL. Alternative parameterization of constraint values.
#' @param critical_value Default: \code{qnorm(1-0.05/2)}. Numeric critical value used for constructing Wald confidence intervals of the form \eqn{\mathrm{estimate} \pm \mathrm{critical\_value} \times (\mathrm{standard\ error})}.
#' @param data Default: NULL. Optional data frame providing context for formula-based model specification.
#' @param weights Default: NULL. Alternative name for \code{observation_weights}, maintained for interface compatibility and R-user friendliness. This is the equivalent to the usual \code{weights} in \code{glm} and \code{lm}.
#' @param no_intercept Default: FALSE. Logical flag to remove intercept, constraining it to 0. The function automatically constructs \code{constraint_vectors} and \code{constraint_values} to achieve this. Calling formulas with a \code{"0+"} in it like \code{y ~ 0 + .} will set this option to TRUE.
#' @param correlation_id,spacetime Default: NULL. N-length vector and N-row matrix of cluster (or subject, group etc.) ids and longitudinal/spatial variables respectively, whereby observations within each grouping of \code{correlation_id} are correlated with respect to the variables submitted to \code{spacetime}.
#' @param correlation_structure Default: NULL. Native implementations of popular variance-covariance structures. Offers options for \code{"exchangeable"}, \code{"spatial-exponential"}, \code{"squared-exponential"}, \code{"ar(1)"}, \code{"spherical"}, \code{"gaussian-cosine"}, \code{"gamma-cosine"}, and \code{"matern"}, along with their aliases. The eponymous correlation structure is fit along with coefficients and dispersion, with correlation estimated using a REML objective. All built-in structures use exponential parameterization. See section "Correlation Structure Estimation" for more details.
#' @param VhalfInv Default: NULL. Matrix representing a fixed, custom square-root-inverse covariance structure for the response variable for longitudinal and spatial modeling. Must be an \eqn{N \times N} matrix where N is number of observations. This matrix \eqn{\mathbf{V}^{-1/2}} serves as a fixed transformation matrix for the response, equivalent to GLS with known covariance \eqn{\mathbf{V}}. When provided, all post-fit inference quantities (\eqn{\mathbf{G}}, \eqn{\mathbf{U}}, trace, variance-covariance matrix, posterior draws) are recomputed from the whitened Gram matrices \eqn{\mathbf{X}^{\top}\mathbf{V}^{-1}\mathbf{X}} for correct calibration under correlated errors. This is known as "whitening" in some literature.
#' @param Vhalf Default: NULL. Matrix representing a fixed, custom square-root covariance structure for the response variable for longitudinal and spatial modeling. Must be an \eqn{N \times N} matrix where N is number of observations. This matrix \eqn{\mathbf{V}^{1/2}} is used when backtransforming coefficients for fitting GLMs with arbitrary correlation structure. When \code{VhalfInv} is supplied but \code{Vhalf} is not, \code{Vhalf} is computed as the inverse of \code{VhalfInv}.
#' @param VhalfInv_fxn Default: NULL. Function for parametric modeling of the covariance structure \eqn{\mathbf{V}^{-1/2}}. Must take a single numeric vector argument \code{"par"} and return an \eqn{N \times N} matrix. When provided with \code{VhalfInv_par_init}, this function is optimized via BFGS to find optimal covariance parameters that minimize the negative REML log-likelihood (or custom loss if \code{custom_VhalfInv_loss} is specified). The function must return a valid square root of the inverse covariance matrix — i.e., if \eqn{\mathbf{V}} is the true covariance, \code{VhalfInv_fxn} should return \eqn{\mathbf{V}^{-1/2}} such that \code{VhalfInv_fxn(par) \%*\% VhalfInv_fxn(par)} = \eqn{\mathbf{V}^{-1}}.
#' @param Vhalf_fxn Default: NULL. Function for efficient computation of \eqn{\mathbf{V}^{1/2}}, used only when optimizing correlation structures with non-canonical-Gaussian response.
#' @param VhalfInv_par_init Default: \code{c()}. Numeric vector of initial parameter values for \code{VhalfInv_fxn} optimization, on unbounded transformed scale (not the natural scale). When provided with \code{VhalfInv_fxn}, triggers optimization of the covariance structure. Length determines the dimension of the parameter space. Defaults to 0 if not provided.
#' @param REML_grad Default: NULL. Function for evaluating the gradient of the objective function (negative REML or custom loss) with respect to the parameters of \code{VhalfInv_fxn}. Must take the same \code{"par"} argument as \code{VhalfInv_fxn}, as well as second argument \code{"model_fit"} for the output of \code{lgspline.fit} and ellipses \code{"..."} as a third argument. It should return a vector of partial derivatives matching the length of par. When provided, enables more efficient optimization via closed-form gradients rather than numerical approximation. Optional — if NULL, BFGS uses numerical gradients.
#' @param custom_VhalfInv_loss Default: NULL. Alternative to negative REML for serving as the objective function for optimizing correlation parameters. Must take the same \code{"par"} argument as \code{VhalfInv_fxn}, as well as second argument \code{"model_fit"} for the output of \code{lgspline.fit} and ellipses \code{"..."} as a third argument. It should return a numeric scalar.
#' @param VhalfInv_logdet Default: NULL. Function for efficient computation of \eqn{\log|\mathbf{V}^{-1/2}|} that bypasses construction of the full \eqn{\mathbf{V}^{-1/2}} matrix. Must take the same parameter vector \code{'par'} as \code{VhalfInv_fxn} and return a scalar value equal to \eqn{\log(\det(\mathbf{V}^{-1/2}))}. When NULL, the determinant is computed directly from \code{VhalfInv}, which can be computationally expensive for large matrices.
#' @param include_warnings Default: TRUE. Logical switch to control display of warning messages during model fitting.
#' @param ... Additional arguments passed to the unconstrained model fitting function.
#'
#' @return A list of class \code{"lgspline"} containing the following components:
#' \describe{
#'   \item{y}{Original response vector.}
#'   \item{ytilde}{Fitted/predicted values on the scale of the response.}
#'   \item{X}{List of design matrices \eqn{\mathbf{X}_{k}} for each partition k, containing basis expansions including intercept, linear, quadratic, cubic, and interaction terms as specified. Returned on the unstandardized scale.}
#'   \item{A}{Constraint matrix \eqn{\mathbf{A}} encoding smoothness constraints at knot points and any user-specified linear constraints. Only a linearly independent subset of columns is retained (via pivoted QR decomposition).}
#'   \item{B}{List of fitted coefficients \eqn{\boldsymbol{\beta}_{k}} for each partition k on the original, unstandardized scale of the predictors and response.}
#'   \item{B_raw}{List of fitted coefficients for each partition on the predictor-and-response standardized scale.}
#'   \item{K}{Number of interior knots with one predictor (number of partitions minus 1 with > 1 predictor).}
#'   \item{p}{Number of basis expansions of predictors per partition.}
#'   \item{q}{Number of predictor variables.}
#'   \item{P}{Total number of coefficients (\eqn{p \times (K+1)}).}
#'   \item{N}{Number of observations.}
#'   \item{penalties}{List containing optimized penalty matrices and components:
#'     \itemize{
#'       \item Lambda: Combined penalty matrix (\eqn{\boldsymbol{\Lambda}}), includes \eqn{\mathbf{L}_{\mathrm{predictor\_list}}} contributions but not \eqn{\mathbf{L}_{\mathrm{partition\_list}}}.
#'       \item L1: Smoothing spline penalty matrix (\eqn{\mathbf{L}_{1}}).
#'       \item L2: Ridge penalty matrix (\eqn{\mathbf{L}_{2}}).
#'       \item L_predictor_list: Predictor-specific penalty matrices (\eqn{\mathbf{L}_{\mathrm{predictor\_list}}}).
#'       \item L_partition_list: Partition-specific penalty matrices (\eqn{\mathbf{L}_{\mathrm{partition\_list}}}).
#'     }
#'   }
#'   \item{knot_scale_transf}{Function for transforming predictors to standardized scale used for knot placement.}
#'   \item{knot_scale_inv_transf}{Function for transforming standardized predictors back to original scale.}
#'   \item{knots}{Matrix of knot locations on original unstandardized predictor scale for one predictor.}
#'   \item{partition_codes}{Vector assigning observations to partitions.}
#'   \item{partition_bounds}{Vector or matrix specifying the boundaries between partitions.}
#'   \item{knot_expand_function}{Internal function for expanding data according to partition structure.}
#'   \item{predict}{Function for generating predictions on new data. For multi-predictor models, \code{take_first_derivatives = TRUE, take_second_derivatives} returns derivatives as a named list of components per predictor variable, rather than a concatenated vector.}
#'   \item{assign_partition}{Function for assigning new observations to partitions.}
#'   \item{family}{GLM family object specifying the error distribution and link function.}
#'   \item{estimate_dispersion}{Logical indicating whether dispersion parameter was estimated.}
#'   \item{unbias_dispersion}{Logical indicating whether dispersion estimates should be unbiased.}
#'   \item{backtransform_coefficients}{Function for converting standardized coefficients to original scale.}
#'   \item{forwtransform_coefficients}{Function for converting coefficients to standardized scale.}
#'   \item{mean_y, sd_y}{Mean and standard deviation of response if standardized.}
#'   \item{og_order}{Original ordering of observations before partitioning.}
#'   \item{order_list}{List containing observation indices for each partition.}
#'   \item{constraint_values, constraint_vectors}{Matrices specifying linear equality constraints if provided.}
#'   \item{make_partition_list}{List containing partition information for > 1-D cases.}
#'   \item{expansion_scales}{Vector of scaling factors used for standardizing basis expansions.}
#'   \item{take_derivative, take_interaction_2ndderivative}{Functions for computing derivatives of basis expansions.}
#'   \item{get_all_derivatives_insample}{Function for computing all derivatives on training data.}
#'   \item{numerics}{Indices of numeric predictors used in basis expansions.}
#'   \item{power1_cols, power2_cols, power3_cols, power4_cols}{Column indices for linear through quartic terms.}
#'   \item{quad_cols}{Column indices for all quadratic terms (including interactions).}
#'   \item{interaction_single_cols, interaction_quad_cols}{Column indices for linear-linear and linear-quadratic interactions.}
#'   \item{triplet_cols}{Column indices for three-way interactions.}
#'   \item{nonspline_cols}{Column indices for terms excluded from spline expansion.}
#'   \item{return_varcovmat}{Logical indicating whether variance-covariance matrix was computed.}
#'   \item{raw_expansion_names}{Names of basis expansion terms.}
#'   \item{std_X, unstd_X}{Functions for standardizing/unstandardizing design matrices.}
#'   \item{parallel_cluster_supplied}{Logical indicating whether a parallel cluster was supplied.}
#'   \item{weights}{List of observation weights per partition.}
#'   \item{G}{List of unscaled unconstrained variance-covariance matrices \eqn{\mathbf{G}_{k}} per partition k if \code{return_G = TRUE}. Computed as \eqn{(\mathbf{X}_{k}^{\top}\mathbf{X}_{k} + \boldsymbol{\Lambda}_{\mathrm{eff}})^{-1}} for partition k. When \code{VhalfInv} is non-\code{NULL}, recomputed from the whitened Gram matrices \eqn{\mathbf{X}_{k}^{\top}\mathbf{V}^{-1}\mathbf{X}_{k}}.}
#'   \item{Ghalf}{List of \eqn{\mathbf{G}_{k}^{1/2}} matrices if \code{return_Ghalf = TRUE}. When \code{VhalfInv} is non-\code{NULL}, the full \eqn{\mathbf{G}_{\mathrm{correct}}^{1/2}} from the penalized GLS information is used for posterior draws and variance-covariance computation.}
#'   \item{U}{Constraint projection matrix \eqn{\mathbf{U}} if \code{return_U = TRUE}. For K=0 and no constraints, returns identity. Otherwise, returns \eqn{\mathbf{U} = \mathbf{I} - \mathbf{G}\mathbf{A}(\mathbf{A}^{\top}\mathbf{G}\mathbf{A})^{-1}\mathbf{A}^{\top}}. Used for computing the variance-covariance matrix \eqn{\sigma^{2}\mathbf{U}\mathbf{G}}.}
#'   \item{sigmasq_tilde}{Estimated (or fixed) dispersion parameter \eqn{\tilde{\sigma}^{2}}. When \code{VhalfInv} is non-\code{NULL} and \code{family} is Gaussian with identity link, computed from whitened residuals as \eqn{\mathrm{mean}((\mathbf{V}^{-1/2}\mathbf{y} - \mathbf{V}^{-1/2}\hat{\mathbf{y}})^{2})} with optional bias correction by multiplying by \eqn{\frac{N}{N-\mathrm{trace}(\mathbf{X}\mathbf{U}\mathbf{G}\mathbf{X}^{\top}}}.}
#'   \item{trace_XUGX}{Effective degrees of freedom (\eqn{\mathrm{trace}(\mathbf{X}\mathbf{U}\mathbf{G}\mathbf{X}^{\top})}), where \eqn{\mathbf{X}\mathbf{U}\mathbf{G}\mathbf{X}^{\top}} serves as the "hat" matrix. When \code{VhalfInv} is non-\code{NULL}, computed as \eqn{\|\mathbf{V}^{-1/2}\mathbf{X}\mathbf{U}\mathbf{G}_{\mathrm{correct}}^{1/2}\|_{F}^{2}} using the full penalized GLS information.}
#'   \item{varcovmat}{Variance-covariance matrix of coefficient estimates if \code{return_varcovmat = TRUE}. Computed as \eqn{\sigma^{2}(\mathbf{U}\mathbf{G}^{1/2})(\mathbf{U}\mathbf{G}^{1/2})^{\top}} for numerical stability. When \code{VhalfInv} is non-\code{NULL}, uses the full \eqn{\mathbf{G}_{\mathrm{correct}}^{1/2}} in place of the block-diagonal \eqn{\mathbf{G}^{1/2}}.}
#'   \item{lagrange_multipliers}{Vector of Lagrangian multipliers. If \code{return_lagrange_multipliers = TRUE} and constraints are present. Each element corresponds to a column of \eqn{\mathbf{A}} and quantifies the sensitivity of the penalized objective to relaxing the associated constraint. \code{NULL} if no constraints are active (\eqn{\mathbf{A}} is \code{NULL} or \code{K == 0}).}
#'   \item{VhalfInv}{The \eqn{\mathbf{V}^{-1/2}} matrix used for implementing correlation structures, if specified.}
#'   \item{VhalfInv_fxn, Vhalf_fxn, VhalfInv_logdet, REML_grad}{Functions for generating \eqn{\mathbf{V}^{-1/2}}, \eqn{\mathbf{V}^{1/2}}, \eqn{\log|\mathbf{V}^{-1/2}|}, and gradient of REML if provided.}
#'   \item{VhalfInv_params_estimates}{Vector of estimated correlation parameters when using \code{VhalfInv_fxn}.}
#'   \item{VhalfInv_params_vcov}{Approximate variance-covariance matrix of estimated correlation parameters from BFGS optimization.}
#'   \item{wald_univariate}{Function for computing univariate Wald statistics and confidence intervals. Returns an S3 object of class \code{wald_lgspline} with dedicated \code{print}, \code{summary}, and \code{plot} methods.}
#'   \item{critical_value}{Critical value used for confidence interval construction.}
#'   \item{generate_posterior}{Function for drawing from the posterior distribution of coefficients. When \code{VhalfInv} is non-\code{NULL}, draws are from the correct joint posterior \eqn{\mathbf{U}\mathbf{G}_{\mathrm{correct}}^{1/2}\mathbf{z}} using the full penalized GLS information, reflecting cross-partition posterior covariance induced by off-diagonal blocks of \eqn{\mathbf{V}^{-1/2}}.}
#'   \item{find_extremum}{Function for optimizing the fitted function. Accepts both numeric column indices and character column names for \code{vars}. When \code{select_vars_fl = TRUE}, L-BFGS-B bounds are correctly subsetted to the optimized variables.}
#'   \item{plot}{Function for visualizing fitted curves.}
#'   \item{quadprog_list}{List containing quadratic programming components if applicable.}
#'   \item{.fit_call_args}{List containing the arguments passed to \eqn{\code{\lgspline}}.}
#' }
#'
#' The returned object has class \code{"lgspline"} and provides comprehensive tools for
#' model interpretation, inference, prediction, and visualization. All
#' coefficients and predictions can be transformed between standardized and
#' original scales using the provided transformation functions. The object includes
#' both frequentist and Bayesian inference capabilities through Wald statistics
#' and posterior sampling. S3 methods \code{\link{logLik.lgspline}} and
#' \code{\link{confint.lgspline}} are available for standard log-likelihood
#' extraction and confidence interval computation, respectively.
#' Advanced customization options are available for
#' analyzing arbitrarily complex study designs.
#'
#' @examples
#'
#' ## ## ## ## Simple Examples ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#' ## Simulate some data, fit using default settings, and plot
#' set.seed(1234)
#' t <- runif(2500, -10, 10)
#' y <- 2*sin(t) + -0.06*t^2 + rnorm(length(t))
#' model_fit <- lgspline(t, y)
#' plot(t, y, main = 'Observed Data vs. Fitted Function, Colored by Partition',
#'      ylim = c(-10, 10))
#' plot(model_fit, add = TRUE)
#' points(t, antiderivative(model_fit)[[1]], cex = 0.5, col = 'grey')
#' legend('topright',
#'        col = c('grey'),
#'        lty = 1,
#'        cex = 0.5,
#'        legend = c('Antiderivative (Up to a Constant of Integration)'))
#'
#' ## Repeat using logistic regression, with univariate inference shown
#' # and alternative function call
#' y <- rbinom(length(y), 1, 1/(1+exp(-std(y))))
#' df <- data.frame(t = t, y = y)
#' model_fit <- lgspline(y ~ spl(t),
#'                       df,
#'                       opt = FALSE, # no tuning penalties
#'                       family = quasibinomial())
#' plot(t, y, main = 'Observed Data vs Fitted Function with Formulas and Derivatives',
#'   ylim = c(-0.5, 1.05), cex.main = 0.8)
#' plot(model_fit,
#'      show_formulas = TRUE,
#'      text_size_formula = 0.65,
#'      legend_pos = 'bottomleft',
#'      legend_args = list(y.intersp = 1.1),
#'      add = TRUE)
#' ## Notice how the coefficients match the formula, and expansions are
#' # homogenous across partitions without reparameterization
#' print(summary(model_fit))
#'
#' ## Overlay first and second derivatives of fitted function respectively
#' derivs <- predict(model_fit,
#'                   new_predictors = sort(t),
#'                   take_first_derivatives = TRUE,
#'                   take_second_derivatives = TRUE)
#' points(sort(t), derivs$first_deriv, col = 'gold', type = 'l')
#' points(sort(t), derivs$second_deriv, col = 'goldenrod', type = 'l')
#' legend('bottomright',
#'        col = c('gold','goldenrod'),
#'        lty = 1,
#'        legend = c('First Derivative', 'Second Derivative'))
#'
#' ## Simple 2D example - including a non-spline effect
#' z <- seq(-2, 2, length.out = length(y))
#' df <- data.frame(Predictor1 = t,
#'                  Predictor2 = z,
#'                  Response = sin(y)+0.1*z)
#' model_fit <- lgspline(Response ~ spl(Predictor1) + Predictor1*Predictor2,
#'                       df)
#'
#' ## Notice, while spline effects change over partitions,
#' # interactions and non-spline effects are constrained to remain the same
#' coefficients <- Reduce('cbind', coef(model_fit))
#' colnames(coefficients) <- paste0('Partition ', 1:(model_fit$K+1))
#' print(coefficients)
#'
#' ## One or two variables can be selected for plotting at a time
#' # even when >= 3 predictors are present
#' plot(model_fit,
#'       custom_title = 'Marginal Relationship of Predictor 1 and Response',
#'       vars = 'Predictor1',
#'       custom_response_lab = 'Response',
#'       show_formulas = TRUE,
#'       legend_pos = 'bottomright',
#'       digits = 4,
#'       text_size_formula = 0.5)
#'
#' ## 3D plots are implemented as well, retaining closed-formulas
#' my_plot <- plot(model_fit,
#'                 show_formulas = TRUE,
#'                 custom_response_lab = 'Response')
#' my_plot
#'
#'
#' ## ## ## ## More Detailed 1D Example ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#' ## 1D data generating functions
#' t <- seq(-9, 9, length.out = 1000)
#' slinky <- function(x) {
#'   (50 * cos(x * 2) -2 * x^2 + (0.25 * x)^4 + 80)
#' }
#' coil <- function(x) {
#'   (100 * cos(x * 2) +-1.5 * x^2 + (0.1 * x)^4 +
#'   (0.05 * x^3) + (-0.01 * x^5) +
#'      (0.00002 * x^6) -(0.000001 * x^7) + 100)
#' }
#' exponential_log <- function(x) {
#'   unlist(c(sapply(x, function(xx) {
#'     if (xx <= 1) {
#'       100 * (exp(xx) - exp(1))
#'     } else {
#'       100 * (log(xx))
#'     }
#'   })))
#' }
#' scaled_abs_gamma <- function(x) {
#'   2*sqrt(gamma(abs(x)))
#' }
#'
#' ## Composite function
#' fxn <- function(x)(slinky(t) +
#'                    coil(t) +
#'                    exponential_log(t) +
#'                    scaled_abs_gamma(t))
#'
#' ## Bind together with random noise
#' dat <- cbind(t, fxn(t) + rnorm(length(t), 0, 50))
#' colnames(dat) <- c('t', 'y')
#' x <- dat[,'t']
#' y <- dat[,'y']
#'
#' ## Fit Model, 4 equivalent ways are shown below
#' model_fit <- lgspline(t, y)
#' model_fit <- lgspline(y ~ spl(t), as.data.frame(dat))
#' model_fit <- lgspline(response = y, predictors = t)
#' model_fit <- lgspline(data = as.data.frame(dat), formula = y ~ .)
#'
#' # This is not valid: lgspline(y ~ ., t)
#' # This is not valid: lgspline(y, data = as.data.frame(dat))
#' # Do not put operations in formulas, not valid: lgspline(y ~ log(t) + spl(t))
#'
#' ## Basic Functionality
#' predict(model_fit, new_predictors = rnorm(1)) # make prediction on new data
#' suppressWarnings(head(leave_one_out(model_fit))) # leave-one-out cross-validated predictions will return NA if leverage is too high
#' coef(model_fit) # extract coefficients
#' summary(model_fit) # model information and Wald inference
#' generate_posterior(model_fit) # generate draws of parameters from posterior distribution
#' find_extremum(model_fit, minimize = TRUE) # find the minimum of the fitted function
#'
#' ## Incorporate range constraints, custom knots, keep penalization identical
#' # across partitions
#' model_fit <- lgspline(y ~ spl(t),
#'                       unique_penalty_per_partition = FALSE,
#'                       custom_knots = cbind(c(-2, -1, 0, 1, 2)),
#'                       data = data.frame(t = t, y = y),
#'                       qp_range_lower = -150,
#'                       qp_range_upper = 150)
#'
#' ## Plotting the constraints and knots
#' plot(model_fit,
#'      custom_title = 'Fitted Function Constrained to Lie Between (-150, 150)',
#'      cex.main = 0.75)
#' # knot locations
#' abline(v = model_fit$knots)
#' # lower bound from quadratic program
#' abline(h = -150, lty = 2)
#' # upper bound from quadratic program
#' abline(h = 150, lty = 2)
#' # observed data
#' points(t, y, cex = 0.24)
#'
#' ## Enforce monotonic increasing constraints on fitted values
#' # K = 4 => 5 partitions
#' t <- seq(-10, 10, length.out = 100)
#' y <- 5*sin(t) + t + 2*rnorm(length(t))
#' model_fit <- lgspline(t,
#'                       y,
#'                       K = 4,
#'                       qp_monotonic_increase = TRUE)
#' plot(t, y, main = 'Monotonic Increasing Function with Respect to Fitted Values')
#' plot(model_fit,
#'      add = TRUE,
#'      show_formulas = TRUE,
#'      legend_pos = 'bottomright',
#'      custom_predictor_lab = 't',
#'      custom_response_lab = 'y')
#'
#' ## ## ## ## 2D Example using Volcano Dataset ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#' ## Prep
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
#' # When fitting with > 1 predictor and large K, including quartic terms
#' # is highly recommended, and/or dropping the second-derivative constraint.
#' # Otherwise, the constraints can impose all partitions to be equal, with one
#' # cubic function fit for all (there is not enough degrees of freedom to fit
#' # unique cubic functions due to the massive amount of constraints).
#' # Below, quartic terms are included and the constraint of second-derivative
#' # smoothness at knots is ignored.
#' model_fit <- lgspline(volcano_long[,c(1, 2)],
#'                       volcano_long[,3],
#'                       include_quadratic_interactions = TRUE,
#'                       K = 49,
#'                       opt = FALSE,
#'                       return_U = FALSE,
#'                       return_varcov = FALSE,
#'                       estimate_variance = TRUE,
#'                       return_Ghalf = FALSE,
#'                       return_G = FALSE,
#'                       include_constrain_second_deriv = FALSE,
#'                       unique_penalty_per_predictor = FALSE,
#'                       unique_penalty_per_partition = FALSE,
#'                       wiggle_penalty = 1e-10, # the fixed wiggle penalty
#'                       flat_ridge_penalty = 0) # the ridge penalty / wiggle penalty
#'
#' ## Plotting on new data with interactive visual + formulas
#' new_input <- expand.grid(seq(min(volcano_long[,1]),
#'                              max(volcano_long[,1]),
#'                              length.out = 250),
#'                          seq(min(volcano_long[,2]),
#'                              max(volcano_long[,2]),
#'                              length.out = 250))
#' plot(model_fit,
#'      new_predictors = new_input,
#'      show_formulas = TRUE,
#'      custom_response_lab = "Height",
#'      custom_title = 'Volcano 3-D Map',
#'      digits = 2)
#'
#' ## ## ## ## Advanced Techniques using Trees Dataset ## ## ## ## ## ## ## ## ## ## ## ## ##
#' ## Goal here is to introduce how lgspline works with non-canonical GLMs and
#' # demonstrate some custom features
#' data('trees')
#'
#' ## L1-regularization constraint function on standardized coefficients
#' # Bound all coefficients to be less than a certain value (l1_bound) in absolute
#' # magnitude such that | B^{(j)}_k | < lambda for all j = 1....p coefficients,
#' # and k = 1...K+1 partitions.
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
#' ## Bounds absolute value of coefficients to be < l1_bound
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
#' # Gamma-distributed response, and custom quadratic-programming constraints,
#' # with qp_score_function/glm_weight_function updated for non-canonical GLMs
#' # as well as quartic terms, keeping the effect of height constant across
#' # partitions, and 3 partitions in total. Hence, this is an advanced-usage
#' # case.
#' # You can modify this code for performing l1-regularization in general.
#' # For canonical GLMs, the default qp_score_function/glm_weight_function are
#' # correct and do not need to be changed
#' # (custom functionality is not needed for canonical GLMs).
#' model_fit <- lgspline(
#'   Volume ~ spl(Girth) + Height*Girth,
#'   data = with(trees, cbind(Girth, Height, Volume)),
#'   family = Gamma(link = 'log'),
#'   keep_weighted_Lambda = TRUE,
#'   glm_weight_function = function(
#'     mu,
#'     y,
#'     order_indices,
#'     family,
#'     dispersion,
#'     observation_weights,
#'    ...){
#'      rep(1/dispersion, length(y))
#'    },
#'    dispersion_function = function(
#'      mu,
#'      y,
#'      order_indices,
#'      family,
#'      observation_weights,
#'      VhalfInv,
#'    ...){
#'     mean(
#'       mu^2/((y-mu)^2)
#'     )
#'   }, # = biased estimate of 1/shape parameter
#'   need_dispersion_for_estimation = TRUE,
#'   unbias_dispersion = TRUE, # multiply dispersion by N/(N-trace(XUGX^{T}))
#'   K = 2, # 3 partitions
#'   opt = FALSE, # keep penalties fixed
#'   unique_penalty_per_partition = FALSE,
#'   unique_penalty_per_predictor = FALSE,
#'   flat_ridge_penalty = 1e-64,
#'   wiggle_penalty = 1e-64,
#'   qp_score_function = function(X, y, mu, order_list, dispersion, VhalfInv,
#'     observation_weights, ...){
#'    t(X) %*% diag(c(1/mu * 1/dispersion)) %*% cbind(y - mu)
#'   }, # updated score for gamma regression with log link
#'   qp_Amat_fxn = function(N, p, K, X, colnm, scales, deriv_fxn, ...) {
#'     l1_constraint_matrix(p, K)
#'   },
#'   qp_bvec_fxn = function(qp_Amat, N, p, K, X, colnm, scales, deriv_fxn, ...) {
#'     l1_bound_vector(qp_Amat, scales, 25)
#'   },
#'   qp_meq_fxn = function(qp_Amat, N, p, K, X, colnm, scales, deriv_fxn, ...) 0
#' )
#'
#' ## Notice, interaction effect is constant across partitions as is the effect
#' # of Height alone
#' Reduce('cbind', coef(model_fit))
#'
#' ## Many constraints, many coefficients, and small sample size makes inference
#' #  using asymptotic variance-covariance matrix untrustworthy.
#' #  Confidence intervals are often too wide or narrow, even for "good" fit.
#' #  Consider bootstrapping or alternative.
#' print(summary(model_fit))
#'
#' ## Plot results
#' plot(model_fit, custom_predictor_lab1 = 'Girth',
#'      custom_predictor_lab2 = 'Height',
#'      custom_response_lab = 'Volume',
#'      custom_title = 'Girth and Height Predicting Volume of Trees',
#'      show_formulas = TRUE)
#'
#' ## Verify magnitude of unstandardized coefficients does not exceed bound (25)
#' print(max(abs(unlist(model_fit$B))))
#'
#' ## Find height and girth where tree volume is closest to 42
#' # Uses custom objective that minimizes MSE discrepancy between predicted
#' # value and 42.
#' # The vanilla find_extremum function can be thought of as
#' # using "function(mu)mu" aka the identity function as the
#' # objective, where mu = "f(t)", our estimated function. The derivative is then
#' # d_mu = "df/dt" with respect to predictors t.
#' # But with more creative objectives, and since we have machinery for
#' # df/dt already available, we can compute gradients for (and optimize)
#' # arbitrary differentiable functions of our predictors too.
#' # For any objective, differentiate w.r.t. to mu, then multiply by d_mu to
#' # satisfy chain rule.
#' # Here, we have objective function: 0.5*(mu-42)^2
#' # and gradient                    : (mu-42)*d_mu
#' # and L-BFGS-B will be used to find the height and girth that most closely
#' # yields a prediction of 42 within the bounds of the observed data.
#' # The d_mu also takes into account link function transforms automatically
#' # for most common link functions, and will return warning + instructions
#' # on how to program the link-function derivatives otherwise.
#'
#' ## Custom acquisition functions for Bayesian optimization could be coded here.
#' find_extremum(
#'   model_fit,
#'   minimize = TRUE,
#'   custom_objective_function = function(mu, sigma, ybest, ...){
#'     0.5*(mu - 42)^2
#'   },
#'   custom_objective_derivative = function(mu, sigma, ybest, d_mu, ...){
#'     (mu - 42) * d_mu
#'   }
#' )
#'
#' ## ## ## ## How to Use Formulas in lgspline ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#' ## Demonstrates splines with multiple mixed predictors and interactions
#'
#' ## Generate data
#' n <- 2500
#' x <- rnorm(n)
#' y <- rnorm(n)
#' z <- sin(x)*mean(abs(y))/2
#'
#' ## Categorical predictors
#' cat1 <- rbinom(n, 1, 0.5)
#' cat2 <- rbinom(n, 1, 0.5)
#' cat3 <- rbinom(n, 1, 0.5)
#'
#' ## Response with mix of effects
#' response <- y + z + 0.1*(2*cat1 - 1)
#'
#' ## Continuous predictors re-named
#' continuous1 <- x
#' continuous2 <- z
#'
#' ## Combine data
#' dat <- data.frame(
#'   response = response,
#'   continuous1 = continuous1,
#'   continuous2 = continuous2,
#'   cat1 = cat1,
#'   cat2 = cat2,
#'   cat3 = cat3
#' )
#'
#' ## Example 1: Basic Model with Default Terms, No Intercept
#' # standardize_response = FALSE often needed when constraining intercepts to 0
#' fit1 <- lgspline(
#'   formula = response ~ 0 + spl(continuous1, continuous2) +
#'     cat1*cat2*continuous1 + cat3,
#'   K = 2,
#'   standardize_response = FALSE,
#'   data = dat
#' )
#' ## Examine coefficients included
#' rownames(fit1$B$partition1)
#' ## Verify intercept term is near 0 up to some numeric tolerance
#' abs(fit1$B[[1]][1]) < 1e-8
#'
#' ## Example 2: Similar Model with Intercept, Other Terms Excluded
#' fit2 <- lgspline(
#'   formula = response ~ spl(continuous1, continuous2) +
#'     cat1*cat2*continuous1 + cat3,
#'   K = 1,
#'   standardize_response = FALSE,
#'   include_cubic_terms = FALSE,
#'   exclude_these_expansions = c( # Not all need to actually be present
#'     '_batman_x_robin_',
#'     '_3_x_4_', # no cat1 x cat2 interaction, coded using column indices
#'     'continuous1xcontinuous2', # no continuous1 x continuous2 interaction
#'     'thejoker'
#'   ),
#'   data = dat
#' )
#' ## Examine coefficients included
#' rownames(Reduce('cbind',coef(fit2)))
#' # Intercept will probably be present and non-0 now
#' abs(fit2$B[[1]][1]) < 1e-8
#'
#' ## ## ## ## Compare Inference to survreg for Weibull AFT Model Validation ##
#' # Only linear predictors, no knots, no penalties, using Weibull AFT Model
#' # The goal here is to ensure that for the special case of no spline effects
#' # and no knots, this implementation will be consistent with other model
#' # implementations.
#' # Also note, that when using models (like Weibull AFT) where dispersion is
#' # being estimated and is required for estimating beta coefficients,
#' # we use a schur complement correction function to adjust (or "correct") our
#' # variance-covariance matrix for both estimation and inference to account for
#' # uncertainty in estimating the dispersion.
#' # Typically the schur_correction_function would return a negative-definite
#' # matrix, as its output is elementwise added to the information matrix prior
#' # to inversion.
#' require(survival)
#' df <- data.frame(na.omit(
#'   pbc[,c('time','trt','stage','hepato','bili','age','status')]
#' ))
#'
#' ## Weibull AFT using lgspline, showing how some custom options can be used to
#' # fit more complicated models
#' model_fit <- lgspline(time ~ trt + stage + hepato + bili + age,
#'                       df,
#'                       family = weibull_family(),
#'                       need_dispersion_for_estimation = TRUE,
#'                       dispersion_function = weibull_dispersion_function,
#'                       glm_weight_function = weibull_glm_weight_function,
#'                       schur_correction_function = weibull_schur_correction,
#'                       unconstrained_fit_fxn = unconstrained_fit_weibull,
#'                       opt = FALSE,
#'                       wiggle_penalty = 0,
#'                       flat_ridge_penalty = 0,
#'                       K = 0,
#'                       status = pbc$status!=0)
#' print(summary(model_fit))
#'
#' ## Survreg results match closely on estimates and inference for coefficients
#' survreg_fit <- survreg(Surv(time, status!=0) ~ trt + stage + hepato + bili + age,
#'                        df)
#' print(summary(survreg_fit))
#'
#' ## sigmasq_tilde = scale^2 of survreg
#' print(c(sqrt(model_fit$sigmasq_tilde), survreg_fit$scale))
#'
#' ## ## ## ## Modelling Correlation Structures ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#' ## Setup
#' n_blocks <- 150 # Number of correlation_ids (subjects)
#' block_size <- 5 # Size of each correlation_ids (number of repeated measures per subj.)
#' N <- n_blocks * block_size # total sample size (balanced here)
#' rho_true <- 0.25  # True correlation
#'
#' ## Generate predictors and mean structure
#' t <- seq(-9, 9, length.out = N)
#' true_mean <- sin(t)
#'
#' ## Create block compound symmetric errors = I(1-p) + Jp
#' errors <- Reduce('rbind',
#'                  lapply(1:n_blocks,
#'                         function(i){
#'                           sigma <- diag(block_size) + rho_true *
#'                             (matrix(1, block_size, block_size) -
#'                                diag(block_size))
#'                           matsqrt(sigma) %*% rnorm(block_size)
#'                         }))
#'
#' ## Generate response with correlated errors
#' y <- true_mean + errors * 0.5
#'
#' ## Fit model with correlation structure
#' # include_warnings = FALSE is a good idea here, since many proposed
#' # correlations will not work
#' model_fit <- lgspline(t,
#'                       y,
#'                       K = 4,
#'                       correlation_id = rep(1:n_blocks, each = block_size),
#'                       correlation_structure = 'exchangeable',
#'                       include_warnings = FALSE
#' )
#'
#' ## Assess overall fit
#' plot(t, y, main = 'Sinosudial Fit Under Correlation Structure')
#' plot(model_fit, add = TRUE, show_formulas = TRUE, custom_predictor_lab = 't')
#'
#' ## Compare estimated vs true correlation
#' # Built-in exchangeable uses rho = exp(-exp(par)), so par in (-Inf, Inf)
#' # maps to rho in (0, 1). Only positive correlation is supported.
#' rho_est <- exp(-exp(model_fit$VhalfInv_params_estimates))
#' print(c("True correlation:" = rho_true,
#'         "Estimated correlation:" = rho_est))
#'
#' ## Quantify uncertainty in correlation estimate with 95% confidence interval
#' # CI is constructed on the working scale and back-transformed
#' se <- c(sqrt(diag(model_fit$VhalfInv_params_vcov))) / sqrt(model_fit$N)
#' ci <- exp(-exp(model_fit$VhalfInv_params_estimates + c(-1.96, 1.96)*se))
#' # Note: the transformation reverses order since exp(-exp(.)) is decreasing
#' ci <- rev(ci)
#' print("95% CI for correlation:")
#' print(ci)
#'
#' ## Also check SD (should be close to 0.5)
#' print(sqrt(model_fit$sigmasq_tilde))
#'
#' ## Toeplitz Simulation Setup, with demonstration of custom functions
#' # and boilerplate. Toep is not implemented by default, because it makes
#' # strong assumptions on the study design and missingness that are rarely met,
#' # with non-obvious workarounds.
#' # If a GLM was to-be-fit, you would also submit a function "Vhalf_fxn" analogous
#' # to VhalfInv_fxn with same argument (par) and an output of an N x N matrix
#' # that yields the inverse of VhalfInv_fxn output.
#' n_blocks <- 150   # Number of correlation_ids
#' block_size <- 4   # Observations per correlation_id
#' N <- n_blocks * block_size # total sample size
#' rho_true <- 0.5    # True correlation within correlation_ids
#' true_intercept <- 2     # True intercept
#' true_slope <- 0.5       # True slope for covariate
#'
#' ## Create design matrix with meaningful predictors
#' Tmat <- matrix(0, N, 2)
#' Tmat[,1] <- 1  # Intercept
#' Tmat[,2] <- cos(rnorm(N))  # Continuous predictor
#'
#' ## True coefficients
#' beta <- c(true_intercept, true_slope)
#'
#' ## Create time and correlation_id variables
#' time_var <- rep(1:block_size, n_blocks)
#' correlation_id_var <- rep(1:n_blocks, each = block_size)
#'
#' ## Create block compound symmetric errors
#' errors <- Reduce('rbind',
#'                  lapply(1:n_blocks, function(i) {
#'                    sigma <- diag(block_size) + rho_true *
#'                      (matrix(1, block_size, block_size) -
#'                         diag(block_size))
#'                    matsqrt(sigma) %*% rnorm(block_size)
#'                  }))
#'
#' ## Generate response with correlated errors and covariate effect
#' y <- Tmat %*% beta + errors * 2
#'
#' ## Toeplitz correlation function
#' # Custom functions can use any parameterization. Here we use
#' # rho_lag = exp(-exp(par[lag])) to keep correlations in (0, 1),
#' # consistent with the built-in structures.
#' VhalfInv_fxn <- function(par) {
#'   # Initialize correlation matrix
#'   corr <- matrix(0, block_size, block_size)
#'
#'   # Construct Toeplitz matrix with correlation by lag
#'   for(i in 1:block_size) {
#'     for(j in 1:block_size) {
#'       lag <- abs(time_var[i] - time_var[j])
#'       if(lag == 0) {
#'         corr[i,j] <- 1
#'       } else if(lag <= length(par)) {
#'         corr[i,j] <- exp(-exp(par[lag]))
#'       }
#'     }
#'   }
#'
#'   ## Matrix square root inverse
#'   corr_inv_sqrt <- matinvsqrt(corr)
#'
#'   ## Expand to full matrix using Kronecker product
#'   kronecker(diag(n_blocks), corr_inv_sqrt)
#' }
#'
#' ## Determinant function (for efficiency)
#' # This avoids taking determinant of N by N matrix
#' VhalfInv_logdet <- function(par) {
#'   # Initialize correlation matrix
#'   corr <- matrix(0, block_size, block_size)
#'
#'   # Construct Toeplitz matrix
#'   for(i in 1:block_size) {
#'     for(j in 1:block_size) {
#'       lag <- abs(time_var[i] - time_var[j])
#'       if(lag == 0) {
#'         corr[i,j] <- 1
#'       } else if(lag <= length(par)) {
#'         corr[i,j] <- exp(-exp(par[lag]))
#'       }
#'     }
#'   }
#'
#'   # Compute log determinant
#'   log_det_invsqrt_corr <- -0.5 * determinant(corr, logarithm=TRUE)$modulus[1]
#'   return(n_blocks * log_det_invsqrt_corr)
#' }
#'
#' ## REML gradient function
#' # The gradient has three terms per parameter (see ?`.reml_grad_from_dV`):
#' #   (i)   0.5 * tr(V^{-1} dV/dpar) => log-determinant of V
#' #   (ii)  -0.5/sigma^2 * t(r) dV/dpar r => residual quadratic form
#' #   (iii) -0.5 * tr(M^+ t(X*)V^{-1}dV V^{-1}X*) => REML correction
#' # where X* = XU is the constrained design and
#' # M = t(X*)V^{-1}X* + t(U) Lambda U is the projected penalized information.
#' REML_grad <- function(par, model_fit, ...) {
#'   ## Initialize gradient vector
#'   n_par <- length(par)
#'   gradient <- numeric(n_par)
#'
#'   ## Process derivatives one parameter at a time
#'   for(p in 1:n_par) {
#'     ## Initialize derivative matrix
#'     dV <- matrix(0, nrow(model_fit$VhalfInv), ncol(model_fit$VhalfInv))
#'     V <- matrix(0, nrow(model_fit$VhalfInv), ncol(model_fit$VhalfInv))
#'
#'     ## Compute full correlation matrix and its derivative for parameter p
#'     for(clust in unique(correlation_id_var)) {
#'       inds <- which(correlation_id_var == clust)
#'       block_size <- length(inds)
#'
#'       ## Initialize block matrices
#'       V_block <- matrix(0, block_size, block_size)
#'       dV_block <- matrix(0, block_size, block_size)
#'
#'       ## Construct Toeplitz matrix and its derivative
#'       for(i in 1:block_size) {
#'         for(j in 1:block_size) {
#'           lag <- abs(time_var[i] - time_var[j])
#'
#'           if(i == j) {
#'             V_block[i,j] <- 1
#'             dV_block[i,j] <- 0
#'           } else {
#'             if(lag <= length(par)) {
#'               ## Correlation via exp(-exp(.)) parameterization
#'               rho_lag <- exp(-exp(par[lag]))
#'               V_block[i,j] <- rho_lag
#'
#'               ## Derivative for the relevant parameter only
#'               ## d/dpar exp(-exp(par)) = -exp(par) * exp(-exp(par))
#'               if(lag == p) {
#'                 dV_block[i,j] <- -exp(par[p]) * rho_lag
#'               }
#'             }
#'           }
#'         }
#'       }
#'
#'       ## Assign blocks to full matrices
#'       V[inds, inds] <- V_block
#'       dV[inds, inds] <- dV_block
#'     }
#'
#'     ## GLM Weights (all 1s for Gaussian identity link)
#'     glm_weights <- rep(1, model_fit$N)
#'
#'     ## (ii) Quadratic form contribution
#'     resid <- model_fit$y - model_fit$ytilde
#'     VinvResid <- model_fit$VhalfInv %*% cbind(resid) / glm_weights
#'     quad_term <- -0.5 * ((t(VinvResid) %*% dV) %*% VinvResid) /
#'       model_fit$sigmasq_tilde
#'
#'     ## (i) Log|V| contribution
#'     trace_term <- 0.5 * sum((model_fit$VhalfInv %*% dV) *
#'                                model_fit$VhalfInv)
#'
#'     ## (iii) REML correction on constrained design X* = XU
#'     U <- t(t(model_fit$U) * rep(c(1, model_fit$expansion_scales),
#'                                 model_fit$K + 1)) /
#'       model_fit$sd_y
#'     VhalfInvX <- model_fit$VhalfInv %*%
#'       collapse_block_diagonal(model_fit$X)[unlist(
#'         model_fit$og_order
#'       ),] %*%
#'       U
#'
#'     ## Build projected penalty t(U) Lambda U
#'     if(length(model_fit$penalties$L_partition_list) != (model_fit$K + 1)){
#'       model_fit$penalties$L_partition_list <- lapply(
#'         1:(model_fit$K + 1), function(k)0
#'       )
#'     }
#'     Lambda_raw <- collapse_block_diagonal(
#'       lapply(1:(model_fit$K + 1),
#'              function(k)
#'                c(1, model_fit$expansion_scales) * (
#'                  model_fit$penalties$L_partition_list[[k]] +
#'                    model_fit$penalties$Lambda) %*%
#'                diag(c(1, model_fit$expansion_scales)) /
#'                model_fit$sd_y^2
#'       )
#'     )
#'     ULambdaU <- t(U) %*% Lambda_raw %*% U
#'
#'     ## Projected penalized information inverse
#'     XVinvX_inv <- invert(crossprod(t(t(VhalfInvX) * c(glm_weights))) +
#'                            ULambdaU)
#'     VInvX <- model_fit$VhalfInv %*% VhalfInvX
#'     sc <- sqrt(norm(VInvX, '2'))
#'     VInvX <- VInvX/sc
#'     dXVinvX <-
#'       (XVinvX_inv %*% t(VInvX)) %*%
#'       (dV %*% VInvX)
#'     XVinvX_term <- -0.5 * colSums(cbind(c(diag(dXVinvX) * sc))) * sc
#'
#'     ## Store gradient component (all three terms)
#'     gradient[p] <- as.numeric(quad_term + trace_term + XVinvX_term)
#'   }
#'
#'   ## Return normalized gradient
#'   return(gradient / model_fit$N)
#' }
#'
#' ## Visualization
#' plot(time_var, y, col = correlation_id_var,
#'      main = "Simulated Data with Toeplitz Correlation")
#'
#' ## Fit model with custom Toeplitz (takes ~ 5-10 minutes on my laptop)
#' # Note, for GLMs, efficiency would be improved by supplying a Vhalf_fxn
#' # although strictly only VhalfInv_fxn and VhalfInv_par_init are needed
#' model_fit <- lgspline(
#'   response = y,
#'   predictors = Tmat[,2],
#'   K = 4,
#'   VhalfInv_fxn = VhalfInv_fxn,
#'   VhalfInv_logdet = VhalfInv_logdet,
#'   REML_grad = REML_grad,
#'   VhalfInv_par_init = log(-log(c(0.5, 0.5, 0.5))),
#'   include_warnings = FALSE
#' )
#'
#' ## Print comparison of true and estimated correlations
#' rho_true <- rep(0.5, 3)
#' rho_est <- exp(-exp(model_fit$VhalfInv_params_estimates))
#' cat("Correlation Estimates:\n")
#' print(data.frame(
#'   "True Correlation" = rho_true,
#'   "Estimated Correlation" = rho_est
#' ))
#'
#' ## Should be ~ 2
#' print(sqrt(model_fit$sigmasq_tilde))
#'
#' ## ## ## ## Parallelism ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#' ## Data generating function
#' a <- runif(500000, -9, 9)
#' b <- runif(500000, -9, 9)
#' c <- rnorm(500000)
#' d <- rpois(500000, 1)
#' y <- sin(a) + cos(b) - 0.2*sqrt(a^2 + b^2) +
#'   abs(a) + b +
#'   0.5*(a^2 + b^2) +
#'   (1/6)*(a^3 + b^3) +
#'   a*b*c -
#'   c +
#'   d +
#'   rnorm(500000, 0, 5)
#'
#' ## Set up cores
#' cl <- parallel::makeCluster(1)
#' on.exit(parallel::stopCluster(cl))
#'
#' ## This example shows some options for what operations can be parallelized
#' # By default, only parallel_eigen and parallel_unconstrained are TRUE
#' # parallel_unconstrained is only for GLMs, for identity link Gaussian
#' # response, use parallel_matmult=TRUE to ensure parallel fitting across
#' # partitions.
#' # G, G^{-1/2}, and G^{1/2} are computed in parallel across each of the
#' # K+1 partitions.
#' # However, parallel_unconstrained only affects GLMs without corr. components
#' # - it does not affect fitting here
#' system.time({
#'   parfit <- lgspline(y ~ spl(a, b) + a*b*c + d,
#'                      data = data.frame(y = y,
#'                                        a = a,
#'                                        b = b,
#'                                        c = c,
#'                                        d = d),
#'                      cl = cl,
#'                      K = 1,
#'                      parallel_eigen = TRUE,
#'                      parallel_unconstrained = TRUE,
#'                      parallel_aga = FALSE,
#'                      parallel_find_neighbors = FALSE,
#'                      parallel_trace = FALSE,
#'                      parallel_matmult = TRUE,
#'                      parallel_make_constraint = FALSE,
#'                      parallel_penalty = FALSE)
#' })
#' parallel::stopCluster(cl)
#' print(summary(parfit))
#'
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{lgspline.fit}} for the low-level fitting interface
#'   \item \code{\link{logLik.lgspline}} for log-likelihood extraction
#'   \item \code{\link{confint.lgspline}} for confidence interval extraction
#'   \item \code{\link{leave_one_out}} for leave-one-out cross-validated predictions
#'   \item \code{\link{blockfit_solve}} for the standalone backfitting solver
#'   \item \code{\link[quadprog]{solve.QP}} for quadratic programming optimization
#'   \item \code{\link[plotly]{plot_ly}} for interactive plotting
#'   \item \code{\link[stats]{kmeans}} for k-means clustering
#'   \item \code{\link[stats]{optim}} for general purpose optimization routines
#' }
#'
#' @export
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
    schur_correction_function = function(X,
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
                                   VhalfInv,
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
    initial_wiggle = exp(c(-25, -14, -3, -0.2, 1)),
    initial_flat = c(0.1, 10),
    wiggle_penalty = 2e-7,
    flat_ridge_penalty = 0.5,
    unique_penalty_per_partition = TRUE,
    unique_penalty_per_predictor = TRUE,
    meta_penalty = 1e-8,
    predictor_penalties = NULL,
    partition_penalties = NULL,
    include_quadratic_terms = TRUE,
    include_cubic_terms = TRUE,
    include_quartic_terms = NULL,
    include_2way_interactions = TRUE,
    include_3way_interactions = TRUE,
    include_quadratic_interactions = FALSE,
    offset = c(),
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
    parallel_unconstrained = TRUE,
    parallel_find_neighbors = FALSE,
    parallel_penalty = FALSE,
    parallel_make_constraint = FALSE,
    unconstrained_fit_fxn = unconstrained_fit_default,
    keep_weighted_Lambda = FALSE,
    iterate_tune = TRUE,
    iterate_final_fit = TRUE,
    blockfit = FALSE,
    qp_score_function = function(X, y, mu, order_list, dispersion, VhalfInv, observation_weights, ...) {
      if(!is.null(observation_weights)) {
        crossprod(X, cbind((y - mu)*observation_weights))
      } else {
        crossprod(X, cbind(y - mu))
      }
    },
    qp_observations = NULL,
    qp_Amat = NULL,
    qp_bvec = NULL,
    qp_meq = 0,
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
    constraint_values = cbind(),
    constraint_vectors = cbind(),
    return_G = TRUE,
    return_Ghalf = TRUE,
    return_U = TRUE,
    estimate_dispersion = TRUE,
    unbias_dispersion = NULL,
    return_varcovmat = TRUE,
    return_lagrange_multipliers = FALSE,
    custom_penalty_mat = NULL,
    cluster_args = c(custom_centers = NA, nstart = 10),
    dummy_dividor = 0.00000000000000000000012345672152894,
    dummy_adder = 0.000000000000000002234567210529,
    verbose = FALSE,
    verbose_tune = FALSE,
    dummy_fit = FALSE,
    auto_encode_factors = TRUE,
    observation_weights = NULL,
    do_not_cluster_on_these = c(),
    neighbor_tolerance = 1 + 1e-8,
    null_constraint = NULL,
    critical_value = qnorm(1-0.05/2),
    data = NULL,
    weights = NULL,
    no_intercept = FALSE,
    correlation_id = NULL,
    spacetime = NULL,
    correlation_structure = NULL,
    VhalfInv = NULL,
    Vhalf = NULL,
    VhalfInv_fxn = NULL,
    Vhalf_fxn = NULL,
    VhalfInv_par_init = c(),
    REML_grad = NULL,
    custom_VhalfInv_loss = NULL,
    VhalfInv_logdet = NULL,
    include_warnings = TRUE,
    ...
  ){

  if(verbose){
    cat('Pre-Processing\n')
  }

  ## [Change 2026-02-12] Backward compat: expansions_only in ... maps
  #  to dummy_fit with deprecation warning
  #  shur -> schur for correct spelling
  dots <- list(...)
  if(!is.null(dots$expansions_only) && dots$expansions_only){
    if(include_warnings){
      warning("'expansions_only' is deprecated; use 'dummy_fit = TRUE' instead.",
              " The full pipeline now executes with coefficients set to 0.")
    }
    dummy_fit <- TRUE
  }
  if(!is.null(dots$shur_correction_function)){
    if(include_warnings){
      warning("'shur_correction_function' option is depreciated, and won't be",
              "applied; repeat with 'schur_correction_function' instead.")
    }
  }

  ## Unbiasing dispersion is automatically done for Gaussian + identity only
  # [Change 2026-02-13] Additional documentation,
  # multiplies the MLE estimate of variance,
  # or exponential dispersion, by N/(N-tr) where
  # tr  = effective degrees of freedom
  if(paste0(family)[1] == 'gaussian' &
     paste0(family)[2] == 'identity' &
     is.null(unbias_dispersion)){
    unbias_dispersion <- TRUE
  } else if(is.null(unbias_dispersion)){
    unbias_dispersion <- FALSE
  }

  ## Update naming conventions, if first argument is a formula and second is a
  # data frame, assumed by R-like interfaces for user convenience.
  if(any(!is.null(predictors)) & any(!is.null(y))){
    if(any(inherits(y,'data.frame') & inherits(predictors, "formula"))){
      data <- y
    }
  }

  ## Update naming conventions, if response supplied in place of "y" for user
  # convenience.
  if(any(is.null(predictors)) & any(!is.null(formula))){
    predictors <- formula
  } else if(any(is.null(predictors))){
    stop('\n \t Predictors argument is NULL without formula supplied.',
    ' Either supply a formula to predictors OR formula argument, or a',
    ' data frame to predictors argument, or a matrix of numeric predictors',
    ' to the predictor argument. \n')
  }

  ## Update naming conventions, if response supplied in place of "y" for user
  # convenience.
  if(any(is.null(y)) & any(!is.null(response))){
    y <- response
    response <- NULL
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
    if(any(inherits(ncluster, 'try-error'))){
      stop('\n \t custom_centers should be a matrix; do not include any other',
      ' arguments within cluster_args if you include custom_centers. \n')
    }
    if(!is.null(K)){
      if(ncluster != (K+1) & include_warnings){
        warning('\n \t K must be equal to number of custom_centers minus 1. ',
        'Updating K for compatibility. \n')
        K <- ncluster - 1
      }
    } else {
      K <- ncluster - 1
    }
  }

  ## Check data and formula argument
  if(!is.null(data) & !inherits(predictors, "formula")) {
    stop("\n \t If submitting data argument, formula must be supplied and',
         ' variables must match. Otherwise, use predictors and y (or response)',
         ' arguments directly.\n",
         "\t Example: lgspline(y ~ spl(x1, x2) + x3 + x4*x5, data = my_data)\n",
         "\t Example: lgspline(y ~ ., data = my_data)\n")
  }

  ## Handle formula interface
  if(inherits(predictors, "formula")) {
    ## Check data argument
    if(is.null(data)) {
      stop("\n \t When using formula interface, data argument must be ",
           "provided.\n",
           "\t Example: lgspline(y ~ spl(x1, x2) + x3 + x4*x5, data = my_data)\n",
           "\t Example: lgspline(y ~ ., data = my_data)\n")
    }

    ## Try to coerce data to data.frame
    tryCatch({
      data <- as.data.frame(data)
    }, error = function(e) {
      stop("\n \t Could not coerce data argument to data.frame. ",
           "Please provide data in a format coercible to data.frame. ",
           "Examples: data.frame, tibble, or matrix.")
    })

    ## [Change 2026-02-12] Auto one-hot encode factor/character columns
    # For reviewers: the decision to include one-hot rather than
    # dummy-intercept is for interpretation as conditional random effects
    # models. We would not like this to drop a random effect level,
    # The flat_ridge_penalty should address lack of identifiability to an extent
    # Users can manually supply the variable names after one-hot encoding
    # if dummy-intercept is truly desired.
    if(auto_encode_factors && !is.null(data)){
      non_numeric_cols <- which(!sapply(data, is.numeric))
      ## Exclude response column from encoding
      resp_name <- NULL
      if(inherits(predictors, "formula")){
        resp_name <- as.character(predictors[[2]])
      }
      non_numeric_cols <- non_numeric_cols[
        !names(non_numeric_cols) %in% resp_name
      ]
      if(length(non_numeric_cols) > 0){
        for(col_idx in rev(non_numeric_cols)){
          col_name <- colnames(data)[col_idx]
          ## One-hot encode this column
          encoded <- create_onehot(data[[col_name]])
          colnames(encoded) <- paste0(col_name, "_", colnames(encoded))
          ## Remove original, append indicators
          data <- cbind(data[, -col_idx, drop = FALSE], encoded)
          ## Track these as linear-without-interactions only if the
          #  original column was not involved in any interaction term
          new_names <- colnames(encoded)
          form_check <- paste(deparse(predictors, width.cutoff = 500L),
                              collapse = "")
          col_in_interaction <- grepl(
            paste0("(\\*|:)\\s*", col_name, "(\\s|\\)|\\+|$)"), form_check
          ) || grepl(
            paste0("(\\+|\\(|~|\\s)", col_name, "\\s*(\\*|:)"), form_check
          )
          if(!col_in_interaction){
            if(is.null(just_linear_without_interactions)){
              just_linear_without_interactions <- new_names
            } else {
              just_linear_without_interactions <- c(
                just_linear_without_interactions, new_names
              )
            }
          }
          ## Update formula if needed: replace col_name with encoded names
          if(inherits(predictors, "formula")){
            form_str <- deparse(predictors)
            replacement <- paste0("(", paste(new_names, collapse = " + "), ")")
            form_str <- gsub(col_name, replacement, form_str, fixed = TRUE)
            predictors <- as.formula(form_str)
          }
        }
      }
    }

    ## Check column names are present
    if(any(is.null(colnames(data)))){
      stop('\n \t Column names of data must be supplied if formula is',
      ' supplied, and column names must match what is provided in formula. \n')
    }

    ## Stringed formula
    form_paste0 <- paste0(predictors)

    ## If formula is y ~ ., replace with y ~ spl(x1, x2, ...)
    if(form_paste0[1] == '~' &
       (gsub(' ', '', form_paste0[3]) %in% c('.', '0+.','.+0')) &
       length(form_paste0) == 3){
         if(form_paste0[2] == '.'){
           stop('\n \t . ~ . is not valid for this function, specify the',
           ' response directly, like "y ~ . \n"')
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

    ## If offset() appears in the formula, remove offset from the formula
    # and set offset = names of terms inside the offset() operator
    if(any(grepl('offset', paste0(predictors)))){
      if(any(grepl("offset\\(", deparse(predictors)))){
        ## Extract the original formula as string
        form_str <- deparse(predictors)

        ## Find all offset terms using regex
        offset_matches <- gregexpr("offset\\(([^)]+)\\)", form_str)
        offset_vars <- regmatches(form_str, offset_matches)[[1]]

        ## Extract variable names from offset terms
        offset_var_names <- gsub("offset\\((.*)\\)", "\\1", offset_vars)

        ## Remove whitespace
        offset_var_names <- trimws(offset_var_names)

        ## Add to offset vector
        offset <- c(offset, offset_var_names)

        ## Replace offset(var) with var in formula
        for(i in seq_along(offset_vars)) {
          form_str <- gsub(offset_vars[i], offset_var_names[i],
                           form_str,
                           fixed=TRUE)
        }

        ## Convert back to formula
        predictors <- as.formula(form_str)
      }
    }

    ## Check that all formula predictors are numeric in data
    # [Change 2026-02-14] Update to also process
    # one-hot factor and characters
    terms <- terms(predictors)
    term_labels <- attr(terms, "term.labels")
    non_numeric <- c(1:ncol(data))[!sapply(data, is.numeric)]
    if(length(non_numeric) > 0) {
      for(v in non_numeric){
        if(any(term_labels == colnames(data)[v])){
          if(auto_encode_factors){
            stop("Auto-encoding failed for column '",
                 colnames(data)[v], "'. Check that it is a factor ",
                 "or character vector.")
          } else {
            stop("Non-numeric columns detected in predictors. ",
                 "Set auto_encode_factors = TRUE to automatically ",
                 "one-hot encode factor/character columns, or ",
                 "convert these to numeric before proceeding. ",
                 "You can use create_onehot() on categorical ",
                 "columns of your dataset to obtain binary ",
                 "indicator variables if you need to. Remember ",
                 "to remove the original categorical variable ",
                 "and append the indicators to your data ",
                 "before adjusting your formula and call ",
                 "to lgspline(...).\n\nSee ?create_onehot for an example.")
          }
        }
      }
    }

    ## Parse rest of formula
    has_3way <- any(attr(terms, "order") == 3)
    if(include_3way_interactions) include_3way_interactions <- has_3way

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

    ## Extract explicit interactions from formula
    formula_interactions <- term_labels[attr(terms, "order") > 1]
    allowed_interactions <- sapply(formula_interactions, function(term) {
      if(grepl(":", term)) {
        vars <- strsplit(term, ":")[[1]]
        # Only keep interactions where NO term is a spline term
        if(!any(vars %in% spline_terms)) {
          return(term)
        }
      }
      return(NULL)
    })

    ## Identify interaction structure using factors matrix
    var_names <- rownames(factors)[-1] # Remove responses

    ## Extract all types of variables in formula
    all_formula_vars <- unique(c(spline_terms, var_names))
    formula_cols <- which(colnames(data) %in% all_formula_vars)

    ## Match custom exclusions nominal to numeric formula columns, if available
    if(inherits(exclude_interactions_for, 'character')){
      exclude_interactions_for <-
        unlist(lapply(exclude_interactions_for,
                                         function(var){
                                            grep(var,
                                                 colnames(data)[formula_cols])
                                         }))

    }
    if(inherits(just_linear_with_interactions,  'character')){
      just_linear_with_interactions <-
        unlist(lapply(just_linear_with_interactions,
                          function(var){
                            grep(var,
                                 colnames(data)[formula_cols])
                          }))

    }
    if(inherits(just_linear_without_interactions, 'character')){
      just_linear_without_interactions <-
        unlist(lapply(just_linear_without_interactions,
                      function(var){
                        grep(var,
                             colnames(data)[formula_cols])
                      }))

    }
    if(length(exclude_these_expansions) > 0){
      for(ii in 1:length(exclude_these_expansions)){
        exps <- exclude_these_expansions[ii]
        if(substr(exps, 1, 1) != "_" |
           substr(exps, length(exps), length(exps)) != "_"){
          for(jj in 1:length(formula_cols)){
            if(grepl(colnames(data)[formula_cols[jj]],
                     exclude_these_expansions[ii])){
              exclude_these_expansions[ii] <- gsub(
                colnames(data)[formula_cols[jj]],
                paste0('_', jj, '_'),
                exclude_these_expansions[ii]
              )
            }
          }
        }
      }
    }
    if(inherits(offset,  'character')){
      offset <-
        unlist(lapply(offset,
                      function(var){
                        grep(var,
                             colnames(data)[formula_cols])
                      }))
    }

    ## [Change 2026-02-17] Resolve character do_not_cluster_on_these to numeric
    #  column indices within the formula path. Without this, make_partitions'
    #  assign_partition closure fails when indexing matrices after transf()
    #  strips colnames. Mirrors resolution logic for offset,
    #  exclude_interactions_for.
    if(length(do_not_cluster_on_these) > 0 &&
       is.character(do_not_cluster_on_these)){
      do_not_cluster_on_these <-
        unlist(lapply(do_not_cluster_on_these,
                      function(var){
                        idx <- which(colnames(data)[formula_cols] == var)
                        if(length(idx) == 0){
                          idx <- grep(var, colnames(data)[formula_cols])
                        }
                        idx
                      }))
      do_not_cluster_on_these <- unique(do_not_cluster_on_these)
    }

    ## Non-spline variables that appear in interactions (from factors matrix)
    interaction_terms <- term_labels[attr(terms, "order") > 1]
    nonspline_interact_vars <- unique(sapply(interaction_terms,
                                                    function(term) {
                                                      if(grepl(":", term)) {
                                                        vars <-
                                                         strsplit(term,":")[[1]]
                                                        vars[!vars %in%
                                                         spline_terms]
                                                      }
                                                    }))

    ## [Change 2026-02-15] prevents errors, since this may not be declared yet
    if(!exists("allowed_interaction_pairs")){
      allowed_interaction_pairs <- list()
    }

    ## Add explicit interactions for spline terms
    # [Change 2026-02-12] 3-way interactions now correctly generated
    # for spl(a, b, c) with 3 continuous variables.
    # This was previously not working in earlier versions < 1.0.
    if(length(spline_terms) > 1) {
      ## Generate all possible interactions between spline terms
      spline_interactions <- utils::combn(spline_terms, 2, simplify=FALSE)
      spline_triplets <- if(length(spline_terms) >= 3 &
                            include_3way_interactions) {
        utils::combn(spline_terms, 3, simplify=FALSE)
      } else {
        list()
      }

      ## Add 2-way interactions
      for(pair in spline_interactions) {
        interaction_terms <- c(interaction_terms, paste(pair, collapse=":"))
      }

      ## Add 3-way interactions
      for(triplet in spline_triplets) {
        interaction_terms <- c(interaction_terms, paste(triplet, collapse=":"))
      }

      ## Also register spline triplets as allowed interaction pairs
      ## so they are not excluded by the grouping logic below
      for(triplet in spline_triplets) {
        allowed_interaction_pairs[[length(allowed_interaction_pairs) + 1]] <-
          list(pair = triplet,
               transforms = c(
                 paste(triplet, collapse = "x"),
                 paste(rev(triplet), collapse = "x")
               )
          )
      }
    }

    ## Get indices of variables in raw expansions (after response removed)
    resp_ind <- which(colnames(data) == paste0(terms[[2]]))
    var_positions <- match(var_names, colnames(data[,formula_cols]))

    ## Get allowed interaction pairs
    allowed_pairs <- lapply(interaction_terms[grepl(":", interaction_terms)],
                            function(term) {
                              vars <- strsplit(term, ":")[[1]]
                              match(vars, colnames(data[,formula_cols]))
                            })

    ## Generate exclusion patterns for non-spline vars
    exclude_patterns <- c()

    ## Get all possible 2-way interactions between ANY variables
    vars <- colnames(data[,formula_cols])
    if(include_2way_interactions){
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
                ## Get interaction terms that could involve these variables
                relevant_terms <- interaction_terms[grepl(paste(triplet_vars,
                                                                collapse="|"),
                                                          interaction_terms)]

                ## Check if this exact triplet exists in any order
                matches_interaction <- any(sapply(strsplit(relevant_terms, ":"),
                                                  function(term) {
                                                    length(term) == 3 &&
                                          all(sort(triplet_vars) == sort(term))
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

    ## Helper-vector for the next step, in isolating which interactions are
    # not part of the same grouping/block.
    # The above only works if we have only 1 additive interaction effect.
    # the code below corrects for when we have multiple (for example, a*b + c*d)
    # and want to prevent unspecified interactions (for example, remove a*d and b*c)
    different_grouping_exclusions <- c()

    ## Explicitly track allowed interactions from * and : terms
    # [Change 2026-02-12] Safe initialization, may already exist from
    # spline triplet registration
    if(!exists("allowed_interaction_pairs")){
      allowed_interaction_pairs <- list()
    }

    ## Track which interactions are explicitly specified
    for(term in term_labels) {
      # Handle * interactions
      if(grepl("\\*", term)) {
        vars <- trimws(strsplit(term, "\\*")[[1]])

        # Ensure exactly 2 variables in the interaction
        if(length(vars) == 2) {
          allowed_interaction_pairs[[length(allowed_interaction_pairs) + 1]] <-
            list(pair = vars,
                 transforms = c(
                   paste(vars, collapse = "x"),
                   paste(rev(vars), collapse = "x")
                 )
            )
        }
      }

      # Handle : interactions
      if(grepl(":", term)) {
        vars <- trimws(strsplit(term, ":")[[1]])

        # Ensure exactly 2 variables in the interaction
        if(length(vars) == 2) {
          allowed_interaction_pairs[[length(allowed_interaction_pairs) + 1]] <-
            list(pair = vars,
                 transforms = c(
                   paste(vars, collapse = "x"),
                   paste(rev(vars), collapse = "x")
                 )
            )
        }
      }
    }

    ## Identify spline block variables
    spline_block_vars <- c()
    for(term in term_labels) {
      if(grepl("^spl\\(.*\\)$", term)) {
        vars <- gsub("^spl\\((.*)\\)$", "\\1", term)
        term_vars <- trimws(strsplit(vars, ",")[[1]])
        spline_block_vars <- c(spline_block_vars, term_vars)
      }
    }

    ## Generate exclusions for any interaction not in the allowed pairs
    vars <- colnames(data[,formula_cols])
    if(include_2way_interactions | include_3way_interactions){
      for(ii in seq_along(vars)) {
        for(jj in seq_along(vars)) {
          if(ii != jj) {

            pair <- c(vars[ii], vars[jj])
            current_transform <- paste(pair, collapse = "x")

            ## Skip if not in formula
            if(any(!(pair %in% all_formula_vars))){
              next
            }

            ## Skip if any variable in exclude_interactions_for
            if(!is.null(exclude_interactions_for)){
              if(any(c(ii, jj) %in% exclude_interactions_for)){
                next
              }
            }

            ## Check if this is an allowed interaction
            is_allowed_interaction <- FALSE
            for(allowed in allowed_interaction_pairs) {
              if(setequal(pair, allowed$pair) ||
                 any(current_transform == allowed$transforms)) {
                is_allowed_interaction <- TRUE
                break
              }
            }

            ## Check if this is a spline block interaction (2-way or 3-way)
            # [Change 2026-02-14] added clarification in comment
            # for both 2-way and 3-way being included here, no other change
            is_spline_block_interaction <- all(pair %in% spline_block_vars)

            ## If not an allowed or spline block interaction,
            # generate exclusion patterns
            if(!is_allowed_interaction && !is_spline_block_interaction) {
              patterns <- get_interaction_patterns(pair)
              different_grouping_exclusions <- c(different_grouping_exclusions,
                                                 patterns)
            }
          }
        }
      }
    }

    ## Remove duplicates
    different_grouping_exclusions <- unique(different_grouping_exclusions)

    ## Remove explicitly allowed interactions from exclusions
    for(term in interaction_terms) {
      vars <- strsplit(term, ":")[[1]]
      allowed <- get_interaction_patterns(vars)
      exclude_patterns <- setdiff(exclude_patterns, allowed)
    }
    exclude_patterns <- unique(c(exclude_patterns,
                                 different_grouping_exclusions))

    ## Convert to positional notation
    vars <- colnames(data)[formula_cols]
    for(ii in seq_along(formula_cols)) {
      exclude_patterns <- gsub(vars[ii], paste0("_", ii, "_"), exclude_patterns)
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
        # its in an interaction
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
    predictors <- data[, formula_cols, drop = FALSE]
    y <- data[, resp_ind]

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
    stop("\n \t Could not coerce predictors to matrix. ",
         "predictors must be either a formula or an object coercible to matrix.",
         "Examples:\n",
         "  Formula: lgspline(y ~ spl(x1, x2) + x3, data = my_data)\n",
         "  Matrix:  lgspline(predictors = Tmat, y = y) \n")
  })

  ## Check numeric type
  if(any(!is.numeric(predictors))){
    stop("\n \t predictors matrix must be numeric. ",
         "Please convert categorical variables to numeric indicators. \n")
  }

  ## Check response for missings
  if(any(is.na(y) | is.nan(y) | !is.finite(y))){
    stop("\n \t NA, NaN, or infinite value detected in response. \n")
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
     length(constraint_values) == 0){
     # constraint_values <-
     #   crossprod(t(constraint_vectors),
     #   crossprod(t(invert(crossprod(cbind(constraint_vectors)))),
     #   cbind(c(null_constraint))))
     ## [Change 2026-02-14] implemented alternative method following
     # suggestions about efficient matrix multiplications
     # Original code from version < 1.0 included
     constraint_values <-
      crossprod(
        t(constraint_vectors),
        crossprod(
          t(invert(crossprod(cbind(constraint_vectors)))),
          cbind(c(null_constraint))
        )
      )
  }

  ## Check nrow of input predictors and matrix coersion
  t <- try({if(nrow(methods::as(predictors,'matrix')) < 3){
    stop('\n \t Need at least 3 observations to fit model \n')
  }}, silent = TRUE)
  if(inherits(t, 'try-error')){
    stop('\n \t Cannot coerce predictors to a matrix \n')
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
    if(any(inherits(custom_knots, 'try-error')) & include_warnings){
      warning('\n \t custom_knots must be a matrix, or should be coercible ',
      'to it. Custom_knots will be ignored. \n')
      custom_knots <- NULL
    }
  }

  ## If ncol(predictors) > 1 and include_quartic_terms is NULL,
  # set include_quartic_terms = TRUE
  # otherwise, if ncol(predictors) == 1 then FALSE
  if(is.null(include_quartic_terms)){
    if(ncol(cbind(predictors)) == 1){
      include_quartic_terms <- FALSE
    } else {
      include_quartic_terms <- TRUE
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
                                 schur_correction_function,
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
                                 initial_wiggle,
                                 initial_flat,
                                 wiggle_penalty,
                                 flat_ridge_penalty,
                                 unique_penalty_per_partition,
                                 unique_penalty_per_predictor,
                                 meta_penalty,
                                 predictor_penalties,
                                 partition_penalties,
                                 include_quadratic_terms,
                                 include_cubic_terms,
                                 include_quartic_terms,
                                 include_2way_interactions,
                                 include_3way_interactions,
                                 include_quadratic_interactions,
                                 offset,
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
                                 blockfit,
                                 qp_score_function,
                                 qp_observations,
                                 qp_Amat,
                                 qp_bvec,
                                 qp_meq,
                                 qp_positive_derivative,
                                 qp_negative_derivative,
                                 qp_positive_2ndderivative,
                                 qp_negative_2ndderivative,
                                 qp_monotonic_increase,
                                 qp_monotonic_decrease,
                                 qp_range_upper,
                                 qp_range_lower,
                                 qp_Amat_fxn,
                                 qp_bvec_fxn,
                                 qp_meq_fxn,
                                 constraint_values,
                                 constraint_vectors,
                                 return_G,
                                 return_Ghalf,
                                 return_U,
                                 estimate_dispersion,
                                 unbias_dispersion,
                                 return_varcovmat,
                                 return_lagrange_multipliers,
                                 custom_penalty_mat,
                                 cluster_args,
                                 dummy_dividor,
                                 dummy_adder,
                                 verbose,
                                 verbose_tune,
                                 dummy_fit | (!is.null(VhalfInv_fxn) &
                                                length(VhalfInv_par_init) > 0),
                                 # expansions_only |
                                 # (!is.null(VhalfInv_fxn) &
                                 #    length(VhalfInv_par_init) > 0),
                                 ## [Change 2026-02-12] dummy_fit parameter
                                 # instead of expansions_only
                                 auto_encode_factors,
                                 observation_weights,
                                 do_not_cluster_on_these,
                                 neighbor_tolerance,
                                 no_intercept,
                                 VhalfInv,
                                 Vhalf,
                                 include_warnings,
                                 og_cols,
                                 ...
  )}, silent = TRUE)

  ## Return try error if model fails to to be fit
  if(any(inherits(model_fit, 'try-error')) & include_warnings){
    warning('\n \t Model fitting error: try verbose = TRUE, checking for NAs,',
            ' adjusting starting tuning grid, or K. If using parallel options,',
            ' check your parallel cluster, submit to cl argument if valid, and make ',
            'sure base R parallel package is loaded. Also make sure your ',
            'formula and overall setup is valid. \n')
    return(model_fit)
  }

  ## [Change 2026-02-12] Return dummy_fit output (replaces expansions_only)
  if(dummy_fit & is.null(VhalfInv_fxn)){
    return(c(model_fit, list(og_cols = colnames(predictors))))
  }

  ## Default correlation structures
  # [Change 2026-02-14] Simplified the code and verified correctness, added
  # some documentation for reviewer ease.
  # See correlation helper functions
  # .compute_dist_block, .rank_dists, reml_grad_from_dV
  if(!is.null(correlation_structure) &
     !is.null(correlation_id)){if(inherits(correlation_structure, 'character') &
                                  length(correlation_structure) == 1 &
                                  length(correlation_id) == nrow(predictors)){
       if(length(spacetime) > 0){
         spacetime <- cbind(spacetime)
         if(nrow(spacetime) != nrow(predictors)){
           stop('\n\t Spacetime must be an N-length vector or N-row matrix ')
         }
       }

       ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
       # Exchangeable (compound symmetric) correlation
       # V_ij = exp(-exp(par)) for i != j within cluster, 1 on diagonal.
       #  Parameterization: rho = exp(-exp(par)), par in (-Inf, Inf),
       #    rho in (0, 1). Only positive correlation is supported.
       ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
       if(correlation_structure %in% c('exchangeable',
                                       'cs',
                                       'CS',
                                       'compoundsymmetric',
                                       'compound-symmetric',
                                       'compound symmetric')){

         ## Construct V^{-1/2} for exchangeable correlation
         VhalfInv_fxn <- function(par) {
           corr <- matrix(0, nrow(predictors), nrow(predictors))
           rho <- exp(-exp(par))
           for(clust in unique(correlation_id)){
             inds <- which(correlation_id == clust)
             block_size <- length(inds)
             V <- (1 - rho) * diag(block_size) +
               rho * matrix(1, block_size, block_size)
             corr[inds, inds] <- matinvsqrt(V)
           }
           corr
         }

         ## V^{1/2} for non-Gaussian responses
         if(paste0(family)[1] != 'gaussian' |
            paste0(family)[2] != 'identity'){
           Vhalf_fxn <- function(par) {
             corr <- matrix(0, nrow(predictors), nrow(predictors))
             rho <- exp(-exp(par))
             for(clust in unique(correlation_id)){
               inds <- which(correlation_id == clust)
               block_size <- length(inds)
               V <- (1 - rho) * diag(block_size) +
                 rho * matrix(1, block_size, block_size)
               corr[inds, inds] <- matsqrt(V)
             }
             corr
           }
         } else {
           Vhalf_fxn <- NULL
         }

         ## log|V^{-1/2}| via block structure
         VhalfInv_logdet <- function(par) {
           log_det <- 0
           rho <- exp(-exp(par))
           for(clust in unique(correlation_id)){
             inds <- which(correlation_id == clust)
             block_size <- length(inds)
             V <- (1 - rho) * diag(block_size) +
               rho * matrix(1, block_size, block_size)
             log_det <- log_det +
               (-0.5 * determinant(V, logarithm = TRUE)$modulus[1])
           }
           log_det
         }

         ## REML gradient
         ## ## ## ## ## ## ## ##
         #  The negative REML has three terms whose derivatives w.r.t. the
         #  correlation parameter par are computed below:
         #
         #    (i)   0.5 * tr(V^{-1} dV/dpar)
         #          sensitivity of log|V| to the correlation parameter
         #
         #    (ii)  -0.5 / sigma^2 * r' dV/dpar r
         #          sensitivity of the whitened residual quadratic form,
         #             where r = V^{-1/2}(y - mu) / sqrt(W)
         #
         #    (iii) -0.5 * tr( ((XU)'V^{-1}(XU) + Lambda)^{-1}
         #                      (XU)'V^{-1} dV/dpar V^{-1}(XU) )
         #          sensitivity of the REML correction log-determinant.
         #             The constraint projection U is absorbed into the
         #             design matrix before this term is evaluated, so
         #             constraints are handled correctly.
         #
         #  For exchangeable correlation with rho = exp(-exp(par)):
         #    dV/dpar = drho * (J - I),   drho = -exp(par) * rho
         ## ## ## ## ## ## ## ##
         REML_grad <- function(par, model_fit, ...) {
           ## Correlation parameter and chain rule factor
           rho <- exp(-exp(par))
           drho <- -exp(par) * rho

           ## Build V and dV/dpar blockwise
           dV <- matrix(0, nrow(predictors), nrow(predictors))
           V <- dV
           for(clust in unique(correlation_id)){
             inds <- which(correlation_id == clust)
             block_size <- length(inds)
             J_minus_I <- matrix(1, block_size, block_size) - diag(block_size)
             V[inds, inds] <- (1 - rho) * diag(block_size) +
               rho * matrix(1, block_size, block_size)
             dV[inds, inds] <- drho * J_minus_I
           }

           ## (ii) Quadratic form term
           glm_weights <- sqrt(c(glm_weight_function(model_fit$ytilde,
                                                     model_fit$y,
                                                     1:model_fit$N,
                                                     model_fit$family,
                                                     model_fit$sigmasq_tilde,
                                                     rep(1, model_fit$N),
                                                     ...))) /
             sqrt(unlist(model_fit$weights)[model_fit$og_order])
           resid <- model_fit$y - model_fit$ytilde
           VinvResid <- model_fit$VhalfInv %**% cbind(resid) / glm_weights
           quad_term <- -0.5 * ((t(VinvResid) %**% dV) %**% VinvResid) /
             model_fit$sigmasq_tilde

           ## (i) Trace term: 0.5 * tr(V^{-1} dV)
           trace_term <- 0.5 * sum(diag(model_fit$VhalfInv %**%
                                          model_fit$VhalfInv %**%
                                          dV))

           ## (iii) Information matrix term on constrained design XU
           U <- t(t(model_fit$U) * rep(c(1, model_fit$expansion_scales),
                                       model_fit$K + 1)) /
             model_fit$sd_y
           VhalfInvX <- model_fit$VhalfInv %**%
             collapse_block_diagonal(model_fit$X)[unlist(
               model_fit$og_order
             ),] %**%
             U

           ## Penalty in the constrained basis
           if(length(model_fit$penalties$L_partition_list) != (model_fit$K + 1)){
             model_fit$penalties$L_partition_list <- lapply(
               1:(model_fit$K + 1), function(k)0
             )
           }
           ULambdaU <- t(U) %**% collapse_block_diagonal(
             lapply(1:(model_fit$K + 1),
                    function(k)
                      c(1, model_fit$expansion_scales) * (
                        model_fit$penalties$L_partition_list[[k]] +
                          model_fit$penalties$Lambda) %**%
                      diag(c(1, model_fit$expansion_scales)) /
                      model_fit$sd_y^2
             )
           ) %**% U

           ## Penalized information on constrained design
           glm_weights <- sqrt(c(glm_weight_function(model_fit$ytilde,
                                                     model_fit$y,
                                                     1:model_fit$N,
                                                     model_fit$family,
                                                     model_fit$sigmasq_tilde,
                                                     rep(1, model_fit$N),
                                                     ...))) *
             sqrt(unlist(model_fit$weights)[model_fit$og_order])
           ## [Change 2026-02-16] Row-wise weighting: scale row i of VhalfInvX by
           #  glm_weights[i]. The previous t(t(M)*v) scales columns, which
           #  silently recycles the N-length vector across P columns (N != P).
           XVinvX_inv <- invert(crossprod(t(t(VhalfInvX) * c(glm_weights))) +
                                  ULambdaU)
           VInvX <- model_fit$VhalfInv %**% VhalfInvX
           sc <- sqrt(norm(VInvX, '2'))
           VInvX <- VInvX / sc
           dXVinvX <-
             (XVinvX_inv %**% t(VInvX)) %**%
             (dV %**% VInvX)
           XVinvX_term <- -0.5 * colSums(cbind(c(diag(dXVinvX) * sc))) * sc

           as.numeric(quad_term + trace_term + XVinvX_term) /
             nrow(predictors)
         }
       } else if(length(spacetime) == 0){
         stop('\n\t "Spacetime" variable must be supplied if correlation ',
              'structure other than exchangeable is selected. Spacetime can be ',
              'spatial coordinates, time coordinates, or any other longitudinal ',
              'vector/matrix of measurements. \n')
       }

       ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
       # Spatial exponential correlation
       # V_ij = exp(-omega * d_ij),  omega = exp(par) > 0.
       ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
       if(correlation_structure %in% c('spatial-exponential',
                                       'spatialexponential',
                                       'exp',
                                       'exponential')){

         VhalfInv_fxn <- function(par) {
           corr <- matrix(0, nrow(predictors), nrow(predictors))
           omega <- exp(par)
           for(clust in unique(correlation_id)){
             inds <- which(correlation_id == clust)
             diffs <- .compute_dist_block(spacetime, inds)
             V <- exp(-omega * diffs)
             corr[inds, inds] <- matinvsqrt(V)
           }
           return(corr)
         }

         if(paste0(family)[1] != 'gaussian' |
            paste0(family)[2] != 'identity'){
           Vhalf_fxn <- function(par) {
             corr <- matrix(0, nrow(predictors), nrow(predictors))
             omega <- exp(par)
             for(clust in unique(correlation_id)){
               inds <- which(correlation_id == clust)
               diffs <- .compute_dist_block(spacetime, inds)
               V <- exp(-omega * diffs)
               corr[inds, inds] <- matsqrt(V)
             }
             return(corr)
           }
         } else {
           Vhalf_fxn <- NULL
         }

         VhalfInv_logdet <- function(par) {
           log_det <- 0
           omega <- exp(par)
           for(clust in unique(correlation_id)){
             inds <- which(correlation_id == clust)
             diffs <- .compute_dist_block(spacetime, inds)
             V <- exp(-omega * diffs)
             log_det <- log_det +
               (-0.5 * determinant(V, logarithm = TRUE)$modulus[1])
           }
           log_det
         }

         ## REML gradient
         #  dV/dpar = -diffs * V_block * omega   (chain rule: d/dpar exp(par) = exp(par))
         REML_grad <- function(par, model_fit, ...) {
           omega <- exp(par)
           dV <- matrix(0, nrow(predictors), nrow(predictors))
           V <- matrix(0, nrow(predictors), nrow(predictors))

           for(clust in unique(correlation_id)){
             inds <- which(correlation_id == clust)
             diffs <- .compute_dist_block(spacetime, inds)
             V_block <- exp(-omega * diffs)
             diag(V_block) <- 1
             dV_block <- -diffs * V_block * omega
             diag(dV_block) <- 0
             V[inds, inds] <- V_block
             dV[inds, inds] <- dV_block
           }

           reml_grad_from_dV(dV, model_fit, glm_weight_function, ...)
         }

         if(length(VhalfInv_par_init) == 0){
           VhalfInv_par_init <- 0
         }
       }

       ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
       # AR(1) correlation
       #  V_ij = rho^{rank_ij},  rho = exp(-exp(par)) in (0,1).
       #  Ranks are computed from pairwise spacetime distances within cluster.
       ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
       if(correlation_structure %in% c('ar1','ar(1)','AR(1)','AR1')){

         VhalfInv_fxn <- function(par) {
           corr <- matrix(0, nrow(predictors), nrow(predictors))
           rho <- exp(-exp(par))
           for(clust in unique(correlation_id)){
             inds <- which(correlation_id == clust)
             ranked <- .rank_dists(spacetime, inds)
             V <- rho^ranked
             corr[inds, inds] <- matinvsqrt(V)
           }
           corr
         }

         if(paste0(family)[1] != 'gaussian' |
            paste0(family)[2] != 'identity'){
           Vhalf_fxn <- function(par) {
             corr <- matrix(0, nrow(predictors), nrow(predictors))
             rho <- exp(-exp(par))
             for(clust in unique(correlation_id)){
               inds <- which(correlation_id == clust)
               ranked <- .rank_dists(spacetime, inds)
               V <- rho^ranked
               corr[inds, inds] <- matsqrt(V)
             }
             corr
           }
         } else {
           Vhalf_fxn <- NULL
         }

         VhalfInv_logdet <- function(par) {
           log_det <- 0
           rho <- exp(-exp(par))
           for(clust in unique(correlation_id)){
             inds <- which(correlation_id == clust)
             ranked <- .rank_dists(spacetime, inds)
             V <- rho^ranked
             log_det <- log_det +
               (-0.5 * determinant(V, logarithm = TRUE)$modulus[1])
           }
           log_det
         }

         ## REML gradient
         #  rho = exp(-exp(par)),  drho/dpar = -exp(par) * rho
         #  V_ij = rho^r_ij,  dV_ij/dpar = r_ij * rho^{r_ij - 1} * drho/dpar
         REML_grad <- function(par, model_fit, ...) {
           rho <- exp(-exp(par))
           drho <- -exp(par) * rho

           dV <- matrix(0, nrow(predictors), nrow(predictors))
           V <- dV
           for(clust in unique(correlation_id)){
             inds <- which(correlation_id == clust)
             ranked <- .rank_dists(spacetime, inds)
             V[inds, inds] <- rho^ranked
             dV[inds, inds] <- ranked * rho^(ranked - 1) * drho
           }

           reml_grad_from_dV(dV, model_fit, glm_weight_function, ...)
         }
       }

       ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
       # Gaussian / squared-exponential / RBF correlation
       #  V_ij = exp(-d_ij^2 / (2 * ell^2)),  ell = exp(par) > 0.
       ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
       if(correlation_structure %in% c('gaussian', 'rbf', 'squared-exponential')){

         VhalfInv_fxn <- function(par) {
           corr <- matrix(0, nrow(predictors), nrow(predictors))
           ell <- exp(par)
           for(clust in unique(correlation_id)){
             inds <- which(correlation_id == clust)
             diffs <- .compute_dist_block(spacetime, inds)
             V <- exp(-diffs^2 / (2 * ell^2))
             diag(V) <- 1
             corr[inds, inds] <- matinvsqrt(V)
           }
           corr
         }

         if(paste0(family)[1] != 'gaussian' ||
            paste0(family)[2] != 'identity'){
           Vhalf_fxn <- function(par) {
             corr <- matrix(0, nrow(predictors), nrow(predictors))
             ell <- exp(par)
             for(clust in unique(correlation_id)){
               inds <- which(correlation_id == clust)
               diffs <- .compute_dist_block(spacetime, inds)
               V <- exp(-diffs^2 / (2 * ell^2))
               diag(V) <- 1
               corr[inds, inds] <- matsqrt(V)
             }
             corr
           }
         } else {
           Vhalf_fxn <- NULL
         }

         VhalfInv_logdet <- function(par) {
           log_det <- 0
           ell <- exp(par)
           for(clust in unique(correlation_id)){
             inds <- which(correlation_id == clust)
             diffs <- .compute_dist_block(spacetime, inds)
             V <- exp(-diffs^2 / (2 * ell^2))
             diag(V) <- 1
             log_det <- log_det +
               (-0.5 * determinant(V, logarithm = TRUE)$modulus[1])
           }
           log_det
         }

         ## REML gradient
         #  ell = exp(par),  dell/dpar = ell
         #  dV_ij/dell = (d^2 / ell^3) * V_ij
         #  dV_ij/dpar = (d^2 / ell^3) * V_ij * ell = (d^2 / ell^2) * V_ij
         REML_grad <- function(par, model_fit, ...) {
           ell <- exp(par)
           dV <- matrix(0, nrow(predictors), nrow(predictors))
           V <- matrix(0, nrow(predictors), nrow(predictors))

           for(clust in unique(correlation_id)){
             inds <- which(correlation_id == clust)
             diffs <- .compute_dist_block(spacetime, inds)
             V_block <- exp(-diffs^2 / (2 * ell^2))
             diag(V_block) <- 1
             ## Chain rule: d/dpar = d/dell * dell/dpar = d/dell * ell
             dV_block <- (diffs^2 / ell^3) * V_block * ell
             diag(dV_block) <- 0
             V[inds, inds] <- V_block
             dV[inds, inds] <- dV_block
           }

           reml_grad_from_dV(dV, model_fit, glm_weight_function, ...)
         }
       }

       ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
       # Spherical correlation
       #  V_ij = 1 - 1.5*h + 0.5*h^3  for h = d/r <= 1, 0 otherwise.
       #  Range r = exp(par) > 0.
       ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
       if(correlation_structure %in% c('spherical', 'cubic',
                                       'Spherical', 'sphere')) {

         VhalfInv_fxn <- function(par) {
           corr <- matrix(0, nrow(predictors), nrow(predictors))
           range_par <- exp(par)
           for(clust in unique(correlation_id)) {
             inds <- which(correlation_id == clust)
             block_size <- length(inds)
             diffs <- .compute_dist_block(spacetime, inds)
             V <- diag(block_size)
             for(i in 1:block_size) {
               for(j in 1:block_size) {
                 if(i != j) {
                   d_val <- diffs[i, j]
                   if(d_val <= range_par) {
                     h <- d_val / range_par
                     V[i, j] <- 1 - 1.5 * h + 0.5 * h^3
                   }
                 }
               }
             }
             corr[inds, inds] <- matinvsqrt(V)
           }
           corr
         }

         if(paste0(family)[1] != 'gaussian' ||
            paste0(family)[2] != 'identity') {
           Vhalf_fxn <- function(par) {
             corr <- matrix(0, nrow(predictors), nrow(predictors))
             range_par <- exp(par)
             for(clust in unique(correlation_id)) {
               inds <- which(correlation_id == clust)
               block_size <- length(inds)
               diffs <- .compute_dist_block(spacetime, inds)
               V <- diag(block_size)
               for(i in 1:block_size) {
                 for(j in 1:block_size) {
                   if(i != j) {
                     d_val <- diffs[i, j]
                     if(d_val <= range_par) {
                       h <- d_val / range_par
                       V[i, j] <- 1 - 1.5 * h + 0.5 * h^3
                     }
                   }
                 }
               }
               corr[inds, inds] <- matsqrt(V)
             }
             corr
           }
         } else {
           Vhalf_fxn <- NULL
         }

         VhalfInv_logdet <- function(par) {
           log_det <- 0
           range_par <- exp(par)
           for(clust in unique(correlation_id)) {
             inds <- which(correlation_id == clust)
             block_size <- length(inds)
             diffs <- .compute_dist_block(spacetime, inds)
             V <- diag(block_size)
             for(i in 1:block_size) {
               for(j in 1:block_size) {
                 if(i != j) {
                   d_val <- diffs[i, j]
                   if(d_val <= range_par) {
                     h <- d_val / range_par
                     V[i, j] <- 1 - 1.5 * h + 0.5 * h^3
                   }
                 }
               }
             }
             log_det <- log_det +
               (-0.5 * determinant(V, logarithm = TRUE)$modulus[1])
           }
           log_det
         }

         ## REML gradient
         #  r = exp(par),  dr/dpar = r
         #  For d <= r:  dV/dr = 1.5*d/r^2 - 1.5*d^3/r^4
         #  dV/dpar = dV/dr * r
         REML_grad <- function(par, model_fit, ...) {
           range_par <- exp(par)
           dV <- matrix(0, nrow(predictors), nrow(predictors))
           V <- matrix(0, nrow(predictors), nrow(predictors))

           for(clust in unique(correlation_id)) {
             inds <- which(correlation_id == clust)
             block_size <- length(inds)
             diffs <- .compute_dist_block(spacetime, inds)

             V_block <- diag(block_size)
             dV_block <- matrix(0, block_size, block_size)

             for(i in 1:block_size) {
               for(j in 1:block_size) {
                 if(i != j) {
                   d_val <- diffs[i, j]
                   if(d_val <= range_par) {
                     h <- d_val / range_par
                     V_block[i, j] <- 1 - 1.5 * h + 0.5 * h^3
                     ## dV/dr * dr/dpar = (1.5*d/r^2 - 1.5*d^3/r^4) * r
                     dV_block[i, j] <- (1.5 * d_val / range_par^2 -
                                          1.5 * d_val^3 / range_par^4) *
                       range_par
                   }
                 }
               }
             }
             V[inds, inds] <- V_block
             dV[inds, inds] <- dV_block
           }

           reml_grad_from_dV(dV, model_fit, glm_weight_function, ...)
         }
       }

       ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
       # Matern correlation
       #  V_ij = (2^{1-nu}/Gamma(nu)) * (sqrt(2*nu)*d/ell)^nu
       #         * K_nu(sqrt(2*nu)*d/ell)
       #  ell = exp(par[1]) > 0,  nu = exp(par[2]) > 0.
       #  No analytic gradient (Bessel function derivative w.r.t. nu).
       ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
       if(correlation_structure %in% c('matern', 'Matern')){

         VhalfInv_fxn <- function(par) {
           corr <- matrix(0, nrow(predictors), nrow(predictors))
           ell <- exp(par[1])
           nu <- exp(par[2])
           for(clust in unique(correlation_id)){
             inds <- which(correlation_id == clust)
             block_size <- length(inds)
             diffs <- .compute_dist_block(spacetime, inds)
             scaled_diffs <- sqrt(2 * nu) * diffs / ell
             V <- matrix(1, block_size, block_size)
             nonzero <- which(scaled_diffs != 0, arr.ind = TRUE)
             if(length(nonzero) > 0){
               V[nonzero] <- (2^(1 - nu) / gamma(nu)) *
                 (scaled_diffs[nonzero])^nu *
                 besselK(scaled_diffs[nonzero], nu)
             }
             corr[inds, inds] <- matinvsqrt(V)
           }
           corr
         }

         if(paste0(family)[1] != 'gaussian' |
            paste0(family)[2] != 'identity'){
           Vhalf_fxn <- function(par) {
             corr <- matrix(0, nrow(predictors), nrow(predictors))
             ell <- exp(par[1])
             nu <- exp(par[2])
             for(clust in unique(correlation_id)){
               inds <- which(correlation_id == clust)
               block_size <- length(inds)
               diffs <- .compute_dist_block(spacetime, inds)
               scaled_diffs <- sqrt(2 * nu) * diffs / ell
               V <- matrix(1, block_size, block_size)
               nonzero <- which(scaled_diffs != 0, arr.ind = TRUE)
               if(length(nonzero) > 0){
                 V[nonzero] <- (2^(1 - nu) / gamma(nu)) *
                   (scaled_diffs[nonzero])^nu *
                   besselK(scaled_diffs[nonzero], nu)
               }
               corr[inds, inds] <- matsqrt(V)
             }
             corr
           }
         } else {
           Vhalf_fxn <- NULL
         }

         VhalfInv_logdet <- function(par) {
           log_det <- 0
           ell <- exp(par[1])
           nu <- exp(par[2])
           for(clust in unique(correlation_id)){
             inds <- which(correlation_id == clust)
             block_size <- length(inds)
             diffs <- .compute_dist_block(spacetime, inds)
             scaled_diffs <- sqrt(2 * nu) * diffs / ell
             V <- matrix(1, block_size, block_size)
             nonzero <- which(scaled_diffs != 0, arr.ind = TRUE)
             if(length(nonzero) > 0){
               V[nonzero] <- (2^(1 - nu) / gamma(nu)) *
                 (scaled_diffs[nonzero])^nu *
                 besselK(scaled_diffs[nonzero], nu)
             }
             log_det <- log_det +
               (-0.5 * determinant(V, logarithm = TRUE)$modulus[1])
           }
           log_det
         }

         ## Finite-difference gradient (no closed form for d/dnu of K_nu)
         REML_grad <- NULL

         if(length(VhalfInv_par_init) == 0){
           VhalfInv_par_init <- c(0, 0)
         }
       }

       ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
       #  Gaussian-cosine correlation
       #  V_ij = exp(-d^2/(2*ell^2)) * cos(omega * d)
       #  ell = exp(par[1]) > 0,  omega = exp(par[2]) > 0.
       ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
       if(correlation_structure %in% c('gaussian-cosine',
                                       'gaussiancosine',
                                       'GaussianCosine')){

         VhalfInv_fxn <- function(par) {
           corr <- matrix(0, nrow(predictors), nrow(predictors))
           ell <- exp(par[1])
           omega <- exp(par[2])
           for(clust in unique(correlation_id)){
             inds <- which(correlation_id == clust)
             block_size <- length(inds)
             diffs <- .compute_dist_block(spacetime, inds)
             V <- matrix(1, block_size, block_size)
             nondiag <- which(diffs > 0, arr.ind = TRUE)
             if(nrow(nondiag) > 0) {
               for(i in 1:nrow(nondiag)) {
                 r <- nondiag[i, 1]; cc <- nondiag[i, 2]
                 d_val <- diffs[r, cc]
                 corr_val <- exp(-d_val^2 / (2 * ell^2)) * cos(omega * d_val)
                 V[r, cc] <- pmin(pmax(corr_val, -1), 1)
               }
             }
             corr[inds, inds] <- matinvsqrt(V)
           }
           corr
         }

         if(paste0(family)[1] != 'gaussian' |
            paste0(family)[2] != 'identity'){
           Vhalf_fxn <- function(par) {
             corr <- matrix(0, nrow(predictors), nrow(predictors))
             ell <- exp(par[1])
             omega <- exp(par[2])
             for(clust in unique(correlation_id)){
               inds <- which(correlation_id == clust)
               block_size <- length(inds)
               diffs <- .compute_dist_block(spacetime, inds)
               V <- matrix(1, block_size, block_size)
               nondiag <- which(diffs > 0, arr.ind = TRUE)
               if(nrow(nondiag) > 0) {
                 for(i in 1:nrow(nondiag)) {
                   r <- nondiag[i, 1]; cc <- nondiag[i, 2]
                   d_val <- diffs[r, cc]
                   corr_val <- exp(-d_val^2 / (2 * ell^2)) * cos(omega * d_val)
                   V[r, cc] <- pmin(pmax(corr_val, -1), 1)
                 }
               }
               corr[inds, inds] <- matsqrt(V)
             }
             corr
           }
         } else {
           Vhalf_fxn <- NULL
         }

         VhalfInv_logdet <- function(par) {
           log_det <- 0
           ell <- exp(par[1])
           omega <- exp(par[2])
           for(clust in unique(correlation_id)){
             inds <- which(correlation_id == clust)
             block_size <- length(inds)
             diffs <- .compute_dist_block(spacetime, inds)
             V <- matrix(1, block_size, block_size)
             nondiag <- which(diffs > 0, arr.ind = TRUE)
             if(nrow(nondiag) > 0) {
               for(i in 1:nrow(nondiag)) {
                 r <- nondiag[i, 1]; cc <- nondiag[i, 2]
                 d_val <- diffs[r, cc]
                 corr_val <- exp(-d_val^2 / (2 * ell^2)) * cos(omega * d_val)
                 V[r, cc] <- pmin(pmax(corr_val, -1), 1)
               }
             }
             log_det <- log_det +
               (-0.5 * determinant(V, logarithm = TRUE)$modulus[1])
           }
           log_det
         }

         ## REML gradient (two parameters)
         #  par[1]: log(ell).  dV/dpar1 = (d^2/ell^2) * gauss * cos * dell/dpar = same * ell
         #  par[2]: log(omega). dV/dpar2 = gauss * (-d*sin(omega*d)) * domega/dpar = same * omega
         REML_grad <- function(par, model_fit, ...) {
           ell <- exp(par[1])
           omega <- exp(par[2])

           dV1 <- matrix(0, nrow(predictors), nrow(predictors))
           dV2 <- matrix(0, nrow(predictors), nrow(predictors))
           V <- matrix(0, nrow(predictors), nrow(predictors))

           for(clust in unique(correlation_id)){
             inds <- which(correlation_id == clust)
             block_size <- length(inds)
             diffs <- .compute_dist_block(spacetime, inds)

             V_block <- matrix(1, block_size, block_size)
             dV1_block <- matrix(0, block_size, block_size)
             dV2_block <- matrix(0, block_size, block_size)

             nondiag <- which(diffs > 0, arr.ind = TRUE)
             if(nrow(nondiag) > 0) {
               for(i in 1:nrow(nondiag)) {
                 r <- nondiag[i, 1]; cc <- nondiag[i, 2]
                 d_val <- diffs[r, cc]
                 gauss <- exp(-d_val^2 / (2 * ell^2))
                 cos_part <- cos(omega * d_val)
                 corr_val <- pmin(pmax(gauss * cos_part, -1), 1)
                 V_block[r, cc] <- corr_val

                 ## d/dpar1: chain rule through ell = exp(par[1])
                 dV1_block[r, cc] <- (d_val^2 / ell^3) * gauss * cos_part * ell
                 ## d/dpar2: chain rule through omega = exp(par[2])
                 dV2_block[r, cc] <- gauss * (-d_val * sin(omega * d_val)) * omega
               }
             }
             V[inds, inds] <- V_block
             dV1[inds, inds] <- dV1_block
             dV2[inds, inds] <- dV2_block
           }

           gradient <- numeric(2)
           gradient[1] <- reml_grad_from_dV(dV1, model_fit,
                                             glm_weight_function, ...)
           gradient[2] <- reml_grad_from_dV(dV2, model_fit,
                                             glm_weight_function, ...)
           return(gradient)
         }

         if(length(VhalfInv_par_init) == 0){
           VhalfInv_par_init <- c(0, 0)
         }
       }

       ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
       # Gamma-cosine correlation
       #  V_ij = [d^{a-1} exp(-b*d) / (Gamma(a)/b^a)] * cos(omega*d)
       #  a = exp(par[1]),  b = exp(par[2]),  omega = exp(par[3]).
       #  All > 0.
       ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
       if(correlation_structure %in% c('gamma-cosine',
                                       'gammacosine',
                                       'GammaCosine')){

         VhalfInv_fxn <- function(par) {
           corr <- matrix(0, nrow(predictors), nrow(predictors))
           shape <- exp(par[1])
           rate <- exp(par[2])
           omega <- exp(par[3])
           norm_const <- gamma(shape) / rate^shape
           for(clust in unique(correlation_id)){
             inds <- which(correlation_id == clust)
             block_size <- length(inds)
             diffs <- .compute_dist_block(spacetime, inds)
             V <- matrix(1, block_size, block_size)
             nondiag <- which(diffs > 0, arr.ind = TRUE)
             if(nrow(nondiag) > 0) {
               for(i in 1:nrow(nondiag)) {
                 r <- nondiag[i, 1]; cc <- nondiag[i, 2]
                 d_val <- diffs[r, cc]
                 gamma_val <- (d_val^(shape - 1) * exp(-rate * d_val)) /
                   norm_const
                 corr_val <- gamma_val * cos(omega * d_val)
                 V[r, cc] <- pmin(pmax(corr_val, -1), 1)
               }
             }
             corr[inds, inds] <- matinvsqrt(V)
           }
           corr
         }

         if(paste0(family)[1] != 'gaussian' |
            paste0(family)[2] != 'identity'){
           Vhalf_fxn <- function(par) {
             corr <- matrix(0, nrow(predictors), nrow(predictors))
             shape <- exp(par[1])
             rate <- exp(par[2])
             omega <- exp(par[3])
             norm_const <- gamma(shape) / rate^shape
             for(clust in unique(correlation_id)){
               inds <- which(correlation_id == clust)
               block_size <- length(inds)
               diffs <- .compute_dist_block(spacetime, inds)
               V <- matrix(1, block_size, block_size)
               nondiag <- which(diffs > 0, arr.ind = TRUE)
               if(nrow(nondiag) > 0) {
                 for(i in 1:nrow(nondiag)) {
                   r <- nondiag[i, 1]; cc <- nondiag[i, 2]
                   d_val <- diffs[r, cc]
                   gamma_val <- (d_val^(shape - 1) * exp(-rate * d_val)) /
                     norm_const
                   corr_val <- gamma_val * cos(omega * d_val)
                   V[r, cc] <- pmin(pmax(corr_val, -1), 1)
                 }
               }
               corr[inds, inds] <- matsqrt(V)
             }
             corr
           }
         } else {
           Vhalf_fxn <- NULL
         }

         VhalfInv_logdet <- function(par) {
           log_det <- 0
           shape <- exp(par[1])
           rate <- exp(par[2])
           omega <- exp(par[3])
           norm_const <- gamma(shape) / rate^shape
           for(clust in unique(correlation_id)){
             inds <- which(correlation_id == clust)
             block_size <- length(inds)
             diffs <- .compute_dist_block(spacetime, inds)
             V <- matrix(1, block_size, block_size)
             nondiag <- which(diffs > 0, arr.ind = TRUE)
             if(nrow(nondiag) > 0) {
               for(i in 1:nrow(nondiag)) {
                 r <- nondiag[i, 1]; cc <- nondiag[i, 2]
                 d_val <- diffs[r, cc]
                 gamma_val <- (d_val^(shape - 1) * exp(-rate * d_val)) /
                   norm_const
                 corr_val <- gamma_val * cos(omega * d_val)
                 V[r, cc] <- pmin(pmax(corr_val, -1), 1)
               }
             }
             log_det <- log_det +
               (-0.5 * determinant(V, logarithm = TRUE)$modulus[1])
           }
           log_det
         }

         ## REML gradient (three parameters, all on log scale)
         REML_grad <- function(par, model_fit, ...) {
           shape <- exp(par[1])
           rate <- exp(par[2])
           omega <- exp(par[3])
           norm_const <- gamma(shape) / rate^shape
           dnorm_dshape <- norm_const * (digamma(shape) - log(rate))
           dnorm_drate <- -shape * norm_const / rate

           dV1 <- matrix(0, nrow(predictors), nrow(predictors))
           dV2 <- matrix(0, nrow(predictors), nrow(predictors))
           dV3 <- matrix(0, nrow(predictors), nrow(predictors))
           V <- matrix(0, nrow(predictors), nrow(predictors))

           for(clust in unique(correlation_id)){
             inds <- which(correlation_id == clust)
             block_size <- length(inds)
             diffs <- .compute_dist_block(spacetime, inds)

             V_block <- matrix(1, block_size, block_size)
             dV1_block <- matrix(0, block_size, block_size)
             dV2_block <- matrix(0, block_size, block_size)
             dV3_block <- matrix(0, block_size, block_size)

             nondiag <- which(diffs > 0, arr.ind = TRUE)
             if(nrow(nondiag) > 0) {
               for(i in 1:nrow(nondiag)) {
                 r <- nondiag[i, 1]; cc <- nondiag[i, 2]
                 d_val <- diffs[r, cc]
                 gamma_part <- (d_val^(shape - 1) * exp(-rate * d_val)) /
                   norm_const
                 cos_part <- cos(omega * d_val)
                 corr_val <- pmin(pmax(gamma_part * cos_part, -1), 1)
                 V_block[r, cc] <- corr_val

                 ## d/dshape of gamma_part, times chain rule shape = exp(par[1])
                 dgamma_shape <- gamma_part *
                   (log(d_val) - dnorm_dshape / norm_const)
                 dV1_block[r, cc] <- shape * dgamma_shape * cos_part

                 ## d/drate of gamma_part, times chain rule rate = exp(par[2])
                 dgamma_rate <- gamma_part *
                   (-d_val - dnorm_drate / norm_const)
                 dV2_block[r, cc] <- rate * dgamma_rate * cos_part

                 ## d/domega of cos_part, times chain rule omega = exp(par[3])
                 dV3_block[r, cc] <- gamma_part *
                   (-d_val * sin(omega * d_val)) * omega
               }
             }
             V[inds, inds] <- V_block
             dV1[inds, inds] <- dV1_block
             dV2[inds, inds] <- dV2_block
             dV3[inds, inds] <- dV3_block
           }

           gradient <- numeric(3)
           gradient[1] <- reml_grad_from_dV(dV1, model_fit,
                                             glm_weight_function, ...)
           gradient[2] <- reml_grad_from_dV(dV2, model_fit,
                                             glm_weight_function, ...)
           gradient[3] <- reml_grad_from_dV(dV3, model_fit,
                                             glm_weight_function, ...)
           return(gradient)
         }

         if(length(VhalfInv_par_init) == 0){
           VhalfInv_par_init <- c(0, 0, 0)
         }
       }
     }
    ## Default starting value
    if(length(VhalfInv_par_init) == 0){
      VhalfInv_par_init <- 0
    }
  }

  ## If fitting correlation-structure components, do so here:
  if(!is.null(VhalfInv_fxn) & length(VhalfInv_par_init) > 0){
      ## Only use custom BFGS if gradient supplied, otherwise stats::optim
      # has a more robust version for finite-difference approximation
      if(is.null(REML_grad)){
        efficient_bfgs <- function(par_in, fn){
          res <- stats::optim(
            par_in,
            function(par)fn(c(par))[[1]],
            method = 'BFGS',
            hessian = TRUE
          )
          res$vcov <- invert(methods::as(res$hessian, 'matrix'))
          return(res)
        }
      }
      if(length(observation_weights) == 0){
        observation_weights <- rep(1, length(y))
      }

      ## Constant used in computing log-likelihood of normal distributions
      neghalf_l2pi <- -0.5*log(2*pi)

      ## Optimize correlation-structure parameters, via VhalfInv
      res <-
        efficient_bfgs(
        c(VhalfInv_par_init),
        fn = function(par){
          tr <- try({
            VhalfInv <- VhalfInv_fxn(par)
            if(!is.null(Vhalf_fxn)){
              Vhalf <- Vhalf_fxn(par)
            } else {
              Vhalf <- invert(VhalfInv)
            }
          })
          if(inherits(tr, 'try-error')){
            return(list(NaN, NaN))
          }
          tr <- try({nrow(VhalfInv) == length(y) &
                     ncol(VhalfInv) == nrow(VhalfInv)}, silent = TRUE)
          if(inherits(tr, 'try-error')){
            stop('\n \t VhalfInv_fxn does not return a matrix.',
                 'Adjust your function. \n')
          } else if(!tr){
            stop('\n \t VhalfInv_fxn does not return an N x N matrix. ',
                 'Adjust your function. \n')
          }
          ## Re-fit
          model_fit <- try({lgspline.fit(predictors,
                                    y,
                                    standardize_response,
                                    standardize_predictors_for_knots,
                                    standardize_expansions_for_fitting,
                                    family,
                                    glm_weight_function,
                                    schur_correction_function,
                                    need_dispersion_for_estimation,
                                    dispersion_function,
                                    K,
                                    model_fit$knots,
                                    cluster_on_indicators,
                                    model_fit$make_partition_list,
                                    previously_tuned_penalties,
                                    smoothing_spline_penalty,
                                    opt,
                                    use_custom_bfgs,
                                    delta,
                                    tol,
                                    initial_wiggle,
                                    initial_flat,
                                    wiggle_penalty,
                                    flat_ridge_penalty,
                                    unique_penalty_per_partition,
                                    unique_penalty_per_predictor,
                                    meta_penalty,
                                    predictor_penalties,
                                    partition_penalties,
                                    include_quadratic_terms,
                                    include_cubic_terms,
                                    include_quartic_terms,
                                    include_2way_interactions,
                                    include_3way_interactions,
                                    include_quadratic_interactions,
                                    offset,
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
                                    blockfit,
                                    qp_score_function,
                                    qp_observations,
                                    qp_Amat,
                                    qp_bvec,
                                    qp_meq,
                                    qp_positive_derivative,
                                    qp_negative_derivative,
                                    qp_positive_2ndderivative,
                                    qp_negative_2ndderivative,
                                    qp_monotonic_increase,
                                    qp_monotonic_decrease,
                                    qp_range_upper,
                                    qp_range_lower,
                                    qp_Amat_fxn,
                                    qp_bvec_fxn,
                                    qp_meq_fxn,
                                    constraint_values,
                                    constraint_vectors,
                                    return_G,
                                    return_Ghalf,
                                    TRUE,#return_U,
                                    TRUE,#estimate dispersion
                                    unbias_dispersion,#unbias dispersion,
                                    TRUE,#return_varcovmat,
                                    FALSE,#return_lagrange_multipliers
                                    custom_penalty_mat,
                                    cluster_args,
                                    dummy_dividor,
                                    dummy_adder,
                                    verbose,
                                    verbose_tune,
                                    dummy_fit,
                                    auto_encode_factors,
                                    observation_weights,
                                    do_not_cluster_on_these,
                                    neighbor_tolerance,
                                    no_intercept,
                                    VhalfInv,
                                    Vhalf,
                                    include_warnings,
                                    og_cols,
                                    ...)}, silent = TRUE)
          if(any(inherits(model_fit, 'try-error'))){
            return(list(NaN, NaN))
          }

          ## Use custom loss if available in family list
          if(!is.null(custom_VhalfInv_loss)){
            return(list(custom_VhalfInv_loss(par,
                                             model_fit,
                                             ...),
                        REML_grad(par,
                        model_fit,
                        ...)))

          ## Otherwise use REML and GLS for non-Gaussian Response
          } else {
            if(verbose) cat('\nREML computation\n')
            VhalfInv <- t(t(VhalfInv) * sqrt(observation_weights))

            ## Crude and conservative estimate if degenerate estimate of sigma^2
            if(model_fit$sigmasq_tilde <= 0){
              degfed <- model_fit$N - model_fit$P + qr(model_fit$A)$rank
              model_fit$sigmasq_tilde <-
                sum((model_fit$y - model_fit$ytilde)^2/degfed)
            }

            ## REML log-likelihood
            # Weighted SSE
            W <- c(glm_weight_function(model_fit$ytilde,
                                       model_fit$y,
                                       unlist(model_fit$order_list),
                                       model_fit$family,
                                       model_fit$sigmasq_tilde,
                                       rep(1, model_fit$N),
                                       ...))
            if(!is.null(model_fit$family$custom_dev.resids)){
              raw <- family$custom_dev.resids(model_fit$y,
                                              model_fit$ytilde,
                                              1:model_fit$N,
                                              model_fit$family,
                                              1+0*model_fit$weights,
                                              ...)
              ## [Change 2026-02-18] Weight residuals by 1/sqrt(W)
              logloss <-  sum((
                VhalfInv %**% cbind(sign(raw)*sqrt(abs(
                  raw
                )) / sqrt(c(W)))
              )^2)
            }
            else if(is.null(model_fit$family$dev.resids)){
              ## [Change 2026-02-18] weight by 1/sqrt(W)
              logloss <-
                sum(c(VhalfInv %**% cbind(
                  sqrt(1/c(W)) * (model_fit$y - model_fit$ytilde)
                ))^2)
            } else {
              ## [Change 2026-02-18] pre-whiten W^{1/2}-weighted with sqrt W
              logloss <- sum(model_fit$family$dev.resids(
                c(VhalfInv %**% cbind(1/sqrt(c(W)) * model_fit$y)),
                c(VhalfInv %**% cbind(1/sqrt(c(W)) * model_fit$ytilde)),
                wt = rep(1, model_fit$N)))
            }
            # -log| V / sigma^2 |
            if(!is.null(VhalfInv_logdet)){
              logdet_VhalfInv <- VhalfInv_logdet(par)
            } else {
              logdet_VhalfInv <- determinant(VhalfInv,
                                             logarithm=TRUE)$modulus[1]
            }

            ## Generalized determinant of inverse information matrix
            eigvals <- eigen(model_fit$varcovmat,
                             symmetric = TRUE,
                             only.values = TRUE)$values
            nonzero_eigvals <- eigvals[eigvals > sqrt(.Machine$double.eps)]
            logdet_varcovB <- -sum(log(nonzero_eigvals))

            ## Full negative-REML, generalized to include
            #  penalties,
            #  constraints,
            #  and glm link functions
            reml_objective <- (
                -logdet_VhalfInv +
                0.5 * model_fit$N * log(model_fit$sigmasq_tilde) +
                0.5 * logloss / model_fit$sigmasq_tilde +
                0.5 * logdet_varcovB
              ) / model_fit$N

            ## Now compute gradient, if function is supplied
            if(!is.null(REML_grad)){
              reml_gradient <- REML_grad(par,
                                         model_fit,
                                         ...)
            } else {
              reml_gradient <- NULL
            }

            if(verbose) cat('\nDone REML computation\n')
            return(list(reml_objective = reml_objective,
                        reml_gradient = reml_gradient))
          }
      }
    )
    if(verbose){
      cat("Done VhalfInv Optimization\n")
    }

    ## Extract variance-covariance matrix of estimates using BFGS approximation
    VhalfInv_params_estimates <- res$par
    VhalfInv <- VhalfInv_fxn(c(VhalfInv_params_estimates))
    if(!is.null(Vhalf_fxn)){
      Vhalf <- Vhalf_fxn(VhalfInv_params_estimates)
    }
    VhalfInv_params_vcov <- res$vcov
    res <- NULL

    ## Re-fit given optimal values for VhalfInv now one final time
    model_fit <- lgspline.fit(predictors,
                              y,
                              standardize_response,
                              standardize_predictors_for_knots,
                              standardize_expansions_for_fitting,
                              family,
                              glm_weight_function,
                              schur_correction_function,
                              need_dispersion_for_estimation,
                              dispersion_function,
                              K,
                              model_fit$knots,
                              cluster_on_indicators,
                              model_fit$make_partition_list,
                              previously_tuned_penalties,
                              smoothing_spline_penalty,
                              opt,
                              use_custom_bfgs,
                              delta,
                              tol,
                              initial_wiggle,
                              initial_flat,
                              wiggle_penalty,
                              flat_ridge_penalty,
                              unique_penalty_per_partition,
                              unique_penalty_per_predictor,
                              meta_penalty,
                              predictor_penalties,
                              partition_penalties,
                              include_quadratic_terms,
                              include_cubic_terms,
                              include_quartic_terms,
                              include_2way_interactions,
                              include_3way_interactions,
                              include_quadratic_interactions,
                              offset,
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
                              blockfit,
                              qp_score_function,
                              qp_observations,
                              qp_Amat,
                              qp_bvec,
                              qp_meq,
                              qp_positive_derivative,
                              qp_negative_derivative,
                              qp_positive_2ndderivative,
                              qp_negative_2ndderivative,
                              qp_monotonic_increase,
                              qp_monotonic_decrease,
                              qp_range_upper,
                              qp_range_lower,
                              qp_Amat_fxn,
                              qp_bvec_fxn,
                              qp_meq_fxn,
                              constraint_values,
                              constraint_vectors,
                              return_G,
                              return_Ghalf,
                              return_U,
                              estimate_dispersion,
                              unbias_dispersion,
                              return_varcovmat,
                              return_lagrange_multipliers,
                              custom_penalty_mat,
                              cluster_args,
                              dummy_dividor,
                              dummy_adder,
                              verbose,
                              verbose_tune,
                              dummy_fit,
                              auto_encode_factors,
                              observation_weights,
                              do_not_cluster_on_these,
                              neighbor_tolerance,
                              no_intercept,
                              VhalfInv,
                              Vhalf,
                              include_warnings,
                              og_cols,
                              ...)
    model_fit$VhalfInv_params_vcov <- abs(VhalfInv_params_vcov)
    model_fit$VhalfInv_params_estimates <- VhalfInv_params_estimates
    model_fit$VhalfInv_fxn <- VhalfInv_fxn
    model_fit$Vhalf_fxn <- Vhalf_fxn
    model_fit$VhalfInv_logdet <- VhalfInv_logdet
    model_fit$REML_grad <- REML_grad
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

  ## Important information thus far missing
  model_fit$critical_value <- critical_value

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
                                           enforce_constraints = FALSE,
                                           max_rejection_draws = 50L,
                                           ...){

    ## Check compatibility, that new_predictors should be a matrix
    if(any(!is.null(new_predictors))){
      new_predictors <- try(methods::as(cbind(new_predictors), 'matrix'), silent = TRUE)
      if(any(inherits(new_predictors, 'try-error'))){
        stop('\n\t new_predictors must be coercible to a matrix.\n')
      }
    }

    ## One approach to avoiding R issues while preserving the rest of
    # lgspline's functionality
    only_1 <- FALSE
    if(nrow(new_predictors) == 1){
      only_1 <- TRUE
      new_predictors <- rbind(new_predictors, new_predictors)
    }

    ## Helpful components
    nc <- model_fit$p
    K  <- model_fit$K
    nr <- model_fit$N

    ## Helper: check whether a single draw satisfies all QP constraints
    #  qp_Amat has dimensions P x n_constraints (solve.QP convention),
    #  so the feasibility check is t(Amat) %*% beta >= bvec.
    has_qp <- !is.null(model_fit$quadprog_list)
    .draw_feasible <- function(beta_vec){
      if(!has_qp || !enforce_constraints) return(TRUE)
      Amat <- model_fit$quadprog_list$qp_Amat
      bvec <- model_fit$quadprog_list$qp_bvec
      meq  <- model_fit$quadprog_list$qp_meq
      lhs  <- c(t(Amat) %**% cbind(beta_vec))
      ## Equality constraints: lhs == bvec (within tolerance)
      eq_ok  <- all(abs(lhs[seq_len(meq)] - bvec[seq_len(meq)]),
                    sqrt(.Machine$double.eps))
      ## Inequality constraints: lhs >= bvec
      ineq_ok <- all(lhs[-seq_len(meq)] >= bvec[-seq_len(meq)] -
                       sqrt(.Machine$double.eps))
      eq_ok & ineq_ok
    }

    ## Single-draw sampler (may loop for accept/reject)
    .one_draw <- function(){

      ## Draw a dispersion parameter, if applicable, from InvG distribution
      if(draw_dispersion){
        shape <- theta_1 +
          0.5*(model_fit$N - model_fit$unbias_dispersion *
                 model_fit$trace_XUGX)
        rate <- theta_2 +
          0.5*(model_fit$N - model_fit$unbias_dispersion *
                 model_fit$trace_XUGX) *
          new_sigmasq_tilde
        if(shape <= 0){
          stop('\n \t Posterior inverse-gamma shape is <= 0, increase ',
               'theta_1 argument to draw a dispersion parameter. \n')
        }
        if(rate <= 0){
          stop('\n \t Posterior inverse-gamma rate is <= 0, increase theta_2 ',
               'argument to draw a dispersion parameter. \n')
        }
        post_draw_sigmasq <- 1/rgamma(1, shape, rate)

        if((is.nan(post_draw_sigmasq) | !is.finite(post_draw_sigmasq)) &
           include_warnings){
          warning("\n\t Infinite/NaN posterior draw of dispersion detected. \n")
          post_draw_sigmasq <- new_sigmasq_tilde
        }
      } else {
        post_draw_sigmasq <- new_sigmasq_tilde
      }

      ## Accept/reject loop (only entered when enforce_constraints = TRUE
      #  and QP constraints are present; otherwise exits after one draw)
      n_reject <- 0L
      repeat {

        ## Draw posterior coefficient errors
        if(!is.null(model_fit$VhalfInv)){

          ## Full N x P unwhitened design (original partition order)
          X_full_post <- collapse_block_diagonal(
            lapply(model_fit$X, model_fit$std)
          )

          W_glm_post <- c(glm_weight_function(
            model_fit$ytilde[unlist(model_fit$order_list)],
            model_fit$y[unlist(model_fit$order_list)],
            1:model_fit$N,
            model_fit$family,
            new_sigmasq_tilde,
            rep(1, model_fit$N), ...)
          )
          W_glm_post <- pmax(W_glm_post, .Machine$double.eps)

          D_post <- model_fit$weights
          if(is.list(D_post)){
            D_post <- unlist(D_post)[unlist(model_fit$order_list)]
          } else {
            D_post <- D_post[unlist(model_fit$order_list)]
          }

          X_full_post <- X_full_post * sqrt(W_glm_post * D_post)

          VinvhalfX_post <- model_fit$VhalfInv[unlist(model_fit$order_list),
                                               unlist(model_fit$order_list)] %**%
            X_full_post

          has_part_pen <-
            length(model_fit$penalties$L_partition_list) == (model_fit$K + 1)
          Lambda_full_post <- collapse_block_diagonal(
            lapply(1:(model_fit$K + 1), function(k){
              if(has_part_pen){
                model_fit$penalties$Lambda +
                  model_fit$penalties$L_partition_list[[k]]
              } else {
                model_fit$penalties$Lambda
              }
            })
          )

          Ghalf_correct_post <- matinvsqrt(gramMatrix(VinvhalfX_post) +
                                             Lambda_full_post)

          z_full <- cbind(rnorm((model_fit$K + 1) * nc))
          post_draw_coefficients_err <-
            (1 / model_fit$sd_y) *
            sqrt(post_draw_sigmasq) *
            (model_fit$U %**% (Ghalf_correct_post %**% z_full))

        } else {

          post_draw_coefficients_err <-
            (1 / model_fit$sd_y) *
            sqrt(post_draw_sigmasq) *
            (model_fit$U %**%
               cbind(Reduce("c", lapply(1:(K + 1), function(k){
                 c(model_fit$Ghalf[[k]] %**%
                     cbind(rnorm(nc)))
               }))))
        }

        ## Assemble raw draw and check feasibility
        post_draw_coefficients <- lapply(1:(K + 1), function(k){
          raw_draw <- post_draw_coefficients_err[1:nc + (k-1)*nc] +
            model_fit$B_raw[[k]]
          raw_draw <- raw_draw * model_fit$sd_y
          raw_draw[1] <- raw_draw[1] + model_fit$mean_y
          return(model_fit$backtransform_coefficients(raw_draw))
        })

        ## Feasibility check on raw (unstandardized) coefficients
        beta_vec_check <- unlist(post_draw_coefficients)
        if(.draw_feasible(beta_vec_check)) break

        n_reject <- n_reject + 1L
        if(n_reject >= max_rejection_draws){
          if(include_warnings){
            warning(
              "\n\t Accept/reject sampling failed to find a feasible draw ",
              "after ", max_rejection_draws, " attempts. ",
              "Returning the MAP estimate for this draw. ",
              "Consider using enforce_constraints = FALSE (the default) ",
              "and treating draws as approximate.\n"
            )
          }
          ## Fall back to MAP
          post_draw_coefficients <- model_fit$B
          break
        }
      }

      ## Return posterior predictive draws
      if(include_posterior_predictive){
        post_pred_mean <- model_fit$predict(
          new_predictors,
          B_predict = post_draw_coefficients)

        post_pred_draw <- posterior_predictive_draw(length(post_pred_mean),
                                                    post_pred_mean,
                                                    sqrt(post_draw_sigmasq),
                                                    ...)

        return(list(post_pred_draw = post_pred_draw,
                    post_draw_coefficients = post_draw_coefficients,
                    post_draw_sigmasq = post_draw_sigmasq))
      } else {
        return(list(
          post_draw_coefficients = post_draw_coefficients,
          post_draw_sigmasq = post_draw_sigmasq
        ))
      }
    }

    res <- lapply(1:num_draws, function(m) .one_draw())

    if(num_draws == 1){
      if(only_1){
        res[[1]][[1]] <- res[[1]][[1]][1]
      }
      return(res[[1]])
    }

    ## Combine results
    post_draw_coefficients <- lapply(res, `[[`, "post_draw_coefficients")
    post_draw_sigmasq      <- lapply(res, `[[`, "post_draw_sigmasq")
    if(include_posterior_predictive){
      post_pred_draw <- Reduce("cbind", lapply(res, `[[`, "post_pred_draw"))

      if(only_1){
        post_pred_draw <- post_pred_draw[1, , drop = FALSE]
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
  ## Find global maximum/minimum
  model_fit$find_extremum <- function(
    vars = NULL,
    quick_heuristic = TRUE,
    initial = NULL,
    B_predict = NULL,
    minimize = FALSE,
    stochastic = FALSE,
    stochastic_draw = function(mu, sigma, ...){
      N <- length(mu)
      rnorm(N, mu, sigma)
    },
    sigmasq_predict = model_fit$sigmasq_tilde,
    custom_objective_function = NULL,
    custom_objective_derivative = NULL,
    ...
  ){
    ## Square-root dispersion is a more convenient parameterization in practice
    sigma_tilde <- sqrt(sigmasq_predict)

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
    if(quick_heuristic | !is.null(initial)){
      best_fitted <- which.max(model_fit$ytilde * (-min_or_max))
      partitions <- which(sapply(1:(model_fit$K+1),function(k){
        best_fitted %in% model_fit$order_list[[k]]
      }))
    }

    ## Go through each partition, optimize the cubic function within
    # remove empty partitions first
    partitions_keep <- c(c(), which(sapply(partitions, function(k){
      nrow(model_fit$X[[k]])
    }) > 0))

    ## Find variables to optimize over
    if(inherits(vars,  'numeric')){
      nms <- paste0('_', vars, '_')
      beta_inds <- which(sapply(model_fit$raw_expansion_names,
                                function(nm)any(grepl(nm, nms))))
      select_vars_fl <- TRUE
    } else if(inherits(vars,'character')){
      if(length(og_cols) == 0){
        stop('\n\t Do not submit character argument to "vars" unless you have',
             ' named columns in the predictors you used to fit the model ',
             ' and the "data" argument was not NULL \n')
      }
      vars <- unlist(sapply(vars, function(v)which(og_cols == v)))
      nms <- paste0('_', vars, '_')
      beta_inds <- which(sapply(model_fit$raw_expansion_names,
                                function(nm)any(grepl(nm, nms))))
      select_vars_fl <- TRUE
    } else {
      beta_inds <- 1:model_fit$p
      select_vars_fl <- FALSE
      vars <- 1:model_fit$q
    }

    ## Loop through partitions (or only the "best" one)
    best_per_partition <- lapply(partitions[partitions_keep], function(k){

      if(any(!is.null(initial))){
        ## [Change 2026-02-17] Coerce initial to a named numeric 1-row matrix
        #  regardless of input type (data.frame, named vector, plain vector).
        #  When initial has names (e.g., from a data.frame), reorder columns
        #  to match the predictor matrix column order so that Time and Dose
        #  don't get swapped.
        init_vec <- as.numeric(unlist(initial))
        init_names <- names(unlist(initial))
        pred_names <- if(!is.null(og_cols)) og_cols else colnames(predictors)
        if(!is.null(init_names) && !is.null(pred_names)){
          reorder <- match(pred_names, init_names)
          if(!any(is.na(reorder))){
            init_vec <- init_vec[reorder]
          }
        }
        predictors_vals <- rbind(init_vec)
        if(ncol(predictors_vals) != ncol(predictors)){
          stop('\n\t initial must have length equal to number of predictors (',
               ncol(predictors), ')\n')
        }
      } else {
        ## Extract best fitted value for initialization
        yk <- model_fit$X[[k]] %**% B_predict[[k]]
        best <- which.max(-yk*min_or_max)
        predictors_vals <- predictors[model_fit$order_list[[k]][best],
                                      , drop=FALSE]
      }

      ## [Change 2026-02-17] Ensure predictors_vals is always a 1-row numeric
      #  matrix with colnames matching the predictor matrix.
      if(!is.matrix(predictors_vals)){
        predictors_vals <- rbind(as.numeric(predictors_vals))
      }
      if(is.null(colnames(predictors_vals))){
        if(!is.null(og_cols)){
          colnames(predictors_vals) <- og_cols
        } else if(!is.null(colnames(predictors))){
          colnames(predictors_vals) <- colnames(predictors)
        }
      }

      ## [Change 2026-02-17] Precompute positional column indices for vars
      #  within predictors_vals. After the earlier resolution block in
      #  find_extremum, `vars` may be a named integer (e.g., c(Time=1) from
      #  character input) or a plain integer vector (from numeric input).
      #  We need the positional index for column subsetting and assignment.
      if(select_vars_fl){
        if(is.character(vars)){
          vars_idx <- match(vars, colnames(predictors_vals))
        } else {
          ## vars is numeric column indices (possibly named);
          ## use the values directly as positional indices
          vars_idx <- as.integer(vars)
        }
        ## Safety check
        if(any(is.na(vars_idx))){
          stop('\n\t Could not resolve vars to column indices of predictors. ',
               'Check that vars matches predictor column names or indices.\n')
        }
      } else {
        vars_idx <- seq_len(ncol(predictors_vals))
      }

      ## [Change 2026-02-17] Compute bounds, subset to vars when optimizing
      #  a subset. Use positional indices for subsetting.
      pred_lower <- apply(predictors, 2, min)
      pred_upper <- apply(predictors, 2, max)
      if(select_vars_fl){
        optim_lower <- pred_lower[vars_idx]
        optim_upper <- pred_upper[vars_idx]
      } else {
        optim_lower <- pred_lower
        optim_upper <- pred_upper
      }

      ## [Change 2026-02-17] Extract starting values using positional indices.
      start_vals <- as.numeric(predictors_vals[, vars_idx, drop = TRUE])

      ## [Change 2026-02-17] Helper to extract first derivatives as a numeric
      #  vector from predict() output. For multi-predictor models, first_deriv
      #  is a named list of per-variable derivative vectors; for single-predictor
      #  models, it is already a numeric vector/scalar.
      .extract_first_deriv <- function(deriv_result){
        fd <- deriv_result$first_deriv
        if(is.list(fd)){
          as.numeric(unlist(fd))
        } else {
          as.numeric(fd)
        }
      }

      ## Quasi-newton optimization
      opt <- stats::optim(
        start_vals,
        fn = function(par){
          ## [Change 2026-02-17] Reconstruct full predictor vector using
          #  positional vars_idx, always as plain numeric
          if(select_vars_fl){
            dummy <- as.numeric(predictors_vals)
            dummy[vars_idx] <- as.numeric(par)
            par <- dummy
          }
          par <- as.numeric(par)
          if(!is.null(custom_objective_function)){
            pred <- model_fit$predict(new_predictors = rbind(par),
                                      parallel = FALSE,
                                      cl = NULL,
                                      chunk_size = NULL,
                                      num_chunks = NULL,
                                      rem_chunks = NULL,
                                      B_predict = B_predict)
            if(stochastic){
              pred <- stochastic_draw(pred, sigma_tilde, ...)
            }
            min_or_max*custom_objective_function(pred,
                                                 sigma_tilde,
                                                 max(-y*min_or_max),
                                                 ...)
          } else {
            pred <- model_fit$predict(new_predictors = rbind(par),
                                      parallel = FALSE,
                                      cl = NULL,
                                      chunk_size = NULL,
                                      num_chunks = NULL,
                                      rem_chunks = NULL,
                                      B_predict = B_predict)
            if(stochastic){
              pred <- stochastic_draw(pred, sigma_tilde, ...)
            }
            min_or_max*pred
          }
        },
        gr = function(par){
          if(select_vars_fl){
            dummy <- as.numeric(predictors_vals)
            dummy[vars_idx] <- as.numeric(par)
            par <- dummy
          }
          par <- as.numeric(par)
          if(!is.null(custom_objective_derivative)) {
            pred <- model_fit$predict(new_predictors = rbind(par),
                                      parallel = FALSE,
                                      cl = NULL,
                                      chunk_size = NULL,
                                      num_chunks = NULL,
                                      rem_chunks = NULL,
                                      B_predict = B_predict)
            deriv_result <- model_fit$predict(
              new_predictors = rbind(par),
              parallel = FALSE,
              cl = NULL,
              chunk_size = NULL,
              num_chunks = NULL,
              rem_chunks = NULL,
              B_predict = B_predict,
              take_first_derivatives = TRUE)
            gr <- .extract_first_deriv(deriv_result)
            gr_par <- rep(0, length(par))
            gr_raw <- min_or_max*custom_objective_derivative(pred,
                                                             sigma_tilde,
                                                             max(-y*min_or_max),
                                                             gr,
                                                             ...)
            gr_raw <- as.numeric(unlist(gr_raw))
            ## [Change 2026-02-22] Assign gradients only to numeric predictor indices
            #  that correspond to variables being optimized
            numeric_vars_idx <- intersect(model_fit$numerics, vars_idx)
            gr_par[numeric_vars_idx] <- gr_raw[match(numeric_vars_idx, vars_idx)]
            if(select_vars_fl){
              gr_par[vars_idx]
            } else {
              gr_par
            }
          } else {
            deriv_result <- model_fit$predict(
              new_predictors = rbind(par),
              parallel = FALSE,
              cl = NULL,
              chunk_size = NULL,
              num_chunks = NULL,
              rem_chunks = NULL,
              B_predict = B_predict,
              take_first_derivatives = TRUE)
            gr <- .extract_first_deriv(deriv_result)
            gr_par <- rep(0, length(par))
            gr_raw <- min_or_max * gr
            ## [Change 2026-02-22] Assign gradients only to numeric predictor indices
            #  that correspond to variables being optimized
            numeric_vars_idx <- intersect(model_fit$numerics, vars_idx)
            gr_par[numeric_vars_idx] <- gr_raw[match(numeric_vars_idx, vars_idx)]
            if(select_vars_fl){
              gr_par[vars_idx]
            } else {
              gr_par
            }
          }
        },
        method = 'L-BFGS-B',
        lower = optim_lower,
        upper = optim_upper
      )

      ## [Change 2026-02-17] Reconstruct full predictor row from optim result
      if(select_vars_fl){
        dummy <- as.numeric(predictors_vals)
        dummy[vars_idx] <- as.numeric(opt$par)
        par <- rbind(dummy)
        colnames(par) <- colnames(predictors_vals)
      } else {
        par <- rbind(as.numeric(opt$par))
        if(!is.null(colnames(predictors_vals))){
          colnames(par) <- colnames(predictors_vals)
        }
      }

      return(par)
    })

    ## Find the global optimum out of all optimal-per-partitions
    best_per_partition <- Reduce("rbind", best_per_partition)
    preds <- model_fit$predict(new_predictors = best_per_partition,
                               B_predict = B_predict)
    global_max <- which.max(-min_or_max*preds)

    ## Return the optimized values
    extr <- best_per_partition[global_max, ,drop=FALSE]
    colnames(extr) <- colnames(predictors)
    return(list(
      t = extr,
      y = preds[global_max]
    ))
  }

  ## One-dimensional plotting function
  # [Change 2026-02-14] Introduced legend_order option
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
                               plot_fxn_1d,
                               legend_args,
                               color_function,
                               legend_order = NULL, # [Change 2026-02-14]
                               ...) {

    ## For preventing stack issues
    model_fit <- modfit
    drop(modfit)

    ## Linear term and name
    if(length(model_fit$power1_cols) > 0){
      xvals <- lapply(model_fit$X, function(x) x[,model_fit$power1_cols[1]])
    } else {
      xvals <- lapply(model_fit$X, function(x) x[,model_fit$nonspline_cols[1]])
    }

    ## For customizing xlab and legend predictor label
    v1 <- colnames(model_fit$X[[1]])[c(model_fit$power1_cols,
                                       model_fit$nonspline_cols)[1]]
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
    cols = color_function(model_fit$K+1)

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
    # plot_fxn_1d can be plot() or points()
    plot_fxn_1d(xvals[[1]],
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
        names(coefs) <- rownames(coefs)
        names(coefs) <- gsub(
          rownames(model_fit$B[[k]])[2],
          xlab,
          names(coefs)
        )
        names(coefs) <- gsub(v1, custom_predictor_lab, names(coefs))
        paste0(custom_formula_lab, " = ", paste(coefs, names(coefs),
                                                collapse = " + "))
      })
      formulas <- gsub('intercept', '', formulas)
      formulas <- gsub('  ', ' ', formulas)

      ## [Change 2026-02-12] Apply custom legend ordering if specified
      if(!is.null(legend_order)){
        if(length(legend_order) == length(formulas)){
          formulas <- formulas[legend_order]
          cols <- cols[legend_order]
        } else if(include_warnings){
          warning("legend_order length does not match number of partitions; ",
                  "ignoring.")
        }
      }

      ## Create base legend arguments
      legend_base_args <- list(
        x = legend_pos,
        legend = formulas,
        col = cols,
        lwd = 2,
        cex = text_size_formula
      )

      ## Merge with user-supplied legend arguments if any
      if(length(legend_args) > 0) {
        ## If legend_args is a named list, use it directly
        if(is.list(legend_args)) {
          legend_final_args <- utils::modifyList(legend_base_args, legend_args)
        }
        ## If it is not a list, try to convert it first
        else {
          legend_args_list <- as.list(legend_args)
          if(!is.null(names(legend_args))) {
            names(legend_args_list) <- names(legend_args)
            legend_final_args <- utils::modifyList(legend_base_args,
                                                   legend_args_list)
          } else {
            legend_final_args <- legend_base_args
          }
        }
      } else {
        legend_final_args <- legend_base_args
      }

      ## Call legend with final arguments
      do.call(graphics::legend, legend_final_args)
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
                               color_function,
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

    ## For swapping out custom labels from formulas
    if(is.null(og_cols)){
      og_cols <- model_fit$raw_expansion_names[c(model_fit$power1_cols,
                                                 model_fit$nonspline_cols)]
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
        names(coefs) <- rownames(coefs)
        names(coefs) <- gsub("\\^2", "<sup>2</sup>", names(coefs))
        names(coefs) <- gsub("\\^3", "<sup>3</sup>", names(coefs))
        names(coefs) <- gsub("\\^4", "<sup>4</sup>", names(coefs))
        names(coefs) <- gsub(og_cols[as.numeric(substr(v1, 2, nchar(v1)-1))],
                             custom_predictor_lab1, names(coefs))
        names(coefs) <- gsub(og_cols[as.numeric(substr(v2, 2, nchar(v2)-1))],
                             custom_predictor_lab2, names(coefs))
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
                      colors = color_function(model_fit$K+1),
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
                             legend_args = list(),
                             new_predictors = NULL,
                             xlim = NULL,
                             ylim = NULL,
                             color_function = NULL,
                             add = FALSE,
                             vars = c(),
                             legend_order = NULL, # [Change 2026-02-14] Include
                             ...){

    ## add = TRUE has the effect of overlaying the plot over an existing one
    # only for 1D
    if(add){
      plot_fxn_1d = graphics::points
    } else {
      plot_fxn_1d = graphics::plot
    }

    ## Check compatibility, that new_predictors should be a matrix
    if(any(!is.null(new_predictors))){
      new_predictors <- try(methods::as(cbind(new_predictors), 'matrix'), silent = TRUE)
      if(any(inherits(new_predictors, 'try-error'))){
        stop('\n \t new_predictors should be coercible to a matrix. \n')
      }
    }

    ## Default text_size_formula depends on q
    if(is.null(text_size_formula)){
      text_size_formula <- ifelse(model_fit_in$q == 1 | length(vars) == 1,
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
                                   expansions_only = TRUE)
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
    if(model_fit_in$q == 1 | length(vars) == 1){
      ## Color function takes in single argument (K+1) and returns colors we use
      if(is.null(color_function)){
        color_function <- grDevices::rainbow
      }
      if(length(vars) == 1){
        ## Isolate variables of interest
        if(inherits(vars, 'numeric')){
          cols <- paste0('_', vars, '_')
        } else if(!any(is.null(og_cols))){
          inds <- which(og_cols %in% vars)
          if(length(inds) != 1){
            stop('\n\t vars is not an original predictor name in the data')
          }
          cols <- paste0('_', inds, '_')
        } else {
          stop('\n\tInput predictors have no names, use column indices for vars')
        }
        keeps <- unlist(c(1, sapply(cols, function(col)grep(col,
                                                            model_fit_in$raw_expansion_names))))
        if(length(keeps) < 2){
          stop('\n\t Column indices provided are not present in data\n')
        }
        for(k in 1:(model_fit_in$K + 1)){
          model_fit_in$B[[k]] <-
            model_fit_in$B[[k]][keeps,,drop=FALSE]
          model_fit_in$B_raw[[k]] <-
            model_fit_in$B_raw[[k]][keeps,,drop=FALSE]
        }
        if(length(model_fit_in$power1_cols) > 0){
          model_fit_in$power1_cols <- model_fit_in$power1_cols[
            model_fit_in$power1_cols %in% keeps
          ]
        }
        if(length(model_fit_in$nonspline_cols) > 0){
          model_fit_in$nonspline_cols <- model_fit_in$nonspline_cols[
            model_fit_in$nonspline_cols %in% keeps
          ]
        }
      }
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
                       plot_fxn_1d,
                       legend_args,
                       color_function,
                       legend_order = legend_order, # [Change 2026-02-14] Incl.
                       ...)
      ## 2-D plotting
    } else if(model_fit_in$q == 2 | length(vars) == 2){
      ## Color function takes in single argument (K+1) and returns colors we use
      if(is.null(color_function)){
        color_function <- grDevices::colorRampPalette(
          RColorBrewer::brewer.pal(8, "Spectral"))
      }
      ## Isolate variables of interest
      if(length(vars) == 2){
        if(inherits(vars, 'numeric')){
          cols <- paste0('_', vars, '_')
        } else if(!is.null(og_cols)){
          inds <- which(og_cols %in% vars)
          if(length(inds) != 2){
            stop('\n\tOne or both vars are not original predictor names of data')
          }
          cols <- paste0('_', inds, '_')
        } else {
          stop('\n\t Original predictor names not present, try numeric indices',
               ' for vars\n')
        }
        keeps <- unlist(c(1, sapply(cols, function(col)grep(col,
                                                            model_fit_in$raw_expansion_names))))
        if(length(keeps) < 3){
          stop('\n\t Column indices provided are not present in data\n')
        }
        for(k in 1:(model_fit_in$K+1)){
          model_fit_in$B[[k]] <-
            model_fit_in$B[[k]][keeps,,drop=FALSE]
          model_fit_in$B_raw[[k]] <-
            model_fit_in$B_raw[[k]][keeps,,drop=FALSE]
        }
        if(length(model_fit_in$power1_cols) > 0){
          model_fit_in$power1_cols <- model_fit_in$power1_cols[
            model_fit_in$power1_cols %in% keeps
          ]
        }
        if(length(model_fit_in$nonspline_cols) > 0){
          model_fit_in$nonspline_cols <- model_fit_in$nonspline_cols[
            model_fit_in$nonspline_cols %in% keeps
          ]
        }
      }
      plot_lgspline_2d(model_fit_in,
                       show_formulas,
                       digits,
                       custom_response_lab,
                       custom_formula_lab,
                       custom_predictor_lab1,
                       custom_predictor_lab2,
                       custom_title,
                       text_size_formula,
                       color_function,
                       ...)
    } else if(include_warnings){
      warning("\n \t No default plotting functions implemented for q > 2 \n")
    }
  }

  ## [Change 2026-02-21] Store arguments
  model_fit$.fit_call_args <- list(
    standardize_response               = standardize_response,
    standardize_predictors_for_knots   = standardize_predictors_for_knots,
    standardize_expansions_for_fitting = standardize_expansions_for_fitting,
    family                             = family,
    glm_weight_function                = glm_weight_function,
    schur_correction_function          = schur_correction_function,
    need_dispersion_for_estimation     = need_dispersion_for_estimation,
    dispersion_function                = dispersion_function,
    K                                  = K,
    cluster_on_indicators              = cluster_on_indicators,
    use_custom_bfgs                    = use_custom_bfgs,
    delta                              = delta,
    tol                                = tol,
    initial_wiggle                     = initial_wiggle,
    initial_flat                       = initial_flat,
    wiggle_penalty                     = wiggle_penalty,
    flat_ridge_penalty                 = flat_ridge_penalty,
    unique_penalty_per_partition        = unique_penalty_per_partition,
    unique_penalty_per_predictor        = unique_penalty_per_predictor,
    meta_penalty                       = meta_penalty,
    predictor_penalties                = predictor_penalties,
    partition_penalties                 = partition_penalties,
    include_quadratic_terms            = include_quadratic_terms,
    include_cubic_terms                = include_cubic_terms,
    include_quartic_terms              = include_quartic_terms,
    include_2way_interactions          = include_2way_interactions,
    include_3way_interactions          = include_3way_interactions,
    include_quadratic_interactions     = include_quadratic_interactions,
    offset                             = offset,
    just_linear_with_interactions      = just_linear_with_interactions,
    just_linear_without_interactions   = just_linear_without_interactions,
    exclude_interactions_for           = exclude_interactions_for,
    exclude_these_expansions           = exclude_these_expansions,
    custom_basis_fxn                   = custom_basis_fxn,
    include_constrain_fitted           = include_constrain_fitted,
    include_constrain_first_deriv      = include_constrain_first_deriv,
    include_constrain_second_deriv     = include_constrain_second_deriv,
    include_constrain_interactions     = include_constrain_interactions,
    cl                                 = cl,
    chunk_size                         = chunk_size,
    parallel_eigen                     = parallel_eigen,
    parallel_trace                     = parallel_trace,
    parallel_aga                       = parallel_aga,
    parallel_matmult                   = parallel_matmult,
    parallel_unconstrained             = parallel_unconstrained,
    parallel_find_neighbors            = parallel_find_neighbors,
    parallel_penalty                   = parallel_penalty,
    parallel_make_constraint           = parallel_make_constraint,
    unconstrained_fit_fxn              = unconstrained_fit_fxn,
    keep_weighted_Lambda               = keep_weighted_Lambda,
    iterate_tune                       = iterate_tune,
    iterate_final_fit                  = iterate_final_fit,
    blockfit                           = blockfit,
    qp_score_function                  = qp_score_function,
    qp_observations                    = qp_observations,
    qp_Amat                           = qp_Amat,
    qp_bvec                           = qp_bvec,
    qp_meq                            = qp_meq,
    qp_positive_derivative             = qp_positive_derivative,
    qp_negative_derivative             = qp_negative_derivative,
    qp_positive_2ndderivative          = qp_positive_2ndderivative,
    qp_negative_2ndderivative          = qp_negative_2ndderivative,
    qp_monotonic_increase              = qp_monotonic_increase,
    qp_monotonic_decrease              = qp_monotonic_decrease,
    qp_range_upper                     = qp_range_upper,
    qp_range_lower                     = qp_range_lower,
    qp_Amat_fxn                       = qp_Amat_fxn,
    qp_bvec_fxn                       = qp_bvec_fxn,
    qp_meq_fxn                        = qp_meq_fxn,
    constraint_values                  = constraint_values,
    constraint_vectors                 = constraint_vectors,
    return_G                           = return_G,
    return_Ghalf                       = return_Ghalf,
    return_U                           = return_U,
    estimate_dispersion                = estimate_dispersion,
    unbias_dispersion                  = unbias_dispersion,
    return_varcovmat                   = return_varcovmat,
    return_lagrange_multipliers        = return_lagrange_multipliers,
    custom_penalty_mat                 = custom_penalty_mat,
    cluster_args                       = cluster_args,
    dummy_dividor                      = dummy_dividor,
    dummy_adder                        = dummy_adder,
    auto_encode_factors                = auto_encode_factors,
    observation_weights                = observation_weights,
    do_not_cluster_on_these            = do_not_cluster_on_these,
    neighbor_tolerance                 = neighbor_tolerance,
    no_intercept                       = no_intercept,
    og_cols                            = og_cols
  )

  ## Set S3 class
  class(model_fit) <- "lgspline"
  return(model_fit)
}

#' Low-Level Fitting for Lagrangian Smoothing Splines
#'
#' @description
#' The core function for fitting Lagrangian smoothing splines with
#' less user-friendliness. Called internally by \code{\link{lgspline}} after
#' formula parsing, factor encoding, and correlation-structure setup.
#'
#' @details
#' \code{lgspline.fit} performs the following steps:
#' \enumerate{
#'   \item Polynomial expansion and predictor standardization.
#'   \item Knot placement and partitioning (k-means or custom).
#'   \item Constraint matrix \eqn{\mathbf{A}} construction. Only a linearly
#'         independent subset of columns is retained via pivoted QR decomposition.
#'   \item Penalty tuning via GCV (exponential parameterization) or
#'         use of previously tuned penalties.
#'   \item Final coefficient estimation via one of three paths:
#'         \itemize{
#'           \item \strong{Blockfit option} (when \code{blockfit = TRUE},
#'                 flat columns are non-empty, \code{K > 0}, and no correlation
#'                 structure): Routes through \code{blockfit_solve} for
#'                 backfitting with mixed spline and non-interactive linear terms.
#'                 Falls back to \code{get_B} on failure.
#'           \item \strong{Standard \code{get_B}} path: Three internal
#'                 computational paths — GEE (damped SQP with correlation
#'                 structures), Gaussian identity (closed-form OLS projection),
#'                 and general GLM (unconstrained fit + Lagrangian projection
#'                 with optional IRLS loop).
#'         }
#'   \item Post-fit inference: \eqn{\mathbf{U}}, trace, dispersion,
#'         variance-covariance matrix, and optionally Lagrange multipliers.
#'         When \code{VhalfInv} is non-\code{NULL}, these are computed from
#'         the whitened Gram matrices
#'         \eqn{\mathbf{X}^{\top}\mathbf{V}^{-1}\mathbf{X}} via the full
#'         penalized GLS information
#'         \eqn{\mathbf{G}_{\mathrm{correct}} =
#'         (\mathbf{X}^{\top}\mathbf{V}^{-1}\mathbf{X} +
#'         \boldsymbol{\Lambda})^{-1}}.
#' }
#'
#' \strong{Dummy fit.} When \code{dummy_fit = TRUE}, an early-return path
#' skips the expensive fitting steps (\code{compute_G_eigen}, \code{get_B},
#' trace computation, variance-covariance matrix) while retaining all
#' penalty, partitioning, and design matrix information. Coefficients are
#' set to zero. This replaces the deprecated \code{expansions_only} argument.
#'
#' @return A list containing the fitted model components, forming the core
#' structure used internally by \code{\link{lgspline}} and its associated methods.
#' This function is primarily intended for internal use or advanced users needing
#' direct access to fitting components. The returned list contains numerous elements,
#' typically including:
#' \describe{
#'   \item{y}{The original response vector provided.}
#'   \item{ytilde}{The fitted values on the original response scale. Set to
#'     \code{rep(0, N)} when \code{dummy_fit = TRUE}.}
#'   \item{X}{A list, with each element the design matrix (\eqn{\mathbf{X}_{k}})
#'     for partition k, on the unstandardized expansion scale.}
#'   \item{A}{The constraint matrix (\eqn{\mathbf{A}}) encoding smoothness and
#'     any other linear equality constraints. Reduced to linearly independent
#'     columns via pivoted QR decomposition.}
#'   \item{B}{A list of the final fitted coefficient vectors
#'     (\eqn{\boldsymbol{\beta}_{k}}) for each partition k, on the original
#'     predictor/response scale.}
#'   \item{B_raw}{A list of fitted coefficient vectors on the internally
#'     standardized scale used during fitting.}
#'   \item{K, p, q, P, N}{Key dimensions: number of internal knots (K), basis
#'     functions per partition (p), original predictors (q), total coefficients
#'     (P), and sample size (N).}
#'   \item{penalties}{A list containing the final penalty components used
#'     (e.g., \code{Lambda}, \code{L1}, \code{L2}, \code{L_predictor_list},
#'     \code{L_partition_list}). See \code{\link{compute_Lambda}}.}
#'   \item{knot_scale_transf, knot_scale_inv_transf}{Functions to transform
#'     predictors to/from the scale used for knot placement.}
#'   \item{knots}{Matrix or vector of knot locations on the original predictor
#'     scale (NULL if K=0 or q > 1).}
#'   \item{partition_codes}{Vector assigning each original observation to a partition.}
#'   \item{partition_bounds}{Internal representation of partition boundaries.}
#'   \item{make_partition_list}{List containing centers, knot midpoints, neighbor
#'     info, and assignment function from partitioning (NULL if K=0 or 1D).
#'     See \code{\link{make_partitions}}.}
#'   \item{knot_expand_function, assign_partition}{Internal functions for
#'     partitioning data. See \code{\link{knot_expand_list}}.}
#'   \item{predict}{The primary function embedded in the object for generating
#'     predictions on new data. For multi-predictor models,
#'     \code{take_first_derivatives = TRUE} returns derivatives as a named list
#'     of per-variable derivative vectors rather than a concatenated vector.
#'     See \code{\link{predict.lgspline}}.}
#'   \item{family}{The \code{\link[stats]{family}} object or custom list used.}
#'   \item{estimate_dispersion, unbias_dispersion}{Logical flags related to
#'     dispersion estimation settings.}
#'   \item{sigmasq_tilde}{The estimated (or fixed, if
#'     \code{estimate_dispersion = FALSE}) dispersion parameter
#'     \eqn{\tilde{\sigma}^{2}}. When \code{VhalfInv} is non-\code{NULL} and
#'     family is Gaussian with identity link, no \code{observation_weights},
#'     computed from whitened residuals
#'     \eqn{\frac{1}{N-edf}(\mathbf{V}^{-1/2}(\mathbf{y} - \hat{\mathbf{y}}))^{2}}.
#'     Set to \code{NA_real_} when \code{dummy_fit = TRUE}. \eqn{N-edf} is
#'     replaced with \eqn{N} when \code{unbias_dispersion = TRUE} which is
#'     the default except for Gaussian response with identity link.}
#'   \item{backtransform_coefficients, forwtransform_coefficients}{Functions to
#'     convert coefficients between standardized and original scales.}
#'   \item{mean_y, sd_y}{Mean and standard deviation used for standardizing
#'     the response.}
#'   \item{og_order, order_list}{Information mapping original data order to
#'     partitioned order.}
#'   \item{constraint_values, constraint_vectors}{User-supplied additional linear
#'     equality constraints.}
#'   \item{expansion_scales}{Scaling factors applied to basis expansions during
#'     fitting (if \code{standardize_expansions_for_fitting = TRUE}).}
#'   \item{take_derivative, take_interaction_2ndderivative,
#'     get_all_derivatives_insample}{Functions related to computing derivatives
#'     of the fitted spline.}
#'   \item{numerics, power1_cols, ..., nonspline_cols}{Integer vectors storing
#'     column indices identifying different types of terms in the basis expansion.}
#'   \item{return_varcovmat}{Logical indicating if variance matrix calculation
#'     was requested.}
#'   \item{raw_expansion_names}{Original generated names for basis expansion
#'     columns (before potential renaming if input predictors had names).}
#'   \item{std_X, unstd_X}{Functions to standardize/unstandardize design matrices
#'     according to \code{expansion_scales}.}
#'   \item{parallel_cluster_supplied}{Logical indicating if a parallel cluster
#'     was used.}
#'   \item{weights}{The original observation weights provided (potentially
#'     reformatted).}
#'   \item{VhalfInv}{The fixed \eqn{\mathbf{V}^{-1/2}} matrix if supplied.}
#'   \item{quadprog_list}{List containing components related to quadratic
#'     programming constraints, if used.}
#'   \item{G}{List of unscaled variance-covariance matrices
#'     \eqn{\mathbf{G}_{k}} per partition, returned if \code{return_G = TRUE}.
#'     When \code{VhalfInv} is non-\code{NULL}, recomputed from whitened Gram
#'     matrices. Set to \code{NULL} when \code{dummy_fit = TRUE}.}
#'   \item{Ghalf}{List of \eqn{\mathbf{G}_{k}^{1/2}} matrices, returned if
#'     \code{return_Ghalf = TRUE}. When \code{VhalfInv} is non-\code{NULL},
#'     the full \eqn{\mathbf{G}_{\mathrm{correct}}^{1/2}} is used for
#'     posterior draws and variance-covariance computation.
#'     Set to \code{NULL} when \code{dummy_fit = TRUE}.}
#'   \item{U}{Constraint projection matrix \eqn{\mathbf{U}}, returned if
#'     \code{return_U = TRUE}. Set to \code{NULL} when
#'     \code{dummy_fit = TRUE}.}
#'   \item{trace_XUGX}{The effective degrees of freedom
#'     \eqn{\mathrm{trace}(\textbf{D}\mathbf{X}\mathbf{U}\mathbf{G}\mathbf{X}^{\top})}.
#'     When \code{VhalfInv} is non-\code{NULL}, computed as
#'     \eqn{\|\textbf{D}^{1/2}\mathbf{V}^{-1/2}\mathbf{X}\mathbf{G}_{\mathrm{correct}}^{1/2}\|_{F}^{2}}. Set to \code{NA_real_} when \code{dummy_fit = TRUE}.}
#'   \item{varcovmat}{The final variance-covariance matrix of the estimated
#'     coefficients. Computed via the outer-product form
#'     \eqn{\sigma^{2}(\mathbf{U}\mathbf{G}^{1/2})(\mathbf{U}\mathbf{G}^{1/2})^{\top}}
#'     for numerical stability. When \code{VhalfInv} is non-\code{NULL}, uses
#'     the full \eqn{\mathbf{G}_{\mathrm{correct}}^{1/2}} in place of
#'     block-diagonal \eqn{\mathbf{G}^{1/2}}. Returned if
#'     \code{return_varcovmat = TRUE}. Set to \code{NULL} when
#'     \code{dummy_fit = TRUE}.}
#'   \item{lagrange_multipliers}{Vector of Lagrangian multipliers if
#'     \code{return_lagrange_multipliers = TRUE} and constraints are active via
#'     the formulation
#'     \eqn{(\mathbf{A}^{\top}\mathbf{G}\mathbf{A})^{-1}\mathbf{A}^{\top}(\hat{\boldsymbol{\beta}} - \boldsymbol{\beta_0})}.
#'     \code{NULL} if no constraints are active. Multipliers for inequality
#'     constraints are not returned}.
#' }
#' Note that the exact components returned depend heavily on the function
#' arguments (e.g., values of \code{return_G}, \code{return_varcovmat}, etc.)
#' and whether \code{dummy_fit = TRUE}.
#'
#' @usage
#' lgspline.fit(predictors, y = NULL, standardize_response = TRUE,
#'              standardize_expansions_for_fitting = TRUE, family = gaussian(),
#'              glm_weight_function, schur_correction_function,
#'              need_dispersion_for_estimation = FALSE,
#'              dispersion_function,
#'              K = NULL, custom_knots = NULL, cluster_on_indicators = FALSE,
#'              make_partition_list = NULL, previously_tuned_penalties = NULL,
#'              smoothing_spline_penalty = NULL, opt = TRUE, use_custom_bfgs = TRUE,
#'              delta = NULL, tol = 10*sqrt(.Machine$double.eps),
#'              initial_wiggle = exp(c(-25, -14, -3, -0.2, 1)),
#'              initial_flat = c(0.1, 10), wiggle_penalty = 2e-07,
#'              flat_ridge_penalty = 0.5, unique_penalty_per_partition = TRUE,
#'              unique_penalty_per_predictor = TRUE, meta_penalty = 1e-08,
#'              predictor_penalties = NULL, partition_penalties = NULL,
#'              include_quadratic_terms = TRUE, include_cubic_terms = TRUE,
#'              include_quartic_terms = FALSE, include_2way_interactions = TRUE,
#'              include_3way_interactions = TRUE,
#'              include_quadratic_interactions = FALSE,
#'              offset = c(), just_linear_with_interactions = NULL,
#'              just_linear_without_interactions = NULL,
#'              exclude_interactions_for = NULL,
#'              exclude_these_expansions = NULL, custom_basis_fxn = NULL,
#'              include_constrain_fitted = TRUE,
#'              include_constrain_first_deriv = TRUE,
#'              include_constrain_second_deriv = TRUE,
#'              include_constrain_interactions = TRUE, cl = NULL, chunk_size = NULL,
#'              parallel_eigen = TRUE, parallel_trace = FALSE, parallel_aga = FALSE,
#'              parallel_matmult = FALSE, parallel_unconstrained = FALSE,
#'              parallel_find_neighbors = FALSE, parallel_penalty = FALSE,
#'              parallel_make_constraint = FALSE,
#'              unconstrained_fit_fxn = unconstrained_fit_default,
#'              keep_weighted_Lambda = FALSE, iterate_tune = TRUE,
#'              iterate_final_fit = TRUE, blockfit = FALSE,
#'              qp_score_function,
#'              qp_observations = NULL, qp_Amat = NULL, qp_bvec = NULL, qp_meq = 0,
#'              qp_positive_derivative = FALSE, qp_negative_derivative = FALSE,
#'              qp_positive_2ndderivative = FALSE, qp_positive_2ndderivative = FALSE,
#'              qp_monotonic_increase = FALSE, qp_monotonic_decrease = FALSE,
#'              qp_range_upper = NULL, qp_range_lower = NULL, qp_Amat_fxn = NULL,
#'              qp_bvec_fxn = NULL, qp_meq_fxn = NULL, constraint_values = cbind(),
#'              constraint_vectors = cbind(), return_G = TRUE, return_Ghalf = TRUE,
#'              return_U = TRUE, estimate_dispersion = TRUE,
#'              unbias_dispersion = TRUE,
#'              return_varcovmat = TRUE, return_lagrange_multipliers = FALSE,
#'              custom_penalty_mat = NULL,
#'              cluster_args = c(custom_centers = NA, nstart = 10),
#'              dummy_dividor = 1.2345672152894e-22,
#'              dummy_adder = 2.234567210529e-18,
#'              verbose = FALSE, verbose_tune = FALSE,
#'              dummy_fit = FALSE, auto_encode_factors = TRUE,
#'              observation_weights = NULL, do_not_cluster_on_these = c(),
#'              neighbor_tolerance = 1 + 1e-16, no_intercept = FALSE,
#'              VhalfInv = NULL, Vhalf = NULL, include_warnings = TRUE,
#'              og_cols = NULL, ...)
#'
#' @inheritParams lgspline
#' @param og_cols Original predictor names
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
                         schur_correction_function = function(X,
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
                                                        VhalfInv,
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
                         initial_wiggle = exp(c(-25, -14, -3, -0.2, 1)),
                         initial_flat = c(0.1, 10),
                         wiggle_penalty = 2e-7,
                         flat_ridge_penalty = 0.5,
                         unique_penalty_per_partition = TRUE,
                         unique_penalty_per_predictor = TRUE,
                         meta_penalty = 1e-8,
                         predictor_penalties = NULL,
                         partition_penalties = NULL,
                         include_quadratic_terms = TRUE,
                         include_cubic_terms = TRUE,
                         include_quartic_terms = FALSE,
                         include_2way_interactions = TRUE,
                         include_3way_interactions = TRUE,
                         include_quadratic_interactions = FALSE,
                         offset = c(),
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
                         blockfit = FALSE,
                         qp_score_function = function(X,
                                                      y,
                                                      mu,
                                                      order_list,
                                                      dispersion,
                                                      VhalfInv,
                                                      observation_weights,
                                                      ...) {
                           if(!is.null(observation_weights)) {
                             crossprod(X, cbind((y - mu)*observation_weights))
                           } else {
                             crossprod(X, cbind(y - mu))
                           }
                         },
                         qp_observations = NULL,
                         qp_Amat = NULL,
                         qp_bvec = NULL,
                         qp_meq = 0,
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
                         constraint_values = cbind(),
                         constraint_vectors = cbind(),
                         return_G = TRUE,
                         return_Ghalf = TRUE,
                         return_U = TRUE,
                         estimate_dispersion = TRUE,
                         unbias_dispersion = TRUE,
                         return_varcovmat = TRUE,
                         return_lagrange_multipliers = FALSE,
                         custom_penalty_mat = NULL,
                         cluster_args = c(custom_centers = NA, nstart = 10),
                         dummy_dividor = 0.00000000000000000000012345672152894,
                         dummy_adder = 0.000000000000000002234567210529,
                         verbose = FALSE,
                         verbose_tune = FALSE,
                         dummy_fit = FALSE,
                         auto_encode_factors = TRUE,
                         observation_weights = NULL,
                         do_not_cluster_on_these = c(),
                         neighbor_tolerance = 1 + 1e-16,
                         no_intercept = FALSE,
                         VhalfInv = NULL,
                         Vhalf = NULL,
                         include_warnings = TRUE,
                         og_cols = NULL,
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

  ## [Change 2026-02-17] Resolve any remaining character entries in
  #  do_not_cluster_on_these to numeric column indices. Handles the
  #  non-formula path where users pass a predictor matrix with colnames
  #  and character do_not_cluster_on_these.
  if(length(do_not_cluster_on_these) > 0 &&
     is.character(do_not_cluster_on_these)){
    pred_colnames <- colnames(predictors)
    if(!is.null(pred_colnames)){
      resolved <- unlist(lapply(
        do_not_cluster_on_these, function(var){
          idx <- which(pred_colnames == var)
          if(length(idx) == 0) idx <- grep(var, pred_colnames)
          idx
        }))
      ## Keep any that were already numeric
      do_not_cluster_on_these <- unique(c(
        do_not_cluster_on_these[is.numeric(do_not_cluster_on_these)],
        resolved
      ))
    } else if(include_warnings){
      warning("Character entries in do_not_cluster_on_these cannot be ",
              "resolved without column names on predictors. ",
              "Use numeric column indices instead.")
      do_not_cluster_on_these <-
        do_not_cluster_on_these[is.numeric(do_not_cluster_on_these)]
    }
  }

  ## Accept raw predictors (the "T" matrix) and get dimensions
  predictors <- methods::as(predictors,'matrix')
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
    stop(
      '\n \t Elements in just_linear_with_interactions, ',
      'just_linear_without_interactions, and/or exclude_interactions_for are',
      ' greater than the number of columns of predictors matrix. \n')
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
    include_2way_interactions <- FALSE
    include_3way_interactions <- FALSE
    include_quadratic_interactions <- FALSE
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
    stop('\n \t Must include at least 1 variable to cluster on if multiple',
         ' variables are present. \n')
  }

  ## K can't be greater than max of number of observations or q
  # for kmeans clustering purposes
  if(K >= max(c(nr, qcols))) {
    if(include_warnings){
      warning('\n \t Max (N, q) too samll for K. K = max(N, q) - 2 will be',
              ' used. \n')
    }
    K <- max(max(c(nr, qcols)) - 2, 0)
  }

  ## Detect if parallel, and K > 0
  if(any(!(is.null(cl))) & K > 0){
    if(inherits(cl, 'cluster')){
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

  ## Standardize outcome for identity link (agnostic to distribution)
  if(paste0(family)[2] == 'identity' &
     #paste0(family)[1] == 'gaussian' &
     length(unique(y)) > 1 &
     standardize_response){

    mean_y <- mean(y)
    sd_y <- try(sd(y),silent = TRUE)
    if(any(inherits(sd_y, 'try-error'))){
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
  colnm_expansions <- colnames(C)
  power1_cols <- 2:(length(numerics) + 1)
  if(length(numerics) == 0){
    power2_cols <- c()
    power1_cols <- c()
    include_constrain_second_deriv <- FALSE
  } else {
    power2_cols <- which(substr(colnm_expansions, nchar(colnm_expansions)-1,
                                nchar(colnm_expansions)) == '^2')
  }
  power3_cols <- which(substr(colnm_expansions, nchar(colnm_expansions)-1,
                              nchar(colnm_expansions)) == '^3')
  power4_cols <- which(substr(colnm_expansions, nchar(colnm_expansions)-1,
                              nchar(colnm_expansions)) == '^4')
  interaction_cols <- grep("_x_", colnm_expansions)
  if(length(numerics) > 2 & length(interaction_cols) > 0){
    triplet_cols <- interaction_cols[
      which(sapply(colnm_expansions[interaction_cols], function(col){
        grepl('_x_',substr(col, regexpr('_x_',col)[[1]]+3,nchar(col)))
      }))]
  } else {
    triplet_cols <- c()
  }
  quad_cols <- which(substr(colnm_expansions, nchar(colnm_expansions)-1,
                            nchar(colnm_expansions)) == "^2")
  interaction_quad_cols <- intersect(
    interaction_cols,quad_cols
  )
  interaction_single_cols <- interaction_cols[!(interaction_cols %in% c(
    triplet_cols, interaction_quad_cols
  ))]

  ## Append non-spline terms
  nonspline_cols <- c(
    which(colnm_expansions %in%
            c(
              paste0("_", just_linear_with_interactions, "_"),
              paste0("_", just_linear_without_interactions, "_")
            )
    )
  )
  nonspline_cols <- nonspline_cols[!(nonspline_cols %in%
                                       c(power1_cols,
                                         interaction_single_cols,
                                         interaction_quad_cols,
                                         triplet_cols))]

  ## Standardize columns of C using expansion/(q0.69 - q0.31)
  # This is a p-1 length vector, it excludes the intercept
  expansion_scales <- apply(C[,-intercept_col,drop=FALSE], 2, function(x){
    1
    if(length(unique(x)) >= 2){
      ## Near sigma for a normal distribution
      # (i.e. this is close to 1 for N(0,1))
      abs(quantile(x, 0.69) - quantile(x, 0.31))
    } else {
      1
    }
  })
  ## For offsets, keep the linear term unscaled
  if(length(offset) > 0){
    offset_ind <- which(colnm_expansions %in% paste0(
      '_', offset, '_'
    ))
    expansion_scales[offset_ind-1] <- 1
  }
  expansion_scales[expansion_scales == 0] <- 1
  names(expansion_scales) <- colnm_expansions[-intercept_col, drop=FALSE]
  # Set back to 1 if not desired
  if(!standardize_expansions_for_fitting){
    expansion_scales <- 0*expansion_scales + 1
  }

  ## Function to un-standardize columns of C
  std_X <- function(unstd_X_in){
    sweep(unstd_X_in, 2, c(1, expansion_scales), "/")
  }
  max_C <- apply(C, 2, max)
  min_C <- apply(C, 2, min)
  max_min_C <- rbind(c(max_C), c(min_C))
  C <- std_X(C)
  unstd_X <- function(std_X_in){
    sweep(std_X_in, 2, c(1, expansion_scales), "*")
  }

  ## If no intercept enforced, include constraint on A indicating this
  if(no_intercept & length(constraint_vectors) < 1){
    constr <- sapply(1:(K+1), function(k){
      vec <- rep(0, nc*(K+1))
      vec[nc*(k-1) + 1] <- 1
      vec
    })
    constraint_vectors <- cbind(constr)
    constraint_values <- 0*constraint_vectors
  } else if(no_intercept){
    constr <- sapply(1:(K+1), function(k){
      vec <- rep(0, nc*(K+1))
      vec[nc*(k-1) + 1] <- 1
      vec
    })
    constraint_vectors <- cbind(constraint_vectors,
                                constr)
    constraint_values <- cbind(rowSums(cbind(constraint_values,
                                             0*constr
    )))
  }

  ## Repeat analogously for offsets if present
  if(length(offset) > 0 & length(constraint_vectors) < 1){
    offset_ind <- which(colnm_expansions %in% paste0(
      '_', offset, '_'
    ))
    constr <- Reduce('cbind', lapply(1:length(offset_ind), function(o){
      rbind(sapply(1:(K+1), function(k){
        vec <- rep(0, nc*(K+1))
        vec[nc*(k-1) + offset_ind[o]] <- 1
        vec
      }))
    }))
    constraint_vectors <- constr
    constraint_values <- cbind(rowSums(constr))
  } else if(length(offset) > 0){
    offset_ind <- which(colnm_expansions %in% paste0(
      '_', offset, '_'
    ))
    constr <- Reduce('cbind', lapply(1:length(offset_ind), function(o){
      rbind(sapply(1:(K+1), function(k){
        vec <- rep(0, nc*(K+1))
        vec[nc*(k-1) + offset_ind[o]] <- 1
        vec
      }))
    }))
    constraint_vectors <- cbind(constraint_vectors,
                                constr)
    constraint_values <- cbind(rowSums(cbind(constraint_values,
                                             constr)))
  }

  ## Adjust coefficients after un-standardizing
  backtransform_coefficients <- function(coef) {
    # Extract intercept and slope coefficients
    intercept <- coef[intercept_col]
    slopes <- coef[-intercept_col]

    # Back-transform slope coefficients
    backtransformed_slopes <- slopes / expansion_scales

    # Combine intercept and back-transformed slopes
    cbind(c(intercept, backtransformed_slopes))
  }

  ## Adjust coefficients for future standardizing
  forwtransform_coefficients <- function(coef) {
    # Extract intercept and slope coefficients
    intercept <- coef[intercept_col]
    slopes <- coef[-intercept_col]

    # Back-transform slope coefficients
    forwtransformed_slopes <- slopes * expansion_scales

    # Combine intercept and back-transformed slopes
    cbind(c(intercept, forwtransformed_slopes))
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

  ## Assign y to their partitions (y_og saves original y unstandardized)
  y <- knot_expand_list(partition_codes,
                        partition_bounds,
                        nr,
                        cbind(y),
                        K)

  ## Get observation weight expansions
  if(any(!is.null(observation_weights))){
    ## Coerce to N x 1 vector if not already
    if(nrow(cbind(observation_weights)) != nr |
       ncol(cbind(observation_weights)) != 1){
      stop('\n \t Observation weights must be an N x 1 vector. \n')
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
    observation_weights <- lapply(1:(K+1), function(k)cbind(
      rep(1, length(y[[k]]))))
    homogenous_weights <- TRUE
  }

  ## Save the original ordering to each partition
  order_list <- knot_expand_list(partition_codes,
                                 partition_bounds,
                                 nr,
                                 cbind(1:nr),
                                 K)
  og_order <- order(unlist(order_list))

  ## If custom variance-covariance structure specified
  if(!is.null(VhalfInv) & !dummy_fit){
    VhalfInv <- try(methods::as(VhalfInv,'matrix'), silent = TRUE)
    if(any(inherits(VhalfInv, 'try-error'))){
      if(include_warnings){
        warning('\n \t VhalfInv cannot be converted to a N by N matrix, it ',
                'will not be considered here. \n')
      }
      VhalfInv <- NULL
    } else if(any(unique(dim(VhalfInv)) != nr)){
      if(include_warnings){
        warning('\n \t VhalfInv is not an N by N matrix; it will not be',
                ' considered here. \n')
      }
      VhalfInv <- NULL
    } else {

      if(verbose){
        cat("Applying Whitening Transform\n")
      }

      ## [Change 2026-02-16] Compute Vhalf unconditionally when
      #  VhalfInv is present. Previously only computed for non-Gaussian
      #  or non-identity link. Needed by get_B and blockfit_solve for
      #  GEE estimation regardless of family.
      if(is.null(Vhalf)){
        Vhalf <- invert(VhalfInv)
      }

      ## Overwrite dev.resids (if present) to match normal approximation
      if(!is.null(family$dev.resids)){
        family$dev.resids <- function(y, mu, wt){
          ((y-mu)^2)*wt
        }
        family$linkfun <- function(mu)mu
      }

      ## [Change 2026-02-16] Do NOT whiten X into per-partition form.
      #  The previous code applied V^{-1/2} to the block-diagonal
      #  design matrix, then extracted only the diagonal blocks back
      #  into the per-partition list. This silently discarded cross-
      #  partition contributions from off-diagonal blocks of V^{-1/2},
      #  corrupting all downstream Gram matrices, cross-products, and
      #  G computations under GEE.
      #
      #  Now, keep X and y un-whitened. The whitening transform is
      #  applied inside get_B and blockfit_solve where the full N x P
      #  matrix structure is available. The unwhitened X_gram used for
      #  penalty tuning is an acceptable approximation (the exact
      #  solve happens in the solver after tuning).
      y_expand_og <- y
      X_expand_og <- X
    }
  }

  ## Return derivatives per-partition of an expanded matrix
  #  [Change 2026-02-16] Swap-out and swap-back removed. X and y are
  #  now always unwhitened, so X_expand_og is identical to X.
  all_derivatives <- function(X,
                              just_first_derivatives = FALSE,
                              just_spline_effects = TRUE){
    lapply(X, function(C){
      make_derivative_matrix(
        nc,
        C,
        power1_cols,
        power2_cols,
        nonspline_cols,
        interaction_single_cols,
        interaction_quad_cols,
        triplet_cols,
        K,
        include_2way_interactions,
        include_3way_interactions,
        include_quadratic_interactions,
        colnm_expansions,
        expansion_scales,
        just_first_derivatives,
        just_spline_effects
      )
    })
  }

  if(verbose){
    cat("2nd Derivative Penalty\n")
  }

  ## Compute integrated squared second derivative of fitted function
  if(!(!(any(is.null(smoothing_spline_penalty))))){
    max_min_C <- std_X(max_min_C)
    smoothing_spline_penalty <-
      get_2ndDerivPenalty_wrapper(K,
                                  colnm_expansions,
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

  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  ## Constraint matrix A construction.
  #  Three branches: (1) 1-D predictor with no nonspline cols,
  #  (2) multi-predictor / nonspline present, (3) K = 0.
  #  All branches standardize A by expansion_scales before returning.
  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

  ## Making a constraint matrix
  A <- 0
  if((K > 0 & length(numerics) == 1 & length(nonspline_cols) == 0)){

    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    ## Branch 1: single continuous predictor, no nonspline columns.
    #  Knot expansions are evaluated at knot_values, then a single
    #  call to make_constraint_matrix builds all smoothness constraints.
    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

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
                                power2_cols,
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
                                colnm_expansions,
                                expansion_scales)

    ## Append user-supplied equality constraints (offset / no_intercept rows
    #  zeroed out so they do not conflict with those special constraints).
    if(length(constraint_vectors) > 0 & length(constraint_values > 0)){
      if(length(offset) > 0){
        offset_ind <- which(colnm_expansions %in% paste0(
          '_', offset, '_'
        ))
        offset_inds <- unlist(lapply(1:(K+1), function(k)nc*(k-1)+offset_ind))
        A[offset_inds,] <- 0
      }
      if(no_intercept){
        A[unlist(lapply(1:(K+1),function(k)nc*(k-1)+1)),] <- 0
      }
      A <- cbind(A, constraint_vectors)
    }
    A <- sweep(A, 1, rep(c(1, expansion_scales), K+1), "/")
    if(any(!is.finite(A))) stop(paste0('\n \t A is not finite \n',
                                       expansion_scales))

  } else if(K > 0){

    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    ## Branch 2: multi-predictor or nonspline columns present.
    #  knot_values may have more rows than K (one row per knot per
    #  predictor dimension), so we process in chunks of K rows and
    #  cbind the resulting constraint columns. A remainder chunk handles
    #  any leftover rows, padded with zeros to keep dimensions consistent.
    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    chunk <- nrow(knot_values) %/% K
    rem <- nrow(knot_values) %% K

    knot_values_perm <- knot_values[1:nrow(knot_values),,drop=FALSE]

    ## Optional parallel construction: each worker handles one chunk
    if(parallel & parallel_make_constraint){
      A <- Reduce("cbind",
                  parallel::parLapply(cl,
                                      1:chunk,
                                      function(i){
                                        knot_values_chunk <-
                                          knot_values_perm[1:K + (i-1)*K,,
                                                           drop=FALSE]
                                        CKnots_chunk <- rbind(
                                          get_polynomial_expansions(
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
                                            ...)
                                        )
                                        rownames(CKnots_chunk) <-
                                          rownames(knot_values_chunk)
                                        make_constraint_matrix(
                                          nc,
                                          CKnots_chunk,
                                          power1_cols,
                                          power2_cols,
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
                                          colnm_expansions,
                                          expansion_scales
                                        )
                                      }))
    } else {

      ## Sequential construction: cbind constraint columns chunk by chunk.
      #  The first iteration drops the placeholder column (A was initialized
      #  to 0, a scalar) before further cbinding.
      for(i in 1:chunk){
        knot_values_chunk <- knot_values_perm[1:K + (i-1)*K,,drop=FALSE]
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
        A <- cbind(A, make_constraint_matrix(nc,
                                             CKnots_chunk,
                                             power1_cols,
                                             power2_cols,
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
                                             colnm_expansions,
                                             expansion_scales))
        if(i == 1){
          A <- A[,-1,drop=FALSE]
        }
      }
    }

    ## Remainder chunk: rows of knot_values not covered by full chunks.
    #  Taken from the tail of knot_values_perm and zero-padded to K rows
    #  so make_constraint_matrix receives a conformable matrix.
    if(rem > 0){
      knot_values_chunk <-
        knot_values_perm[rev(c(nrow(knot_values_perm):1)[1:rem]),,drop=FALSE]
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

      A <- cbind(A, make_constraint_matrix(nc,
                                           CKnots_chunk,
                                           power1_cols,
                                           power2_cols,
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
                                           colnm_expansions,
                                           expansion_scales))
    }

    ## Remove all 0 columns
    A <- A[,which(apply(abs(A), 2, sum) > 1e-16),drop=FALSE]

    ## Append user-supplied equality constraints and standardize,
    #  zeroing offset / no_intercept rows as in Branch 1.
    if(length(constraint_vectors) > 0){
      if(length(offset) > 0){
        offset_ind <- which(colnm_expansions %in% paste0(
          '_', offset, '_'
        ))
        offset_inds <- unlist(lapply(1:(K+1), function(k)nc*(k-1)+offset_ind))
        A[offset_inds,] <- 0
      }
      if(no_intercept){
        A[unlist(lapply(1:(K+1),function(k)nc*(k-1)+1)),] <- 0
      }
      A <- cbind(A, constraint_vectors)
    }
    A <- sweep(A, 1, rep(c(1, expansion_scales), K+1), "/")
    if(any(!is.finite(A))) stop(paste0('\n \t A is not finite \n',
                                       expansion_scales))
    if(any(is.na(A))) stop(paste0('\n \t A is NA somewhere ',
                                  '(any(is.na(A)) == TRUE) \n',
                                  expansion_scales))

  } else {

    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    ## Branch 3: K = 0, no knots.
    #  No smoothness constraints exist; A is either the user-supplied
    #  constraint_vectors alone or NULL (unconstrained).
    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    if(length(constraint_vectors) > 0){
      A <- cbind(constraint_vectors)
      A <- sweep(A, 1, rep(c(1, expansion_scales), K+1), "/")
    } else {
      A <- NULL
    }
  }
  ## [Change 2026-02-15] select only a linearly-independent subset of columns
  if(!(any(is.null(A)))){
    qr_A <- qr(A)
    A_rank <- qr_A$rank
    if(A_rank < ncol(A)){
      A <- qr.Q(qr_A)[, 1:A_rank, drop = FALSE]
    }
    nca <- ncol(A)
  }

  ## With only one predictor, we really only need one penalty
  if(ncol(predictors) == 1){
    unique_penalty_per_predictor <- FALSE
  }

  ## Convert constraint_values, penalty setup,
  #  parallel export, X^TWX gram, SQP setup (quadprog).
  if(verbose){
    cat("Predictor-and-Partiiton Penalty Setup\n")
  }

  ## Getting unique initial penalties for predictors/partitions,
  #  if not specified
  penalty_vec <- c()
  if(unique_penalty_per_predictor & any(is.null(predictor_penalties))){
    predictor_penalties <- sapply(colnm_expansions[c(power1_cols,
                                                     nonspline_cols)],
                                  function(j)exp(rnorm(1, 0, 0.00001)))
    names(predictor_penalties) <- paste0('predictor',
                                         colnm_expansions[c(power1_cols,
                                                            nonspline_cols)])
    penalty_vec <- c(penalty_vec, predictor_penalties)
  } else if(unique_penalty_per_predictor){
    if(length(predictor_penalties) !=
       length(c(power1_cols, nonspline_cols))){
      stop(
        '\n \t Custom predictor_penalties is not the same length as number ',
        'of predictors in model. The number of penalties should coincide with ',
        'the number of predictors, if supplied. \n')
    }
    if(any(predictor_penalties <= 0)){
      stop(
        '\n \t All predictor_penalties must be > 0 if supplied. You can set',
        ' unique_penalty_per_predictor = FALSE to remove predictor penalties. \n')
    }
    names(predictor_penalties) <- paste0('predictor',
                                         colnm_expansions[c(power1_cols,
                                                            nonspline_cols)])
    penalty_vec <- c(penalty_vec, predictor_penalties)
  }
  if(unique_penalty_per_partition & any(is.null(partition_penalties))){
    partition_penalties <- sapply(1:(K+1),
                                  function(j)exp(rnorm(1, 0, 0.00001)))
    names(partition_penalties) <- paste0('partition', 1:(K+1))
    penalty_vec <- c(penalty_vec, partition_penalties)
  } else if(unique_penalty_per_partition){
    if(length(partition_penalties) !=
       (K+1)){
      stop('\n \t Custom partition_penalties is not the same length as number ',
           'of partitions in model. Try setting K manually, and ensuring that ',
           'the length of partition_penalties = K + 1. \n')
    }
    if(any(partition_penalties <= 0)){
      stop(
        '\n \t All partition_penalties must be > 0 if supplied. You can set',
        ' unique_penalty_per_partition = FALSE to remove partition penalties.\n')
    }
    names(partition_penalties) <- paste0('partition', 1:(K+1))
    penalty_vec <- c(penalty_vec, partition_penalties)
  }

  if(verbose){
    cat("Parallel and Weighting Setup\n")
  }

  ## Export components for parallel processing
  if(parallel && !is.null(cl)) {

    shared_env <- new.env(parent = emptyenv())

    shared_vars <- list(
      A = A,
      nca = ncol(A),
      K = K,
      nc = nc,
      nr = nr,
      chunk_size = chunk_size,
      num_chunks = num_chunks,
      rem_chunks = rem_chunks,
      penalty_vec = penalty_vec,
      unique_penalty_per_partition = unique_penalty_per_partition,
      keep_weighted_Lambda = keep_weighted_Lambda,
      custom_penalty_mat = custom_penalty_mat,
      glm_weight_function = glm_weight_function,
      schur_correction_function = schur_correction_function,
      unconstrained_fit_fxn = unconstrained_fit_fxn,
      observation_weights = observation_weights
    )

    for(nm in names(shared_vars)) {
      assign(nm, get(nm, envir = environment()), envir = shared_env)
    }

    tryCatch({
      parallel::clusterExport(cl, "shared_env", envir = environment())
    }, error = function(e) {
      stop("Failed to export 'shared_env' to cluster: ",
           e$message, call. = FALSE)
    })

    parallel::clusterExport(cl, "efficient_matrix_mult", envir = environment())

    parallel::clusterEvalQ(cl, {
      `%**%` <- efficient_matrix_mult
      if (exists("shared_env")) {
        list2env(as.list(shared_env), envir = .GlobalEnv)
      } else {
        warning("shared_env not found on worker node.")
      }
    })

  }

  ## X^{\top}WX
  if(( (paste0(family)[1] == 'gaussian' &
        paste0(family)[2] == 'identity')) &
     !homogenous_weights){
    X <- lapply(1:(K+1), function(k){
      X[[k]] * c(sqrt(observation_weights[[k]]))
    })
  }
  X_gram <- compute_gram_block_diagonal(X,
                                        parallel & parallel_matmult,
                                        cl,
                                        chunk_size,
                                        num_chunks,
                                        rem_chunks)
  if(( (paste0(family)[1] == 'gaussian' &
        paste0(family)[2] == 'identity')) &
     !homogenous_weights){
    X <- lapply(1:(K+1), function(k){
      X[[k]] / c(sqrt(observation_weights[[k]]))
    })
  }

  if(verbose){
    cat("SQP Setup\n")
  }

  ## Update quadprog variable, if the correct arguments are made
  if(qp_negative_derivative | qp_monotonic_decrease | qp_negative_2ndderivative |
     qp_positive_derivative | qp_monotonic_increase | qp_positive_2ndderivative |
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

  ## Quadprog setup
  if(quadprog){

    ## Save results, when multiple inequality constraints present
    #  collect as the script moves through
    qp_Amat_list <- list()
    qp_bvec_list <- list()
    qp_meq_list <- list()

    ## Block-diagonal collapse
    X_block <- Reduce("rbind", lapply(1:(K+1), function(k){
      dummy <- 0*X[[k]]
      Reduce("cbind",lapply(1:(K+1),function(j){
        if(nrow(X[[k]]) == 0){
          return(X[[k]])
        } else if(j == k) X[[k]] else 0*X[[k]]
      }))
    }))

    if(any(!is.null(qp_observations))){
      qp_observations <- try(c(qp_observations), silent = TRUE)
      if(any(inherits(qp_observations, 'try-error'))){
        stop('\n \t qp_observations must coercible to a numeric vector \n')
      }
      X_block <- X_block[c(unlist(order_list)) %in% qp_observations,,drop=FALSE]
    }

    ## Range constraints
    if(!(any(is.null(qp_range_upper))) | !any(is.null(qp_range_lower))){
      if(paste0(family)[2] != 'identity'){
        if(any(!is.null(qp_range_upper))){
          qp_range_upper <- family$linkfun(qp_range_upper)
        }
        if(any(!is.null(qp_range_lower))){
          qp_range_lower <- family$linkfun(qp_range_lower)
        }
        if(any(is.na(c(qp_range_upper, qp_range_lower))) |
           any(!is.finite(c(qp_range_upper,
                            qp_range_lower)))){
          stop('\n\tQuadratic programming upper/lower range constraints ',
               'are NA or not finite after link function transformation. ',
               'Supply bounds on raw-response scale, which are in-range of ',
               'link function transformations.\n')
        }
      }

      if(!(any(is.null(qp_range_upper))) & !any(is.null(qp_range_lower))){
        qp_Amat <- cbind(t(X_block), -t(X_block))
        if(length(qp_range_lower) == 1){
          if(length(qp_range_upper == 1)){
            qp_Amat <- t(unique(t(qp_Amat)))
          }
          qp_bvec <- rep(qp_range_lower, ncol(qp_Amat)/2)
        } else {
          qp_bvec <- qp_range_lower
        }
        qp_bvec_lower <- (qp_bvec - mean_y) / sd_y
        if(length(qp_range_upper) == 1){
          qp_bvec <- rep(qp_range_upper, ncol(qp_Amat)/2)
        } else {
          qp_bvec <- qp_range_upper
        }

        ## Save results, if multiple constraints present
        qp_bvec_upper <- -(qp_bvec - mean_y) / sd_y
        qp_bvec <- c(qp_bvec_lower, qp_bvec_upper)
        qp_Amat_list[[length(qp_Amat_list) + 1]] <- qp_Amat
        qp_bvec_list[[length(qp_bvec_list) + 1]] <- qp_bvec
        qp_meq_list[[length(qp_meq_list) + 1]] <- 0

      } else if(!(any(is.null(qp_range_upper)))){
        qp_Amat <- -t(X_block)
        if(length(qp_range_upper) == 1){
          qp_Amat <- t(unique(t(qp_Amat)))
          qp_bvec <- rep(qp_range_upper, ncol(qp_Amat))
        } else {
          qp_bvec <- qp_range_upper
        }
        qp_bvec <- -(qp_bvec - mean_y) / sd_y

        ## Keep appending the constraints
        qp_Amat_list[[length(qp_Amat_list) + 1]] <- qp_Amat
        qp_bvec_list[[length(qp_bvec_list) + 1]] <- qp_bvec
        qp_meq_list[[length(qp_meq_list) + 1]] <- 0

      } else if(!any(is.null(qp_range_lower))){
        qp_Amat <- t(X_block)
        if(length(qp_range_lower) == 1){
          qp_Amat <- t(unique(t(qp_Amat)))
          qp_bvec <- rep(qp_range_lower, ncol(qp_Amat))
        } else {
          qp_bvec <- qp_range_lower
        }
        qp_bvec <- (qp_bvec - mean_y) / sd_y
        qp_Amat_list[[length(qp_Amat_list) + 1]] <- qp_Amat
        qp_bvec_list[[length(qp_bvec_list) + 1]] <- qp_bvec
        qp_meq_list[[length(qp_meq_list) + 1]] <- 0
      }
    }

    ## Positive first derivative
    if(qp_positive_derivative){
      derivs <- make_derivative_matrix(
        nc, X_block, power1_cols, power2_cols, nonspline_cols,
        interaction_single_cols, interaction_quad_cols, triplet_cols,
        K, include_2way_interactions, include_3way_interactions,
        include_quadratic_interactions, colnm_expansions, expansion_scales,
        just_first_derivatives = TRUE, just_spline_effects = FALSE)
      first_derivative_constraints <- Reduce(
        "rbind",
        lapply(
          derivs$first_derivative, function(deriv_matrix) {
            t(Reduce('rbind', lapply(1:nrow(deriv_matrix), function(i) {
              matrix(c(deriv_matrix[i,]), nrow = 1)
            })))
          }
        )
      )
      qp_Amat <- cbind(first_derivative_constraints)
      qp_Amat <- t(unique(t(qp_Amat)))
      qp_bvec <- rep(0, ncol(qp_Amat))
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- qp_Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- qp_bvec
      qp_meq_list[[length(qp_meq_list) + 1]] <- 0

      ## Negative first derivative
    } else if(qp_negative_derivative){
      derivs <- make_derivative_matrix(
        nc, X_block, power1_cols, power2_cols, nonspline_cols,
        interaction_single_cols, interaction_quad_cols, triplet_cols,
        K, include_2way_interactions, include_3way_interactions,
        include_quadratic_interactions, colnm_expansions, expansion_scales,
        just_first_derivatives = TRUE, just_spline_effects = FALSE)
      first_derivative_constraints <- Reduce(
        "rbind",
        lapply(derivs$first_derivative, function(deriv_matrix) {
          t(Reduce('rbind', lapply(1:nrow(deriv_matrix), function(i) {
            -matrix(c(deriv_matrix[i,]), nrow = 1)
          })))
        }
        )
      )
      qp_Amat <- cbind(first_derivative_constraints)
      qp_Amat <- t(unique(t(qp_Amat)))
      qp_bvec <- rep(0, ncol(qp_Amat))
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- qp_Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- qp_bvec
      qp_meq_list[[length(qp_meq_list) + 1]] <- 0
    }

    ## Positive second derivative
    if(qp_positive_2ndderivative){
      derivs <- make_derivative_matrix(
        nc, X_block, power1_cols, power2_cols, nonspline_cols,
        interaction_single_cols, interaction_quad_cols, triplet_cols,
        K, include_2way_interactions, include_3way_interactions,
        include_quadratic_interactions, colnm_expansions, expansion_scales,
        just_first_derivatives = FALSE, just_spline_effects = FALSE)
      second_derivative_constraints <- Reduce(
        "rbind",
        lapply(
          derivs$second_derivative, function(deriv_matrix) {
            t(Reduce('rbind', lapply(1:nrow(deriv_matrix), function(i) {
              matrix(c(deriv_matrix[i,]), nrow = 1)
            })))
          }
        )
      )
      qp_Amat <- cbind(second_derivative_constraints)
      qp_Amat <- t(unique(t(qp_Amat)))
      qp_bvec <- rep(0, ncol(qp_Amat))
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- qp_Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- qp_bvec
      qp_meq_list[[length(qp_meq_list) + 1]] <- 0

      ## Negative second derivative
    } else if(qp_negative_2ndderivative){
      derivs <- make_derivative_matrix(
        nc, X_block, power1_cols, power2_cols, nonspline_cols,
        interaction_single_cols, interaction_quad_cols, triplet_cols,
        K, include_2way_interactions, include_3way_interactions,
        include_quadratic_interactions, colnm_expansions, expansion_scales,
        just_first_derivatives = FALSE, just_spline_effects = FALSE)
      second_derivative_constraints <- Reduce(
        "rbind",
        lapply(derivs$second_derivative, function(deriv_matrix) {
          t(Reduce('rbind', lapply(1:nrow(deriv_matrix), function(i) {
            -matrix(c(deriv_matrix[i,]), nrow = 1)
          })))
        }
        )
      )
      qp_Amat <- cbind(second_derivative_constraints)
      qp_Amat <- t(unique(t(qp_Amat)))
      qp_bvec <- rep(0, ncol(qp_Amat))
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- qp_Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- qp_bvec
      qp_meq_list[[length(qp_meq_list) + 1]] <- 0
    }

    ## Always increasing with respect order of observations in data
    if(qp_monotonic_increase){
      value_constraints <- t(Reduce(
        'rbind',
        lapply(2:nrow(X_block),
               function(i) {
                 matrix(c(X_block[i,] - X_block[i-1,]), nrow = 1)
               })
      )
      )
      qp_Amat <- cbind(value_constraints)
      qp_Amat <- t(unique(t(qp_Amat)))
      qp_bvec <- rep(0, ncol(qp_Amat))
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- qp_Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- qp_bvec
      qp_meq_list[[length(qp_meq_list) + 1]] <- 0

      ## Always decreasing with respect to order of observations in data
    } else if(qp_monotonic_decrease){
      value_constraints <- -t(Reduce(
        'rbind',
        lapply(2:nrow(X_block),
               function(i) {
                 matrix(c(X_block[i,] - X_block[i-1,]), nrow = 1)
               }
        )
      )
      )
      qp_Amat <- cbind(value_constraints)
      qp_Amat <- t(unique(t(qp_Amat)))
      qp_bvec <- rep(0, ncol(qp_Amat))
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- qp_Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- qp_bvec
      qp_meq_list[[length(qp_meq_list) + 1]] <- 0
    }

    ## Bind all inequality constraints together
    if(!is.null(qp_Amat_fxn) &
       !is.null(qp_bvec_fxn) &
       !is.null(qp_meq_fxn)){
      qp_Amat <- qp_Amat_fxn(nr, nc, K, X_block, colnm_expansions,
                             expansion_scales, all_derivatives, ...)
      qp_bvec <- qp_bvec_fxn(qp_Amat, nr, nc, K, X_block, colnm_expansions,
                             expansion_scales, all_derivatives, ...)
      qp_meq <- qp_meq_fxn(qp_Amat, nr, nc, K, X_block, colnm_expansions,
                           expansion_scales, all_derivatives, ...)
      qp_Amat_list[[length(qp_Amat_list) + 1]] <- qp_Amat
      qp_bvec_list[[length(qp_bvec_list) + 1]] <- qp_bvec
      qp_meq_list[[length(qp_meq_list) + 1]] <- qp_meq
    }

    qp_Amat <- do.call(cbind, qp_Amat_list)
    qp_bvec <- do.call(c, qp_bvec_list)
    qp_meq <- sum(unlist(qp_meq_list))

    X_block <- NULL
  }

  ## Convert constraint_values to per-partition list
  #  Convert non-0 null vectors to (K+1) list of corresponding partitions
  if(length(constraint_values) > 0){

    constraint_values <- lapply(1:(K+1),function(k){
      vec <- cbind(constraint_values)[1:nc + (k-1)*nc,,drop=FALSE]
      vec[1,] <- (vec[1,] - mean_y)
      vec * c(1, expansion_scales) / sd_y
    })

  }

  if(verbose){
    cat("Tune Smoothing Spline Penalty\n")
  }

  ## Weighting for Gaussian response, then tune_Lambda call
  #  This is to incorporate weights efficiently for linear regression outcomes
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
  if(!(!any(is.null(previously_tuned_penalties))) & !dummy_fit){
    tL <- try({
      tune_Lambda(
        y, X, X_gram, smoothing_spline_penalty, A, K, nc, nr,
        opt, use_custom_bfgs, C, colnm_expansions,
        wiggle_penalty, flat_ridge_penalty,
        initial_wiggle, initial_flat,
        unique_penalty_per_predictor, unique_penalty_per_partition,
        penalty_vec, meta_penalty, family, unconstrained_fit_fxn,
        keep_weighted_Lambda, iterate_tune, qp_score_function,
        quadprog, qp_Amat, qp_bvec, qp_meq, tol, sd_y, delta,
        constraint_values,
        parallel, parallel_eigen, parallel_trace, parallel_aga,
        parallel_matmult, parallel_unconstrained,
        cl, chunk_size, num_chunks, rem_chunks, shared_env,
        custom_penalty_mat, order_list,
        glm_weight_function, schur_correction_function,
        need_dispersion_for_estimation, dispersion_function,
        observation_weights, homogenous_weights,
        blockfit, just_linear_without_interactions,
        Vhalf, VhalfInv,
        verbose_tune, include_warnings, ...)}, silent = TRUE)
    if(inherits(tL, 'try-error')){
      if(include_warnings) print(tL)
      return(tL)
    }
  } else {
    ## Penalties already tuned and supplied by the user
    tL <- previously_tuned_penalties
    previously_tuned_penalties <- NULL
  }
  flat_ridge_penalty <- tL$flat_ridge_penalty
  wiggle_penalty <- tL$wiggle_penalty

  ##  [Change 2026-02-16] dummy_fit early return
  if(dummy_fit){
    if(verbose) cat('Dummy Fit Early Return\n')

    if(length(c(numerics, nonspline_cols)) > 1 & K > 0){
      assign_partition <- partitions$assign_partition
    } else if(length(numerics) == 1 & K > 0){
      assign_partition <- function(x)rowMeans(cbind(x))
    } else {
      assign_partition <- function(x)0.5
    }

    B <- lapply(1:(K+1), function(k){
      b <- cbind(rep(0, nc))
      rownames(b) <- colnm_expansions
      names(b) <- colnm_expansions
      b
    })
    names(B) <- paste0('partition', 1:(K+1))
    B_raw <- B
    X_out <- lapply(X, unstd_X)

    if(K == 0){
      knots <- NULL
    } else {
      knots <- inv_transf(knot_values)
      if(length(numerics) == 1 & length(nonspline_cols) == 0){
        rownames(knots) <- paste0(1:K, '_', 2:(K+1))
      }
    }

    if(quadprog){
      quadprog_list <- list(
        qp_Amat = qp_Amat,
        qp_bvec = qp_bvec,
        qp_meq = qp_meq
      )
    } else {
      quadprog_list <- list(NA)
    }

    return(list(
      y = y_og, ytilde = rep(0, nr), X = X_out, A = A, B = B,
      B_raw = B_raw, K = K, p = nc, q = ncol(predictors),
      P = (K+1)*nc, N = nr, penalties = tL,
      knot_scale_transf = transf, knot_scale_inv_transf = inv_transf,
      knots = knots, partition_codes = partition_codes,
      partition_bounds = partition_bounds,
      knot_expand_function = knot_expand_list,
      assign_partition = assign_partition,
      make_partition_list = partitions,
      order_list = order_list, og_order = og_order,
      expansion_scales = expansion_scales,
      raw_expansion_names = colnm_expansions,
      family = family, mean_y = mean_y, sd_y = sd_y,
      numerics = numerics, power1_cols = power1_cols,
      power2_cols = power2_cols, power3_cols = power3_cols,
      power4_cols = power4_cols, quad_cols = quad_cols,
      interaction_single_cols = interaction_single_cols,
      interaction_quad_cols = interaction_quad_cols,
      triplet_cols = triplet_cols, nonspline_cols = nonspline_cols,
      constraint_values = constraint_values,
      constraint_vectors = constraint_vectors,
      backtransform_coefficients = backtransform_coefficients,
      forwtransform_coefficients = forwtransform_coefficients,
      std_X = std_X, unstd_X = unstd_X,
      weights = observation_weights_og,
      quadprog_list = quadprog_list,
      G = NULL, Ghalf = NULL, U = NULL,
      sigmasq_tilde = NA_real_, trace_XUGX = NA_real_,
      varcovmat = NULL, VhalfInv = VhalfInv
    ))
  }

  ## This replaces the original get_B call with conditional cal
  #  to blockfit_solve or get_B depending on model configuration
  if(verbose){
    cat("Prep for final fitting\n")
  }

  ## Final fit
  if(K == 0){
    ## ensuring compatibility with no A
    if(any(is.null(A))){
      A <- cbind(rep(0, (K+1)*nc))
      A <- cbind(A, A)
      nca <- 2
    }
  }
  Xy <- vectorproduct_block_diagonal(X, y, K)
  schur_corrections <- lapply(1:(K+1), function(k)0)
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
                            schur_corrections)

  if(verbose){
    cat('Last fit\n')
  }

  ## Determine which columns of X are "flat" (non-interactive linear terms
  #  that should share a single coefficient across all partitions).
  #  These correspond to just_linear_without_interactions in the expansion.
  flat_cols <- c()
  if(length(just_linear_without_interactions) > 0){
    flat_cols <- which(colnm_expansions %in%
                         paste0("_", just_linear_without_interactions, "_"))
    ## Exclude any that ended up in other expansion categories
    flat_cols <- flat_cols[!(flat_cols %in% c(power1_cols,
                                              interaction_single_cols,
                                              interaction_quad_cols,
                                              triplet_cols))]
  }

  ## [Change 2026-02-16] Use backfitting for blockfit,
  #  Divide the design matrix into "blocks" for smooth + interaction terms and
  #  linear terms
  #  otherwise call get_B directly.
  has_qp_ineq <- quadprog

  use_blockfit <- blockfit &&
    length(flat_cols) > 0 &&
    K > 0

  ## Get coefficient and correlation matrix estimates
  return_G_getB <- TRUE

  ## If using the blockfit, backfitting routine
  if(use_blockfit){

    if(verbose){
      cat("Blockfit backfit\n")
    }

    B_list <- try({
      blockfit_solve(
        X = X,
        y = y,
        flat_cols = flat_cols,
        K = K,
        nc = nc,
        Lambda = tL$Lambda,
        L_partition_list = tL$L_partition_list,
        unique_penalty_per_partition = unique_penalty_per_partition,
        A = A,
        nca = nca,
        constraint_values = constraint_values,
        X_gram = X_gram,
        Ghalf_full = G_list$Ghalf,
        GhalfInv_full = G_list$GhalfInv,
        family = family,
        order_list = order_list,
        glm_weight_function = glm_weight_function,
        schur_correction_function = schur_correction_function,
        need_dispersion_for_estimation = need_dispersion_for_estimation,
        dispersion_function = dispersion_function,
        observation_weights = observation_weights,
        homogenous_weights  = homogenous_weights,
        iterate = iterate_final_fit,
        tol = tol,
        parallel_eigen = parallel & parallel_eigen,
        cl = cl,
        chunk_size = chunk_size,
        num_chunks = num_chunks,
        rem_chunks = rem_chunks,
        return_G_getB = return_G_getB,
        quadprog = quadprog,
        qp_Amat = qp_Amat,
        qp_bvec = qp_bvec,
        qp_meq = qp_meq,
        qp_score_function = qp_score_function,
        keep_weighted_Lambda = keep_weighted_Lambda,
        max_backfit_iter = 100,
        Vhalf = Vhalf,
        VhalfInv = VhalfInv,
        include_warnings = include_warnings,
        verbose         = verbose,
        ...
      )
    }, silent = TRUE
    )

    if(any(inherits(B_list, 'try-error'))){
      if(include_warnings){
        print(B_list)
        warning('\n \t blockfit_solve failed, falling back to get_B \n')
      }
      ## Fall back to standard get_B if blockfit_solve fails
      use_blockfit <- FALSE
    }
  }

  ## Standard path: call get_B (includes when blockfit fallback triggered)
  if(!use_blockfit){
    B_list <-
      try({get_B(
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
        constraint_values,
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
        ...)}, silent = TRUE)
  }

  if(any(inherits(B_list, 'try-error'))){
    if(include_warnings) print(B_list)
    stop('\n \t Failure in fitting final model \n')
  }
  B <- B_list$B
  G_list <- B_list$G_list

  ## [Change 2026-02-16] The post-fit code for computing and formatting terms
  #  (scale backtransform, weights, backtransform X/y after VhalfInv transforms,
  #  assign_partition, B_raw, predict_function, plotting, U, trace, dispersion,
  #  varcov, Lagrange multipliers, return_list)

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
  # with VhalfInv NOT involved
  if(!is.null(VhalfInv)){
    X <- X_expand_og
    y <- y_expand_og
    X_expand_og <- NULL
    y_expand_og <- NULL
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
  B <- lapply(B, function(b)b * sd_y)

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
    names(B[[k]]) <- colnm_expansions
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
                               expansions_only = FALSE){

    ## [Change 2026-02-22] SINGLE merged encoding block with guard.
    #  Previous code had two duplicate encoding blocks, both of which
    #  could fire on plain numeric matrices that acquired spurious
    #  colnames from cbind() in predict.lgspline. When colnames like
    #  "predictors_val" didn't match og_cols, the block would:
    #   (a) fail to find any one-hot prefix match,
    #    (b) add zero-filled columns for every "missing" og_col,
    #    (c) select only the zero columns, discarding actual data.
    #
    ## Fix: only perform zero-fill + reorder when at least one column
    #  genuinely needs one-hot encoding (i.e., its name is a prefix of
    #  some og_col). Otherwise strip the spurious colnames and proceed.
    if(any(!is.null(new_predictors))){

      if(!is.null(og_cols) &&
         (inherits(new_predictors, "data.frame") ||
          (is.matrix(new_predictors) &&
           !is.null(colnames(new_predictors))))){

        np_df <- as.data.frame(new_predictors)
        col_names <- colnames(np_df)

        ## Determine if any column actually requires one-hot expansion:
        ## a column needs encoding only when (a) it is NOT already in
        ## og_cols, AND (b) og_cols contains entries matching its
        ## one-hot prefix pattern (e.g., "Subject" -> "Subject_1", ...).
        cols_needing_encoding <- character(0)
        for(cn in col_names){
          if(cn %in% og_cols) next
          onehot_pattern <- paste0("^", cn, "_")
          if(any(grepl(onehot_pattern, og_cols))){
            cols_needing_encoding <- c(cols_needing_encoding, cn)
          }
        }

        if(length(cols_needing_encoding) > 0){
          ## At least one column genuinely needs one-hot expansion.
          for(col_name in cols_needing_encoding){
            onehot_pattern <- paste0("^", col_name, "_")
            candidate_og   <- og_cols[grepl(onehot_pattern, og_cols)]
            prefix    <- paste0(col_name, "_")
            og_levels <- substr(candidate_og,
                                nchar(prefix) + 1L,
                                nchar(candidate_og))
            user_vals <- as.character(np_df[[col_name]])

            onehot_mat <- matrix(0L,
                                 nrow = nrow(np_df),
                                 ncol = length(og_levels))
            colnames(onehot_mat) <- candidate_og
            for(lv_idx in seq_along(og_levels)){
              onehot_mat[user_vals == og_levels[lv_idx], lv_idx] <- 1L
            }

            unrecognized <- unique(user_vals[!user_vals %in% og_levels])
            if(length(unrecognized) > 0 && include_warnings){
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

          ## Add zero columns for any og_cols not yet present
          ## (e.g., one-hot levels that were in training but not
          ## represented in new_predictors)
          missing_cols <- setdiff(og_cols, colnames(np_df))
          if(length(missing_cols) > 0){
            np_df[missing_cols] <- 0L
          }
          ## Reorder columns to match training order
          if(all(og_cols %in% colnames(np_df))){
            np_df <- np_df[, og_cols, drop = FALSE]
          }

          new_predictors <- as.matrix(np_df)

        } else if(all(col_names %in% og_cols)){
          ## Column names already match og_cols exactly — just
          ## reorder to match training column order, then strip names.
          np_df <- np_df[, og_cols, drop = FALSE]
          new_predictors <- as.matrix(np_df)

        } else {
          ## Column names don't match og_cols and don't match any
          ## one-hot prefix. This is a spurious colname (e.g., from
          ## cbind() in predict.lgspline). Strip names and proceed
          ## with the raw numeric matrix as-is.
          ## (Do NOT zero-fill or reorder.)
        }
      }

      ## Final coercion to plain numeric matrix
      new_predictors <- try(
        methods::as(cbind(new_predictors), 'matrix'),
        silent = TRUE)
      if(any(inherits(new_predictors, 'try-error'))){
        stop('\n \t New predictors should be able to be coerced into ',
             'matrix form. If you are passing factor/character columns, ',
             'make sure the column names match the original variable ',
             'names used at fit time. \n')
      }
      ## Strip colnames so downstream code uses positional indexing
      ## only — prevents any further name-based misinterpretation.
      colnames(new_predictors) <- NULL
    }

    ## Avoid R rbind issue with 1 row only for certain internal functions
    if(nrow(new_predictors) == 1){
      new_predictors <- rbind(new_predictors, new_predictors)
      only_1 <- TRUE
    } else {
      only_1 <- FALSE
    }

    ## Accept predictors as matrix
    new_predictors <- transf(methods::as(new_predictors, 'matrix'))

    ## Needed for determining knot locations
    partition_codes_new <- assign_partition(new_predictors)

    ## Back transform, now that knots have been established
    new_predictors <- inv_transf(new_predictors)

    ## Cubic/polynomial expansions
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
    if(expansions_only){
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
        nc, C_new, power1_cols, power2_cols, nonspline_cols,
        interaction_single_cols, interaction_quad_cols, triplet_cols,
        K, include_2way_interactions, include_3way_interactions,
        include_quadratic_interactions, colnm_expansions, expansion_scales,
        !take_second_derivatives)

      if(only_1){
        partition_codes_new <- partition_codes_new[1]
      }

      ## Account for derivatives of link transform
      if (is.null(family$linkinvderiv) || is.null(family$linkinvderiv2)) {
        if (family$link == 'inverse') {
          family$linkinvderiv <- function(mu) 1 / mu
          family$linkinvderiv2 <- function(mu) 2 / mu^3
        } else if (family$link == 'logit') {
          family$linkinvderiv <- function(mu) mu * (1 - mu)
          family$linkinvderiv2 <- function(mu) mu * (1 - mu) * (1 - 2 * mu)
        } else if (family$link == 'log') {
          family$linkinvderiv <- function(mu) mu
          family$linkinvderiv2 <- function(mu) mu
        } else if (family$link == 'identity') {
          family$linkinvderiv <- function(mu) 1
          family$linkinvderiv2 <- function(mu) 0
        } else if (family$link == 'probit') {
          family$linkinvderiv <- function(mu) dnorm(qnorm(mu))
          family$linkinvderiv2 <- function(mu) -dnorm(qnorm(mu)) * qnorm(mu)
        } else if (family$link == 'sqrt') {
          family$linkinvderiv <- function(mu) 2 * sqrt(mu)
          family$linkinvderiv2 <- function(mu) mu^(-1/2)
        } else if (family$link == 'inverse.sqrt') {
          family$linkinvderiv <- function(mu) 2 * mu^(3/2)
          family$linkinvderiv2 <- function(mu) 3 * mu^(1/2)
        } else if (family$link == 'cloglog') {
          family$linkinvderiv <- function(mu) -1 / log(1 - mu) * (1 - mu)
          family$linkinvderiv2 <- function(mu) 2 / (log(1 - mu))^2 * (1 - mu)
        } else if (family$link == 'cauchit') {
          family$linkinvderiv <- function(mu) pi * (1 + (qcauchy(mu))^2)
          family$linkinvderiv2 <- function(mu) -2 * pi *
            qcauchy(mu) * dcauchy(qcauchy(mu))
        } else if (family$link == 'log1p') {
          family$linkinvderiv <- function(mu) 1 / (1 + mu)
          family$linkinvderiv2 <- function(mu) -1 / (1 + mu)^2
        } else {
          if (include_warnings) {
            warning(
              '\n\t',
              'Link function not recognized: supply a custom "linkinvderiv"',
              ' function to your custom "family" object to properly compute ',
              'derivatives accounting for link function transforms for GLMs.',
              ' The derivatives will be returned on the link-transformed ',
              'scale.\n'
            )
          }
          family$linkinvderiv <- function(mu) 1
          family$linkinvderiv2 <- function(mu) 0
        }
      }

      n_deriv_vars <- length(derivs$first_derivative)
      deriv_names <- names(derivs$first_derivative)
      if(is.null(deriv_names)){
        deriv_names <- paste0("var_", seq_len(n_deriv_vars))
      }

      Cprime_new <- Reduce("rbind",
                           lapply(1:n_deriv_vars,
                                  function(var){
                                    d <- derivs$first_derivative[[var]]
                                    if(only_1) return(d[1,,drop=FALSE])
                                    else return(d)
                                  }))

      Xprime_new <- knot_expand_list(
        partition_codes_new, partition_bounds,
        length(partition_codes_new), Cprime_new, K)

      ## [MD 2026-02-22] Reorder preds_prime to match final_preds order.
      #  Previously missing; causes wrong element-wise pairing with
      #  linkinvderiv(final_preds) for non-identity links when
      # partition order differs from observation order.
      preds_prime <-
        unlist(
          matmult_block_diagonal(
            Xprime_new[keep_blocks], B_predict[keep_blocks],
            length(keep_blocks) - 1, parallel, cl,
            chunk_size, num_chunks, rem_chunks))[order(unlist(order_list))]

      if(n_deriv_vars == 1){
        final_preds_prime <- family$linkinvderiv(final_preds) * preds_prime
      } else {
        n_obs <- if(only_1) 1L else length(partition_codes_new)
        final_preds_prime <- lapply(seq_len(n_deriv_vars), function(v){
          rows <- ((v - 1) * n_obs + 1):(v * n_obs)
          family$linkinvderiv(final_preds) * preds_prime[rows]
        })
        names(final_preds_prime) <- deriv_names
      }

      if(take_second_derivatives){
        Cdprime_new <- Reduce("rbind",
                              lapply(1:length(derivs$second_derivative),
                                     function(var){
                                       d <- derivs$second_derivative[[var]]
                                       if(only_1) return(d[1,,drop=FALSE])
                                       else return(d)
                                     }))

        Xdprime_new <- knot_expand_list(
          partition_codes_new, partition_bounds,
          length(partition_codes_new), Cdprime_new, K)

        ## [Change 2026-02-22, see above]
        #  Same reorder for second derivatives.
        preds_dprime <-
          unlist(
            matmult_block_diagonal(
              Xdprime_new[keep_blocks], B_predict[keep_blocks],
              length(keep_blocks) - 1, parallel, cl,
              chunk_size, num_chunks, rem_chunks))[order(unlist(order_list))]

        if(n_deriv_vars == 1){
          final_preds_dprime <-
            family$linkinvderiv2(final_preds)*preds_prime^2 +
            family$linkinvderiv(final_preds)*preds_dprime
        } else {
          n_obs <- if(only_1) 1L else length(partition_codes_new)
          final_preds_dprime <- lapply(seq_len(n_deriv_vars), function(v){
            rows <- ((v - 1) * n_obs + 1):(v * n_obs)
            pp <- preds_prime[rows]
            dp <- preds_dprime[rows]
            family$linkinvderiv2(final_preds) * pp^2 +
              family$linkinvderiv(final_preds) * dp
          })
          names(final_preds_dprime) <- deriv_names
        }
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

  ## For compatibility without knots
  if(K == 0){
    knots <- NULL
  } else if(qcols == 1) {
    knots <- inv_transf(knot_values)
  }

  ## For saving quadratic programming components
  if(quadprog){
    quadprog_list <- list(
      qp_Amat = qp_Amat,
      qp_bvec = qp_bvec,
      qp_meq = qp_meq,
      qp_info = B_list$qp_info
    )
  } else {
    quadprog_list <- list(NA)
  }
  qp_Amat <- NULL
  qp_bvec <- NULL
  qp_meq <- NULL

  ## return_list construction,optional components
  #  (G, Ghalf, U, Lagrange multipliers, dispersion, trace, varcov),
  #  and final return.')

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
                      "unbias_dispersion" = unbias_dispersion,
                      "backtransform_coefficients" = backtransform_coefficients,
                      "forwtransform_coefficients" = forwtransform_coefficients,
                      "mean_y" = mean_y,
                      "sd_y" = sd_y,
                      "og_order" = og_order,
                      "order_list" = order_list,
                      "constraint_values" = constraint_values,
                      "constraint_vectors" = constraint_vectors,
                      "make_partition_list" = partitions,
                      "expansion_scales" = expansion_scales,
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
                      "raw_expansion_names" = colnm_expansions,
                      "std_X" = std_X,
                      "unstd_X" = unstd_X,
                      "parallel_cluster_supplied" = parallel,
                      "weights" = observation_weights_og,
                      "VhalfInv" = VhalfInv,
                      "quadprog_list" = quadprog_list)

  if(verbose){
    cat("Optional Components\n")
  }

  ## We need U and sigma^2 to compute sigma^2*UG
  if(return_varcovmat){
    return_U <- TRUE
    estimate_dispersion <- TRUE
  }

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

    ## When quadprog was used, form the effective constraint matrix from
    #  active constraints (equalities + binding inequalities) stored in
    #  qp_info$Amat_active. Fall back to A if qp_info is unavailable.
    has_active_qp <- quadprog &&
      !is.null(B_list$qp_info) &&
      !is.null(B_list$qp_info$Amat_active) &&
      ncol(B_list$qp_info$Amat_active) > 0

    if(has_active_qp){
      ## Amat_active is already in the block-coefficient space (P rows),
      ## one column per active constraint (equality or binding inequality).
      A_for_U   <- B_list$qp_info$Amat_active
      nca_for_U <- ncol(A_for_U)
    } else {
      A_for_U   <- A
      nca_for_U <- nca
    }

    if(K == 0 & length(constraint_values) == 0 & !quadprog){
      return_list$U <- diag(nc*(K+1))
      if(any(is.null(A))){
        A <- cbind(rep(0, (K+1)*nc))
        A <- cbind(A, A)
        nca <- 2
      }
    } else if(is.null(return_list$VhalfInv)){
      ## Standard block-diagonal path: use active constraint matrix
      #  so U projects onto the null space of all binding constraints.
      return_list$U <- get_U(
        G_list$G,
        A_for_U,
        K,
        nc,
        nca_for_U
      )
    } else {
      ## GEE / VhalfInv path: build full whitened G_correct then project
      #  using the active constraint matrix.
      X_full_tr <- collapse_block_diagonal(X)

      ## W~ without obs weights to avoid double-counting D
      prelim_disp_tr <- if(need_dispersion_for_estimation){
        dispersion_function(mu = ytilde,
                            y = y_og,
                            order_indices = 1:nr,
                            family = family,
                            observation_weights = observation_weights_og,
                            VhalfInv = VhalfInv,
                            ...)
      } else {
        1
      }
      W_glm_tr <- c(glm_weight_function(
        ytilde[unlist(order_list)], y_og[unlist(order_list)],
        1:nr, family, prelim_disp_tr, rep(1, nr), ...))
      W_glm_tr <- pmax(W_glm_tr, .Machine$double.eps)

      ## Observation weights in partition order
      D_tr <- observation_weights_og[unlist(order_list)]

      ## Weight X by sqrt(W~ * D)
      X_full_tr <- t(t(X_full_tr) * sqrt(W_glm_tr * D_tr))

      VinvhalfX_tr <- VhalfInv[unlist(order_list),
                               unlist(order_list)] %**%
        X_full_tr

      has_part_pen_tr <- length(tL$L_partition_list) == (K + 1)
      Lambda_full_tr <- collapse_block_diagonal(
        lapply(1:(K + 1), function(k){
          if(has_part_pen_tr){
            tL$Lambda + tL$L_partition_list[[k]]
          } else {
            tL$Lambda
          }
        })
      )

      ## G_correct = (X^T W~ D V^{-1} X + Lambda)^{-1}
      #  G_correct is the "correct" G because it isn't just block diagonal
      G_correct_tr <- invert(crossprod(VinvhalfX_tr) + Lambda_full_tr)

      ## Recompute U with the correct dense G and active constraint matrix.
      #  When quadprog is active, A_for_U includes columns for binding
      #  inequality constraints so the projection accounts for all active
      #  constraints, not just the smoothness equalities.
      GA <- G_correct_tr %**% A_for_U
      return_list$U <- diag(nrow(G_correct_tr)) -
        GA %**% tcrossprod(invert(crossprod(A_for_U, GA)), A_for_U)
    }
  }

  ## [Change 2026-02-16] Optionally return Lagrangian multipliers
  if(return_lagrange_multipliers && !is.null(A) && K > 0){
    if(verbose) cat("Lagrange Multipliers\n")

    ## When blockfit_solve or get_B ran a QP solve, Lagrange multipliers
    #  are returned directly from solve.QP via qp_info$lagrangian.
    #  These are on the combined constraint space (equalities + inequalities),
    #  one multiplier per column of Amat_combined at convergence.
    #  For the non-QP path, multipliers are recovered analytically from
    #  the unconstrained estimate via (A^T G A)^{-1} A^T (beta_hat - beta_0).
    has_qp_lagrangian <- !is.null(B_list$qp_info) &&
      !is.null(B_list$qp_info$lagrangian)

    if(has_qp_lagrangian){
      ## Multipliers come directly from the last successful solve.QP call.
      #  The full vector covers all columns of Amat_combined; multipliers
      #  for inactive inequality constraints will be at or near zero.
      lagrange_multipliers <- B_list$qp_info$lagrangian

    } else {
      ## Non-QP path: recover equality-constraint multipliers analytically.
      #  Uses the identity lambda = (A^T G A)^{-1} A^T (G X^T y - beta_0)
      #  where G X^T y is the unconstrained penalized estimate.
      Bhat_unc <- unlist(
        matmult_block_diagonal(
          G_list$G,
          lapply(1:(K+1), function(k) Xy[[k]]),
          K,
          parallel = parallel & parallel_matmult,
          cl = cl,
          chunk_size, num_chunks, rem_chunks
        )
      )
      AGAinv <- invert(
        AGAmult_wrapper(
          G_list$G, A, K, nc, nca,
          parallel = parallel & parallel_aga,
          cl = cl,
          chunk_size, num_chunks, rem_chunks
        )
      )
      AtBhat <- crossprod(A, cbind(Bhat_unc))
      if(length(constraint_values) > 0){
        c_vec <- crossprod(A, unlist(constraint_values))
        lagrange_multipliers <- AGAinv %**% cbind(AtBhat - c_vec)
      } else {
        lagrange_multipliers <- AGAinv %**% AtBhat
      }
    }

    return_list$lagrange_multipliers <- lagrange_multipliers

  } else if(return_lagrange_multipliers && (is.null(A) || K == 0)){
    return_list$lagrange_multipliers <- NULL
  }

  ## Estimate sigma^2
  if(estimate_dispersion){
    if(verbose){
      cat("Variance Est \n")
    }

    if(K == 0 && is.null(VhalfInv)){
      trace_XUGX <- sum(unlist(sapply(
        matmult_block_diagonal(
          G_list$G, X_gram, K,
          parallel = FALSE, cl = NULL,
          chunk_size, num_chunks, rem_chunks),
        diag)))
    } else if(is.null(VhalfInv)){
      ## [Change 2026-02-15] use more computationally stable version
      trace_XUGX <- compute_trace_H(
        G_list$G,
        tL$Lambda,
        A,
        invert(AGAmult_wrapper(G_list$G,
                               A, K, nc, nca,
                               parallel = parallel & parallel_aga,
                               cl = cl,
                               chunk_size, num_chunks, rem_chunks)),
        nc, nca, K,
        parallel & parallel_trace,
        cl = cl,
        chunk_size, num_chunks, rem_chunks,
        unique_penalty_per_partition,
        tL$L_partition_list)

    } else if(!is.null(VhalfInv)){
      if(!return_U){
        warning("\n \t return_U is False, cannot compute trace - using P.\n")
        trace_XUGX <- nc*(K+1)
      } else {

        ## tr(U G_correct * X^T W~ D V^{-1} X)
        #  assumes U has been constructed, thus, these components are in memory
        M1 <- return_list$U %**% tcrossprod(G_correct_tr, return_list$U)
        M2 <- crossprod(VinvhalfX_tr)
        trace_XUGX <- sum(M1 * M2)
      }
    }
    if(trace_XUGX < 0 & include_warnings){
      warning('\n \t Trace of hat matrix is < 0, which most often indicates a',
              ' failure of convergence when fitting (i.e. the constrained ',
              'maximum likelihood estimate was not found). Try re-fitting, ',
              'different knot locations, greater penalties, or a less ',
              'complicated model.\n')
    }

    ## Scale factor for unbiased dispersion: N / (N - edf).
    #  When unbias_dispersion = FALSE, scale_by = 1 and the raw mean
    #  squared residual is returned without degrees-of-freedom correction.
    if(unbias_dispersion){
      scale_by <- nr/(nr - trace_XUGX)
    } else {
      scale_by <- 1
    }

    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    ## Dispersion estimation: four sub-paths by family x VhalfInv.
    #  Gaussian identity: closed-form weighted MSE on raw residuals.
    #    With VhalfInv: whiten residuals first so that sigma^2 is on
    #    the decorrelated scale (V^{-1/2}(y - yhat)).
    #  All other families: delegate to dispersion_function, which
    #    handles deviance-based or custom dispersion calculations.
    #    VhalfInv is passed through so the custom function can account for
    #    the correlation structure if needed.
    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    if(paste0(family)[1] == 'gaussian' & paste0(family)[2] == 'identity'){
      if(!is.null(VhalfInv)){
        ## Whiten residuals: r_w = V^{-1/2} (y - yhat)
        #  sigma^2 = mean(D * r_w^2) * scale_by
        resid_w <- c(VhalfInv %**% cbind(y_og - ytilde))
        return_list$sigmasq_tilde <-
          mean(observation_weights_og * resid_w^2) * scale_by
      } else {
        ## Standard weighted MSE on original-scale residuals
        return_list$sigmasq_tilde <-
          mean(observation_weights_og * (y_og - ytilde)^2) * scale_by
      }
    } else {
      if(!is.null(VhalfInv)){
        ## GLM with correlation: pass VhalfInv so dispersion_function
        #  can whiten internally if required by the custom estimator.
        return_list$sigmasq_tilde <- dispersion_function(
          mu = cbind(ytilde),
          y = cbind(y_og),
          order_indices = 1:length(y_og),
          family = family,
          observation_weights = observation_weights_og,
          VhalfInv = VhalfInv,
          ...
        ) * scale_by
      } else {
        ## GLM without correlation: VhalfInv = NULL signals independence.
        return_list$sigmasq_tilde <- dispersion_function(
          mu = ytilde,
          y = y_og,
          order_indices = 1:length(y_og),
          family = family,
          observation_weights = observation_weights_og,
          VhalfInv = NULL,
          ...
        ) * scale_by
      }
    }

    return_list$trace_XUGX <- trace_XUGX

  } else {
    return_list$sigmasq_tilde <- 1
  }

  if(return_varcovmat){
    if(verbose){
      cat("VarCov Mat \n")
    }

    ## [Change 2026-02-17] When VhalfInv is present, the block-diagonal G
    #  from compute_G_eigen ignores cross-partition contributions from
    #  off-diagonal blocks of V^{-1/2}, so G_correct =
    #  (X^T V^{-1} X + Lambda)^{-1} must be computed as a full P x P matrix.
    #  Without VhalfInv the block-diagonal path is exact and unchanged.
    if(!is.null(VhalfInv)){

      ## Full N x P block-diagonal design matrix, original (unwhitened) scale
      X_full <- collapse_block_diagonal(X)[unlist(og_order), , drop = FALSE]

      ## GLM weights W~ at fitted values, without observation
      #  weights (pass rep(1,N)) to avoid double-counting D.
      W_glm_vc <- c(glm_weight_function(ytilde, y_og, 1:nr, family,
                                        return_list$sigmasq_tilde,
                                        rep(1, nr), ...))
      W_glm_vc <- pmax(W_glm_vc, .Machine$double.eps)

      ## Combined weight: sqrt(1/W~ * D) row-wise on X.
      #  For Gaussian identity with homogeneous weights: reduces to X unchanged.
      combined_wt_vc <- sqrt(W_glm_vc * observation_weights_og)
      X_weighted_vc <- t(t(X_full) * combined_wt_vc)

      ## Whitened weighted design: V^{-1/2} (W~ D)^{1/2} X  (N x P)
      VinvhalfX <- VhalfInv %**% X_weighted_vc

      ## Full penalized GLS Gram: X^T W~ D V^{-1} X + Lambda  (P x P)
      Lambda_full <- collapse_block_diagonal(
        lapply(1:(K + 1), function(k){
          if(unique_penalty_per_partition &&
             length(tL$L_partition_list) == (K + 1)){
            tL$Lambda + tL$L_partition_list[[k]]
          } else {
            tL$Lambda
          }
        })
      )
      gram_gls <- crossprod(VinvhalfX) + Lambda_full

      ## G_correct = (X^T W~ D V^{-1} X + Lambda)^{-1}
      Ghalf_correct <- matinvsqrt(gram_gls)

      ## varcovmat = sigma^2 * (U Ghalf_correct)(U Ghalf_correct)^T
      UGhalf <- return_list$U %**% Ghalf_correct
      return_list$varcovmat <-
        tcrossprod(UGhalf) *
        return_list$sigmasq_tilde # contains sd_y^2 already

    } else {

      ## Standard block-diagonal path (exact when no VhalfInv)
      #  [Change 2026-02-14] (UG^{1/2})(UG^{1/2})^{\top} parameterization
      UGhalf <- matmult_U(return_list$U, G_list$Ghalf, nc, K)
      return_list$varcovmat <-
        tcrossprod(UGhalf) *
        return_list$sigmasq_tilde # contains sd_y^2 already
    }

    ## Un-standardize (applies to both paths)
    #  [Change 2026-02-15] corrected "each = K + 1" to "times = K + 1"
    d <- rep(c(1, 1/expansion_scales), times = K + 1)
    return_list$varcovmat <- t(t(return_list$varcovmat * d) * d)

    if(any(diag(return_list$varcovmat) < 0) & include_warnings){
      warning("\n \t Variance-covariance matrix has diagonal elements < 0,",
              " model most likely did not converge when fitting. Try ",
              "re-fitting, a simpler model, changing knot locations, or ",
              "increasing the penalties. \n")
      for(ij in 1:nrow(return_list$varcovmat)){
        return_list$varcovmat[ij, ij] <- max(0,
                                             return_list$varcovmat[ij, ij])
      }
    }
  }


  ## Afterwards, update X to be unstandardized
  return_list$X <- lapply(return_list$X,
                          unstd_X)

  return(return_list)
}
