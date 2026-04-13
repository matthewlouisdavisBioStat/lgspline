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
#' @keywords monomial smoothing regression parametric constrained lagrangian multiplier
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
#' estimated function with respect to predictors, using a monomial basis.
#'
#' The method of Lagrangian multipliers is used to derive a polynomial regression spline
#' that enforces the following smoothing constraints:
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
#' This package supports model fitting for multiple spline and non-spline effects,
#' GLM families, Weibull accelerated failure time (AFT) models,
#' Cox proportional-Hazards models, negative-binomial regression,
#' arbitrary correlation structures, shape constraints, and extensive
#' customization for user-defined models and constraints.
#'
#' In addition, parallel processing capabilities and comprehensive
#' tools for visualization, frequentist, and Bayesian inference are provided.
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
#' @section Response and Predictor Setup:
#' These arguments control the primary data inputs and the initial
#' standardization steps applied before knot placement and fitting.
#' @param predictors Default: NULL. Numeric matrix or data frame of predictor
#'   variables, or a formula when using the formula interface.
#' @param y Default: NULL. Numeric response variable vector.
#' @param formula Default: NULL. Optional statistical formula for model specification,
#'   supporting \code{spl()} (and the alias \code{s()}) for spline terms.
#' @param response Default: NULL. Alternative name for response variable.
#' @param standardize_response Default: TRUE. Logical indicator controlling whether
#'   the response variable should be centered and scaled before model fitting. Only
#'   offered for identity link functions.
#' @param standardize_predictors_for_knots Default: TRUE. Logical flag controlling
#'   whether predictors are internally standardized for partitioning / knot
#'   placement. The exact transformation is handled inside
#'   \code{\link{make_partitions}} and depends on the effective clustering
#'   dimension.
#' @param standardize_expansions_for_fitting Default: TRUE. Logical switch to
#'   standardize polynomial basis expansions during model fitting. Design matrices,
#'   variance-covariance matrices, and coefficients are backtransformed after fitting.
#'   \eqn{\mathbf{U}} and \eqn{\mathbf{G}} remain on the transformed scale;
#'   \code{B_raw} corresponds to coefficients on the expansion-standardized scale.
#' @param family Default: \code{gaussian()}. GLM family specifying the error distribution
#'   and link function. Minimally requires: family name, link name, linkfun, linkinv,
#'   variance.
#' @param data Default: NULL. Optional data frame for formula-based model specification.
#' @param weights Default: NULL. Alias for \code{observation_weights}.
#' @param observation_weights Default: NULL. Numeric vector of observation-specific
#'   weights for generalized least squares estimation.
#' @param no_intercept Default: FALSE. Logical flag to constrain intercept to 0.
#'   Formulas with \code{"0+"} set this to TRUE automatically.
#'
#' @section GLM Customization:
#' These options let you override the default GLM working-weight,
#' dispersion, and partition-wise unconstrained fitting behavior.
#' @param glm_weight_function Default: function returning \code{family$variance(mu)},
#'   optionally multiplied by \code{observation_weights}. Codes the mean-variance
#'   relationship (the diagonal \eqn{\mathbf{W}} matrix) used for updating
#'   \eqn{\mathbf{G}} after obtaining constrained estimates.
#' @param schur_correction_function Default: function returning list of zeros.
#'   Computes Schur complements \eqn{\mathbf{S}} added to \eqn{\mathbf{G}}:
#'   \eqn{\mathbf{G}^{*} = (\mathbf{G}^{-1} + \mathbf{S})^{-1}}.
#' @param need_dispersion_for_estimation Default: FALSE. Logical indicator specifying
#'   whether a dispersion parameter is required for coefficient estimation (e.g.
#'   Weibull AFT).
#' @param dispersion_function Default: function returning mean squared residuals.
#'   Custom function for estimating the exponential dispersion parameter.
#' @param unconstrained_fit_fxn Default: \code{\link{unconstrained_fit_default}}.
#'   Custom function for fitting unconstrained models per partition.
#' @param keep_weighted_Lambda Default: FALSE. Logical flag to retain GLM weights in
#'   penalty constraints using Tikhonov parameterization. Advised for non-canonical GLMs.
#'
#' @section Knots and Partitioning:
#' These arguments determine how the predictor space is partitioned and how
#' knot locations are chosen or reused.
#' @param K Default: NULL. Integer specifying the number of knot locations. Intuitively,
#'   total partitions minus 1.
#' @param custom_knots Default: NULL. Optional matrix providing user-specified knot
#'   locations in 1-D.
#' @param cluster_on_indicators Default: FALSE. Logical flag for whether indicator
#'   variables should be used for clustering knot locations.
#' @param make_partition_list Default: NULL. Optional list allowing direct specification
#'   of custom partition assignments. The \code{make_partition_list} returned by one
#'   model can be supplied here to reuse knot locations.
#' @param cluster_args Default: \code{c(custom_centers = NA, nstart = 10)}. Named
#'   vector of arguments controlling clustering. If the first argument is not
#'   \code{NA}, it is treated as custom cluster centers (typically an
#'   \eqn{(K+1) \times q} matrix). Otherwise, default k-means is used.
#' @param do_not_cluster_on_these Default: \code{c()}. Predictor columns to exclude from
#'   clustering. Accepts numeric column indices or character column names.
#' @param neighbor_tolerance Default: \code{1 + 1e-8}. Numeric tolerance for determining
#'   neighboring partitions using k-means clustering. Intended for internal use.
#'
#' @section Penalty:
#' These arguments configure the smoothing penalty itself and the optional
#' tuning procedure, using exact leave-one-out by default and GCV optionally.
#' @param previously_tuned_penalties Default: NULL. Optional list of pre-computed penalty
#'   components from a previous model fit.
#' @param smoothing_spline_penalty Default: NULL. Optional custom smoothing spline penalty
#'   matrix.
#' @param opt Default: TRUE. Logical switch controlling automatic penalty optimization.
#' @param tuning_criterion Default: \code{"loo"}. Character scalar selecting the
#'   tuning criterion. Use \code{"loo"} for exact leave-one-out on the transformed
#'   tuning problem, or \code{"gcv"} for the generalized cross-validation
#'   criterion. The LOO path computes the needed hat-matrix diagonal exactly
#'   from blockwise constrained-\eqn{\mathbf{G}} quantities, without explicitly
#'   forming the full projection matrix or full hat matrix. In empirical
#'   diagnostics, the observation-wise derivative of the LOO leverage term can
#'   be numerically delicate even when the overall tuning criterion and fitted
#'   penalties remain well behaved; users who prefer a more conservative
#'   optimization path can set \code{use_custom_bfgs = FALSE}. Shared wiggle
#'   and flat tuning directions are differentiated directly, while optional
#'   predictor- and partition-specific penalties continue to use a lower-cost
#'   ratio approximation. For very large samples, generalized cross-validation
#'   is often the more practical choice; as a rough guideline, \code{"gcv"}
#'   is recommended once the sample size is above about 250,000.
#' @param use_custom_bfgs Default: TRUE. Selects between a native damped-BFGS
#'   implementation with closed-form gradients or base R's BFGS with finite-difference
#'   gradients. The native path is usually faster, while the finite-difference
#'   fallback can be preferable when tuning under exact LOO and the leverage
#'   derivative is numerically noisy.
#' @param delta Default: NULL. Numeric pseudocount for stabilizing optimization in
#'   non-identity link function scenarios.
#' @param tol Default: \code{10*sqrt(.Machine$double.eps)}. Numeric convergence tolerance.
#' @param gcv_gamma Default: 1.4. Numeric scalar, at least 1, used only when
#'   \code{tuning_criterion = "gcv"}. Multiplies the effective degrees of freedom
#'   in the GCV denominator during automatic penalty tuning. It is accepted but
#'   ignored when \code{tuning_criterion = "loo"}.
#' @param initial_wiggle Default: \code{c(2e-12, 2e-7, 2e-4, 0.2)}. Numeric vector of
#'   initial grid points for wiggle penalty optimization, on the raw (non-negative) scale.
#' @param initial_flat Default: \code{c(0.5, 5)}. Numeric vector of initial grid points
#'   for ridge penalty optimization, on the raw scale (ratio of ridge to wiggle).
#' @param wiggle_penalty Default: 2e-7. Numeric penalty on the integrated squared second
#'   derivative, governing function smoothness.
#' @param flat_ridge_penalty Default: 0.5. Numeric flat ridge penalty for intercepts and
#'   linear terms only. Multiplied by \code{wiggle_penalty} to obtain total ridge penalty.
#' @param unique_penalty_per_partition Default: TRUE. Logical flag allowing penalty
#'   magnitude to differ across partitions.
#' @param unique_penalty_per_predictor Default: TRUE. Logical flag allowing penalty
#'   magnitude to differ between predictors.
#' @param meta_penalty Default: 1e-8. Numeric regularization coefficient for
#'   predictor- and partition-specific penalties during tuning. On the raw scale,
#'   the implemented meta-penalty shrinks these penalty multipliers toward 1; the
#'   wiggle penalty receives only a tiny stabilizing penalty by default.
#' @param predictor_penalties Default: NULL. Optional vector of custom penalties per
#'   predictor, on the raw (positive) scale.
#' @param partition_penalties Default: NULL. Optional vector of custom penalties per
#'   partition, on the raw (positive) scale.
#' @param custom_penalty_mat Default: NULL. Optional \eqn{p \times p} custom penalty
#'   matrix for individual partitions, replacing the default ridge on linear/intercept
#'   terms. Run with \code{dummy_fit = TRUE} first to inspect expansion structure.
#'
#' @section Basis Expansions:
#' These arguments control which polynomial and interaction terms are included
#' in the partition-specific design matrices.
#' @param include_quadratic_terms Default: TRUE. Logical switch to include squared
#'   predictor terms.
#' @param include_cubic_terms Default: TRUE. Logical switch to include cubic predictor
#'   terms.
#' @param include_quartic_terms Default: NULL. Includes quartic terms; when NULL, set to
#'   FALSE for single predictor and TRUE otherwise. Highly recommended for multi-predictor
#'   models to avoid over-specified constraints.
#' @param include_2way_interactions Default: TRUE. Logical switch for linear two-way
#'   interactions.
#' @param include_3way_interactions Default: TRUE. Logical switch for three-way
#'   interactions.
#' @param include_quadratic_interactions Default: FALSE. Logical switch for
#'   linear-quadratic interaction terms.
#' @param offset Default: Empty vector. Column indices/names to include as offsets.
#'   Coefficients for offset terms are automatically constrained to 1.
#' @param just_linear_with_interactions Default: NULL. Integer or character vector
#'   specifying predictors to retain as linear terms while still allowing
#'   interactions.
#' @param just_linear_without_interactions Default: NULL. Integer or character
#'   vector specifying predictors to retain only as linear terms without
#'   interactions. Eligible for blockfitting.
#' @param exclude_interactions_for Default: NULL. Integer or character vector of
#'   predictors to exclude from all interaction terms.
#' @param exclude_these_expansions Default: NULL. Character vector of basis expansions to
#'   exclude. Named columns of data, or in the form \code{"_1_"}, \code{"_2_"},
#'   \code{"_1_x_2_"}, \code{"_2_^2"} etc.
#' @param custom_basis_fxn Default: NULL. Optional user-defined function for custom basis
#'   expansions. See \code{\link{get_polynomial_expansions}}.
#' @param auto_encode_factors Default: TRUE. Logical switch to automatically one-hot encode
#'   factor or character variables when using the formula interface.
#'
#' @section Constraints:
#' These arguments govern the smoothness equalities and any additional user
#' supplied linear equality constraints.
#' @param include_constrain_fitted Default: TRUE. Logical switch to constrain fitted
#'   values at knot points.
#' @param include_constrain_first_deriv Default: TRUE. Logical switch to constrain first
#'   derivatives at knot points.
#' @param include_constrain_second_deriv Default: TRUE. Logical switch to constrain second
#'   derivatives at knot points.
#' @param include_constrain_interactions Default: TRUE. Logical switch to constrain
#'   interaction terms at knot points.
#' @param constraint_values Default: \code{cbind()}. Optional matrix encoding
#'   nonzero equality targets paired with \code{constraint_vectors}. When left
#'   empty, added equality constraints are treated as homogeneous.
#' @param constraint_vectors Default: \code{cbind()}. Optional matrix of
#'   user-supplied equality-constraint vectors, appended to the internally
#'   generated smoothness constraints.
#' @param null_constraint Default: NULL. Alternative parameterization for a
#'   nonzero equality target when \code{constraint_vectors} is supplied and
#'   \code{constraint_values} is left empty.
#'
#' @section Quadratic Programming:
#' These arguments activate built-in or custom inequality constraints handled
#' through quadratic programming.
#' @param qp_score_function Default: \eqn{\mathbf{X}^{\top}(\mathbf{y} - \boldsymbol{\mu})}.
#'   Score function for quadratic programming, blockfit, and GEE formulations.
#'   Accepts arguments \code{"X, y, mu, order_list, dispersion, VhalfInv,
#'   observation_weights, ..."}.
#' @param qp_observations Default: NULL. Numeric vector of observation indices
#'   at which built-in QP constraints are evaluated. Useful for reducing the
#'   size of the constrained system.
#' @param qp_Amat Default: NULL. Optional pre-built QP constraint matrix.
#'   In the current pipeline its presence marks QP handling as active, but the
#'   built-in constructor does not merge it into the assembled constraint set;
#'   use \code{qp_Amat_fxn} for custom assembled constraints.
#' @param qp_bvec Default: NULL. Optional pre-built QP right-hand side paired
#'   with \code{qp_Amat}. Like \code{qp_Amat}, it is currently treated as an
#'   advanced placeholder rather than merged into the built-in constructor.
#' @param qp_meq Default: 0. Optional number of equality constraints paired
#'   with \code{qp_Amat}. Like \code{qp_Amat}, it is currently treated as an
#'   advanced placeholder rather than merged into the built-in constructor.
#' @param qp_positive_derivative Default: FALSE. Constrain function to have positive first derivatives. Accepts: \code{FALSE} (no constraint), \code{TRUE} (all predictors), or a character/integer vector naming specific predictor variables to constrain. For example, \code{qp_positive_derivative = "Dose"} constrains only the Dose variable, while \code{qp_positive_derivative = c(1, 3)} constrains columns 1 and 3 of the predictor matrix.
#' @param qp_negative_derivative Default: FALSE. Constrain function to have negative first derivatives. Same input types as \code{qp_positive_derivative}. Can be used simultaneously with \code{qp_positive_derivative} on different variables.
#' @param qp_positive_2ndderivative Default: FALSE. Constrain function to have positive (convex) second derivatives. Same input types as \code{qp_positive_derivative}.
#' @param qp_negative_2ndderivative Default: FALSE. Constrain function to have negative (concave) second derivatives. Same input types as \code{qp_positive_derivative}.
#' @param qp_monotonic_increase Default: FALSE. Logical only. Constrain fitted values to be monotonically increasing in observation order.
#' @param qp_monotonic_decrease Default: FALSE. Logical only. Constrain fitted values to be monotonically decreasing in observation order.
#' @param qp_range_upper Default: NULL. Numeric upper bound for constrained fitted values.
#' @param qp_range_lower Default: NULL. Numeric lower bound for constrained fitted values.
#' @param qp_Amat_fxn Default: NULL. Custom function generating Amat.
#' @param qp_bvec_fxn Default: NULL. Custom function generating bvec.
#' @param qp_meq_fxn Default: NULL. Custom function generating meq.
#'
#' @section Parallel Processing:
#' These arguments control which computational subroutines may run in parallel
#' and how work is chunked across cluster workers.
#' @param cl Default: NULL. Parallel processing cluster object
#'   (use \code{parallel::makeCluster()}).
#' @param chunk_size Default: NULL. Integer specifying custom chunk size for parallel
#'   processing.
#' @param parallel_eigen Default: TRUE. Logical flag for parallel eigenvalue decomposition.
#' @param parallel_trace Default: FALSE. Logical flag for parallel trace computation.
#' @param parallel_aga Default: FALSE. Logical flag for parallel \eqn{\mathbf{G}} and
#'   \eqn{\mathbf{A}} matrix operations.
#' @param parallel_matmult Default: FALSE. Logical flag for parallel block-diagonal matrix
#'   multiplication.
#' @param parallel_unconstrained Default: TRUE. Logical flag for parallel unconstrained
#'   MLE for non-identity-link-Gaussian models.
#' @param parallel_find_neighbors Default: FALSE. Logical flag for parallel neighbor
#'   identification.
#' @param parallel_penalty Default: FALSE. Logical flag for parallel penalty matrix
#'   construction.
#' @param parallel_make_constraint Default: FALSE. Logical flag for parallel constraint
#'   matrix generation.
#'
#' @section Tuning Control:
#' These options control iterative updates during penalty tuning and the final
#' constrained fit.
#' @param iterate_tune Default: TRUE. Logical switch for iterative optimization during
#'   penalty tuning.
#' @param iterate_final_fit Default: TRUE. Logical switch for iterative optimization in
#'   final model fitting.
#' @param blockfit Default: TRUE. Logical switch for backfitting with mixed spline and
#'   non-interactive linear terms. Requires flat columns, \code{K > 0}, and no active
#'   correlation structure. Falls back to \code{get_B} on failure.
#'
#' @section Return Control:
#' These arguments determine which intermediate matrices and inferential
#' quantities are retained in the returned fit object.
#' @param return_G Default: TRUE. Logical switch to return the unscaled unconstrained
#'   variance-covariance matrix \eqn{\mathbf{G}}.
#' @param return_Ghalf Default: TRUE. Logical switch to return
#'   \eqn{\mathbf{G}^{1/2}}.
#' @param return_U Default: TRUE. Logical switch to return the constraint projection
#'   matrix \eqn{\mathbf{U}}.
#' @param estimate_dispersion Default: TRUE. Logical flag to estimate dispersion after
#'   fitting.
#' @param unbias_dispersion Default: NULL. Logical switch to multiply dispersion by
#'   \eqn{N/(N - \mathrm{trace}(\mathbf{H}))}. When NULL, set to TRUE for Gaussian
#'   identity link and FALSE otherwise.
#' @param return_varcovmat Default: TRUE. Logical switch to return the variance-covariance
#'   matrix of estimated coefficients. Needed for Wald inference.
#' @param exact_varcovmat Default: FALSE. Logical switch to replace the default
#'   asymptotic (Bayesian posterior) variance-covariance matrix with the exact
#'   frequentist variance-covariance matrix of the constrained estimator. The
#'   asymptotic version uses the Hessian of
#'   the penalized log-likelihood:
#'   \eqn{\tilde{\sigma}^{2}\mathbf{U}\mathbf{G}\mathbf{U}^{\top}}.
#'   The exact version additionally corrects for the penalty's contribution as a shrinkage
#'   prior, giving:
#'   \deqn{\tilde{\sigma}^{2}\mathbf{U}\mathbf{G}\mathbf{U}^{\top}
#'         - \tilde{\sigma}^{2}\mathbf{U}\mathbf{G}\boldsymbol{\Lambda}\mathbf{G}\mathbf{U}^{\top}}
#'   When a correlation structure is present (\code{VhalfInv} non-NULL),
#'   \eqn{\mathbf{G}_{\mathrm{correct}}} replaces the block-diagonal \eqn{\mathbf{G}}.
#'   For Gaussian identity link (with or without correlation structure), the result is
#'   the exact variance-covariance matrix of the constrained estimate.
#'   The returned object still stores the result in \code{varcovmat}. Requires
#'   \code{return_varcovmat = TRUE}.
#' @param return_lagrange_multipliers Default: FALSE. Logical switch to return the
#'   Lagrangian multiplier vector.
#'
#' @section Correlation Structures:
#' These arguments enable built-in or custom working-correlation structures for
#' longitudinal, clustered, or spatially indexed responses. In the notation
#' used throughout \code{\link{Details}}, the correlated penalized
#' information is written as
#' \eqn{\mathbf{G}^{-1} =
#' \mathbf{G}_{\mathrm{on}}^{-1} + \mathbf{G}_{\mathrm{off}}^{-1}}.
#' When the cross-partition part is low rank, the internal Woodbury helpers
#' factor
#' \eqn{\mathbf{G}_{\mathrm{off}}^{-1} =
#' \mathbf{E}\mathbf{J}\mathbf{E}^{\top}}
#' and accelerate the correlated solve without changing the final estimator.
#' @param correlation_id,spacetime Default: NULL. N-length vector and N-row matrix of
#'   cluster ids and longitudinal/spatial variables, respectively.
#' @param correlation_structure Default: NULL. Native implementations: \code{"exchangeable"},
#'   \code{"spatial-exponential"}, \code{"squared-exponential"}, \code{"ar(1)"},
#'   \code{"spherical"}, \code{"gaussian-cosine"}, \code{"gamma-cosine"},
#'   \code{"matern"}, and aliases. Estimated via REML.
#' @param VhalfInv Default: NULL. Fixed custom \eqn{N \times N} square-root-inverse
#'   covariance matrix \eqn{\mathbf{V}^{-1/2}}. Triggers GLS with known covariance.
#'   Post-fit inference recomputed from whitened Gram matrices.
#' @param Vhalf Default: NULL. Fixed custom \eqn{N \times N} square-root covariance
#'   \eqn{\mathbf{V}^{1/2}}. Computed as inverse of \code{VhalfInv} if not supplied.
#' @param VhalfInv_fxn Default: NULL. Parametric function for \eqn{\mathbf{V}^{-1/2}};
#'   takes single numeric vector \code{"par"}, returns \eqn{N \times N} matrix.
#'   Optimized via BFGS when \code{VhalfInv_par_init} is provided.
#' @param Vhalf_fxn Default: NULL. Optional function for efficient computation of
#'   \eqn{\mathbf{V}^{1/2}} from the same parameter vector used by
#'   \code{VhalfInv_fxn}. When omitted, \code{Vhalf} is obtained by explicit
#'   matrix inversion of \code{VhalfInv}.
#' @param VhalfInv_par_init Default: \code{c()}. Initial parameter values for
#'   \code{VhalfInv_fxn} optimization, on unbounded transformed scale.
#' @param REML_grad Default: NULL. Function for the gradient of the negative REML (or
#'   custom loss) with respect to the parameters of \code{VhalfInv_fxn}. Takes
#'   \code{"par"}, \code{"model_fit"}, and \code{"..."}.
#' @param custom_VhalfInv_loss Default: NULL. Alternative to negative REML for the
#'   correlation parameter objective function. Takes \code{"par"}, \code{"model_fit"},
#'   and \code{"..."}.
#' @param VhalfInv_logdet Default: NULL. Function for efficient \eqn{\log|\mathbf{V}^{-1/2}|}
#'   computation. Takes same \code{"par"} as \code{VhalfInv_fxn}.
#'
#' @section Grouped Argument Lists:
#' For convenience, related arguments can be bundled into named lists.
#' When a grouped argument is non-NULL, its entries overwrite the corresponding
#' individual arguments. Individual arguments remain available for backward compatibility.
#'
#' \describe{
#'   \item{\code{penalty_args}}{Groups: \code{wiggle_penalty},
#'     \code{flat_ridge_penalty}, \code{unique_penalty_per_partition},
#'     \code{unique_penalty_per_predictor}, \code{meta_penalty},
#'     \code{predictor_penalties}, \code{partition_penalties},
#'     \code{custom_penalty_mat}, \code{previously_tuned_penalties},
#'     \code{smoothing_spline_penalty}.}
#'   \item{\code{tuning_args}}{Groups: \code{opt},
#'     \code{use_custom_bfgs}, \code{delta}, \code{tol},
#'     \code{tuning_criterion},
#'     \code{gcv_gamma},
#'     \code{initial_wiggle}, \code{initial_flat},
#'     \code{iterate_tune}, \code{iterate_final_fit}.}
#'   \item{\code{expansion_args}}{Groups: \code{include_quadratic_terms},
#'     \code{include_cubic_terms}, \code{include_quartic_terms},
#'     \code{include_2way_interactions}, \code{include_3way_interactions},
#'     \code{include_quadratic_interactions},
#'     \code{just_linear_with_interactions},
#'     \code{just_linear_without_interactions},
#'     \code{exclude_interactions_for}, \code{exclude_these_expansions},
#'     \code{custom_basis_fxn}, \code{offset}.}
#'   \item{\code{constraint_args}}{Groups:
#'     \code{include_constrain_fitted},
#'     \code{include_constrain_first_deriv},
#'     \code{include_constrain_second_deriv},
#'     \code{include_constrain_interactions},
#'     \code{constraint_values}, \code{constraint_vectors},
#'     \code{no_intercept}.}
#'   \item{\code{qp_args}}{Groups all \code{qp_*} arguments.}
#'   \item{\code{parallel_args}}{Groups: \code{cl},
#'     \code{chunk_size}, and all \code{parallel_*} flags.}
#'   \item{\code{covariance_args}}{Groups: \code{correlation_id},
#'     \code{spacetime}, \code{correlation_structure},
#'     \code{VhalfInv}, \code{Vhalf}, \code{VhalfInv_fxn},
#'     \code{Vhalf_fxn}, \code{VhalfInv_par_init},
#'     \code{REML_grad}, \code{custom_VhalfInv_loss},
#'     \code{VhalfInv_logdet}.}
#'   \item{\code{return_args}}{Groups: \code{return_G},
#'     \code{return_Ghalf}, \code{return_U},
#'     \code{estimate_dispersion}, \code{unbias_dispersion},
#'     \code{return_varcovmat}, \code{exact_varcovmat},
#'     \code{return_lagrange_multipliers}.}
#'   \item{\code{glm_args}}{Groups: \code{glm_weight_function},
#'     \code{schur_correction_function},
#'     \code{need_dispersion_for_estimation},
#'     \code{dispersion_function}, \code{unconstrained_fit_fxn},
#'     \code{keep_weighted_Lambda}.}
#' }
#'
#' @section Miscellaneous:
#' These remaining arguments affect inference defaults, numerical safeguards,
#' verbosity, and developer-oriented diagnostics.
#' @param critical_value Default: \code{qnorm(1-0.05/2)}. Numeric critical value for
#'   Wald confidence intervals.
#' @param dummy_dividor Default: 0.00000000000000000000012345672152894. Small numeric
#'   constant to prevent division by zero.
#' @param dummy_adder Default: 0.000000000000000002234567210529. Small numeric constant
#'   to prevent division by zero.
#' @param verbose Default: FALSE. Logical flag to print general progress messages.
#' @param verbose_tune Default: FALSE. Logical flag to print detailed progress during
#'   penalty tuning.
#' @param dummy_fit Default: FALSE. Runs the full pipeline but sets coefficients to zero,
#'   allowing inspection of design matrix structure, penalty matrices, and partitioning.
#'   Replaces the deprecated \code{expansions_only} argument.
#' @param include_warnings Default: TRUE. Logical switch to control display of warnings.
#' @param penalty_args Default: NULL. Optional named list grouping penalty-related
#'   arguments. See section "Grouped Argument Lists".
#' @param tuning_args Default: NULL. Optional named list grouping tuning-related
#'   arguments.
#' @param expansion_args Default: NULL. Optional named list grouping basis expansion
#'   arguments.
#' @param constraint_args Default: NULL. Optional named list grouping constraint arguments.
#' @param qp_args Default: NULL. Optional named list grouping quadratic programming
#'   arguments.
#' @param parallel_args Default: NULL. Optional named list grouping parallel processing
#'   arguments.
#' @param covariance_args Default: NULL. Optional named list grouping correlation
#'   structure arguments.
#' @param return_args Default: NULL. Optional named list grouping return-control
#'   arguments.
#' @param glm_args Default: NULL. Optional named list grouping GLM customization
#'   arguments.
#' @param ... Additional arguments passed to the unconstrained model fitting function.
#'
#' @return A list of class \code{"lgspline"} containing model components:
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
#'   \item{predict}{Function for generating predictions on new data. For multi-predictor models, \code{take_first_derivatives = TRUE, take_second_derivatives} returns derivatives as a named list of components per predictor variable, rather than a concatenated vector. When \code{new_predictors} contains columns not present in the data, extraneous columns are silently dropped before prediction.}
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
#'   \item{weights}{Original observation weights on the data scale. When no
#'     weights were supplied, this is a vector of ones.}
#'   \item{G}{List of unscaled partition-wise information inverses
#'     \eqn{\mathbf{G}_{k}} if \code{return_G = TRUE}. These are the blockwise
#'     quantities stored on the fitting scale; correlation-aware trace,
#'     posterior, and variance calculations additionally use dense GLS analogues
#'     internally when needed.}
#'   \item{Ghalf}{List of \eqn{\mathbf{G}_{k}^{1/2}} matrices if
#'     \code{return_Ghalf = TRUE}. As with \code{G}, dense GLS square-root
#'     factors may also be constructed internally for correlation-aware
#'     post-fit calculations.}
#'   \item{U}{Constraint projection matrix \eqn{\mathbf{U}} if \code{return_U = TRUE}. For K=0 and no constraints, returns identity. Otherwise, returns \eqn{\mathbf{U} = \mathbf{I} - \mathbf{G}\mathbf{A}(\mathbf{A}^{\top}\mathbf{G}\mathbf{A})^{-1}\mathbf{A}^{\top}}. Used for computing the variance-covariance matrix \eqn{\sigma^{2}\mathbf{U}\mathbf{G}}.}
#'   \item{sigmasq_tilde}{Estimated (or fixed) dispersion parameter
#'     \eqn{\tilde{\sigma}^{2}}. For Gaussian identity fits without correlation,
#'     this is the weighted mean squared residual with optional bias correction.
#'     When \code{VhalfInv} is non-\code{NULL}, Gaussian-identity residuals are
#'     whitened before this calculation.}
#'   \item{trace_XUGX}{Effective degrees of freedom (\eqn{\mathrm{trace}(\mathbf{X}\mathbf{U}\mathbf{G}\mathbf{X}^{\top})}), where \eqn{\mathbf{X}\mathbf{U}\mathbf{G}\mathbf{X}^{\top}} serves as the "hat" matrix. When \code{VhalfInv} is non-\code{NULL}, computed as \eqn{\|\mathbf{V}^{-1/2}\mathbf{X}\mathbf{U}\mathbf{G}_{\mathrm{correct}}^{1/2}\|_{F}^{2}} using the full penalized GLS information.}
#'   \item{varcovmat}{Variance-covariance matrix of coefficient estimates if \code{return_varcovmat = TRUE}. Computed as \eqn{\sigma^{2}(\mathbf{U}\mathbf{G}^{1/2})(\mathbf{U}\mathbf{G}^{1/2})^{\top}} for numerical stability. When \code{VhalfInv} is non-\code{NULL}, uses the full \eqn{\mathbf{G}_{\mathrm{correct}}^{1/2}} in place of the block-diagonal \eqn{\mathbf{G}^{1/2}}.}
#'   \item{lagrange_multipliers}{Vector of Lagrangian multipliers if \code{return_lagrange_multipliers = TRUE}. For equality-only fits these correspond to the active columns of \eqn{\mathbf{A}}; when quadratic-programming constraints are active they are taken directly from \code{solve.QP} and therefore refer to the combined equality/inequality constraint system. \code{NULL} if no constraints are active (\eqn{\mathbf{A}} is \code{NULL} or \code{K == 0}).}
#'   \item{VhalfInv}{The \eqn{\mathbf{V}^{-1/2}} matrix used for implementing correlation structures, if specified.}
#'   \item{VhalfInv_fxn, Vhalf_fxn, VhalfInv_logdet, REML_grad}{Functions for generating \eqn{\mathbf{V}^{-1/2}}, \eqn{\mathbf{V}^{1/2}}, \eqn{\log|\mathbf{V}^{-1/2}|}, and gradient of REML if provided.}
#'   \item{VhalfInv_params_estimates}{Vector of estimated correlation parameters when using \code{VhalfInv_fxn}.}
#'   \item{VhalfInv_params_vcov}{Approximate variance-covariance matrix of estimated correlation parameters from BFGS optimization.}
#'   \item{wald_univariate}{Function for computing univariate Wald statistics and confidence intervals. Returns an S3 object of class \code{"wald_lgspline"} with dedicated \code{print}, \code{summary}, \code{plot}, \code{coef}, and \code{confint} methods. The \code{print} method uses \code{printCoefmat()} for standard R coefficient table formatting with significance stars.}
#'   \item{critical_value}{Critical value used for confidence interval construction.}
#'   \item{generate_posterior}{Function for drawing from the posterior distribution of coefficients. When \code{VhalfInv} is non-\code{NULL}, draws are from the correct joint posterior \eqn{\mathbf{U}\mathbf{G}_{\mathrm{correct}}^{1/2}\mathbf{z}} using the full penalized GLS information, reflecting cross-partition posterior covariance induced by off-diagonal blocks of \eqn{\mathbf{V}^{-1/2}}.}
#'   \item{find_extremum}{Function for optimizing the fitted function. Accepts both numeric column indices and character column names for \code{vars}. When \code{select_vars_fl = TRUE}, L-BFGS-B bounds are correctly subsetted to the optimized variables.}
#'   \item{plot}{Function for visualizing fitted curves.}
#'   \item{quadprog_list}{List containing quadratic programming components if applicable.}
#'   \item{.fit_call_args}{List containing the arguments passed to
#'     \code{\link{lgspline}}.}
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
#' ## Simulate some data, fit using default settings without tuning, and plot
#' set.seed(1234)
#' t <- runif(2500, -10, 10)
#' y <- 2*sin(t) + -0.06*t^2 + rnorm(length(t))
#' model_fit <- lgspline(t, y, opt = FALSE)
#' plot(t, y, main = 'Observed Data vs. Fitted Function, Colored by Partition',
#'      ylim = c(-10, 10))
#' plot(model_fit, add = TRUE)
#'
#' \donttest{
#' ## Repeat using logistic regression, with univariate inference shown
#' # and alternative function call
#' y <- rbinom(length(y), 1, 1/(1+exp(-std(y))))
#' df <- data.frame(t = t, y = y)
#' model_fit <- lgspline(y ~ spl(t),
#'                       df,
#'                       family = binomial())
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
#' model_fit <- lgspline(t, y, opt = FALSE)
#' model_fit <- lgspline(y ~ spl(t), as.data.frame(dat), opt = FALSE)
#' model_fit <- lgspline(response = y, predictors = t, opt = FALSE)
#' model_fit <- lgspline(data = as.data.frame(dat), formula = y ~ ., opt = FALSE)
#'
#' # This is not valid: lgspline(y ~ ., t)
#' # This is not valid: lgspline(y, data = as.data.frame(dat))
#' # Do not put operations in formulas, not valid: lgspline(y ~ log(t) + spl(t))
#'
#' ## Basic Functionality
#' predict(model_fit, new_predictors = rnorm(1)) # make prediction on new data
#' loo_vals <- suppressWarnings(head(leave_one_out(model_fit)))
#' loo_vals # may contain NA when leverage is too high
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
#' ## Posterior draws under constraint
#' draw <- generate_posterior(model_fit, enforce_qp_constraints = TRUE)
#' pr <- predict(model_fit, B_predict = draw$post_draw_coefficients)
#' points(t, pr, col = 'grey')
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
#'                       return_varcovmat = FALSE,
#'                       estimate_dispersion = TRUE,
#'                       return_Ghalf = FALSE,
#'                       return_G = FALSE,
#'                       include_constrain_second_deriv = FALSE,
#'                       unique_penalty_per_predictor = FALSE,
#'                       unique_penalty_per_partition = FALSE,
#'                       wiggle_penalty = 1e-10, # the fixed wiggle penalty
#'                       flat_ridge_penalty = 1e-2) # the ridge penalty / wiggle penalty
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
#' ## Get AUC
#' area_under_volcano <- integrate(model_fit,
#'                                lower = apply(volcano_long, 2, min)[1:2],
#'                                upper = apply(volcano_long, 2, max)[1:2])
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
#'      ...){
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
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   data("pbc", package = "survival")
#'   df <- data.frame(na.omit(
#'     pbc[, c("time", "trt", "stage", "hepato", "bili", "age", "status")]
#'   ))
#'
#'   ## Weibull AFT using lgspline, showing how some custom options can be used to
#'   # fit more complicated models
#'   model_fit <- lgspline(time ~ trt + stage + hepato + bili + age,
#'                         df,
#'                         family = weibull_family(),
#'                         need_dispersion_for_estimation = TRUE,
#'                         dispersion_function = weibull_dispersion_function,
#'                         glm_weight_function = weibull_glm_weight_function,
#'                         schur_correction_function = weibull_schur_correction,
#'                         unconstrained_fit_fxn = unconstrained_fit_weibull,
#'                         opt = FALSE,
#'                         wiggle_penalty = 0,
#'                         flat_ridge_penalty = 0,
#'                         K = 0,
#'                         status = df$status != 0)
#'   print(summary(model_fit))
#'
#'   ## Survreg results match closely on estimates and inference for coefficients
#'   survreg_fit <- survival::survreg(
#'     survival::Surv(time, status != 0) ~ trt + stage + hepato + bili + age,
#'     df
#'   )
#'   print(summary(survreg_fit))
#'
#'   ## sigmasq_tilde = scale^2 of survreg
#'   print(c(sqrt(model_fit$sigmasq_tilde), survreg_fit$scale))
#' }
#'
#' ## ## ## ## Modelling Correlation Structures ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#' ## Setup
#' n_blocks <- 200 # Number of correlation_ids (subjects)
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
#' #  CI is constructed on the working scale and back-transformed
#' ci_transformed <- confint(model_fit)['Correlation parameter 1',]
#' ci_natural <- sort(exp(-exp(ci_transformed)))
#' print("95% CI for correlation:")
#' print(ci_natural)
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
#' n_blocks <- 250   # Number of correlation_ids
#' block_size <- 8   # Observations per correlation_id
#' N <- n_blocks * block_size # total sample size
#' sigma_true <- 2   # Marginal standard deviation
#'
#' ## True Toeplitz components
#' # This example uses a convex combination of two geometric lag kernels:
#' # corr(h) = mix * rho_fast^h + (1 - mix) * rho_slow^h
#' # which is Toeplitz and positive definite for mix in (0, 1) and
#' # rho_fast, rho_slow in (0, 1).
#' rho_fast_true <- 0.25
#' rho_slow_true <- 0.75
#' mix_true <- 0.40
#'
#' ## Create time and correlation_id variables
#' time_var <- rep(1:block_size, n_blocks)
#' correlation_id_var <- rep(1:n_blocks, each = block_size)
#'
#' ## Create nonlinear predictor-response relationship
#' # Not sinusoidal and not polynomial.
#' t_base <- seq(-2, 2, length.out = block_size)
#' t <- rep(t_base, n_blocks) + rnorm(N, sd = 0.10)
#' f_true <- function(t) {
#'   1.4 + 0.9 * atan(1.8 * t) + 0.8 * exp(-1.2 * (t - 0.4)^2)
#' }
#'
#' ## Generate mean structure
#' mu_true <- f_true(t)
#'
#' ## Toeplitz correlation helper
#' corr_from_components <- function(rho_fast, rho_slow, mix) {
#'   corr <- matrix(0, block_size, block_size)
#'   for(i in 1:block_size) {
#'     for(j in 1:block_size) {
#'       lag <- abs(i - j)
#'       if(lag == 0) {
#'         corr[i, j] <- 1
#'       } else {
#'         corr[i, j] <- mix * rho_fast^lag + (1 - mix) * rho_slow^lag
#'       }
#'     }
#'   }
#'   corr
#' }
#'
#' ## Toeplitz correlation function
#' # Custom functions can use any parameterization. Here we map:
#' #   par[1] -> rho_fast = exp(-exp(par[1]))
#' #   par[2] -> rho_slow = exp(-exp(par[2]))
#' #   par[3] -> mix      = plogis(par[3])
#' # so the parameter space is unconstrained, while the resulting Toeplitz
#' # correlation matrix remains valid.
#' corr_from_par <- function(par) {
#'   rho_fast <- exp(-exp(par[1]))
#'   rho_slow <- exp(-exp(par[2]))
#'   mix <- plogis(par[3])
#'   corr_from_components(rho_fast, rho_slow, mix)
#' }
#'
#' ## Create block Toeplitz errors from the same family we will fit
#' corr_true <- corr_from_components(rho_fast_true, rho_slow_true, mix_true)
#' errors <- Reduce('c',
#'                  lapply(1:n_blocks, function(i) {
#'                    c(matsqrt(corr_true) %*% rnorm(block_size))
#'                  }))
#'
#' ## Generate response with correlated errors and nonlinear covariate effect
#' y <- mu_true + sigma_true * errors
#'
#' VhalfInv_fxn <- function(par) {
#'   corr <- corr_from_par(par)
#'   kronecker(diag(n_blocks), matinvsqrt(corr))
#' }
#'
#' Vhalf_fxn <- function(par) {
#'   corr <- corr_from_par(par)
#'   kronecker(diag(n_blocks), matsqrt(corr))
#' }
#'
#' ## Determinant function (for efficiency)
#' # This avoids taking determinant of N by N matrix
#' VhalfInv_logdet <- function(par) {
#'   corr <- corr_from_par(par)
#'   log_det_invsqrt_corr <- -0.5 * determinant(corr, logarithm = TRUE)$modulus[1]
#'   n_blocks * log_det_invsqrt_corr
#' }
#'
#' ## GLM weights for REML gradient helper
#' # For Gaussian identity, these are all 1.
#' glm_weight_function <- function(mu, y, order_indices, family,
#'                                 dispersion, observation_weights, ...) {
#'   rep(1, length(mu))
#' }
#'
#' ## REML gradient function
#' # The helper reml_grad_from_dV computes the three REML terms once dV / dpar
#' # is supplied. For this parameterization, dV / dpar has closed form.
#' REML_grad <- function(par, model_fit, ...) {
#'   rho_fast <- exp(-exp(par[1]))
#'   rho_slow <- exp(-exp(par[2]))
#'   mix <- plogis(par[3])
#'
#'   dV1_block <- matrix(0, block_size, block_size)
#'   dV2_block <- matrix(0, block_size, block_size)
#'   dV3_block <- matrix(0, block_size, block_size)
#'
#'   for(i in 1:block_size) {
#'     for(j in 1:block_size) {
#'       lag <- abs(i - j)
#'       if(lag > 0) {
#'         ## d/dpar[1] through rho_fast = exp(-exp(par[1]))
#'         dV1_block[i, j] <- -mix * lag * exp(par[1]) * rho_fast^lag
#'         ## d/dpar[2] through rho_slow = exp(-exp(par[2]))
#'         dV2_block[i, j] <- -(1 - mix) * lag * exp(par[2]) * rho_slow^lag
#'         ## d/dpar[3] through mix = plogis(par[3])
#'         dV3_block[i, j] <- mix * (1 - mix) * (rho_fast^lag - rho_slow^lag)
#'       }
#'     }
#'   }
#'
#'   dV1 <- kronecker(diag(n_blocks), dV1_block)
#'   dV2 <- kronecker(diag(n_blocks), dV2_block)
#'   dV3 <- kronecker(diag(n_blocks), dV3_block)
#'
#'   gradient <- numeric(3)
#'   gradient[1] <- lgspline::reml_grad_from_dV(dV1, model_fit,
#'                                    glm_weight_function, ...)
#'   gradient[2] <- reml_grad_from_dV(dV2, model_fit,
#'                                    glm_weight_function, ...)
#'   gradient[3] <- reml_grad_from_dV(dV3, model_fit,
#'                                    glm_weight_function, ...)
#'   gradient
#' }
#'
#' ## Visualization
#' plot(t, y, col = correlation_id_var,
#'      main = 'Simulated Data with Toeplitz Correlation')
#'
#' ## Fit model with custom Toeplitz
#' model_fit <- lgspline(
#'   response = y,
#'   predictors = t,
#'   K = 4,
#'   standardize_response = FALSE,
#'   VhalfInv_fxn = VhalfInv_fxn,
#'   Vhalf_fxn = Vhalf_fxn,
#'   VhalfInv_logdet = VhalfInv_logdet,
#'   REML_grad = REML_grad,
#'   VhalfInv_par_init = c(0, -1, 0),
#'   include_warnings = FALSE
#' )
#'
#' ## Print comparison of true and estimated correlations
#' lag_values <- 1:(block_size - 1)
#' corr_true_by_lag <- sapply(lag_values, function(h) {
#'   mix_true * rho_fast_true^h + (1 - mix_true) * rho_slow_true^h
#' })
#' rho_fast_est <- exp(-exp(model_fit$VhalfInv_params_estimates[1]))
#' rho_slow_est <- exp(-exp(model_fit$VhalfInv_params_estimates[2]))
#' mix_est <- plogis(model_fit$VhalfInv_params_estimates[3])
#' corr_est_by_lag <- sapply(lag_values, function(h) {
#'   mix_est * rho_fast_est^h + (1 - mix_est) * rho_slow_est^h
#' })
#' cat('Toeplitz Correlation Estimates by Lag:\n')
#' print(data.frame(
#'   Lag = lag_values,
#'   True.Correlation = round(corr_true_by_lag, 4),
#'   Estimated.Correlation = round(corr_est_by_lag, 4)
#' ))
#'
#' ## Should be ~ 2
#' print(sqrt(model_fit$sigmasq_tilde))
#'
#' ## ## ## ## Parallelism ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#' if (requireNamespace("parallel", quietly = TRUE)) {
#'   ## Data generating function
#'   a <- runif(500000, -9, 9)
#'   b <- runif(500000, -9, 9)
#'   c <- rnorm(500000)
#'   d <- rpois(500000, 1)
#'   y <- sin(a) + cos(b) - 0.2*sqrt(a^2 + b^2) +
#'     abs(a) + b +
#'     0.5*(a^2 + b^2) +
#'     (1/6)*(a^3 + b^3) +
#'     a*b*c -
#'     c +
#'     d +
#'     rnorm(500000, 0, 5)
#'
#'   ## Set up cores
#'   cl <- parallel::makeCluster(1)
#'   on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
#'
#'   ## This example shows some options for what operations can be parallelized
#'   # By default, only parallel_eigen and parallel_unconstrained are TRUE
#'   # parallel_unconstrained is only for GLMs, for identity link Gaussian
#'   # response, use parallel_matmult=TRUE to ensure parallel fitting across
#'   # partitions.
#'   # G, G^{-1/2}, and G^{1/2} are computed in parallel across each of the
#'   # K+1 partitions.
#'   # However, parallel_unconstrained only affects GLMs without corr. components
#'   # - it does not affect fitting here
#'   system.time({
#'     parfit <- lgspline(y ~ spl(a, b) + a*b*c + d,
#'                        data = data.frame(y = y,
#'                                          a = a,
#'                                          b = b,
#'                                          c = c,
#'                                          d = d),
#'                        cl = cl,
#'                        K = 1,
#'                        parallel_eigen = TRUE,
#'                        parallel_unconstrained = TRUE,
#'                        parallel_aga = FALSE,
#'                        parallel_find_neighbors = FALSE,
#'                        parallel_trace = FALSE,
#'                        parallel_matmult = TRUE,
#'                        parallel_make_constraint = FALSE,
#'                        parallel_penalty = FALSE)
#'   })
#'   print(summary(parfit))
#' }
#' }
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
    dispersion_function = function(mu, y, order_indices, family,
                                   observation_weights, VhalfInv,
                                   ...) {

      ## If covariance present
      if(!is.null(VhalfInv)){
        VhalfInv <- VhalfInv[order_indices, order_indices]
        c(mean(
          (
            tcrossprod(VhalfInv, t(y-mu))
          )^2 /
            family$variance(mu)
        ))
        ## If no covariance present
      } else{
        c(mean(
          (
            y - mu
          )^2 /
            family$variance(mu)
        ))
      }
    },
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
    tuning_criterion = "loo",
    gcv_gamma = 1.4,
    initial_wiggle = c(2e-12, 2e-7, 2e-4, 0.2),
    initial_flat = c(0.5, 5),
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
    blockfit = TRUE,
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
    unbias_dispersion = NULL,
    return_varcovmat = TRUE,
    exact_varcovmat = FALSE,
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
    ## [Change 2026-03-01] Include list arguments for organization
    penalty_args = NULL,
    tuning_args = NULL,
    expansion_args = NULL,
    constraint_args = NULL,
    qp_args = NULL,
    parallel_args = NULL,
    covariance_args = NULL,
    return_args = NULL,
    glm_args = NULL,
    ...
){

  if(verbose){
    cat('Pre-Processing\n')
  }

  # [Change 2026-03-01] Unpack grouped argument lists.
  #  Each list argument, when non-NULL, overwrites the corresponding
  #  raw argument. This allows users to specify:
  #    lgspline(t, y, penalty_args = list(wiggle_penalty = 1e-4))
  #  instead of:
  #    lgspline(t, y, wiggle_penalty = 1e-4)
  #  Both forms are valid. The grouped form takes precedence.
  .unpack_group <- function(group_list, env){
    if(!is.null(group_list)){
      if(!is.list(group_list)){
        stop('\n \t Grouped argument (e.g., penalty_args) must be a list.\n')
      }
      for(nm in names(group_list)){
        if(exists(nm, envir = env, inherits = FALSE)){
          assign(nm, group_list[[nm]], envir = env)
        } else if(include_warnings){
          warning('Grouped argument "', nm,
                  '" is not a recognized lgspline parameter; ignoring.')
        }
      }
    }
  }
  local_env <- environment()
  .unpack_group(penalty_args, local_env)
  .unpack_group(tuning_args, local_env)
  .unpack_group(expansion_args, local_env)
  .unpack_group(constraint_args, local_env)
  .unpack_group(qp_args, local_env)
  .unpack_group(parallel_args, local_env)
  .unpack_group(covariance_args, local_env)
  .unpack_group(return_args, local_env)
  .unpack_group(glm_args, local_env)

  ## [Change 2026-03-02] Delegate input processing to process_input()
  processed <- process_input(
    predictors = predictors,
    y = y,
    formula = formula,
    response = response,
    data = data,
    weights = weights,
    observation_weights = observation_weights,
    family = family,
    K = K,
    custom_knots = custom_knots,
    auto_encode_factors = auto_encode_factors,
    include_2way_interactions = include_2way_interactions,
    include_3way_interactions = include_3way_interactions,
    just_linear_with_interactions = just_linear_with_interactions,
    just_linear_without_interactions = just_linear_without_interactions,
    exclude_interactions_for = exclude_interactions_for,
    exclude_these_expansions = exclude_these_expansions,
    offset = offset,
    no_intercept = no_intercept,
    do_not_cluster_on_these = do_not_cluster_on_these,
    include_quartic_terms = include_quartic_terms,
    cluster_args = cluster_args,
    include_warnings = include_warnings,
    dummy_fit = dummy_fit,
    include_constrain_second_deriv = include_constrain_second_deriv,
    ...
  )

  ## Unpack processed results
  predictors                       <- processed$predictors
  y                                <- processed$y
  og_cols                          <- processed$og_cols
  replace_colnames                 <- processed$replace_colnames
  just_linear_with_interactions    <- processed$just_linear_with_interactions
  just_linear_without_interactions <- processed$just_linear_without_interactions
  exclude_interactions_for         <- processed$exclude_interactions_for
  exclude_these_expansions         <- processed$exclude_these_expansions
  offset                           <- processed$offset
  no_intercept                     <- processed$no_intercept
  do_not_cluster_on_these          <- processed$do_not_cluster_on_these
  observation_weights              <- processed$observation_weights
  K                                <- processed$K
  include_3way_interactions        <- processed$include_3way_interactions
  include_quartic_terms            <- processed$include_quartic_terms
  data                             <- processed$data
  custom_knots                     <- processed$custom_knots
  dummy_fit                        <- processed$dummy_fit
  include_constrain_second_deriv   <- processed$include_constrain_second_deriv

  ## Alternative parameterization of null constraint for ease of use
  if(any(!is.null(null_constraint)) &
     length(constraint_vectors) > 0 &
     length(constraint_values) == 0){
    constraint_values <-
      crossprod(
        t(constraint_vectors),
        crossprod(
          t(invert(crossprod(cbind(constraint_vectors)))),
          cbind(c(null_constraint))
        )
      )
  }

  ## Unbias dispersion default
  if(is.null(unbias_dispersion)){
    if(!(any(!(paste0(family)[1:4] == paste0(gaussian())[1:4])))){
      unbias_dispersion <- TRUE
    } else {
      unbias_dispersion <- FALSE
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
                                 tuning_criterion,
                                 gcv_gamma,
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
                                 exact_varcovmat,
                                 return_lagrange_multipliers,
                                 custom_penalty_mat,
                                 cluster_args,
                                 dummy_dividor,
                                 dummy_adder,
                                 verbose,
                                 verbose_tune,
                                 dummy_fit | (!is.null(VhalfInv_fxn) &
                                                length(VhalfInv_par_init) > 0),
                                 auto_encode_factors,
                                 observation_weights,
                                 do_not_cluster_on_these,
                                 neighbor_tolerance,
                                 no_intercept,
                                 VhalfInv,
                                 Vhalf,
                                 include_warnings,
                                 og_cols,
                                 NULL, # factor_groups
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
    if(replace_colnames){
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

         ## V^{1/2} used directly whenever correlation-aware code paths
         #  need it, including Gaussian identity fits.
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
                                         tuning_criterion,
                                         gcv_gamma,
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
                                         FALSE, # exact varcovmat
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
                                         NULL, # factor_groups
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
                              tuning_criterion,
                              gcv_gamma,
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
                              exact_varcovmat,
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
                              NULL, # factor_groups
                              ...)
    model_fit$VhalfInv_params_vcov <- VhalfInv_params_vcov
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

      est <- unlist(model_fit$B)
      se <- sqrt(scale_vcovmat_by * diag(model_fit$varcovmat))
      stat <- est / se
      interval_lb <- se * (stat - cv)
      interval_ub <- se * (stat + cv)

      ## If normal errors, use exact t-test
      if(!(any(!(paste0(family)[1:4] == paste0(gaussian())[1:4])))){
        df_resid <- model_fit$N - model_fit$trace_XUGX
        pval <- 2 * (1 - pt(abs(stat), df = df_resid))
        test_type <- "t"
      } else {
        df_resid <- Inf
        pval <- 2 * (1 - pnorm(abs(stat)))
        test_type <- "z"
      }

      ## Coefficient table for printCoefmat
      coef_table <- cbind(
        Estimate   = est,
        Std.Error  = se,
        Statistic  = stat,
        Lower      = interval_lb,
        Upper      = interval_ub,
        p.value    = pval
      )
      if(test_type == "t"){
        colnames(coef_table)[3] <- "t value"
      } else {
        colnames(coef_table)[3] <- "z value"
      }
      colnames(coef_table)[4] <- paste0("Lower ", round(100*(1 -
                                                               2*(1 - pnorm(cv))), 1), "%")
      colnames(coef_table)[5] <- paste0("Upper ", round(100*(1 -
                                                               2*(1 - pnorm(cv))), 1), "%")
      colnames(coef_table)[6] <- paste0("Pr(>|",test_type,"|)")

      ## Build return object
      result <- list(
        coefficients  = coef_table,
        est           = est,
        se            = se,
        stat          = stat,
        interval_lb   = interval_lb,
        interval_ub   = interval_ub,
        pval          = pval,
        df.residual   = df_resid,
        test_type     = test_type,
        critical_value = cv,
        scale_vcovmat_by = scale_vcovmat_by,
        K             = model_fit$K,
        p             = model_fit$p,
        N             = model_fit$N,
        family        = model_fit$family
      )
      class(result) <- "wald_lgspline"
      return(result)
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
                                             function(N_obs, mean, sqrt_dispersion, ...)
                                               rnorm(N_obs, mean, sqrt_dispersion),
                                           draw_dispersion = TRUE,
                                           include_posterior_predictive = FALSE,
                                           num_draws = 1,
                                           enforce_constraints = FALSE,
                                           enforce_qp_constraints = NULL,
                                           max_slice_iterations = 1000,
                                           override_B_raw = NULL,
                                           override_U = NULL,
                                           override_Ghalf_correct = NULL,
                                           override_Ghalf = NULL,
                                           override_VhalfInv = NULL,
                                           override_sigmasq_tilde = NULL,
                                           override_trace_XUGX = NULL,
                                           ...){

    if(!is.null(enforce_qp_constraints)){
      enforce_constraints <- enforce_qp_constraints
    }

    if(any(!is.null(new_predictors))){
      new_predictors <- try(methods::as(cbind(new_predictors), 'matrix'),
                            silent = TRUE)
      if(any(inherits(new_predictors, 'try-error'))){
        stop('\n\t new_predictors must be coercible to a matrix.\n')
      }
    }

    only_1 <- FALSE
    if(nrow(new_predictors) == 1){
      only_1 <- TRUE
      new_predictors <- rbind(new_predictors, new_predictors)
    }

    p_expansions <- model_fit$p
    K  <- model_fit$K
    N_obs <- model_fit$N

    ## Resolve overrides: use injected values when provided, else
    #  fall back to model_fit's captured environment.
    use_B_raw      <- if(!is.null(override_B_raw)) override_B_raw
    else model_fit$B_raw
    use_U          <- if(!is.null(override_U)) override_U
    else model_fit$U
    use_VhalfInv   <- if(!is.null(override_VhalfInv)) override_VhalfInv
    else model_fit$VhalfInv
    use_Ghalf      <- if(!is.null(override_Ghalf)) override_Ghalf
    else model_fit$Ghalf
    use_trace_XUGX <- if(!is.null(override_trace_XUGX)) override_trace_XUGX
    else model_fit$trace_XUGX

    ## When draw_dispersion is TRUE the InvGamma rate uses
    #  new_sigmasq_tilde, which the caller already sets to the
    #  re-estimated value. override_sigmasq_tilde only matters when
    #  draw_dispersion is FALSE and the caller wants the re-estimated
    #  point estimate used as-is.
    if(!is.null(override_sigmasq_tilde) && !draw_dispersion){
      new_sigmasq_tilde <- override_sigmasq_tilde
    }

    ## Extract inequality constraints from QP. First qp_meq columns
    #  are equalities already handled by U; the rest are inequalities
    #  that need slice sampling.
    has_ineq <- FALSE
    has_qp   <- !is.null(model_fit$quadprog_list) &&
      !is.null(model_fit$quadprog_list$qp_Amat)

    if(has_qp && enforce_constraints){
      qp_Amat_full <- model_fit$quadprog_list$qp_Amat
      qp_bvec_full <- model_fit$quadprog_list$qp_bvec
      qp_meq_val   <- model_fit$quadprog_list$qp_meq
      n_total <- ncol(qp_Amat_full)
      if(qp_meq_val < n_total){
        ineq_cols <- (qp_meq_val + 1):n_total
        ineq_Amat <- qp_Amat_full[, ineq_cols, drop = FALSE]
        ineq_bvec <- qp_bvec_full[ineq_cols]
        has_ineq  <- TRUE
      }
    }

    ## t(Amat) %*% beta >= bvec on standardized (B_raw) scale
    .check_feasible <- function(b){
      if(!has_ineq) return(TRUE)
      all(c(crossprod(ineq_Amat, cbind(b))) >=
            ineq_bvec - sqrt(.Machine$double.eps))
    }

    ## Elliptical slice step (Murray et al 2010). Proposes on
    #    f(th) = beta_cur * cos(th) + nu * sin(th)
    #  with nu ~ N(0, sigma^2 L L'), shrinks the bracket until
    #  a feasible point is found.
    .ess_step <- function(beta_cur, L_post, sigma){
      nu <- sigma * c(L_post %**% cbind(rnorm(ncol(L_post))))
      theta <- runif(1, 0, 2 * pi)
      theta_min <- theta - 2 * pi
      theta_max <- theta
      for(iter in 1:max_slice_iterations){
        bp <- beta_cur * cos(theta) + nu * sin(theta)
        if(.check_feasible(bp)) return(bp)
        if(theta < 0) theta_min <- theta else theta_max <- theta
        theta <- runif(1, theta_min, theta_max)
      }
      if(include_warnings){
        warning('\n\t Slice sampler: no feasible point in ',
                max_slice_iterations, ' iterations, returning ',
                'current state.\n')
      }
      beta_cur
    }

    ## InvGamma draw for dispersion
    .draw_sigmasq <- function(){
      if(!draw_dispersion) return(new_sigmasq_tilde)
      half_df <- 0.5 * (N_obs - model_fit$unbias_dispersion *
                          use_trace_XUGX)
      shape <- theta_1 + half_df
      rate  <- theta_2 + half_df * new_sigmasq_tilde
      if(shape <= 0){
        stop('\n\t Posterior inverse-gamma shape <= 0, ',
             'increase theta_1.\n')
      }
      if(rate <= 0){
        stop('\n\t Posterior inverse-gamma rate <= 0, ',
             'increase theta_2.\n')
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

    ## Back-transform a draw on B_raw scale to per-partition
    #  original-scale coefficients.
    .backtransform <- function(beta_raw){
      lapply(1:(K + 1), function(k){
        bk <- beta_raw[1:p_expansions + (k - 1) * p_expansions] * model_fit$sd_y
        bk[1] <- bk[1] + model_fit$mean_y
        model_fit$backtransform_coefficients(bk)
      })
    }

    ## Precompute L_post = (1/sd_y) * U %*% Ghalf. Does not depend
    #  on sigmasq (sigmasq only scales it), so compute once.
    if(!is.null(use_VhalfInv)){

      ## When Ghalf_correct was provided by the correlation re-estimation
      #  path we can skip rebuilding the whitened Gram entirely.
      if(!is.null(override_Ghalf_correct)){
        L_post <- (1 / model_fit$sd_y) *
          (use_U %**% override_Ghalf_correct)
      } else {
        ## GLS path: L_post via Ghalf_correct from whitened Gram
        ord <- unlist(model_fit$order_list)
        X_full <- collapse_block_diagonal(
          lapply(model_fit$X, model_fit$std)
        )
        W_glm <- c(glm_weight_function(
          model_fit$ytilde[ord], model_fit$y[ord],
          1:N_obs, model_fit$family, new_sigmasq_tilde,
          rep(1, N_obs), ...))
        W_glm <- pmax(W_glm, .Machine$double.eps)
        D <- model_fit$weights
        if(is.list(D)) D <- unlist(D)[model_fit$og_order] else D <- D[model_fit$og_order]
        X_full <- X_full * sqrt(W_glm * D)
        VinvhalfX <- use_VhalfInv[ord, ord] %**% X_full
        has_part_pen <- length(model_fit$penalties$L_partition_list) ==
          (K + 1)
        Lambda_full <- collapse_block_diagonal(
          lapply(1:(K + 1), function(k){
            if(has_part_pen){
              model_fit$penalties$Lambda +
                model_fit$penalties$L_partition_list[[k]]
            } else {
              model_fit$penalties$Lambda
            }
          })
        )
        L_post <- (1 / model_fit$sd_y) *
          (use_U %**% matinvsqrt(crossprod(VinvhalfX) + Lambda_full))
      }

    } else {
      ## Block-diagonal path
      P_total  <- p_expansions * (K + 1)
      Ghalf_bd <- matrix(0, P_total, P_total)
      for(k in 1:(K + 1)){
        rows <- ((k - 1) * p_expansions + 1):(k * p_expansions)
        Ghalf_bd[rows, rows] <- use_Ghalf[[k]]
      }
      L_post <- (1 / model_fit$sd_y) * (use_U %**% Ghalf_bd)
    }

    beta_mode <- unlist(use_B_raw)

    ## Package one draw into the return list
    .package <- function(beta_raw, s2){
      coefs <- .backtransform(beta_raw)
      if(include_posterior_predictive){
        pm <- model_fit$predict(new_predictors, B_predict = coefs)
        pp <- posterior_predictive_draw(length(pm), pm, sqrt(s2), ...)
        list(post_pred_draw = pp,
             post_draw_coefficients = coefs,
             post_draw_sigmasq = s2)
      } else {
        list(post_draw_coefficients = coefs,
             post_draw_sigmasq = s2)
      }
    }

    ## Draw. Constrained path is MCMC (chain state across draws);
    #  unconstrained path is iid.
    if(has_ineq){
      beta_cur <- beta_mode
      res <- vector("list", num_draws)
      for(m in 1:num_draws){
        s2 <- .draw_sigmasq()
        beta_cur <- .ess_step(beta_cur, L_post, sqrt(s2))
        res[[m]] <- .package(beta_cur, s2)
      }
    } else {
      res <- lapply(1:num_draws, function(m){
        s2 <- .draw_sigmasq()
        z  <- rnorm(ncol(L_post))
        b  <- beta_mode + sqrt(s2) * c(L_post %**% cbind(z))
        .package(b, s2)
      })
    }

    ## Format output
    if(num_draws == 1){
      out <- res[[1]]
      if(only_1 && !is.null(out$post_pred_draw)){
        out$post_pred_draw <- out$post_pred_draw[1]
      }
      return(out)
    }

    out <- list(
      post_draw_coefficients = lapply(res, `[[`, "post_draw_coefficients"),
      post_draw_sigmasq      = lapply(res, `[[`, "post_draw_sigmasq"))
    if(include_posterior_predictive){
      out$post_pred_draw <- Reduce("cbind",
                                   lapply(res, `[[`, "post_pred_draw"))
      if(only_1){
        out$post_pred_draw <- out$post_pred_draw[1, , drop = FALSE]
      }
    }
    out
  }

  ## Find global maximum/minimum
  model_fit$find_extremum <- function(
    vars = NULL,
    quick_heuristic = TRUE,
    initial = NULL,
    B_predict = NULL,
    minimize = FALSE,
    stochastic = FALSE,
    stochastic_draw = function(mu, sigma, ...){
      N_obs <- length(mu)
      rnorm(N_obs, mu, sigma)
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
    if(inherits(extr, 'matrix')){
      rownames(extr) <- NULL
    }
    return(list(
      t = extr,
      y = preds[global_max]
    ))
  }

  ## One-dimensional plotting function
  # [Change 2026-02-14] Introduced legend_order option
  plot_lgspline_1d <- function(modfit,
                               show_formulas,
                               formula_B = NULL,
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
      formula_terms <- if(is.null(formula_B)) model_fit$B else formula_B
      formulas <- sapply(1:(model_fit$K+1), function(k) {
        coefs <- round(formula_terms[[k]], digits)
        names(coefs) <- rownames(coefs)
        term_names <- rownames(formula_terms[[k]])
        plotted_term_name <- rownames(model_fit$B[[k]])[2]
        term_names <- safe_replace_var(term_names, plotted_term_name, xlab)
        names(coefs) <- term_names
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
                               formula_B = NULL,
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
      formula_terms <- if(is.null(formula_B)) model_fit$B else formula_B
      formulas <- sapply(1:(model_fit$K+1), function(k) {
        coefs <- round(formula_terms[[k]], digits)
        names(coefs) <- rownames(coefs)
        names(coefs) <- gsub("\\^2", "<sup>2</sup>", names(coefs))
        names(coefs) <- gsub("\\^3", "<sup>3</sup>", names(coefs))
        names(coefs) <- gsub("\\^4", "<sup>4</sup>", names(coefs))
        term_names <- rownames(formula_terms[[k]])
        term_names <- safe_replace_var(
          term_names,
          og_cols[as.numeric(substr(v1, 2, nchar(v1)-1))],
          custom_predictor_lab1
        )
        term_names <- safe_replace_var(
          term_names,
          og_cols[as.numeric(substr(v2, 2, nchar(v2)-1))],
          custom_predictor_lab2
        )
        names(coefs) <- term_names
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
                             include_all_terms_in_formulas = FALSE,
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

    ## Default text_size_formula depends on q_predictors
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

    ## Keep the full partition equations available for legend / hover text
    #  when plotting a marginal relationship but requesting the complete
    #  formula instead of only the displayed predictor terms.
    formula_B <- NULL
    if(show_formulas && include_all_terms_in_formulas){
      formula_B <- model_fit_in$B
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
                       formula_B,
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
                       formula_B,
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
    tuning_criterion                   = tuning_criterion,
    gcv_gamma                          = gcv_gamma,
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
    exact_varcovmat                    = exact_varcovmat,
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
#'   \item Penalty tuning via exact leave-one-out by default, or generalized
#'         cross-validation when \code{tuning_criterion = "gcv"}, or use of
#'         previously tuned penalties.
#'   \item Final coefficient estimation via one of three paths:
#'         \itemize{
#'           \item \strong{Blockfit option} (when \code{blockfit = TRUE},
#'                 flat columns are non-empty, \code{K > 0}, and no correlation
#'                 structure): Routes through \code{blockfit_solve} for
#'                 backfitting with mixed spline and non-interactive linear terms.
#'                 Falls back to \code{get_B} on failure.
#'           \item \strong{Standard \code{get_B}} path: Three internal
#'                 computational paths: GEE (damped SQP with correlation
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
#'   \item{sigmasq_tilde}{The estimated (or fixed) dispersion parameter
#'     \eqn{\tilde{\sigma}^{2}}. For Gaussian identity fits with
#'     \code{VhalfInv} non-\code{NULL}, this is computed from whitened
#'     residuals \eqn{\mathbf{V}^{-1/2}(\mathbf{y} - \hat{\mathbf{y}})},
#'     multiplied by the observation weights and the optional
#'     bias-correction factor. When \code{estimate_dispersion = FALSE},
#'     set to 1. Omitted when \code{dummy_fit = TRUE}.}
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
#'   \item{exact_varcovmat}{Not returned as a standalone component; this
#'     argument only controls whether \code{varcovmat}, when requested, is
#'     left as the default asymptotic/Laplace version or replaced by the
#'     exact frequentist correction available for Gaussian identity fits.}
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
#'     matrices. Omitted when \code{dummy_fit = TRUE}.}
#'   \item{Ghalf}{List of \eqn{\mathbf{G}_{k}^{1/2}} matrices, returned if
#'     \code{return_Ghalf = TRUE}. When \code{VhalfInv} is non-\code{NULL},
#'     the full \eqn{\mathbf{G}_{\mathrm{correct}}^{1/2}} is used for
#'     posterior draws and variance-covariance computation.
#'     Omitted when \code{dummy_fit = TRUE}.}
#'   \item{U}{Constraint projection matrix \eqn{\mathbf{U}}, returned if
#'     \code{return_U = TRUE}. Omitted when \code{dummy_fit = TRUE}.}
#'   \item{trace_XUGX}{The effective degrees-of-freedom trace term.
#'     When \code{VhalfInv} is non-\code{NULL}, it is computed from the
#'     full penalized GLS information rather than the block-diagonal
#'     approximation. Omitted when \code{dummy_fit = TRUE}.}
#'   \item{varcovmat}{The final variance-covariance matrix of the estimated
#'     coefficients. Computed via the outer-product form
#'     \eqn{\sigma^{2}(\mathbf{U}\mathbf{G}^{1/2})(\mathbf{U}\mathbf{G}^{1/2})^{\top}}
#'     for numerical stability. When \code{VhalfInv} is non-\code{NULL}, uses
#'     the full \eqn{\mathbf{G}_{\mathrm{correct}}^{1/2}} in place of
#'     block-diagonal \eqn{\mathbf{G}^{1/2}}. Returned if
#'     \code{return_varcovmat = TRUE}. By default this is the asymptotic
#'     (Laplace/posterior) variance-covariance matrix; when
#'     \code{exact_varcovmat = TRUE}, it is replaced in-place by the exact
#'     frequentist correction available for Gaussian identity fits.
#'     Omitted when \code{dummy_fit = TRUE}.}
#'   \item{lagrange_multipliers}{Vector of Lagrangian multipliers if
#'     \code{return_lagrange_multipliers = TRUE}. For equality-only fits these
#'     follow the formulation
#'     \eqn{(\mathbf{A}^{\top}\mathbf{G}\mathbf{A})^{-1}\mathbf{A}^{\top}(\hat{\boldsymbol{\beta}} - \boldsymbol{\beta_0})}.
#'     When quadratic-programming constraints are active they are taken
#'     directly from \code{solve.QP} and therefore refer to the combined
#'     equality/inequality constraint system. \code{NULL} if no constraints
#'     are active}.
#' }
#' Note that the exact components returned depend heavily on the function
#' arguments (e.g., values of \code{return_G}, \code{return_varcovmat}, etc.)
#' and whether \code{dummy_fit = TRUE}.
#'
#' @usage
#' lgspline.fit(predictors, y = NULL, standardize_response = TRUE,
#'              standardize_predictors_for_knots = TRUE,
#'              standardize_expansions_for_fitting = TRUE, family = gaussian(),
#'              glm_weight_function, schur_correction_function,
#'              need_dispersion_for_estimation = FALSE,
#'              dispersion_function,
#'              K = NULL, custom_knots = NULL, cluster_on_indicators = FALSE,
#'              make_partition_list = NULL, previously_tuned_penalties = NULL,
#'              smoothing_spline_penalty = NULL, opt = TRUE, use_custom_bfgs = TRUE,
#'              delta = NULL, tol = 10*sqrt(.Machine$double.eps),
#'              tuning_criterion = "loo", gcv_gamma = 1.4,
#'              initial_wiggle = c(2e-12, 2e-7, 2e-4, 0.2),
#'              initial_flat = c(0.5, 5), wiggle_penalty = 2e-07,
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
#'              iterate_final_fit = TRUE, blockfit = TRUE,
#'              qp_score_function,
#'              qp_observations = NULL, qp_Amat = NULL, qp_bvec = NULL, qp_meq = 0,
#'              qp_positive_derivative = FALSE, qp_negative_derivative = FALSE,
#'              qp_positive_2ndderivative = FALSE, qp_negative_2ndderivative = FALSE,
#'              qp_monotonic_increase = FALSE, qp_monotonic_decrease = FALSE,
#'              qp_range_upper = NULL, qp_range_lower = NULL, qp_Amat_fxn = NULL,
#'              qp_bvec_fxn = NULL, qp_meq_fxn = NULL, constraint_values = cbind(),
#'              constraint_vectors = cbind(), return_G = TRUE, return_Ghalf = TRUE,
#'              return_U = TRUE, estimate_dispersion = TRUE,
#'              unbias_dispersion = TRUE,
#'              return_varcovmat = TRUE, exact_varcovmat = FALSE,
#'              return_lagrange_multipliers = FALSE,
#'              custom_penalty_mat = NULL,
#'              cluster_args = c(custom_centers = NA, nstart = 10),
#'              dummy_dividor = 1.2345672152894e-22,
#'              dummy_adder = 2.234567210529e-18,
#'              verbose = FALSE, verbose_tune = FALSE,
#'              dummy_fit = FALSE, auto_encode_factors = TRUE,
#'              observation_weights = NULL, do_not_cluster_on_these = c(),
#'              neighbor_tolerance = 1 + 1e-16, no_intercept = FALSE,
#'              VhalfInv = NULL, Vhalf = NULL, include_warnings = TRUE,
#'              og_cols = NULL,
#'              factor_groups = NULL, ...)
#'
#' @inheritParams lgspline
#' @param predictors Numeric matrix or data frame of predictor variables on the
#'   low-level input scale expected by \code{lgspline.fit}. Unlike
#'   \code{\link{lgspline}}, this interface does not parse formulas or a
#'   separate \code{data} argument.
#' @param include_quartic_terms Default: FALSE. Logical switch to include quartic
#'   predictor terms at this low-level interface.
#' @param parallel_unconstrained Default: FALSE. Logical flag for parallel
#'   unconstrained MLE for non-identity-link-Gaussian models.
#' @param unbias_dispersion Default: TRUE. Logical switch to multiply
#'   dispersion by \eqn{N/(N - \mathrm{trace}(\mathbf{H}))}. Unlike
#'   \code{\link{lgspline}}, no wrapper-level auto-resolution is performed here.
#' @param auto_encode_factors Default: TRUE. Compatibility flag carried through
#'   from higher-level preprocessing. Direct calls to \code{lgspline.fit}
#'   should usually supply already encoded predictors and use
#'   \code{factor_groups} when sum-to-zero constraints are needed.
#' @param neighbor_tolerance Default: \code{1 + 1e-16}. Numeric tolerance for
#'   determining neighboring partitions using k-means clustering. Intended for
#'   internal use.
#' @param og_cols Default: NULL. Original predictor names
#' @param factor_groups Named list mapping original factor variable names to
#' integer vectors of their corresponding one-hot indicator column positions
#' within the predictor matrix. Each element enforces a sum-to-zero equality
#' constraint on the linear-term coefficients of its indicator columns within
#' every partition, ensuring identifiability when all factor levels are
#' included without a reference/dropped level. For a group with indicator
#' columns at positions \code{j1, j2, ..., jm}, the constraint is
#' \eqn{\sum_{i=1}^{m} \beta_{ji,k} = 0} for each partition \eqn{k}.
#' Groups with fewer than two resolved positions are silently ignored.
#' Populated automatically by \code{\link{process_input}} when
#' \code{auto_encode_factors = TRUE}; users calling \code{lgspline.fit}
#' directly should construct this list manually when passing one-hot encoded
#' predictors without a reference level. Default \code{NULL} (no
#' sum-to-zero constraints imposed).
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
                         dispersion_function = function(mu, y, order_indices, family,
                                                        observation_weights, VhalfInv,
                                                        ...) {

                           ## If covariance present
                           if(!is.null(VhalfInv)){
                             VhalfInv <- VhalfInv[order_indices, order_indices]
                             c(mean(
                               (
                                 tcrossprod(VhalfInv, t(y-mu))
                               )^2 /
                                 family$variance(mu)
                             ))
                             ## If no covariance present
                           } else{
                             c(mean(
                               (
                                 y - mu
                               )^2 /
                                 family$variance(mu)
                             ))
                           }
                         },
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
                         tuning_criterion = "loo",
                         gcv_gamma = 1.4,
                         initial_wiggle = c(2e-12, 2e-7, 2e-4, 0.2),
                         initial_flat = c(0.5, 5),
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
                         blockfit = TRUE,
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
                         exact_varcovmat = FALSE,
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
                         factor_groups = NULL,
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
  q_predictors <- ncol(predictors)
  N_obs <- nrow(predictors)

  ## [Change 2026-03-05] Resolve character just_linear_* and
  #  exclude_interactions_for arguments to integer column indices.
  #  process_input handles this on the formula path; this block covers
  #  the direct-call path where the user passes a predictor matrix or
  #  data.frame with named columns and character selector arguments.
  #  Resolution uses colnames(predictors); if no names are available and
  #  character entries remain, a warning is issued and they are dropped
  #  (integer entries are always kept).
  pred_colnames <- colnames(predictors)

  .resolve_to_indices <- function(arg, pred_colnames, arg_name,
                                  include_warnings) {
    if (is.null(arg) || length(arg) == 0)        return(arg)
    if (!any(is.character(arg)))                  return(arg)
    if (is.null(pred_colnames)) {
      if (include_warnings) {
        warning("\n \t Character entries in ", arg_name,
                " cannot be resolved: predictors matrix has no column names.",
                " Use integer column indices instead.",
                " Character entries will be dropped.\n")
      }
      ## Keep any integer entries that were already present
      return(arg[!is.character(arg)])
    }
    resolved <- unlist(lapply(arg, function(v) {
      if (is.character(v)) {
        idx <- which(pred_colnames == v)
        if (length(idx) == 0L) idx <- grep(v, pred_colnames, fixed = TRUE)
        idx
      } else {
        v
      }
    }))
    unique(as.integer(resolved))
  }

  just_linear_with_interactions <-
    .resolve_to_indices(just_linear_with_interactions,
                        pred_colnames,
                        "just_linear_with_interactions",
                        include_warnings)

  just_linear_without_interactions <-
    .resolve_to_indices(just_linear_without_interactions,
                        pred_colnames,
                        "just_linear_without_interactions",
                        include_warnings)

  exclude_interactions_for <-
    .resolve_to_indices(exclude_interactions_for,
                        pred_colnames,
                        "exclude_interactions_for",
                        include_warnings)

  ## Return error message if any terms are > q_predictors
  vecdummy <- c(1,
                just_linear_with_interactions,
                just_linear_without_interactions,
                exclude_interactions_for)
  if(any(
    c(1,
      just_linear_with_interactions,
      just_linear_without_interactions,
      exclude_interactions_for) > q_predictors
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
  numerics <- 1:q_predictors

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
  if(length(exclude_interactions_for) >= (q_predictors - 1)){
    include_2way_interactions <- FALSE
    include_3way_interactions <- FALSE
    include_quadratic_interactions <- FALSE
  }
  ## With only two terms available for interactions, exclude 3-ways
  if(length(exclude_interactions_for) >= (q_predictors - 2)){
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

  ## Number of basis expansions per-partition (p_expansions)
  p_expansions <- ncol(C)

  ## In 1-D, set K = number of constraints given by custom knots
  # if not null
  if(any(!(is.null(custom_knots))) & q_predictors == 1){
    if(is.null(K)){
      K <- nrow(custom_knots)
    }
  }

  ## Default K
  orig_null <- FALSE
  if(is.null(K)){
    orig_null <- TRUE
    K <- round(max(min(24/(1 +
                             1*(q_predictors > 1) +
                             1*(paste0(family)[1] != 'gaussian' |
                                  paste0(family)[2] != 'identity')),
                       N_obs/p_expansions),
                   0)/(1 +
                         1*(q_predictors > 1) +
                         1*(paste0(family)[1] != 'gaussian' |
                              paste0(family)[2] != 'identity')))
  }
  if(K == 0){
    unique_penalty_per_partition <- FALSE
  }

  ## Catch error where we need to cluster on some variables, but we have none
  # allowed
  if(length(do_not_cluster_on_these) >= q_predictors &
     length(q_predictors) > 1 &
     K > 0){
    stop('\n \t Must include at least 1 variable to cluster on if multiple',
         ' variables are present. \n')
  }

  ## K can't be greater than max of number of observations or q_predictors
  # for kmeans clustering purposes
  if(K >= max(c(N_obs, q_predictors))) {
    if(include_warnings){
      warning('\n \t Max (N, q_predictors) too small for K. K = max(N, q_predictors) - 2 will be',
              ' used. \n')
    }
    K <- max(max(c(N_obs, q_predictors)) - 2, 0)
  }

  ## Default chunking parameters (overridden if parallel cluster detected)
  chunk_size <- K + 1L
  num_chunks <- 1L
  rem_chunks <- 0L

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

  ## Keep knot-scale transforms for backward compatibility
  #  These are still exposed in the return list as knot_scale_transf /
  #  knot_scale_inv_transf, but predictors are no longer transformed
  #  in-place here. Standardization for clustering now happens inside
  #  make_partitions().
  if(standardize_predictors_for_knots){
    if(q_predictors == 1){
      minns <- apply(predictors, 2, min)
      maxxs <- apply(predictors, 2, max)
    } else {
      means <- apply(predictors, 2, mean)
      sds   <- apply(predictors, 2, function(x) tryCatch(sd(x), error = function(err) 1))
    }
  } else {
    if(q_predictors == 1){
      maxxs       <- rep(1, q_predictors)
      minns       <- rep(0, q_predictors)
      dummy_adder  <- 0
      dummy_dividor <- 0
    } else {
      means        <- rep(0, q_predictors)
      sds          <- rep(1, q_predictors)
      dummy_adder  <- 0
      dummy_dividor <- 0
    }
  }

  ## Transform function for cardinal knot placement (kept for return list / compat)
  transf <- function(X) {
    if(q_predictors == 1){
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

  ## Inverse transform function for cardinal knot placement (kept for return list / compat)
  inv_transf <- function(Xsc) {
    if(q_predictors == 1){
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

  ## Unified partitioning (1-D and multi-D) --------------------------------
  #  Previously, 1-D used equally-spaced quantile knots on the standardized
  #  scale while multi-D used kmeans.  Now both paths use make_partitions()
  #  with internal standardization, which:
  #    - clusters on the standardized scale for numerical stability
  #    - returns centers, knots, and assign_partition on the RAW scale
  #    - uses the integer partition-code scheme (k - 0.5) throughout
  #  This eliminates the need to call inv_transf on knot_values and to call
  #  transf / inv_transf inside the predict function.

  if(!any(is.null(make_partition_list))){
    ## User supplied a pre-computed partition list (e.g., from a dummy_fit).
    #  Assumed to already have raw-scale centers/knots and a raw-accepting
    #  assign_partition (true for objects returned by the new make_partitions).
    partitions <- make_partition_list

  } else if(K > 0 && !is.null(custom_knots) && q_predictors == 1){
    ## 1-D custom knots: user specifies exact breakpoint locations on the raw
    #  scale.  Build a lightweight partition list using those breakpoints
    #  directly, bypassing kmeans.
    knot_mat <- as.matrix(custom_knots)
    breaks   <- sort(knot_mat[, 1])

    ## Synthetic centers: midpoints of each interval (including ?????????Inf tails)
    all_bounds  <- c(-Inf, breaks, Inf)
    center_vals <- sapply(seq_len(length(all_bounds) - 1), function(i){
      lo <- all_bounds[i];   hi <- all_bounds[i + 1]
      if(is.infinite(lo))  return(hi - abs(hi) - 1)
      if(is.infinite(hi))  return(lo + abs(lo) + 1)
      (lo + hi) / 2
    })
    custom_centers_mat <- matrix(center_vals, ncol = q_predictors)

    ## Simple interval-based assign_partition (faster than kNN for 1-D)
    partitions <- list(
      centers = custom_centers_mat,
      knots   = knot_mat,
      assign_partition = function(new_data){
        vals        <- rowMeans(cbind(new_data))
        assignments <- findInterval(vals, breaks) + 1L
        assignments <- pmin(assignments, K + 1L)
        assignments - 0.5
      },
      neighbors = lapply(seq_len(K + 1), function(i){
        nb <- integer(0)
        if(i > 1)     nb <- c(nb, i - 1L)
        if(i < K + 1) nb <- c(nb, i + 1L)
        nb
      }),
      standardize_transf     = function(X) X,
      standardize_inv_transf = function(X) X,
      centers_std            = custom_centers_mat
    )
    rownames(partitions$centers) <- paste0("center_", seq_len(nrow(partitions$centers)))

  } else if(K > 0){
    ## General path: kmeans via make_partitions (handles both 1-D and multi-D).
    #  Predictors are passed on the RAW scale; make_partitions standardizes
    #  internally and returns raw-scale centers / knots.
    partitions <- make_partitions(
      predictors,
      cluster_args,
      cluster_on_indicators,
      K,
      parallel & parallel_find_neighbors,
      cl,
      do_not_cluster_on_these,
      neighbor_tolerance,
      standardize      = standardize_predictors_for_knots,
      standardize_mode = "auto",
      dummy_adder      = dummy_adder,
      dummy_dividor    = dummy_dividor
    )
  } else {
    ## K == 0: no knots, single partition.
    partitions <- NULL
  }

  ## Extract knot values (RAW scale) and establish the integer partition-code
  #  scheme used by knot_expand_list.
  if(K > 0){
    knot_values      <- partitions$knots                           # raw scale
    partition_codes  <- partitions$assign_partition(predictors)    # raw predictors accepted
    partition_bounds <- seq_len(nrow(partitions$centers))          # integer: 1, 2, ..., K+1
  } else {
    knot_values      <- matrix(nrow = 0, ncol = q_predictors)
    partition_codes  <- rep(0.5, N_obs)
    partition_bounds <- c()
  }

  ## NOTE: predictors remain on the raw scale from this point on.
  #  The inv_transf() call that previously appeared here is no longer needed.

  if(verbose){
    cat("Expansion Standardize\n")
  }

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
  # This is a p_expansions-1 length vector, it excludes the intercept
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

  ## Sum-to-zero constraints for one-hot encoded factors
  #  Each factor group contributes one equality constraint per partition,
  #  enforcing that the indicator coefficients sum to 0 within that
  #  partition. This keeps the fully encoded factor effects identifiable.
  if(!is.null(factor_groups) && length(factor_groups) > 0){
    colnm_expansions_temp <- colnames(C)
    factor_constraint_list <- list()

    for(fg_name in names(factor_groups)){
      fg_cols <- factor_groups[[fg_name]]

      ## Find the linear expansion positions for these indicator columns.
      fg_expansion_positions <- which(
        colnm_expansions_temp %in% paste0("_", fg_cols, "_")
      )

      ## Need at least 2 indicator columns for a meaningful sum-to-zero constraint.
      if(length(fg_expansion_positions) > 1){
        ## Build one constraint column per partition.
        fg_constr <- sapply(1:(K+1), function(k){
          vec <- rep(0, p_expansions * (K + 1))
          vec[p_expansions * (k - 1) + fg_expansion_positions] <- 1
          vec
        })
        factor_constraint_list[[fg_name]] <- fg_constr

        if(verbose){
          cat(sprintf(
            "Sum-to-zero constraint for factor '%s': %d indicators\n",
            fg_name, length(fg_expansion_positions)
          ))
        }
      }
    }

    ## Combine all factor constraints and append to constraint_vectors.
    #  constraint_values for these new constraints are all 0 (sum = 0).
    if(length(factor_constraint_list) > 0){
      all_factor_constr <- do.call(cbind, factor_constraint_list)
      zero_vals <- matrix(0,
                          nrow = nrow(all_factor_constr),
                          ncol = ncol(all_factor_constr))

      if(length(constraint_vectors) > 0 &&
         is.matrix(constraint_vectors) &&
         nrow(constraint_vectors) == nrow(all_factor_constr)){
        ## Existing user constraints present: append columns
        constraint_vectors <- cbind(constraint_vectors, all_factor_constr)
        constraint_values  <- cbind(constraint_values,  zero_vals)
      } else if(length(constraint_vectors) == 0 ||
                !is.matrix(constraint_vectors)){
        ## No prior constraints: initialize from factor constraints only
        constraint_vectors <- all_factor_constr
        constraint_values  <- zero_vals
      }
      ## If dimensions do not match (unusual), skip silently to avoid error
    }
  }

  ## If no intercept enforced, include constraint on A indicating this
  if(no_intercept & length(constraint_vectors) < 1){
    constr <- sapply(1:(K+1), function(k){
      vec <- rep(0, p_expansions*(K+1))
      vec[p_expansions*(k-1) + 1] <- 1
      vec
    })
    constraint_vectors <- cbind(constr)
    constraint_values <- 0*constraint_vectors
  } else if(no_intercept){
    constr <- sapply(1:(K+1), function(k){
      vec <- rep(0, p_expansions*(K+1))
      vec[p_expansions*(k-1) + 1] <- 1
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
        vec <- rep(0, p_expansions*(K+1))
        vec[p_expansions*(k-1) + offset_ind[o]] <- 1
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
        vec <- rep(0, p_expansions*(K+1))
        vec[p_expansions*(k-1) + offset_ind[o]] <- 1
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
                        N_obs,
                        C,
                        K)

  ## Assign y to their partitions (y_og saves original y unstandardized)
  y <- knot_expand_list(partition_codes,
                        partition_bounds,
                        N_obs,
                        cbind(y),
                        K)

  ## Get observation weight expansions
  if(any(!is.null(observation_weights))){
    ## Coerce to N x 1 vector if not already
    if(nrow(cbind(observation_weights)) != N_obs |
       ncol(cbind(observation_weights)) != 1){
      stop('\n \t Observation weights must be an N x 1 vector. \n')
    }
    observation_weights_og <- observation_weights
    homogenous_weights <- (length(unique(observation_weights_og)) == 1)
    observation_weights <-
      knot_expand_list(partition_codes,
                       partition_bounds,
                       N_obs,
                       cbind(observation_weights),
                       K)
  } else {
    observation_weights_og <- rep(1, N_obs)
    observation_weights <- lapply(1:(K+1), function(k)cbind(
      rep(1, length(y[[k]]))))
    homogenous_weights <- TRUE
  }

  ## Save the original ordering to each partition
  order_list <- knot_expand_list(partition_codes,
                                 partition_bounds,
                                 N_obs,
                                 cbind(1:N_obs),
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
    } else if(any(unique(dim(VhalfInv)) != N_obs)){
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
        p_expansions,
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
                                  p_expansions,
                                  parallel & parallel_penalty,
                                  cl)
    colnames(smoothing_spline_penalty) <- colnames(C)
  }

  if(verbose){
    cat("Constraint Matrix\n")
  }

  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  ## Constraint matrix A construction.
  #  Two branches:
  #  (1) multi-predictor / nonspline present, (2) K = 0.
  #  All branches standardize A by expansion_scales before returning.
  #
  #  knot_values are now on the RAW scale (returned directly by
  #  make_partitions or constructed from raw-scale custom_knots).
  #  All inv_transf() calls that previously converted standardized knots
  #  back to raw scale have been removed; knot_values_chunk and knot_values
  #  are passed directly to get_polynomial_expansions().
  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

  ## Making a constraint matrix
  A <- 0
  if(K > 0){

    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    ## Branch 1: multi-predictor or nonspline columns present.
    #  knot_values may have more rows than K (one row per knot per
    #  predictor dimension), so we process in chunks of K rows and
    #  cbind the resulting constraint columns. A remainder chunk handles
    #  any leftover rows, padded with zeros to keep dimensions consistent.
    #  knot_values are already on raw scale; no inv_transf() needed

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
                                        ## knot_values_chunk already on raw scale
                                        CKnots_chunk <- rbind(
                                          get_polynomial_expansions(
                                            knot_values_chunk,
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
                                          p_expansions,
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
        ## knot_values_chunk already on raw scale
        CKnots_chunk <- rbind(
          get_polynomial_expansions(knot_values_chunk,
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
        A <- cbind(A, make_constraint_matrix(p_expansions,
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

    ## Remove all 0 columns
    A <- A[,which(apply(abs(A), 2, sum) > 1e-16),drop=FALSE]

    ## Append user-supplied equality constraints and standardize,
    #  zeroing offset / no_intercept rows as in Branch 1.
    if(length(constraint_vectors) > 0){
      if(length(offset) > 0){
        offset_ind <- which(colnm_expansions %in% paste0(
          '_', offset, '_'
        ))
        offset_inds <- unlist(lapply(1:(K+1), function(k)p_expansions*(k-1)+offset_ind))
        A[offset_inds,] <- 0
      }
      if(no_intercept){
        A[unlist(lapply(1:(K+1),function(k)p_expansions*(k-1)+1)),] <- 0
      }
      A <- cbind(A, constraint_vectors)
    }
    A <- sweep(A, 1, rep(c(1, expansion_scales), K+1), "/")
    if(any(!is.finite(A))) stop(paste0('\n \t A is not finite \n',
                                       expansion_scales))
    if(any(is.na(A))) stop(paste0('\n \t A is NA somewhere ',
                                  '(any(is.na(A)) == TRUE) \n',
                                  expansion_scales))

    ## Remainder chunk: rows of knot_values not covered by full chunks.
    #  Taken from the tail of knot_values_perm and zero-padded to K rows
    #  so make_constraint_matrix receives a conformable matrix.
    if(rem > 0){
      knot_values_chunk <-
        knot_values_perm[rev(c(nrow(knot_values_perm):1)[1:rem]),,drop=FALSE]
      ## knot_values_chunk already on raw scale
      temp_dat <- knot_values_chunk
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

      A <- cbind(A, make_constraint_matrix(p_expansions,
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

  } else {

    ## Branch 2: K = 0, no knots.
    #  No smoothness constraints exist; A is either the user-supplied
    #  constraint_vectors alone or NULL (unconstrained).

    if(length(constraint_vectors) > 0){
      A <- cbind(constraint_vectors)
      A <- sweep(A, 1, rep(c(1, expansion_scales), K+1), "/")
    } else {
      A <- NULL
    }
  }
  ## Select only a linearly-independent subset of columns
  if(!(any(is.null(A)))){
    qr_A <- qr(A)
    A_rank <- qr_A$rank
    if(A_rank < ncol(A)){
      A <- qr.Q(qr_A)[, 1:A_rank, drop = FALSE]
    }
    R_constraints <- ncol(A)
  }

  ## With only one predictor, we really only need one penalty
  if(ncol(predictors) == 1){
    unique_penalty_per_predictor <- FALSE
  }

  tuning_criterion <- match.arg(tuning_criterion, c("loo", "gcv"))

  ## Modified GCV multiplier for tuning
  if(tuning_criterion == "gcv" &&
     (!is.numeric(gcv_gamma) || length(gcv_gamma) != 1 ||
      !is.finite(gcv_gamma) || gcv_gamma < 1)){
    stop("\n \t gcv_gamma must be a finite numeric scalar >= 1. ",
         "Set gcv_gamma = 1 to recover ordinary GCV. \n")
  }

  ## Convert constraint_values, penalty setup,
  #  parallel export, X^TWX gram, SQP setup (quadprog).
  if(verbose){
    cat("Predictor-and-Partition Penalty Setup\n")
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
  shared_env <- NULL
  if(parallel && !is.null(cl)) {

    shared_vars <- list(
      A = A,
      R_constraints = if(is.null(A)) 0L else ncol(A),
      K = K,
      p_expansions = p_expansions,
      N_obs = N_obs,
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
      observation_weights = observation_weights,
      efficient_matrix_mult = efficient_matrix_mult
    )

    export_env <- list2env(shared_vars, parent = emptyenv())

    tryCatch({
      parallel::clusterExport(cl, names(shared_vars), envir = export_env)
      parallel::clusterEvalQ(cl, {
        `%**%` <- efficient_matrix_mult
        NULL
      })
    }, error = function(e) {
      stop("Failed to export parallel worker state to cluster: ",
           e$message, call. = FALSE)
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

  ## [Change 2026-03-06] Delegate QP setup to process_qp()
  qp_result <- process_qp(
    X = X,
    K = K,
    p_expansions = p_expansions,
    order_list = order_list,
    colnm_expansions = colnm_expansions,
    expansion_scales = expansion_scales,
    power1_cols = power1_cols,
    power2_cols = power2_cols,
    nonspline_cols = nonspline_cols,
    interaction_single_cols = interaction_single_cols,
    interaction_quad_cols = interaction_quad_cols,
    triplet_cols = triplet_cols,
    include_2way_interactions = include_2way_interactions,
    include_3way_interactions = include_3way_interactions,
    include_quadratic_interactions = include_quadratic_interactions,
    family = family,
    mean_y = mean_y,
    sd_y = sd_y,
    N_obs = N_obs,
    qp_observations = qp_observations,
    qp_positive_derivative = qp_positive_derivative,
    qp_negative_derivative = qp_negative_derivative,
    qp_positive_2ndderivative = qp_positive_2ndderivative,
    qp_negative_2ndderivative = qp_negative_2ndderivative,
    qp_monotonic_increase = qp_monotonic_increase,
    qp_monotonic_decrease = qp_monotonic_decrease,
    qp_range_upper = qp_range_upper,
    qp_range_lower = qp_range_lower,
    qp_Amat_fxn = qp_Amat_fxn,
    qp_bvec_fxn = qp_bvec_fxn,
    qp_meq_fxn = qp_meq_fxn,
    qp_Amat = qp_Amat,
    qp_bvec = qp_bvec,
    qp_meq = qp_meq,
    all_derivatives_fxn = all_derivatives,
    og_cols = og_cols,
    include_warnings = include_warnings,
    ...
  )
  quadprog  <- qp_result$quadprog
  qp_Amat  <- qp_result$qp_Amat
  qp_bvec  <- qp_result$qp_bvec
  qp_meq   <- qp_result$qp_meq

  ## Convert constraint_values to per-partition list
  #  Convert non-0 null vectors to (K+1) list of corresponding partitions
  if(length(constraint_values) > 0){

    constraint_values <- lapply(1:(K+1),function(k){
      vec <- cbind(constraint_values)[1:p_expansions + (k-1)*p_expansions,,drop=FALSE]
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
        y = y,
        X = X,
        X_gram = X_gram,
        smoothing_spline_penalty = smoothing_spline_penalty,
        A = A,
        K = K,
        p_expansions = p_expansions,
        N_obs = N_obs,
        opt = opt,
        use_custom_bfgs = use_custom_bfgs,
        C = C,
        colnm_expansions = colnm_expansions,
        wiggle_penalty = wiggle_penalty,
        flat_ridge_penalty = flat_ridge_penalty,
        initial_wiggle = initial_wiggle,
        initial_flat = initial_flat,
        unique_penalty_per_predictor = unique_penalty_per_predictor,
        unique_penalty_per_partition = unique_penalty_per_partition,
        penalty_vec = penalty_vec,
        meta_penalty = meta_penalty,
        family = family,
        unconstrained_fit_fxn = unconstrained_fit_fxn,
        keep_weighted_Lambda = keep_weighted_Lambda,
        iterate = iterate_tune,
        qp_score_function = qp_score_function,
        quadprog = quadprog,
        qp_Amat = qp_Amat,
        qp_bvec = qp_bvec,
        qp_meq = qp_meq,
        tol = tol,
        sd_y = sd_y,
        delta = delta,
        tuning_criterion = tuning_criterion,
        gcv_gamma = gcv_gamma,
        constraint_value_vectors = constraint_values,
        parallel = parallel,
        parallel_eigen = parallel_eigen,
        parallel_trace = parallel_trace,
        parallel_aga = parallel_aga,
        parallel_matmult = parallel_matmult,
        parallel_unconstrained = parallel_unconstrained,
        cl = cl,
        chunk_size = chunk_size,
        num_chunks = num_chunks,
        rem_chunks = rem_chunks,
        shared_env = shared_env,
        custom_penalty_mat = custom_penalty_mat,
        order_list = order_list,
        glm_weight_function = glm_weight_function,
        schur_correction_function = schur_correction_function,
        need_dispersion_for_estimation = need_dispersion_for_estimation,
        dispersion_function = dispersion_function,
        observation_weights = observation_weights,
        homogenous_weights = homogenous_weights,
        blockfit = blockfit,
        just_linear_without_interactions = just_linear_without_interactions,
        Vhalf = Vhalf,
        VhalfInv = VhalfInv,
        verbose = verbose_tune,
        include_warnings = include_warnings,
        ...)}, silent = TRUE)
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

    ## Unified assign_partition: always from partitions object (raw scale).
    ## K = 0 gets the constant function; K > 0 gets the kmeans closure.
    if(K > 0 && !is.null(partitions)){
      assign_partition <- partitions$assign_partition
    } else {
      assign_partition <- function(x) 0.5
    }

    B <- lapply(1:(K+1), function(k){
      b <- cbind(rep(0, p_expansions))
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
      ## knot_values already on raw scale
      knots <- knot_values
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
      y = y_og, ytilde = rep(0, N_obs), X = X_out, A = A, B = B,
      B_raw = B_raw, K = K, p = p_expansions, q = ncol(predictors),
      P = (K+1)*p_expansions, N = N_obs, penalties = tL,
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
      og_cols = og_cols,
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

  ## This replaces the original get_B call with conditional call
  #  to blockfit_solve or get_B depending on model configuration
  if(verbose){
    cat("Prep for final fitting\n")
  }

  ## Final fit
  if(K == 0){
    ## ensuring compatibility with no A
    if(any(is.null(A))){
      A <- cbind(rep(0, (K+1)*p_expansions))
      A <- cbind(A, A)
      R_constraints <- 2
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
        p_expansions = p_expansions,
        Lambda = tL$Lambda,
        L_partition_list = tL$L_partition_list,
        unique_penalty_per_partition = unique_penalty_per_partition,
        A = A,
        R_constraints = R_constraints,
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
        p_expansions,
        R_constraints,
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

  ## Unified assign_partition: always from partitions object (raw scale).
  ## K = 0 gets the constant function; K > 0 gets the kmeans closure from
  ## make_partitions, which accepts raw-scale data and standardizes internally.
  if(K > 0 && !is.null(partitions)){
    assign_partition <- partitions$assign_partition
  } else {
    assign_partition <- function(x) 0.5
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
  .shared_predict <- new.env(parent = emptyenv())
  .shared_predict$varcovmat <- NULL
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

    ## [Change 2026-03-04] se.fit: if TRUE, return standard errors and
    #  confidence intervals for predictions.
    #
    #  Derivation (identity link, Gaussian case):
    #    eta_i = x_i' beta   (linear predictor at new point i)
    #    Var(eta_i) = x_i' Var(beta) x_i = x_i' varcovmat x_i
    #
    #  varcovmat is on the unstandardized coefficient / original-y scale
    #  (sigmasq_tilde * UGU', un-standardized by expansion_scales).
    #
    #  Confidence intervals are constructed on the link scale:
    #    eta_i +/- cv * se(eta_i)
    #  then back-transformed via linkinv() to respect response constraints.
    #
    #  The x_i vector is the UNSTANDARDIZED expansion row for observation i,
    #  placed into its partition's block of the full P-dimensional coefficient
    #  vector.  Row assignment uses the NEW data's partition codes, not
    #  the training data's.
    #
    #  Requirements: return_varcovmat = TRUE at fit time.
    #
    #  The varcovmat is accessed via .shared_predict$varcovmat, a shared
    #  environment that is populated after the varcovmat is computed in
    #  lgspline.fit's return path.

    ## Validate se.fit requirements
    if(se.fit){
      if(!return_varcovmat || is.null(.shared_predict$varcovmat)){
        if(include_warnings){
          warning('\n\t se.fit = TRUE requires return_varcovmat = TRUE ',
                  'at model fitting time. Standard errors will not be ',
                  'computed.\n')
        }
        se.fit <- FALSE
      }
    }

    ## Encoding / coercion block
    if(any(!is.null(new_predictors))){

      ## Filter extraneous columns from new_predictors
      if(!is.null(og_cols) &&
         (inherits(new_predictors, "data.frame") ||
          (is.matrix(new_predictors) &&
           !is.null(colnames(new_predictors))))){

        np_colnames <- colnames(as.data.frame(new_predictors))
        relevant_cols <- np_colnames[sapply(np_colnames, function(cn){
          if(cn %in% og_cols) return(TRUE)
          onehot_pattern <- paste0("^", cn, "_")
          if(any(grepl(onehot_pattern, og_cols))) return(TRUE)
          return(FALSE)
        })]

        if(length(relevant_cols) > 0 &&
           length(relevant_cols) < length(np_colnames)){
          new_predictors <- as.data.frame(new_predictors)[,
                                                          relevant_cols, drop = FALSE]
        }
      }

      if(!is.null(og_cols) &&
         (inherits(new_predictors, "data.frame") ||
          (is.matrix(new_predictors) &&
           !is.null(colnames(new_predictors))))){

        np_df <- as.data.frame(new_predictors)
        col_names <- colnames(np_df)

        cols_needing_encoding <- character(0)
        for(cn in col_names){
          if(cn %in% og_cols) next
          onehot_pattern <- paste0("^", cn, "_")
          if(any(grepl(onehot_pattern, og_cols))){
            cols_needing_encoding <- c(cols_needing_encoding, cn)
          }
        }

        if(length(cols_needing_encoding) > 0){
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
              if(unrecognized!=0){
                warning("\n\t predict: value(s) [",
                        paste(unrecognized, collapse = ", "),
                        "] in column '", col_name,
                        "' were not seen during fitting and will be treated as ",
                        "all-zero (reference level).\n")
              }
            }

            col_idx <- which(colnames(np_df) == col_name)
            np_df <- cbind(np_df[, -col_idx, drop = FALSE],
                           as.data.frame(onehot_mat))
          }

          missing_cols <- setdiff(og_cols, colnames(np_df))
          if(length(missing_cols) > 0){
            np_df[missing_cols] <- 0L
          }
          if(all(og_cols %in% colnames(np_df))){
            np_df <- np_df[, og_cols, drop = FALSE]
          }

          new_predictors <- as.matrix(np_df)

        } else if(all(col_names %in% og_cols)){
          np_df <- np_df[, og_cols, drop = FALSE]
          new_predictors <- as.matrix(np_df)

        } else {
          ## Spurious colnames from cbind(); strip and proceed
        }
      }

      new_predictors <- try(
        methods::as(cbind(new_predictors), 'matrix'),
        silent = TRUE)
      if(any(inherits(new_predictors, 'try-error'))){
        stop('\n \t New predictors should be able to be coerced into ',
             'matrix form. If you are passing factor/character columns, ',
             'make sure the column names match the original variable ',
             'names used at fit time. \n')
      }
      colnames(new_predictors) <- NULL
    }

    ## Avoid R rbind issue with 1 row
    if(nrow(new_predictors) == 1){
      new_predictors <- rbind(new_predictors, new_predictors)
      only_1 <- TRUE
    } else {
      only_1 <- FALSE
    }

    ## Partition assignment now works on raw-scale predictors.
    #  assign_partition handles any needed standardization internally.
    partition_codes_new <- assign_partition(new_predictors)

    ## Polynomial expansions (raw / unstandardized)
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

    ## Knot expansions ?????????????????? C_new is UNSTANDARDIZED (matches original code).
    ## B_predict is fully backtransformed to the same scale.
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

    ## Re-order predictions: map from partition-blocked to original order
    order_list_new <- knot_expand_list(
      partition_codes_new,
      partition_bounds,
      length(partition_codes_new),
      cbind(1:nrow(C_new)),
      K)

    ## Only use non-empty partitions
    keep_blocks <- which(sapply(1:(K+1),function(k){
      nrow(X_new[[k]]) > 0
    }))
    order_list_keep <- order_list_new[keep_blocks]

    ## Predictions on the link scale
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
          rem_chunks))[order(unlist(order_list_keep))]

    if(only_1){
      preds <- preds[1]
    }
    final_preds <- family$linkinv(preds)

    ## [Change 2026-03-04] Compute standard errors and confidence intervals.
    #
    #  For each new observation i in partition k:
    #    x_i^{full} = (0,...,0, x_i^{(k)}, 0,...,0)  (P-dimensional)
    #    Var(eta_i) = (x_i^{full})' varcovmat (x_i^{full})
    #
    #  CRITICAL: use order_list_new (new-data partition assignments),
    #  NOT the closure's order_list (training data).
    se_link <- NULL
    lower <- NULL
    upper <- NULL

    if(se.fit){
      vcov <- .shared_predict$varcovmat
      P_total <- p_expansions * (K + 1)
      N_new <- nrow(C_new)

      ## Build full N_new x P block-diagonal design matrix for new data.
      X_full_new <- matrix(0, nrow = N_new, ncol = P_total)
      for(k in keep_blocks){
        if(nrow(X_new[[k]]) == 0) next
        ## order_list_new[[k]]: original row indices (in 1:N_new) for partition k
        row_inds <- order_list_new[[k]]
        col_start <- (k - 1) * p_expansions + 1
        col_end   <- k * p_expansions
        X_full_new[row_inds, col_start:col_end] <- X_new[[k]]
      }

      ## Quadratic form: Var(eta_i) = x_i' V x_i
      XV <- X_full_new %**% vcov
      var_eta <- rowSums(XV * X_full_new)

      ## Guard against numerical negatives
      var_eta <- pmax(var_eta, 0)
      se_link_vec <- sqrt(var_eta)

      if(only_1){
        se_link_vec <- se_link_vec[1]
        preds_for_ci <- preds[1]
      } else {
        preds_for_ci <- preds
      }

      ## CI on link scale, back-transformed to response scale
      eta_lower <- preds_for_ci - cv * se_link_vec
      eta_upper <- preds_for_ci + cv * se_link_vec
      lower <- family$linkinv(eta_lower)
      upper <- family$linkinv(eta_upper)
      se_link <- se_link_vec
    }

    ## If returning derivatives
    if(take_first_derivatives | take_second_derivatives){
      derivs <- make_derivative_matrix(
        p_expansions, C_new, power1_cols, power2_cols, nonspline_cols,
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

      preds_prime <-
        unlist(
          matmult_block_diagonal(
            Xprime_new[keep_blocks], B_predict[keep_blocks],
            length(keep_blocks) - 1, parallel, cl,
            chunk_size, num_chunks, rem_chunks))[order(unlist(order_list_keep))]

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

        preds_dprime <-
          unlist(
            matmult_block_diagonal(
              Xdprime_new[keep_blocks], B_predict[keep_blocks],
              length(keep_blocks) - 1, parallel, cl,
              chunk_size, num_chunks, rem_chunks))[order(unlist(order_list_keep))]

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

      ## [Change 2026-03-05] Remove the last NA dummy observation
      if(only_1 & n_deriv_vars > 1){
        if(take_first_derivatives){
          final_preds_prime <- lapply(final_preds_prime, function(x){
            x[1]
          })
        }
        if(take_second_derivatives){
          final_preds_dprime <- lapply(final_preds_dprime, function(x){
            x[1]
          })
        }
      } else if(only_1){
        if(take_first_derivatives){
          final_preds_prime <- final_preds_prime[1]
        }
        if(take_second_derivatives){
          final_preds_dprime <- final_preds_dprime[1]
        }
      }

      out <- list(
        preds = final_preds,
        first_deriv = final_preds_prime,
        second_deriv = final_preds_dprime
      )

      if(se.fit){
        out$se.fit <- se_link
        out$lower <- lower
        out$upper <- upper
      }

      return(out)

    } else if(se.fit){
      return(list(
        fit = final_preds,
        se.fit = se_link,
        lower = lower,
        upper = upper,
        cv = cv
      ))
    } else {
      return(final_preds)
    }
  }

  ## Get fitted values
  ytilde <- predict_function()

  ## Clean knots on the raw scale.
  if(K == 0){
    knots <- NULL
  } else {
    knots <- knot_values
    if(length(numerics) == 1 & length(nonspline_cols) == 0){
      rownames(knots) <- paste0(1:K, '_', 2:(K+1))
    }
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

  ## return_list construction, optional components
  #  (G, Ghalf, U, Lagrange multipliers, dispersion, trace, varcov),
  #  and final return.

  ## List of items to return
  return_list <- list("y" = y_og,
                      "ytilde" = ytilde,
                      "X" = X,
                      "A" = A,
                      "B" = B,
                      "B_raw" = B_raw,
                      "K" = K,
                      "p" = p_expansions,
                      "q" = ncol(predictors),
                      "P" = (K+1)*p_expansions,
                      "N" = N_obs,
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

  ## We need U to compute sigma^2*UG
  #  Change 2026-02-26: dispersion does not need to be separately estimated
  #  for several GLMs, including Poisson and logistic regression
  if(return_varcovmat){
    return_U <- TRUE
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
      A_for_U          <- B_list$qp_info$Amat_active
      R_constraints_for_U <- ncol(A_for_U)
    } else {
      A_for_U          <- A
      R_constraints_for_U <- R_constraints
    }

    no_constraint_identity_U <- K == 0 &
      length(constraint_values) == 0 &
      !quadprog

    if(no_constraint_identity_U && any(is.null(A))){
      A <- cbind(rep(0, (K+1)*p_expansions))
      A <- cbind(A, A)
      R_constraints <- 2
      A_for_U <- A
      R_constraints_for_U <- R_constraints
    }

    if(no_constraint_identity_U && is.null(return_list$VhalfInv)){
      return_list$U <- diag(p_expansions*(K+1))
    } else if(is.null(return_list$VhalfInv)){
      ## Standard block-diagonal path: use active constraint matrix
      #  so U projects onto the null space of all binding constraints.
      return_list$U <- get_U(
        G_list$G,
        A_for_U,
        K,
        p_expansions,
        R_constraints_for_U
      )
    } else {
      ## GEE / VhalfInv path: build full whitened G_correct then project
      #  using the active constraint matrix.
      X_full_tr <- collapse_block_diagonal(X)

      ## W~ without obs weights to avoid double-counting D
      prelim_disp_tr <- if(need_dispersion_for_estimation){
        dispersion_function(mu = ytilde,
                            y = y_og,
                            order_indices = 1:N_obs,
                            family = family,
                            observation_weights = observation_weights_og,
                            VhalfInv = VhalfInv,
                            ...)
      } else {
        1
      }
      W_glm_tr <- c(glm_weight_function(
        ytilde[unlist(order_list)], y_og[unlist(order_list)],
        1:N_obs, family, prelim_disp_tr, rep(1, N_obs), ...))
      W_glm_tr <- pmax(W_glm_tr, .Machine$double.eps)

      ## Observation weights in partition order
      D_tr <- observation_weights_og[unlist(order_list)]

      ## Apply observation / working weights in the whitened system to
      #  match the fitting paths used by get_B() and blockfit_solve().
      VinvhalfX_tr <- VhalfInv[unlist(order_list),
                               unlist(order_list)] %**%
        X_full_tr
      VinvhalfX_tr <- t(t(VinvhalfX_tr) * sqrt(W_glm_tr * D_tr))

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

      ## G_correct built from the whitened-system penalized Gram.
      #  G_correct is the "correct" G because it isn't just block diagonal
      G_correct_tr <- invert(crossprod(VinvhalfX_tr) + Lambda_full_tr)

      if(no_constraint_identity_U){
        ## Correlated K = 0 fits still need dense GLS quantities for the
        #  downstream trace and REML calculations, even though U is identity.
        return_list$U <- diag(nrow(G_correct_tr))
      } else {
        ## Recompute U with the correct dense G and active constraint matrix.
        #  When quadprog is active, A_for_U includes columns for binding
        #  inequality constraints so the projection accounts for all active
        #  constraints, not just the smoothness equalities.
        GA <- G_correct_tr %**% A_for_U
        return_list$U <- diag(nrow(G_correct_tr)) -
          GA %**% tcrossprod(invert(crossprod(A_for_U, GA)), A_for_U)
      }
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
          G_list$G, A, K, p_expansions, R_constraints,
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
                               A, K, p_expansions, R_constraints,
                               parallel = parallel & parallel_aga,
                               cl = cl,
                               chunk_size, num_chunks, rem_chunks)),
        p_expansions, R_constraints, K,
        parallel & parallel_trace,
        cl = cl,
        chunk_size, num_chunks, rem_chunks,
        unique_penalty_per_partition,
        tL$L_partition_list)

    } else if(!is.null(VhalfInv)){
      if(!return_U){
        warning("\n \t return_U is False, cannot compute trace - using P.\n")
        trace_XUGX <- p_expansions*(K+1)
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
      scale_by <- N_obs/(N_obs - trace_XUGX)
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
    if(estimate_dispersion){
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
    } else {
      return_list$sigmasq_tilde <- 1
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
      W_glm_vc <- c(glm_weight_function(ytilde, y_og, 1:N_obs, family,
                                        return_list$sigmasq_tilde,
                                        rep(1, N_obs), ...))
      W_glm_vc <- pmax(W_glm_vc, .Machine$double.eps)

      ## Combined weight in the whitened system.
      combined_wt_vc <- sqrt(W_glm_vc * observation_weights_og)

      ## Whitened weighted design: D^{1/2} W~^{1/2} V^{-1/2} X  (N x P)
      VinvhalfX <- VhalfInv %**% X_full
      VinvhalfX <- t(t(VinvhalfX) * combined_wt_vc)

      ## Full penalized GLS Gram in the whitened system + Lambda  (P x P)
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
      UGhalf <- matmult_U(return_list$U, G_list$Ghalf, p_expansions, K)
      return_list$varcovmat <-
        tcrossprod(UGhalf) *
        return_list$sigmasq_tilde # contains sd_y^2 already
    }

    ## Un-standardize (both paths)
    d <- rep(c(1, 1/expansion_scales), times = K + 1)
    return_list$varcovmat <- t(t(return_list$varcovmat * d) * d)

    if(any(diag(return_list$varcovmat) < 0) & include_warnings){
      warning("\n \t Variance-covariance matrix has diagonal elements < 0,",
              " model most likely did not converge. Try re-fitting, a simpler",
              " model, different knot locations, or larger penalties. \n")
      for(ij in 1:nrow(return_list$varcovmat)){
        return_list$varcovmat[ij, ij] <- max(0, return_list$varcovmat[ij, ij])
      }
    }

    ## [Change 2026-03-02] Exact frequentist variance-covariance matrix.
    #
    #  The asymptotic (Bayesian posterior) version:
    #    Varcov_asymptotic = sigma~^2 * U G U^T
    #
    #  The exact frequentist version (from derivation in notes):
    #    Var(B_hat) = sigma^2 * U G (X^T X) G U^T
    #    Since G^{-1} = X^T X + Lambda:  X^T X = G^{-1} - Lambda
    #    => Var(B_hat) = sigma^2 * U G (G^{-1} - Lambda) G U^T
    #                  = sigma^2 * U G U^T - sigma^2 * U G Lambda G U^T
    #
    #  Both terms carry sigma~^2. varcovmat already = sigma~^2 * UGU^T,
    #  so:
    #    varcovmat_exact = varcovmat - sigma~^2 * U G Lambda G U^T
    #
    #  For Gaussian identity (with or without correlation), exact.
    #  For other families, asymptotically correct with estimated dispersion.
    #
    #  [Change 2026-03-02] This REPLACES varcovmat in-place (varcovmat_exact
    #  was the prior naming; now the single returned varcovmat is the exact
    #  version when exact_varcovmat = TRUE).
    if(exact_varcovmat){
      if(verbose) cat("Exact Frequentist VarCov Mat \n")

      ## Full block-diagonal penalty (including partition-specific terms)
      has_part_pen_ex <- length(tL$L_partition_list) == (K + 1)
      Lambda_full_ex <- collapse_block_diagonal(
        lapply(1:(K + 1), function(k){
          if(has_part_pen_ex){
            tL$Lambda + tL$L_partition_list[[k]]
          } else {
            tL$Lambda
          }
        })
      )

      ## Un-standardization diagonal: same as used for varcovmat above
      d_ex <- rep(c(1, 1/expansion_scales), times = K + 1)

      if(!is.null(VhalfInv)){
        ## WITH correlation structure:
        #  G_correct_tr and Ghalf_correct are in scope from the varcovmat block.
        #  correction = sigma~^2 * U G_correct Lambda G_correct U^T
        #  Factored stably as: tcrossprod(U %**% G_correct %**% Lambda %**% Ghalf_correct)
        GLG_half_ex   <- G_correct_tr %**% Lambda_full_ex %**% Ghalf_correct
        UGLGhalf_ex   <- return_list$U %**% GLG_half_ex
        correction_ex <- tcrossprod(UGLGhalf_ex)

        ## Apply sigma~^2 to the correction (both terms now carry it)
        correction_ex <- return_list$sigmasq_tilde * correction_ex

        ## Un-standardize
        correction_ex_unstd <- t(t(correction_ex * d_ex) * d_ex)

      } else {
        ## WITHOUT correlation structure:
        #  G is block-diagonal; Ghalf[[k]] available from G_list.
        #  Ghalf_k Lambda_k Ghalf_k^T is the within-partition contribution.
        #  correction = sigma~^2 * U blockdiag(Ghalf_k Lambda_k Ghalf_k^T) U^T
        GhalfLG_block_ex <- lapply(1:(K + 1), function(k){
          Lk <- if(has_part_pen_ex){
            tL$Lambda + tL$L_partition_list[[k]]
          } else {
            tL$Lambda
          }
          ## Numerically stable symmetric form: G^{1/2} Lambda (G^{1/2})^T
          G_list$Ghalf[[k]] %**% Lk %**% t(G_list$Ghalf[[k]])
        })

        GhalfLG_bd_ex <- collapse_block_diagonal(GhalfLG_block_ex)
        correction_ex <- return_list$U %**% GhalfLG_bd_ex %**% t(return_list$U)

        ## Apply sigma~^2 to the correction
        correction_ex <- return_list$sigmasq_tilde * correction_ex

        ## Un-standardize
        correction_ex_unstd <- t(t(correction_ex * d_ex) * d_ex)
      }

      ## Replace varcovmat in-place with the exact version
      return_list$varcovmat <- return_list$varcovmat - correction_ex_unstd

      ## Diagnose: negative diagonal can occur under heavy penalization.
      if(any(diag(return_list$varcovmat) < 0) & include_warnings){
        warning(
          "\n \t Exact frequentist varcovmat has diagonal elements < 0.",
          " This can occur under heavy penalization where the penalty",
          " correction exceeds the asymptotic variance. Consider using",
          " the asymptotic version (exact_varcovmat = FALSE) or reducing",
          " the wiggle_penalty. \n"
        )
      }
    }
  }
  ## Save for predicting confidence bands
  .shared_predict$varcovmat <- return_list$varcovmat

  ## Afterwards, update X to unstandardized
  return_list$X <- lapply(return_list$X, unstd_X)

  return(return_list)
}
