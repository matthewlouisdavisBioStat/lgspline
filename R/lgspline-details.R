#' Lagrangian Multiplier Smoothing Splines: Mathematical Details
#'
#' @description
#' This document provides the mathematical and implementation details for Lagrangian Multiplier
#' Smoothing Splines.
#'
#' @section Statistical Problem Formulation:
#'
#' Consider a dataset with observed predictors \eqn{\mathbf{T}} (an \eqn{N \times q} matrix) and
#' corresponding response values \eqn{\mathbf{y}} (an \eqn{N \times 1} vector). We assume the
#' relationship follows a generalized linear model with an unknown smooth function \eqn{f}:
#'
#' \deqn{y_i \sim \mathcal{D}(g^{-1}(f(\mathbf{t}_i)), \sigma^2)}
#'
#' where \eqn{\mathcal{D}} is a distribution, \eqn{g} is a link function,
#' \eqn{\mathbf{t}_i} is the \eqn{i}th row of \eqn{\mathbf{T}}, and \eqn{\sigma^2} is a dispersion parameter.
#'
#' The objective is to estimate the unknown function \eqn{f} that balances:
#' \itemize{
#'   \item Goodness of fit: How well the model fits the observed data
#'   \item Smoothness: Avoiding excessive fluctuations in the estimated function
#'   \item Interpretability: Understanding the relationship between predictors and response
#' }
#'
#' Lagrangian Multiplier Smoothing Splines address this by:
#' \enumerate{
#'   \item Partitioning the predictor space into \eqn{K+1} regions
#'   \item Fitting local polynomial models within each partition
#'   \item Explicitly enforcing smoothness where partitions meet using Lagrangian multipliers
#'   \item Penalizing the integrated squared second derivative of the estimated function
#' }
#'
#' Unlike other smoothing spline formulations, this technique ensures that no
#' post-fitting algebraic rearrangement or disentangelement of a spline basis is needed
#' to obtain interpretable models. The relationship between predictor and response is explicit,
#' and the basis expansions for each partition are homogeneous.
#'
#' @section Overview:
#' Lagrangian Multiplier Smoothing Splines fit piecewise polynomial regression models to
#' partitioned data, with smoothness at partition boundaries enforced through Lagrangian
#' multipliers. The approach is penalized by the integrated squared second derivative of the
#' estimated function.
#'
#' Unlike traditional smoothing splines that implicitly derive piecewise polynomials through
#' optimization, or regression splines using specialized bases (e.g., B-splines), this method:
#'
#' \itemize{
#'   \item Explicitly represents polynomial basis functions in a natural form
#'   \item Uses Lagrangian multipliers to enforce smoothness constraints
#'   \item Maintains interpretability of coefficient estimates
#' }
#'
#' The fitted model is directly interpretable without the need for reparameterization following fitting.
#' This implementation accommodates non-spline terms and interactions, GLMs, correlation structures, and inequality constraints in addition
#' to linear regression assuming Gaussian response. Extensive customization is offered for users to adapt lgspline for
#' their own modelling frameworks.
#'
#' Core notation:
#' \itemize{
#'   \item \eqn{\mathbf{y}_{(N \times 1)}}: Response vector
#'   \item \eqn{\mathbf{T}_{(N \times q)}}: Matrix of predictors
#'   \item \eqn{\mathbf{X}_{(N \times P)}}: Block-diagonal matrix of polynomial expansions
#'   \item \eqn{\boldsymbol{\Lambda}_{(P \times P)}}: Penalty matrix
#'   \item \eqn{\tilde{\boldsymbol{\beta}}_{(P \times 1)}}: Constrained coefficient estimates
#'   \item \eqn{\mathbf{G}_{(P \times P)}}: \eqn{(\mathbf{X}^T\mathbf{W}\mathbf{X} + \boldsymbol{\Lambda})^{-1}}
#'   \item \eqn{\mathbf{A}_{(P \times r)}}: Constraint matrix ensuring smoothness
#'   \item \eqn{\mathbf{U}_{(P \times P)}}: \eqn{\mathbf{I} - \mathbf{G}\mathbf{A}(\mathbf{A}^T\mathbf{G}\mathbf{A})^{-1}\mathbf{A}^T}
#' }
#'
#' @section Model Formulation and Estimation:
#'
#' \subsection{Model Structure}{
#' The method decomposes the predictor space into K+1 partitions and fits polynomial
#' regression models within each partition, constraining the fitted function to be smooth at
#' partition boundaries.
#'
#' For a single predictor, the function within each partition k is represented as:
#' \deqn{f_k(t) = \beta_k^{(0)} + \beta_k^{(1)}t + \beta_k^{(2)}t^2 + \beta_k^{(3)}t^3 + ...}
#'
#' More generally, for each partition k, the model takes the form:
#' \deqn{f_k(\mathbf{t}) = \mathbf{x}^T\boldsymbol{\beta}_k}
#'
#' Where \eqn{\mathbf{x}} contains polynomial basis functions (intercept, linear, quadratic,
#' cubic terms, and their interactions) and \eqn{\boldsymbol{\beta}_k} are the corresponding coefficients.
#'
#' Smoothness constraints enforce that the function value, first and second derivatives match at
#' adjacent partition boundaries:
#' \deqn{f_k(t_k) = f_{k+1}(t_k)}
#' \deqn{f'_k(t_k) = f'_{k+1}(t_k)}
#' \deqn{f''_k(t_k) = f''_{k+1}(t_k)}
#'
#' These constraints are expressed as linear equations in the \eqn{\mathbf{A}} matrix such that
#' \eqn{\mathbf{A}^T\boldsymbol{\beta} = \mathbf{0}} implies the smoothness conditions are satisfied.
#' }
#'
#' \subsection{Estimation Procedure}{
#' The estimation procedure follows these key steps:
#'
#' 1. Unconstrained estimation:
#'    \deqn{\hat{\boldsymbol{\beta}} = \mathbf{G}\mathbf{X}^T\mathbf{y}}
#'    where \eqn{\mathbf{G} = (\mathbf{X}^T\mathbf{W}\mathbf{X} + \boldsymbol{\Lambda})^{-1}} for weighted design matrix \eqn{\mathbf{W}} and penalty matrix \eqn{\boldsymbol{\Lambda}}.
#'
#' 2. Apply smoothness constraints:
#'    \deqn{\tilde{\boldsymbol{\beta}} = \hat{\boldsymbol{\beta}} - \mathbf{G}\mathbf{A}(\mathbf{A}^T\mathbf{G}\mathbf{A})^{-1}\mathbf{A}^T\hat{\boldsymbol{\beta}}}
#'    This can be computed efficiently as:
#'    \deqn{\tilde{\boldsymbol{\beta}} = \mathbf{U}\hat{\boldsymbol{\beta}}}
#'    where \eqn{\mathbf{U} = \mathbf{I} - \mathbf{G}\mathbf{A}(\mathbf{A}^T\mathbf{G}\mathbf{A})^{-1}\mathbf{A}^T}.
#'
#' 3. For GLMs, iterative refinement:
#'    \itemize{
#'      \item Update \eqn{\mathbf{G}} based on current estimates
#'      \item Recompute \eqn{\tilde{\boldsymbol{\beta}}} using the new \eqn{\mathbf{G}}
#'      \item Continue until convergence
#'    }
#'
#' 4. For inequality constraints (shape or range constraints):
#'    \itemize{
#'      \item Sequential quadratic programming using \code{\link[quadprog]{solve.QP}}
#'      \item Initial values from the equality-constrained estimates
#'    }
#'
#' 5. Variance estimation:
#'    \deqn{\tilde{\sigma}^2 = \frac{\|\mathbf{y} - \mathbf{X}\tilde{\boldsymbol{\beta}}\|^2}{N - \text{trace}(\mathbf{X}\mathbf{U}\mathbf{G}\mathbf{X}^T)}}
#'
#' The variance-covariance matrix of \eqn{\tilde{\boldsymbol{\beta}}} is estimated as:
#' \deqn{\text{Var}(\tilde{\boldsymbol{\beta}}) = \tilde{\sigma}^2\mathbf{U}\mathbf{G}}
#'
#' This approach offers computational efficiency through:
#' \itemize{
#'   \item Parallel computation for each partition
#'   \item Inversion of only small matrices (one per partition)
#'   \item Efficient calculation of \eqn{\mathbf{X}\mathbf{U}\mathbf{G}\mathbf{X}^T} trace
#' }
#' }
#'
#' @section Knot Selection and Partitioning:
#'
#' \subsection{Univariate Case}{
#' For a single predictor, knots are placed at evenly-spaced quantiles of the predictor variable:
#' \itemize{
#'   \item The predictor range is divided into K+1 regions using K interior knots
#'   \item Each observation is assigned to a partition based on these knot boundaries
#'   \item Custom knots can be specified through the \code{custom_knots} parameter
#' }
#' }
#'
#' \subsection{Multivariate Case}{
#' For multiple predictors, a k-means clustering approach is used:
#' \enumerate{
#'   \item K+1 cluster centers are identified using k-means on standardized predictors
#'   \item Neighboring centers are determined using a midpoint criterion:
#'     \itemize{
#'       \item Centers i and j are neighbors if the midpoint between them is closer to both i and j than to any other center
#'     }
#'   \item Observations are assigned to the nearest cluster center
#' }
#'
#' The default number of knots (K) is determined adaptively based on:
#' \itemize{
#'   \item Sample size (N)
#'   \item Number of basis terms per partition (p)
#'   \item Number of predictors (q)
#'   \item Model complexity (GLM family)
#' }
#' }
#'
#' @section Custom Model Specification:
#' The package provides interfaces for extending to custom distributions and estimation procedures.
#'
#' \subsection{Family Structure}{
#' A list containing GLM components:
#' \describe{
#'   \item{family}{Character name of distribution}
#'   \item{link}{Character name of link function}
#'   \item{linkfun}{Function transforming response to linear predictor scale}
#'   \item{linkinv}{Function transforming linear predictor to response scale}
#'   \item{custom_dev.resids}{Optional function for GCV optimization criterion}
#' }
#' }
#'
#' \subsection{Unconstrained Fitting}{
#' Core function for unconstrained coefficient estimation per partition:
#' \preformatted{
#' function(X, y, LambdaHalf, Lambda, keep_weighted_Lambda, family,
#'          tol, K, parallel, cl, chunk_size, num_chunks, rem_chunks,
#'          order_indices, weights, status, ...) {
#'   # Returns p-length coefficient vector
#' }
#' }
#' }
#'
#' \subsection{Dispersion Estimation}{
#' Computes scale/dispersion parameter:
#' \preformatted{
#' function(mu, y, order_indices, family, observation_weights, ...) {
#'   # Returns scalar dispersion estimate
#'   # Defaults to 1 (no dispersion)
#' }
#' }
#' }
#'
#' \subsection{GLM Weights}{
#' Computes observation weights:
#' \preformatted{
#' function(mu, y, order_indices, family, dispersion,
#'          observation_weights, ...) {
#'   # Returns weights vector
#'   # Defaults to family variance with optional weights
#' }
#' }
#' }
#'
#' \subsection{Schur Corrections}{
#' Adjusts covariance matrix for parameter uncertainty:
#' \preformatted{
#' function(X, y, B, dispersion, order_list, K, family,
#'          observation_weights, ...) {
#'   # Required for valid standard errors when dispersion affects coefficient estimation
#' }
#' }
#' }
#'
#' @section Constraint Systems:
#'
#' \subsection{Linear Equality Constraints}{
#' Linear constraints enforce exact relationships between coefficients, implemented through
#' the Lagrangian multiplier approach and integrated with smoothness constraints.
#'
#' Constraints are specified as \eqn{\mathbf{h}^T\boldsymbol{\beta} = \mathbf{h}^T\boldsymbol{\beta}_0} via:
#' \preformatted{
#' lgspline(
#'   y ~ spl(x),
#'   constraint_vectors = h,     # Matrix of constraint vectors
#'   constraint_values = beta0   # Target values
#' )
#' }
#'
#' Common applications include:
#' \itemize{
#'   \item No-intercept models (via \code{no_intercept = TRUE})
#'   \item Offset terms with coefficients constrained to 1
#'   \item Hypothesis testing with nested models
#' }
#' }
#'
#' \subsection{Inequality Constraints}{
#' Inequality constraints enforce bounds or directional relationships on the fitted function
#' through sequential quadratic programming.
#'
#' Built-in constraints include:
#' \itemize{
#'   \item Range constraints: \code{qp_range_lower} and \code{qp_range_upper}
#'   \item Monotonicity: \code{qp_monotonic_increase} or \code{qp_monotonic_decrease}
#'   \item Derivative constraints: \code{qp_positive_derivative} or \code{qp_negative_derivative}
#' }
#'
#' Custom inequality constraints can be defined through:
#' \itemize{
#'   \item \code{qp_Amat_fxn}: Function returning constraint matrix
#'   \item \code{qp_bvec_fxn}: Function returning constraint vector
#'   \item \code{qp_meq_fxn}: Function returning number of equality constraints
#' }
#' }
#'
#' @section Correlation Structure Estimation:
#' Correlation patterns in clustered data are modeled using a matrix whitening approach based on
#' restricted maximum likelihood (REML). The correlation structure is parameterized through the
#' square root inverse matrix \eqn{\mathbf{V}^{-1/2}}.
#' \subsection{Default Correlation Structures}{
#'
#' The package provides several built-in correlation structures for modeling spatial and temporal dependence.
#' These structures can be specified using the "correlation_structure" parameter.
#'
#' The "cluster" parameter identifies groups of observations that share a correlation structure.
#'
#' With exception of exchangeable correlation, all structures use the "spacetime" parameter (an N-row matrix)
#' to determine distances or dissimilarities between observations.
#'
#' \describe{
#'   \item{Exchangeable}{
#'     'exchangeable', 'cs', 'CS', 'compoundsymmetric', 'compound-symmetric', 'compound symmetric'
#'
#'     A single constant correlation between any observations \eqn{\nu} within the same cluster. Uses tanh parameterization \eqn{\nu = tanh^{-1}(\rho)}.
#'   }
#'   \item{Spatial Exponential}{
#'    'spatial-exponential', 'spatialexponential', 'exp', 'exponential'
#'
#'    Correlation decays exponentially with distance: \eqn{e^{-\omega d}} where \eqn{d} is the Euclidean distance
#'    and \eqn{\omega > 0} is the rate parameter (parameterized using softplus: \eqn{\omega = \log(1+e^{\rho})}).
#'    This is mathematically equivalent to the spatial-power correlation \eqn{\theta^d} where \eqn{\theta = e^{-\omega}}
#'    but provides better numerical properties during optimization.
#'   }
#'   \item{AR(1)}{
#'     'ar1', 'ar(1)', 'AR(1)', 'AR1'
#'
#'     Correlation depends on the rank difference between observations: \eqn{\nu^r} where \eqn{r} is the rank difference.
#'     Uses intra-cluster rankings of spacetime distances, with parameterization of \eqn{\nu = tanh^{-1}(\rho)}.
#'   }
#'   \item{Gaussian/Squared Exponential}{
#'     'gaussian', 'rbf', 'squared-exponential'
#'
#'     Smooth decay with squared distance: \eqn{exp(-d^2/(2l^2))} where \eqn{l} is the length scale
#'     (parameterized using softplus: \eqn{l = \log(1+e^{\rho})}).
#'   }
#'   \item{Spherical}{
#'     'spherical', 'Spherical', 'cubic', 'sphere'
#'
#'     Polynomial decay with exact cutoff at range \eqn{r}: \eqn{1 - 1.5(d/r) + 0.5(d/r)^3} for \eqn{d \le r}, 0 otherwise.
#'     Range parameter \eqn{r} is parameterized using softplus: \eqn{r = \log(1+e^{\rho})}.
#'   }
#'   \item{Matérn}{
#'     'matern', 'Matern'
#'
#'     Flexible correlation with adjustable smoothness: \eqn{(2^{1-\nu}/\Gamma(\nu))(\sqrt{2\nu}d/l)^\nu K_\nu(\sqrt{2\nu}d/l)}.
#'     Has two parameters: length scale \eqn{l} and smoothness \eqn{\nu}, both parameterized using softplus.
#'
#'     While theoretically advantageous, there does not exist an analytical gradient for this correlation structure,
#'     so it will be slower and potentially more error-prone to fit than the others.
#'   }
#'   \item{Gamma-Cosine}{
#'     'gamma-cosine', 'gammacosine', 'GammaCosine'
#'
#'     Flexible correlation that can model both positive and negative correlations:
#'     \eqn{(d^{\alpha-1}e^{-\beta d})/(\Gamma(\alpha)/\beta^{\alpha}) \cdot \cos(\omega d)}.
#'     Has three parameters: shape \eqn{\alpha}, rate \eqn{\beta} (both parameterized using softplus), and frequency \eqn{\omega}.
#'     When \eqn{\alpha=1} and \eqn{\omega=0}, reduces to exponential correlation.
#'     When \eqn{\beta} is small and \eqn{\omega=0}, approximates power-law.
#'     When \eqn{\omega > 0}, allows for oscillation between positive and negative correlations.
#'   }
#'   \item{Gaussian-Cosine}{
#'     'gaussian-cosine', 'gaussiancosine', 'GaussianCosine'
#'
#'     Smooth correlation that allows for oscillations: \eqn{exp(-d^2/(2l^2)) \cdot \cos(\omega d)}.
#'     Has two parameters: length scale \eqn{l} (parameterized using softplus) and frequency \eqn{\omega}.
#'     When \eqn{\omega=0}, reduces to standard Gaussian correlation.
#'     When \eqn{\omega > 0}, allows for oscillation between positive and negative correlations.
#'   }
#' }
#'
#' \strong{Parameter Interpretation:}
#'
#' Correlation parameters are estimated on transformed scales to ensure valid ranges.
#' For proper interpretation, these estimates must be converted back to their natural scales.
#'
#' For correlation structures using the \code{tanh} transformation (Exchangeable, AR(1)),
#' the correlation parameter \eqn{\rho} is bounded between -1 and 1:
#'
#' \preformatted{
#' # Correlation estimate
#' cor_est <- tanh(model_fit$VhalfInv_params_estimates)
#'
#' # 95% CI
#' se <- sqrt(diag(model_fit$VhalfInv_params_vcov)) / sqrt(model_fit$N)
#' ci <- tanh(model_fit$VhalfInv_params_estimates + c(-1.96, 1.96) * se)
#' }
#'
#' For correlation structures using the softplus transformation (Gaussian/Squared Exponential, Spherical, Matérn):
#'
#' \preformatted{
#' # Length scale estimate
#' length_scale <- log(1 + exp(model_fit$VhalfInv_params_estimates))
#'
#' # 95% CI
#' se <- sqrt(diag(model_fit$VhalfInv_params_vcov)) / sqrt(model_fit$N)
#' ci <- log(1 + exp(model_fit$VhalfInv_params_estimates + c(-1.96, 1.96) * se))
#' }
#'
#' For Spatial-Exponential correlation with softplus parameterization:
#'
#' \preformatted{
#' # Rate parameter omega estimate
#' omega_est <- log(1 + exp(model_fit$VhalfInv_params_estimates))
#'
#' # Equivalent power-law parameter theta estimate
#' theta_est <- exp(-omega_est)
#'
#' # 95% CI for the working parameter
#' se <- sqrt(diag(model_fit$VhalfInv_params_vcov)) / sqrt(model_fit$N)
#' rho_ci <- model_fit$VhalfInv_params_estimates + c(-1.96, 1.96) * se
#'
#' # Transform CI to omega scale
#' omega_ci <- log(1 + exp(rho_ci))
#'
#' # Transform CI to theta scale (note the reversed order due to negative exponent)
#' theta_ci <- exp(-omega_ci[c(2,1)])
#' }
#'
#' For Matérn correlation (two parameters):
#'
#' \preformatted{
#' # Length scale and smoothness estimates
#' length_scale <- log(1 + exp(model_fit$VhalfInv_params_estimates[1]))
#' nu <- log(1 + exp(model_fit$VhalfInv_params_estimates[2]))
#'
#' # 95% CIs
#' se <- sqrt(diag(model_fit$VhalfInv_params_vcov)) / sqrt(model_fit$N)
#' length_scale_ci <- log(1 + exp(model_fit$VhalfInv_params_estimates[1] + c(-1.96, 1.96) * se[1]))
#' nu_ci <- log(1 + exp(model_fit$VhalfInv_params_estimates[2] + c(-1.96, 1.96) * se[2]))
#' }
#'
#' For Gamma-Cosine correlation (three parameters):
#'
#' \preformatted{
#' # Shape (phi_1), rate (phi_2), and omega (frequency) estimates
#' shape <- log(1 + exp(model_fit$VhalfInv_params_estimates[1]))  # shape = softplus(phi_1)
#' rate <- log(1 + exp(model_fit$VhalfInv_params_estimates[2]))   # rate = softplus(phi_2)
#' omega <- model_fit$VhalfInv_params_estimates[3]                # No transformation
#'
#' # 95% CIs
#' se <- sqrt(diag(model_fit$VhalfInv_params_vcov)) / sqrt(model_fit$N)
#' shape_ci <- log(1 + exp(model_fit$VhalfInv_params_estimates[1] + c(-1.96, 1.96) * se[1]))
#' rate_ci <- log(1 + exp(model_fit$VhalfInv_params_estimates[2] + c(-1.96, 1.96) * se[2]))
#' omega_ci <- model_fit$VhalfInv_params_estimates[3] + c(-1.96, 1.96) * se[3]
#' }
#'
#' For Gaussian-Cosine correlation (two parameters):
#'
#' \preformatted{
#' # Length scale and omega (frequency) estimates
#' length_scale <- log(1 + exp(model_fit$VhalfInv_params_estimates[1]))
#' omega <- model_fit$VhalfInv_params_estimates[2]  # No transformation
#'
#' # 95% CIs
#' se <- sqrt(diag(model_fit$VhalfInv_params_vcov)) / sqrt(model_fit$N)
#' length_scale_ci <- log(1 + exp(model_fit$VhalfInv_params_estimates[1] + c(-1.96, 1.96) * se[1]))
#' omega_ci <- model_fit$VhalfInv_params_estimates[2] + c(-1.96, 1.96) * se[2]
#' }
#'
#' The variance-covariance matrix \code{model_fit$VhalfInv_params_vcov} contains parameter
#' uncertainty on the transformed scale. When reporting results, both point estimates and
#' confidence intervals should be transformed to the appropriate scale.
#' }
#'
#' \subsection{Correlation Objective}{
#' The default correlation structures include the following, as covered above:
#' \itemize{
#'   \item \strong{Exchangeable}: Constant correlation within clusters
#'   \item \strong{Spatial Exponential}: Exponential decay with distance
#'   \item \strong{AR(1)}: Correlation based on rank differences
#'   \item \strong{Gaussian/Squared Exponential}: Smooth decay with squared distance
#'   \item \strong{Spherical}: Polynomial decay with exact cutoff
#'   \item \strong{Matern}: Flexible correlation with adjustable smoothness
#'   \item \strong{Gamma-Cosine}: Combined gamma decay with oscillation
#'   \item \strong{Gaussian-Cosine}: Combined Gaussian decay with oscillation
#' }
#'
#' The REML objective function combines correlation structure, parameter estimates, and
#' smoothness constraints:
#'
#' \deqn{\frac{1}{N}\big(-\log|\mathbf{V}^{-1/2}| + \frac{N}{2}\log(\sigma^2) + \frac{1}{2\sigma^2}(\mathbf{y} - \mathbf{X}\tilde{\boldsymbol{\beta}})^{T}\mathbf{V}^{-1}(\mathbf{y} - \mathbf{X}\tilde{\boldsymbol{\beta}}) + \frac{1}{2}\log|\sigma^2\mathbf{U}\mathbf{G} | \big)}
#'
#' Analytical gradients are provided for efficient optimization of correlation parameters.
#' }
#'
#' \subsection{Custom Correlation Functions}{
#' Custom correlation structures can be specified through:
#' \itemize{
#'   \item \code{VhalfInv_fxn}: Creates \eqn{\mathbf{V}^{-1/2}}
#'   \item \code{Vhalf_fxn}: Creates \eqn{\mathbf{V}^{1/2}} (for non-Gaussian responses)
#'   \item \code{REML_grad}: Provides analytical gradient
#'   \item \code{VhalfInv_logdet}: Efficient determinant computation
#'   \item \code{custom_VhalfInv_loss}: Alternative optimization objective
#' }
#' }
#'
#' @section Penalty Construction and Optimization:
#'
#' \subsection{Smoothing Spline Penalty}{
#' The penalty matrix \eqn{\boldsymbol{\Lambda}} is constructed as the integrated squared second
#' derivative of the estimated function over the support of the predictors:
#'
#' \deqn{\int_a^b \|f''(t)\|^2 dt = \int_a^b \|\mathbf{x}''^\top \boldsymbol{\beta}_k\|^2 dt = \boldsymbol{\beta}_k^{T}\left\{\int_a^b \mathbf{x}''\mathbf{x}''^{\top} dt\right\} \boldsymbol{\beta}_k = \boldsymbol{\beta}_k^{T}\boldsymbol{\Lambda}_s\boldsymbol{\beta}_k}
#'
#' For a one-dimensional cubic function, the second derivative basis is \eqn{\mathbf{x}'' = (0, 0, 2, 6t)^{T}} and the penalty matrix has a closed-form expression:
#'
#' \deqn{\boldsymbol{\Lambda}_s = \begin{pmatrix}
#' 0 & 0 & 0 & 0 \\
#' 0 & 0 & 0 & 0 \\
#' 0 & 0 & 4(b-a) & 6(b^2-a^2) \\
#' 0 & 0 & 6(b^2-a^2) & 12(b^3-a^3)
#' \end{pmatrix}}
#'
#' The full penalty matrix combines the smoothing spline penalty with a ridge penalty on
#' non-spline terms and optional predictor-specific or partition-specific penalties:
#'
#' \deqn{\boldsymbol{\Lambda} = \lambda_s (\boldsymbol{\Lambda}_s + \lambda_r\boldsymbol{\Lambda}_r + \sum_{l=1}^L \lambda_l \boldsymbol{\Lambda}_l)}
#'
#' where:
#' \itemize{
#'   \item \eqn{\lambda_s} is the global smoothing parameter (\code{wiggle_penalty})
#'   \item \eqn{\boldsymbol{\Lambda}_s} is the smoothing spline penalty
#'   \item \eqn{\boldsymbol{\Lambda}_r} is a ridge penalty on lower-order terms (\code{flat_ridge_penalty})
#'   \item \eqn{\lambda_l} are optional predictor/partition-specific penalties
#' }
#'
#' This penalty structure defines a meaningful metric for function smoothness, pulling estimates toward linear functions rather than simply shrinking coefficients toward zero.
#' }
#'
#' \subsection{Penalty Optimization}{
#' A penalized variant of Generalized Cross Validation (GCV) is used to find optimal penalties:
#'
#' \deqn{\text{GCV} = \frac{\sum_{i=1}^N r_i^2}{N \big(1-\frac{1}{N}\text{trace}(\mathbf{X}\mathbf{U}\mathbf{G}\mathbf{X}^T) \big)^2} + \frac{mp}{2} \sum_{l=1}^{L}(\log (1+\lambda_{l})-1)^2}
#'
#' where:
#' \itemize{
#'   \item \eqn{r_i = y_i - \tilde{y}_i} are residuals (or replaced with custom alternative for GLMs)
#'   \item \eqn{\text{trace}(\mathbf{X}\mathbf{U}\mathbf{G}\mathbf{X}^T)} represents effective degrees of freedom
#'   \item \eqn{mp} is the "meta-penalty" term that regularizes predictor/partition-specific penalties
#' }
#'
#' For GLMs, a pseudo-count approach is used to ensure valid link transformations, or
#' \code{custom_dev.resids} can be provided to replace the sum-of-squared errors.
#'
#' Optimization employs:
#' \itemize{
#'   \item Grid search for initial values
#'   \item Damped BFGS with analytical gradients
#'   \item Automated restarts and error handling
#'   \item Inflation factor \eqn{((N+2)/(N-2))^2} for final penalties
#' }
#' }
#' @references
#' Buse, A., & Lim, L. (1977). Cubic Splines as a Special Case of Restricted Least Squares. Journal of the American Statistical Association, 72, 64-68.
#'
#' Craven, P., & Wahba, G. (1978). Smoothing noisy data with spline functions. Numerische Mathematik, 31, 377-403.
#'
#' Eilers, P. H. C., & Marx, B. D. (1996). Flexible smoothing with B-splines and penalties. Statistical Science, 11, 89-121.
#'
#' Ezhov, N., Neitzel, F., & Petrovic, S. (2018). Spline approximation, Part 1: Basic methodology. Journal of Applied Geodesy, 12, 139-155.
#'
#' Kisi, O., Heddam, S., Parmar, K. S., Petroselli, A., Külls, C., & Zounemat-Kermani, M. (2025). Integration of Gaussian process regression and K means clustering for enhanced short term rainfall runoff modeling. Scientific Reports, 15, 7444.
#'
#' Nocedal, J., & Wright, S. J. (2006). Numerical Optimization (2nd ed.). Springer.
#'
#' Reinsch, C. H. (1967). Smoothing by spline functions. Numerische Mathematik, 10, 177-183.
#'
#' Searle, S. R., Casella, G., & McCulloch, C. E. (2009). Variance Components. Wiley Series in Probability and Statistics. Wiley.
#'
#' Silverman, B. W. (1984). Spline Smoothing: The Equivalent Variable Kernel Method. The Annals of Statistics, 12, 898-916.
#'
#' Wahba, G. (1990). Spline Models for Observational Data. Society for Industrial and Applied Mathematics.
#'
#' Wood, S. N. (2006). Generalized Additive Models: An Introduction with R. Chapman & Hall/CRC, Boca Raton.
#'
#' @name Details
NULL
