#' Lagrangian Multiplier Smoothing Splines: Mathematical Details
#'
#' @description
#' This document provides the mathematical and implementation details for
#' Lagrangian Multiplier Smoothing Splines as implemented in \pkg{lgspline}.
#'
#' The material is presented such that a programmer or statistician of
#' reasonable experience and background can understand and implement the
#' procedure from scratch, and also potentially critique some of the modelling
#' choices that went into designing this package.
#'
#' Informally, \pkg{lgspline} answers the following question:
#' How can we best adapt a useful functionality of basis splines, under the
#' alternative interpretation of smoothing as explicit external constraints instead?
#'
#' The obvious benefit is a much more flexible and interpretable final model that for
#' non-experienced users is simply easier to understand without post-hoc processing, and
#' for experienced users can be used to customize models more easily.
#'
#' The drawback is that the interpretation of constraints as external adds a new layer of complexity
#' to each step of the model fitting process, whereas for implicit design matrix
#' construction these complications are largely bypassed, and computational stability
#' is compromised.
#'
#' While it is true a B-spline can always be converted back into monomial form,
#' tensor-product splines that generalize this to multiple dimensions often
#' explodes the number and degree of interaction terms, the conversion may not
#' be computationally stable, and it is not available in standard software.
#'
#' @section Statistical Problem Formulation:
#'
#' Consider an \eqn{N \times q} matrix of predictors
#' \eqn{\mathbf{T} = (\mathbf{t}_1, \dots, \mathbf{t}_N)^{\top}} and an
#' \eqn{N \times 1} response vector \eqn{\mathbf{y} = (y_1, \dots, y_N)^{\top}}.
#' We assume the relationship follows a generalized linear model with unknown
#' smooth function \eqn{f}:
#' \deqn{y_i \sim \mathcal{D}(g^{-1}(f(\mathbf{t}_i)),\, \sigma^2)}
#' where \eqn{\mathcal{D}} is a distribution (e.g. exponential family or related) with mean
#' \eqn{\mu_i = g^{-1}(f(\mathbf{t}_i))}, link function \eqn{g(\cdot)}, and
#' dispersion parameter \eqn{\sigma^2}. For Gaussian response with identity
#' link, observations are independently distributed as
#' \eqn{y_i \mid \mathbf{t}_i, \sigma^2 \sim \mathcal{N}(f(\mathbf{t}_i), \sigma^2)}.
#'
#' The objective is to estimate \eqn{f} by:
#' \enumerate{
#'   \item Partitioning the predictor space into \eqn{K+1} mutually exclusive regions.
#'   \item Fitting local polynomial models within each partition.
#'   \item Enforcing smoothness at partition boundaries via Lagrangian multipliers.
#'   \item Penalizing the integrated squared second derivative to discourage roughness.
#' }
#'
#' Unlike other smoothing spline formulations, no post-fitting algebraic
#' rearrangement or disentanglement of a spline basis is needed to obtain
#' interpretable models. The polynomial expansions are homogeneous across
#' partitions, and the relationship between predictor and response is explicit
#' at the coefficient level.
#'
#' To anchor the notation, in the single-predictor cubic case one would write
#' \deqn{\hat{f}(t_i) = \hat{\beta}_{(0)} + \hat{\beta}_{(1)}t_i +
#'   \hat{\beta}_{(2)}t_i^2 + \hat{\beta}_{(3)}t_i^3 =
#'   \mathbf{x}_i^{\top}\hat{\boldsymbol{\beta}},}
#' where \eqn{\mathbf{x}_i = (1, t_i, t_i^2, t_i^3)^{\top}}. The LMSS
#' formulation preserves exactly this kind of polynomial representation, but now
#' does so within each partition and then forces neighboring pieces to agree in
#' the smoothness conditions described below.
#'
#' Core notation used throughout:
#' \itemize{
#'   \item \eqn{\mathbf{y}_{(N \times 1)}}: Response vector.
#'   \item \eqn{\mathbf{T}_{(N \times q)}}: Matrix of predictors.
#'   \item \eqn{\mathbf{X}_{(N \times P)}}: Block-diagonal matrix of polynomial
#'     expansions, with diagonal blocks \eqn{\mathbf{X}_k} of dimension
#'     \eqn{n_k \times p}.
#'   \item \eqn{\boldsymbol{\Lambda}_{(P \times P)}}: Block-diagonal penalty matrix,
#'     with blocks \eqn{\boldsymbol{\Lambda}_k} of dimension \eqn{p \times p}.
#'   \item \eqn{\hat{\boldsymbol{\beta}}_{(P \times 1)}}: Unconstrained penalized
#'     estimate.
#'   \item \eqn{\tilde{\boldsymbol{\beta}}_{(P \times 1)}}: Constrained coefficient
#'     estimates.
#'   \item \eqn{\mathbf{G}_{(P \times P)}}: Block-diagonal matrix with blocks
#'     \eqn{\mathbf{G}_k = (\mathbf{X}_k^{\top}\mathbf{W}_k\mathbf{D}_k\mathbf{X}_k + \boldsymbol{\Lambda}_k)^{-1}},
#'     where \eqn{\mathbf{W}_k} and \eqn{\mathbf{D}_k} are defined below.
#'   \item \eqn{\mathbf{A}_{(P \times r)}}: Constraint matrix encoding smoothness
#'     conditions. Reduced to linearly independent columns via pivoted QR
#'     decomposition.
#'   \item \eqn{\mathbf{U}_{(P \times P)}}:
#'     \eqn{\mathbf{I} - \mathbf{G}\mathbf{A}(\mathbf{A}^{\top}\mathbf{G}\mathbf{A})^{-1}\mathbf{A}^{\top}}.
#'   \item \eqn{\mathbf{D}_{(N \times N)}}: Diagonal matrix of user-supplied
#'     observation weights (\code{observation_weights} or \code{weights}).
#'     Defaults to the identity. These play the role of prior precision on
#'     individual observations: a weight of 2 is equivalent to seeing that
#'     observation twice.
#'   \item \eqn{\mathbf{W}_{(N \times N)}}: Diagonal matrix of GLM working
#'     weights. In the implementation these diagonal entries are whatever is
#'     returned by \code{glm_weight_function}; by default this is
#'     \code{family$variance(mu)}, optionally multiplied by user-supplied
#'     observation weights. For Gaussian response with identity link,
#'     \eqn{\mathbf{W} = \mathbf{I}}. For other families,
#'     \eqn{\mathbf{W}} depends on the current fitted values and is updated at
#'     each Newton--Raphson iteration. For the common canonical families used
#'     by default, this matches the familiar Fisher-scoring weighting role.
#'   \item \eqn{\mathbf{V}_{(N \times N)}}: Correlation matrix of errors.
#'     When no correlation structure is specified, \eqn{\mathbf{V} = \mathbf{I}}.
#'     Otherwise supplied via \code{VhalfInv} or estimated through
#'     \code{VhalfInv_fxn}.
#' }
#'
#' In the Gaussian identity case with unit weights and no correlation,
#' \eqn{\mathbf{G}_k = (\mathbf{X}_k^{\top}\mathbf{X}_k + \boldsymbol{\Lambda}_k)^{-1}}
#' and most formulas simplify accordingly. When \eqn{\mathbf{D}} or
#' \eqn{\mathbf{W}} appear in a formula, the product \eqn{\mathbf{W}\mathbf{D}}
#' means ``GLM working weights times observation weights''; whenever one of
#' them is the identity it drops out.
#'
#' Before these quantities reach the main fitting stage, the user-facing inputs
#' are parsed, standardized, and organized by \code{\link{process_input}}. When
#' the formula interface is used and \code{auto_encode_factors = TRUE}, that
#' preprocessing step also relies on helpers such as \code{\link{create_onehot}}
#' to encode factor levels before the design reaches \code{\link{lgspline.fit}()}.
#' The notation in the remainder of this document therefore refers to the internal
#' objects that actually enter \code{\link{lgspline.fit}()}, not necessarily the raw
#' objects originally supplied by the user.
#'
#' From the user side, many of the arguments that control these internal objects
#' can be supplied either individually or through the grouped lists
#' \code{penalty_args}, \code{tuning_args}, \code{expansion_args},
#' \code{constraint_args}, \code{qp_args}, \code{parallel_args},
#' \code{covariance_args}, \code{return_args}, and \code{glm_args}, as
#' documented in \code{\link{lgspline}}. These grouped lists are unpacked before
#' dispatch into the same fitting pipeline, so they are a convenience layer
#' rather than a separate modeling abstraction. A closely related exploratory
#' mode is \code{dummy_fit = TRUE} in \code{\link{lgspline}} or
#' \code{\link{lgspline.fit}}, which runs the preprocessing, partition
#' construction, expansion building, and penalty setup without solving for
#' nonzero coefficients, making it a practical way to inspect objects such as
#' \code{X}, \code{A}, the returned \code{make_partition_list} from
#' \code{\link{make_partitions}}, and the assembled \code{penalties} from
#' \code{\link{compute_Lambda}} before a full fit.
#'
#'
#' @section Model Formulation and Estimation:
#'
#' \subsection{Piecewise Polynomial Structure}{
#' For \eqn{K} knots (one predictor) or \eqn{K+1} partitions (multiple
#' predictors) there are \eqn{K+1} mutually exclusive partitions
#' \eqn{\mathcal{P}_0, \dots, \mathcal{P}_{K}}. Each observation \eqn{i}
#' belongs to exactly one partition. Within partition \eqn{k}, the function is
#' represented as a polynomial of degree \eqn{p-1} in each predictor:
#' \deqn{\hat{f}_k(\mathbf{t}) = \mathbf{x}^{\top}\tilde{\boldsymbol{\beta}}_k}
#' where \eqn{\mathbf{x}} collects the polynomial basis terms (intercept,
#' linear, quadratic, cubic, and optionally quartic and interaction terms) and
#' \eqn{\tilde{\boldsymbol{\beta}}_k} are the corresponding coefficients. In
#' one predictor, the same idea can be written more explicitly as
#' \deqn{\hat{f}(t_i) = \sum_{k=0}^{K}
#'   \mathbf{x}_{ik}^{\top}\hat{\boldsymbol{\beta}}_k
#'   \mathbf{1}(t_i \in \mathcal{P}_k),}
#' which highlights that the unconstrained problem is just a collection of
#' local polynomial regressions. The expansions are homogeneous across
#' partitions, so coefficients are directly comparable. This is implemented via
#' \code{\link{get_polynomial_expansions}}.
#'
#' The exact contents of \eqn{\mathbf{x}} are controlled by the basis-expansion
#' arguments documented in \code{\link{lgspline}}: \code{include_quadratic_terms},
#' \code{include_cubic_terms}, \code{include_quartic_terms},
#' \code{include_2way_interactions}, \code{include_3way_interactions},
#' \code{include_quadratic_interactions}, \code{exclude_interactions_for},
#' \code{exclude_these_expansions}, and \code{custom_basis_fxn}. Likewise,
#' \code{just_linear_with_interactions} and
#' \code{just_linear_without_interactions} determine which predictors remain
#' structurally linear even though they still participate in the same
#' partition-wise polynomial bookkeeping described here.
#'
#' Letting \eqn{p} denote the number of basis terms per partition,
#' \eqn{P = p(K+1)} is the total number of coefficients. The full design matrix
#' \eqn{\mathbf{X}} and penalty matrix \eqn{\boldsymbol{\Lambda}} are
#' block-diagonal with \eqn{K+1} blocks, so unconstrained estimation reduces to
#' \eqn{K+1} independent penalized regressions, which appears as follows for the identity link case:
#' \deqn{\hat{\boldsymbol{\beta}}_k = \mathbf{G}_k \mathbf{X}_k^{\top} \mathbf{W}_k\mathbf{D}_k\mathbf{y}_k, \quad
#'   \mathbf{G}_k = (\mathbf{X}_k^{\top}\mathbf{W}_k\mathbf{D}_k\mathbf{X}_k + \boldsymbol{\Lambda}_k)^{-1}.}
#' For Gaussian identity with unit weights this reduces to the familiar
#' \eqn{\mathbf{G}_k = (\mathbf{X}_k^{\top}\mathbf{X}_k + \boldsymbol{\Lambda}_k)^{-1}}.
#' The block-diagonal structure means these can be computed in parallel across
#' partitions. In the user-facing interface this is realized by supplying a
#' cluster through \code{cl}, optionally controlling work splitting with
#' \code{chunk_size}, and enabling stages such as \code{parallel_eigen}
#' for the eigendecompositions and, in non-Gaussian Path 3, \code{parallel_unconstrained}
#' for the partition-wise unconstrained fits; nearby stages can likewise use
#' \code{parallel_penalty} and \code{parallel_make_constraint}. The eigenvalue
#' decomposition and matrix square roots of each \eqn{\mathbf{G}_k} are
#' computed by \code{\link{compute_G_eigen}}, and can be returned in the fitted
#' object as \code{G} and \code{Ghalf} when \code{return_G = TRUE} and
#' \code{return_Ghalf = TRUE}.
#'
#' Fitted values for the canonical Gaussian case appear as
#' \eqn{\tilde{\mathbf{y}} = \mathbf{X}\tilde{\boldsymbol{\beta}} = \mathbf{H}\mathbf{y}}
#' for \eqn{\mathbf{H} = \mathbf{X}\mathbf{U}\mathbf{G}\mathbf{X}^{\top}}.
#' }
#'
#' \subsection{Smoothing Constraints and the Constraint Matrix}{
#' Without further intervention the piecewise polynomial will be discontinuous.
#' The central idea of LMSS is that smoothness is not hidden inside a special
#' basis, but instead imposed directly where neighboring partitions meet.
#' At each knot \eqn{t_{k,k+1}} between neighboring partitions \eqn{k} and
#' \eqn{k+1}, up to three smoothing constraints are imposed:
#' \enumerate{
#'   \item Continuity: \eqn{\mathbf{x}_{k,k+1}^{\top}\boldsymbol{\beta}_k = \mathbf{x}_{k,k+1}^{\top}\boldsymbol{\beta}_{k+1}}.
#'   \item First-derivative continuity: \eqn{\mathbf{x}_{k,k+1}^{\prime\top}\boldsymbol{\beta}_k = \mathbf{x}_{k,k+1}^{\prime\top}\boldsymbol{\beta}_{k+1}}.
#'   \item Second-derivative continuity: \eqn{\mathbf{x}_{k,k+1}^{\prime\prime\top}\boldsymbol{\beta}_k = \mathbf{x}_{k,k+1}^{\prime\prime\top}\boldsymbol{\beta}_{k+1}}.
#' }
#' where \eqn{\mathbf{x}^{\prime}} and \eqn{\mathbf{x}^{\prime\prime}} are
#' elementwise first and second derivatives of the basis with respect to
#' \eqn{\mathbf{t}}. For the familiar cubic single-predictor basis
#' \eqn{\mathbf{x} = (1, t, t^2, t^3)^{\top}}, these derivative vectors are
#' \deqn{\mathbf{x}' = (0, 1, 2t, 3t^2)^{\top}, \qquad
#'   \mathbf{x}'' = (0, 0, 2, 6t)^{\top}.}
#' With \eqn{K} knots this yields up to \eqn{3K} scalar
#' constraints (for a single predictor; more for multiple predictors with
#' interactions), collected as linear equations
#' \eqn{\mathbf{A}^{\top}\boldsymbol{\beta} = \mathbf{0}}
#' in a \eqn{P \times r} matrix \eqn{\mathbf{A}}. The constraint matrix is
#' built by \code{\link{make_constraint_matrix}} and returned in the fitted
#' object as \code{A}.
#'
#' In higher dimensions or with many partitions, the constraints can become
#' over-specified and force the model toward a single global polynomial. In these
#' cases it is recommended to drop second-derivative constraints or include quartic
#' terms, allowing the model to fit a richer surface while maintaining
#' perceived smoothness at knots. The appropriate constraint level can be
#' controlled via \code{include_constrain_fitted},
#' \code{include_constrain_first_deriv}, and
#' \code{include_constrain_second_deriv}.
#' The companion flag \code{include_constrain_interactions} determines whether
#' the analogous mixed-partial constraints are imposed for interaction terms,
#' and \code{no_intercept} adds the special homogeneous equality constraint that
#' fixes the intercept at zero (the same behavior triggered by using
#' \code{0 +} in the formula interface).
#'
#' Before computing the projection \eqn{\mathbf{U}}, the constraint matrix is
#' reduced to a linearly independent subset of columns via pivoted QR
#' decomposition. This avoids numerical instability from redundant constraints
#' and ensures \eqn{\mathbf{A}^{\top}\mathbf{G}\mathbf{A}} is invertible.
#' }
#'
#' \subsection{Lagrangian Projection}{
#' The constrained estimate is derived via Lagrangian multipliers. Define the
#' \eqn{P \times P} projection matrix:
#' \deqn{\mathbf{U} = \mathbf{I} - \mathbf{G}\mathbf{A}(\mathbf{A}^{\top}\mathbf{G}\mathbf{A})^{-1}\mathbf{A}^{\top}.}
#' Then the constrained estimate is:
#' \deqn{\tilde{\boldsymbol{\beta}} = \mathbf{U}\hat{\boldsymbol{\beta}}.}
#' The matrix \eqn{\mathbf{U}} has the property that
#' \eqn{\mathbf{U}\mathbf{G}\mathbf{U}^{\top} = \mathbf{U}\mathbf{G}}, which is
#' used extensively in variance estimation and posterior draws. In words, the
#' unconstrained penalized estimate is projected back into the coefficient space
#' that satisfies the smoothness restrictions, and all subsequent uncertainty
#' calculations inherit that same projected geometry.
#' The projection is computed via \code{\link{get_U}} and, when requested,
#' returned in the fitted object as \code{U} through \code{return_U = TRUE}.
#'
#' When the constraints are inhomogeneous
#' (\eqn{\mathbf{A}^{\top}\boldsymbol{\beta} = \mathbf{c}} with
#' \eqn{\mathbf{c} \neq \mathbf{0}}), a particular solution
#' \eqn{\boldsymbol{\beta}_0} satisfying
#' \eqn{\mathbf{A}^{\top}\boldsymbol{\beta}_0 = \mathbf{c}} is added back
#' after projection, yielding the full Lagrangian solution
#' \eqn{\mathbf{U}\hat{\boldsymbol{\beta}} + (\mathbf{I} - \mathbf{U})\boldsymbol{\beta}_0}.
#' In \code{\link{lgspline}} and \code{\link{lgspline.fit}}, users realize
#' this by supplying extra equality columns in \code{constraint_vectors}
#' together with matching right-hand sides in \code{constraint_values};
#' \code{null_constraint} provides the alternate shorthand documented in
#' \code{\link{lgspline}} when \code{constraint_vectors} is supplied and
#' \code{constraint_values} is left empty.
#'
#' In practice \eqn{\mathbf{U}} is never explicitly formed during fitting.
#' The constrained estimate is obtained from a transformed OLS residual problem (the
#' \eqn{\mathbf{G}^{1/2}\mathbf{r}^{*}} trick) in four steps:
#' \enumerate{
#'   \item Obtain the unconstrained partition-wise
#'     unconstrained estimate \eqn{\hat{\boldsymbol{\beta}}}.
#'   \item Set \eqn{\mathbf{y}^{*} = \mathbf{G}^{-1/2}\hat{\boldsymbol{\beta}}} and
#'     \eqn{\mathbf{X}^{*} = \mathbf{G}^{1/2}\mathbf{A}}.
#'   \item Fit the linear model \eqn{\mathbb{E}[\mathbf{y}^{*}] = \mathbf{X}^{*}\boldsymbol{\gamma}}
#'     by OLS using QR decomposition.
#'   \item Compute the residuals \eqn{\mathbf{r}^{*} = \mathbf{y}^{*} - \mathbf{X}^{*}(\mathbf{X}^{*\top}\mathbf{X}^{*})^{-1}\mathbf{X}^{*\top}\mathbf{y}^{*}} from that transformed OLS
#'     fit and recover the constrained estimate by
#'     \eqn{\tilde{\boldsymbol{\beta}} = \mathbf{G}^{1/2}\mathbf{r}^{*}}.
#' }
#'
#' A scaling factor \eqn{1/\sqrt{K+1}} is applied to both
#' \eqn{\mathbf{X}^{*}} and \eqn{\mathbf{y}^{*}} prior to the OLS call
#' and divided out afterward, improving numerical conditioning when the
#' constraint matrix has many rows.
#'
#' The most expensive operation in this approach is the QR decomposition of the
#' \eqn{P \times r} matrix \eqn{\mathbf{X}^{*} = \mathbf{G}^{1/2}\mathbf{A}},
#' which is far cheaper than working with the full \eqn{P \times P} system directly.
#' Without correlation or SQP constraints, \eqn{\mathbf{G}} is stored and operated upon as a
#' list of \eqn{K+1} small \eqn{p \times p} matrices rather than the full
#' \eqn{P \times P} block-diagonal, saving substantial memory when \eqn{K}
#' is large and allowing for parallelism.
#'
#' When correlation is present, \eqn{\mathbf{X}^{\top}\mathbf{V}^{-1}\mathbf{X}}
#' is no longer block-diagonal, so the full-dimensional system must be handled
#' directly unless the Woodbury acceleration
#' (see \code{\link{.woodbury_decompose_V}}) applies.
#' When additional inequality constraints are present, the code either augments the equality system with a
#' partition-wise active-set refinement (block-separable case) or falls back to
#' dense SQP via the Goldfarb--Idnani dual active-set method implemented in
#' \code{\link[quadprog]{solve.QP}}.
#' }
#'
#' @section GLM Extension and Iterative Updates:
#'
#' \subsection{Working Quantities}{
#' The most direct implementation entry point for GLMs is
#' \code{\link{get_B}}. In \code{\link{.get_B_glm_nocorr}}, which is the
#' no-correlation non-Gaussian path used when \code{correlation_structure},
#' \code{VhalfInv}, and \code{VhalfInv_fxn} are all absent, the algorithm starts
#' with partition-wise unconstrained Newton--Raphson fits obtained through
#' \code{unconstrained_fit_fxn} (by default
#' \code{\link{unconstrained_fit_default}} with helpers
#' \code{\link{damped_newton_r}} and \code{\link{nr_iterate}}). This is the
#' path that benefits from the block-diagonal structure, so the partitions can
#' first be fit separately and only then be projected back into the smoothness
#' constraint space.
#'
#' The Gaussian identity model gives the cleanest closed-form derivation, but
#' once the response mean is related to the linear predictor through a link
#' function, the code works with the same quantities that appear in the working-response
#' updates used by Fisher scoring. For GLMs with mean \eqn{\mu_i = g^{-1}(\eta_i)} and linear
#' predictor \eqn{\eta_i = \mathbf{x}_i^{\top}\boldsymbol{\beta}}, the usual
#' working quantities are
#' \deqn{z_i = \eta_i + \frac{y_i - \mu_i}{g'(\mu_i)},}
#' and a diagonal working-weight matrix \eqn{\mathbf{W}} whose entries are
#' supplied by \code{glm_weight_function} at the current fitted values.
#' Here \eqn{g'(\mu) = \partial \mu / \partial \eta} is the derivative of the
#' inverse link function. Under the default \code{glm_weight_function}, the
#' diagonal entries reduce to \code{family$variance(mu)} (optionally including
#' observation weights), which for the common canonical families plays the same
#' role as the familiar Fisher-scoring weights. These quantities appear
#' explicitly in the GEE path and in the blockfit GLM solver, where the code
#' recomputes \eqn{\mathbf{z}} and \eqn{\mathbf{W}} from the current linear
#' predictor.
#'
#' The key distinction is therefore structural. In \code{\link{.get_B_glm_nocorr}},
#' the package does \emph{not} run a single global working-response solve on the constrained system.
#' Instead it (1) obtains partition-wise unconstrained penalized estimates,
#' (2) projects those estimates onto the equality-constraint space, and
#' (3) when iteration is requested through \code{iterate_final_fit} in
#' \code{\link{lgspline}} (mapped internally to \code{iterate}), recomputes the
#' weighted Gram matrices and corresponding \eqn{\mathbf{G}^{1/2}} and
#' \eqn{\mathbf{G}^{-1/2}} at the current constrained estimate, then re-projects.
#' The shared recomputation step is handled numerically by
#' \code{.solver_recompute_G_at_estimate} in \code{solver_utils.R}.
#'
#' By contrast, correlation and block-sharing change the algebra enough that an
#' working-response/Fisher-scoring-style outer loop is more natural. In the non-Gaussian GEE
#' path, \code{\link{.get_B_gee_glm}} and \code{\link{.get_B_gee_glm_woodbury}}
#' work in the whitened full system because \eqn{\mathbf{V}^{-1/2}} couples the
#' partitions; these paths do use damped Newton updates, with the first iterate a
#' constrained Newton step and later iterates solved by \code{solve.QP}. In the
#' blockfit solver, \code{\link{blockfit_solve}} delegates to
#' \code{.bf_case_glm_no_corr} or \code{.bf_case_glm_gee}, both of which update
#' working responses and weights in damped Newton--Raphson outer loops while the
#' inner problem is solved by weighted backfitting. That is why the same
#' quantities \eqn{\mathbf{z}} and \eqn{\mathbf{W}} appear across GLM paths,
#' even though the no-correlation \code{get_B} path is best understood as
#' ``Newton--Raphson first, then project'', whereas the GEE and blockfit paths are
#' better understood as damped working-response iterations on a coupled system.
#'
#' It is useful to visualize the no-correlation extension as an
#' iterate-dependent projection built around a fixed unconstrained anchor. Let
#' \eqn{\hat{\boldsymbol{\beta}}} denote the partition-wise unconstrained
#' penalized estimate obtained before smoothness is imposed. In this path,
#' \eqn{\hat{\boldsymbol{\beta}}} is held fixed while the working information is
#' updated from the current constrained fit. Writing \eqn{\mathbf{G}^{(s)}} for
#' the recomputed information inverse, define
#' \deqn{\mathbf{U}^{(s)} =
#'   \mathbf{I} - \mathbf{G}^{(s)}\mathbf{A}
#'   \left\{\mathbf{A}^{\top}\mathbf{G}^{(s)}\mathbf{A}\right\}^{-1}
#'   \mathbf{A}^{\top}.}
#' The corresponding iterate can then be written as
#' \deqn{\tilde{\boldsymbol{\beta}}^{(s+1)} =
#'   \mathbf{U}^{(s)}\hat{\boldsymbol{\beta}}.}
#' This makes the non-identity-link extension easy to picture:
#' \eqn{\hat{\boldsymbol{\beta}}} does not move, but
#' \eqn{\mathbf{U}^{(s)}} and therefore
#' \eqn{\tilde{\boldsymbol{\beta}}^{(s)}} do, because the local curvature of the
#' penalized log-likelihood changes with the current mean function. So the
#' mathematics are closest to a fixed-point sequence of ``recompute weighted
#' information, then re-project'', rather than a single monolithic working-response fit of
#' the constrained problem.
#' }
#'
#' \subsection{Step Control}{
#' The code uses different step-control rules in different paths.
#'
#' In the dense non-Gaussian GEE path, the coefficient update is explicitly
#' damped:
#' \deqn{\boldsymbol{\beta}^{(s+1)} = (1 - \alpha_s)\boldsymbol{\beta}^{(s)} + \alpha_s\,\boldsymbol{\beta}_{\mathrm{cand}}^{(s)},}
#' with \eqn{\alpha_s = 2^{-d_s}} and \eqn{d_s} increased whenever the proposed
#' step gives non-finite or larger deviance. The loop stops after 10 consecutive
#' rejections or 100 total iterations.
#'
#' In the no-correlation non-Gaussian path, there is no analogous line search
#' on the outer projection loop. Instead, the algorithm repeatedly recomputes
#' \eqn{\mathbf{G}} at the current constrained estimate and applies a fresh
#' Lagrangian projection. The partition-wise unconstrained anchor itself is still
#' obtained by Newton--Raphson through \code{unconstrained_fit_fxn} (default
#' \code{\link{unconstrained_fit_default}} with helpers
#' \code{\link{damped_newton_r}} and \code{\link{nr_iterate}}). Outer
#' iteration stops when the mean absolute coefficient change falls below
#' \code{tol}, or when the update starts getting worse, in which case the
#' previous iterate is restored.
#' }
#'
#' \subsection{Convergence}{
#' Convergence monitoring also depends on the path.
#'
#' For dense non-Gaussian GEE, the accepted iterate must reduce the monitored
#' deviance, and the loop terminates once both the maximum coefficient change
#' and the deviance improvement fall below \code{tol} after a short burn-in.
#'
#' For the non-Gaussian no-correlation path, convergence is based on the mean
#' absolute coefficient change between successive constrained projections. If the
#' change begins to increase, the algorithm reverts to the previous iterate and
#' stops.
#'
#' For the blockfit GLM solvers, deviance is monitored across outer iterations,
#' but the coefficient update itself is produced by the weighted backfitting
#' solve rather than by an explicit convex combination with the previous iterate.
#' }
#'
#' \subsection{Estimation Paths}{
#' The implementation dispatches to one of three main paths depending on the
#' model configuration, primarily through \code{correlation_structure},
#' \code{VhalfInv}, \code{VhalfInv_fxn}, and \code{blockfit} in
#' \code{\link{lgspline}}.
#'
#' \strong{Path 1: Correlation structure present (GEE).}
#' The design matrix \eqn{\mathbf{X}} and response \eqn{\mathbf{y}} arrive
#' unwhitened; whitening is applied internally to the full \eqn{N \times P}
#' block-diagonal design to preserve cross-partition correlation.
#'
#' \emph{Path 1a} (Gaussian identity + GEE) solves the whitened penalized GLS
#' system directly. In the dense version,
#' \eqn{\tilde{\mathbf{G}} = (\tilde{\mathbf{X}}^{\top}\tilde{\mathbf{X}} +
#' \boldsymbol{\Lambda}_{\mathrm{block}})^{-1}} is formed explicitly and the full
#' Lagrangian projection is carried out in \eqn{P}-space. When
#' \eqn{\mathbf{V}^{-1} - \mathbf{I}} has low off-diagonal rank, the Woodbury
#' variant replaces the dense solve with a block-diagonal-plus-low-rank
#' representation.
#'
#' \emph{Path 1b} (non-Gaussian GEE) runs a damped SQP-style iteration on the
#' whitened system through \code{\link{.get_B_gee_glm}} or
#' \code{\link{.get_B_gee_glm_woodbury}}. The first iteration uses a
#' constrained Newton step obtained from the projection matrix \eqn{\mathbf{U}}.
#' Later iterations solve the quadratic approximation with
#' \code{\link[quadprog]{solve.QP}}. The Woodbury version precomputes
#' \eqn{(\mathbf{V}^{-1} - \mathbf{I})\mathbf{X}} once, recomputes only the
#' weighted low-rank correction at each damped Newton iteration, and falls
#' back to the dense path if that representation becomes invalid.
#'
#'
#' \strong{Path 2: Gaussian identity, no correlation.}
#' The constrained estimate is obtained by a single Lagrangian projection from
#' the per-partition penalized least-squares cross-products. No outer iteration
#' is needed. When \eqn{K = 0} and there are no additional constraints, this is
#' just the ordinary penalized closed form.
#'
#' \strong{Path 3: Non-Gaussian GLM, no correlation.}
#' Unconstrained estimates are obtained separately within partitions by
#' \code{unconstrained_fit_fxn}, by default
#' \code{\link{unconstrained_fit_default}} with Newton--Raphson helpers
#' \code{\link{damped_newton_r}} and \code{\link{nr_iterate}}. These are then
#' projected onto the equality constraint space by \code{\link{get_B}}. If
#' \code{iterate = TRUE}, the code updates the weighted
#' Gram matrices, Schur corrections, and square-root information factors at the
#' current constrained estimate and re-applies the projection until the fixed-
#' point iteration stabilizes. In the notation above, this is the path where
#' \eqn{\hat{\boldsymbol{\beta}}} is treated as fixed while
#' \eqn{\mathbf{G}^{(s)}}, \eqn{\mathbf{U}^{(s)}}, and hence
#' \eqn{\tilde{\boldsymbol{\beta}}^{(s)}} are updated until the projected
#' sequence settles.
#' }
#'
#' \subsection{Dispersion Estimation}{
#' Dispersion is handled through the user-visible functions
#' \code{need_dispersion_for_estimation} and \code{dispersion_function}.
#' Whenever \code{need_dispersion_for_estimation = TRUE}, the code recomputes a
#' current dispersion estimate from the present fitted means and passes that
#' value into the working-weight and Schur-correction calculations. Hence the
#' fitting code does not assume a single hard-coded Pearson estimator; it uses
#' whatever \code{dispersion_function} returns.
#'
#' For non-Gaussian fits, optional Schur corrections are computed from the
#' current iterate and added to the information matrix before eigendecomposition
#' or QP construction. Partition-wise no-correlation paths add these corrections
#' block by block; dense GEE paths add the collapsed correction to the full
#' whitened information matrix.
#' }
#'
#' @section Accommodating Correlation Structures:
#'
#' Up to this point, the presentation mirrors the independent-error version of
#' LMSS. Correlation enters as an extension of exactly the same constrained and
#' penalized framework rather than as a separate model class. The main change is
#' that whitening couples the partitions, so several of the computational
#' simplifications from the block-diagonal case no longer apply.
#'
#' \subsection{Parametric Correlation Structures}{
#' Suppose \eqn{\mathrm{Cov}(\mathbf{y}) = \sigma^2\,\mathbf{V}(\boldsymbol{\theta})}
#' for a known parametric family indexed by \eqn{\boldsymbol{\theta}} (e.g.,
#' AR(1) with \eqn{\theta = \rho}, Matern with
#' \eqn{\boldsymbol{\theta} = (\ell, \nu)}, exchangeable with
#' \eqn{\theta = \rho}). The penalized generalized least-squares problem
#' becomes
#' \deqn{\min_{\boldsymbol{\beta}}\;
#'   (\mathbf{y} - \mathbf{X}\boldsymbol{\beta})^{\top}\mathbf{V}^{-1}
#'   (\mathbf{y} - \mathbf{X}\boldsymbol{\beta})
#'   + \boldsymbol{\beta}^{\top}\boldsymbol{\Lambda}\boldsymbol{\beta}
#'   \quad \text{s.t.}\; \mathbf{A}^{\top}\boldsymbol{\beta} = \mathbf{0}.}
#' }
#'
#' \subsection{Loss of Block-Diagonal Structure}{
#' Unlike the independent-errors case,
#' \eqn{\mathbf{X}^{\top}\mathbf{V}^{-1}\mathbf{X}} is generally
#' \emph{not} block-diagonal because \eqn{\mathbf{V}^{-1}} introduces
#' cross-partition covariance. The unconstrained estimator
#' \deqn{\hat{\boldsymbol{\beta}} = (\mathbf{X}^{\top}\mathbf{V}^{-1}\mathbf{X} + \boldsymbol{\Lambda})^{-1}\mathbf{X}^{\top}\mathbf{V}^{-1}\mathbf{y}}
#' requires a full \eqn{P \times P} solve. However, when \eqn{\mathbf{V}^{-1}}
#' has banded or structured sparsity (as for AR(1) or compactly supported
#' correlation functions), the fill-in is limited and sparse matrix methods
#' remain efficient. For dense \eqn{\mathbf{V}^{-1}}, the cost is
#' \eqn{O(P^3)}, which remains tractable for the moderate values of \eqn{P}
#' typical in smoothing spline applications.
#'
#' The constraint projection proceeds identically, with
#' \eqn{\mathbf{G} = (\mathbf{X}^{\top}\mathbf{V}^{-1}\mathbf{X} + \boldsymbol{\Lambda})^{-1}}.
#' }
#'
#' \subsection{Whitening and Permutation}{
#' In the implementation, the correlation matrix \eqn{\mathbf{V}} is supplied
#' through the fitted-object components \code{Vhalf} and \code{VhalfInv},
#' either directly or from user functions \code{Vhalf_fxn} and
#' \code{VhalfInv_fxn}. Because the data are stored in partition ordering
#' (all observations from partition 0, then partition 1, etc.), while
#' \eqn{\mathbf{V}} is in the original observation ordering, a permutation is
#' applied:
#' \eqn{\mathbf{V}_{\mathrm{perm}}^{-1/2} = \mathbf{V}^{-1/2}[\boldsymbol{\pi}, \boldsymbol{\pi}]},
#' where \eqn{\boldsymbol{\pi}} maps original indices to partition-ordered
#' indices. The whitened design and response are then
#' \eqn{\tilde{\mathbf{X}} = \mathbf{V}_{\mathrm{perm}}^{-1/2}\mathbf{X}_{\mathrm{block}}}
#' and
#' \eqn{\tilde{\mathbf{y}} = \mathbf{V}_{\mathrm{perm}}^{-1/2}\mathbf{y}}.
#'
#' The \eqn{\mathbf{X}} and \eqn{\mathbf{y}} inputs to \code{\link{lgspline.fit}}
#' are preserved in their unwhitened form. Whitening is applied internally
#' within \code{\link{get_B}} and \code{\link{blockfit_solve}} where the full
#' \eqn{N \times P} structure is available, since applying
#' \eqn{\mathbf{V}^{-1/2}} to only the diagonal blocks of the partitioned
#' design matrix would silently discard cross-partition contributions and
#' corrupt the Gram matrix. This was a bug present in versions < 1.0. For
#' built-in correlation structures, the required square-root matrices are
#' assembled numerically with helper linear-algebra routines such as
#' \code{\link{matsqrt}}, \code{\link{matinvsqrt}}, and \code{\link{invert}};
#' these same helpers are also used when \code{Vhalf} must be recovered from a
#' user-supplied \code{VhalfInv}, or vice versa.
#' }
#'
#' \subsection{GEE Deviance Monitoring}{
#' For non-Gaussian models with correlation, the deviance used for convergence
#' monitoring is computed in the whitened (decorrelated) space. The raw deviance
#' residuals \eqn{r_i = \mathrm{sign}(d_i)\sqrt{|d_i|}} are divided by
#' \eqn{\sqrt{w_i}} and premultiplied by \eqn{\mathbf{V}_{\mathrm{perm}}^{-1/2}}
#' before squaring and averaging:
#' \deqn{D_{\mathrm{GEE}} = \frac{1}{N}\left\|\mathbf{V}_{\mathrm{perm}}^{-1/2}\,
#'   \mathrm{diag}(\mathbf{w})^{-1/2}\mathbf{r}\right\|^{2},}
#' where \eqn{\mathbf{w}} is the vector of Newton--Raphson working weights at
#' the current iterate, clamped below at \eqn{\sqrt{\varepsilon_{\mathrm{mach}}}}.
#' }
#'
#' \subsection{REML Estimation of Correlation Parameters}{
#' Correlation parameters \eqn{\boldsymbol{\theta}} are estimated by minimizing
#' a negative restricted log-likelihood (REML) objective. The criterion
#' implemented in \pkg{lgspline} is a central-limit-theorem-based working
#' approximation to a Laplace-style marginal likelihood criterion,
#' applied here solely to correlation structure estimation rather than
#' penalty parameter selection.
#' Let \eqn{\mathbf{D} = \mathrm{diag}(d_i)} be the observation weight matrix,
#' \eqn{\tilde{\mathbf{W}} = \mathrm{diag}(\tilde{w}_i)} the GLM working weight
#' matrix at the current fitted values, \eqn{\mathbf{V}} the correlation matrix
#' parameterized by \eqn{\boldsymbol{\rho}} (a vector on the unconstrained
#' real line), and \eqn{\tilde{\sigma}^{2}} the dispersion profiled at its
#' restricted maximum likelihood estimate. The negative REML objective
#' implemented in \pkg{lgspline}, scaled by \eqn{1/N}, is
#' \deqn{-\ell_{R}(\boldsymbol{\rho}) = \frac{1}{N}\!\left[
#'   -\log|\mathbf{V}^{-1/2}|
#'   + \frac{N}{2}\log\tilde{\sigma}^{2}
#'   + \frac{1}{2\tilde{\sigma}^{2}}
#'     (\mathbf{y} - \boldsymbol{\mu})^{\top}
#'     \mathbf{D}\tilde{\mathbf{W}}^{-1}\mathbf{V}^{-1}
#'     (\mathbf{y} - \boldsymbol{\mu})
#'   + \frac{1}{2}\log|(\tilde{\sigma}^{2}\mathbf{U}\mathbf{G})^{-1}|^{+}
#' \right],}
#' where \eqn{|\cdot|^{+}} denotes the generalized determinant (product of
#' nonzero eigenvalues), and
#' \eqn{\boldsymbol{\mu} = g^{-1}(\mathbf{X}\tilde{\boldsymbol{\beta}})} are
#' the fitted values on the response scale.
#'
#' Gradients with respect to correlation parameters are available in
#' closed form for all built-in structures except Matern, which uses
#' finite-difference approximation due to the complexity of differentiating
#' the modified Bessel function \eqn{K_{\nu}} with respect to \eqn{\nu}. See
#' \code{\link{reml_grad_from_dV}} for the full gradient derivation and
#' notation. Custom analytic gradients can be supplied through \code{REML_grad},
#' and a fully custom criterion can replace REML through
#' \code{custom_VhalfInv_loss}. The Toeplitz example in \code{\link{lgspline}}
#' demonstrates how to supply custom correlation structures with user-defined
#' gradient functions. Optimization over these working correlation parameters is
#' then carried out by the same quasi-Newton engine used elsewhere in the
#' package, namely \code{\link{efficient_bfgs}} with fallback to
#' \code{\link{approx_grad}} when needed. In the user-facing interface, this
#' machinery is activated through \code{correlation_structure},
#' \code{correlation_id}, \code{spacetime}, \code{VhalfInv}, \code{Vhalf},
#' \code{VhalfInv_fxn}, \code{Vhalf_fxn}, \code{VhalfInv_par_init},
#' \code{REML_grad}, \code{custom_VhalfInv_loss}, and
#' \code{VhalfInv_logdet}.
#'
#' The gradient of the negative REML has three terms per parameter:
#' \enumerate{
#'   \item \eqn{\frac{1}{2}\mathrm{tr}(\mathbf{V}^{-1}\partial\mathbf{V}/\partial\theta_j)}:
#'     the log-determinant contribution.
#'   \item \eqn{-\frac{1}{2\tilde{\sigma}^{2}}\mathbf{r}^{\top}
#'     (\partial\mathbf{V}/\partial\theta_j)\mathbf{r}}:
#'     the residual quadratic form contribution, where
#'     \eqn{\mathbf{r} = \mathrm{diag}(\sqrt{d_i}/\sqrt{\tilde{w}_i})
#'     \mathbf{V}^{-1/2}(\mathbf{y} - \boldsymbol{\mu})}.
#'   \item \eqn{-\frac{1}{2}\mathrm{tr}\!\left(\mathbf{M}^{+}\mathbf{X}_{*}^{\top}
#'     \mathbf{V}^{-1}(\partial\mathbf{V}/\partial\theta_j)\mathbf{V}^{-1}
#'     \tilde{\mathbf{W}}\mathbf{D}\mathbf{X}_{*}\right)}:
#'     the REML correction, where
#'     \eqn{\mathbf{X}_{*} = \mathbf{X}\mathbf{U}} is the constrained design and
#'     \eqn{\mathbf{M} = \mathbf{X}_{*}^{\top}\mathbf{V}^{-1}
#'     \tilde{\mathbf{W}}\mathbf{D}\mathbf{X}_{*}
#'     + \mathbf{U}^{\top}\boldsymbol{\Lambda}\mathbf{U}} is the projected
#'     penalized information.
#' }
#'
#' For each supported correlation family, the derivatives
#' \eqn{\partial\mathbf{V}/\partial\theta_j} are available in closed form,
#' enabling analytic gradient computation for use with the quasi-Newton
#' optimizer.
#' }
#'
#' \subsection{Connection to Standard Mixed Model REML}{
#' The penalized spline model admits a standard mixed model conditional
#' representation. The quadratic penalty \eqn{\boldsymbol{\Lambda}} acts as the
#' inverse prior covariance of a Gaussian random effect on the spline
#' coefficients, with higher-order basis components treated as random effects
#' with variance \eqn{\tau^{2}} and unpenalized components (intercept, low-order
#' polynomial terms) acting as fixed effects. The smoothing parameter satisfies
#' \eqn{\lambda = \tau^{2}/\sigma^{2}}. This is the same approach used by
#' \pkg{mgcv} for fitting random effects; from the Bayesian perspective, the
#' difference between conditional random effects models and fixed effects with
#' normal-mean-0 priors is just how the ``penalty'' is obtained (i.e., REML, ML,
#' or GCV). Users can obtain Monte Carlo draws under this Laplace-style
#' approximation with the \code{\link{generate_posterior}} function. In
#' practice, with a non-0 \code{flat_ridge_penalty}, even intercepts and
#' linear terms can be interpreted as random effects.
#'
#' For non-Gaussian responses the criterion is motivated by a working quadratic
#' approximation to the log-likelihood. Under standard regularity conditions,
#' \eqn{\tilde{\mathbf{W}}^{-1/2}(\mathbf{y} - \boldsymbol{\mu})} is
#' approximately Gaussian by the central limit theorem, yielding the quadratic
#' term
#' \eqn{(\mathbf{y} - \boldsymbol{\mu})^{\top}\mathbf{D}\tilde{\mathbf{W}}^{-1}\mathbf{V}^{-1}(\mathbf{y} - \boldsymbol{\mu})}.
#' This is the sense in which the REML criterion is a CLT-based working
#' approximation: rather than the exact Laplace approximation to the marginal
#' likelihood, which would require the full penalized log-likelihood Hessian,
#' the implemented criterion substitutes the Fisher information evaluated at
#' the current working estimates. This can be viewed as a profile
#' quasi-likelihood or method-of-moments estimator for the correlation
#' parameters. When the response is Gaussian with identity link, the
#' approximation is exact and coincides with the classical REML criterion.
#'
#' The REML correction term
#' \eqn{\log|(\tilde{\sigma}^{2}\mathbf{U}\mathbf{G})^{-1}|^{+}} is the
#' generalized log-determinant of the precision matrix associated with the
#' penalized coefficients. It plays the same role as the classical mixed model
#' term
#' \eqn{\log|\mathbf{X}^{\top}\mathbf{V}^{-1}\mathbf{X} + \boldsymbol{\Lambda}|},
#' with the generalized determinant ensuring that only nonzero eigenvalues
#' contribute when rank deficiency arises from smoothness constraints or
#' identifiability conditions encoded in \eqn{\mathbf{A}}.
#'
#' When \eqn{\boldsymbol{\Lambda}} is full rank, the criterion coincides with
#' the exact marginal likelihood obtained by integrating out
#' \eqn{\boldsymbol{\beta}} under its Gaussian prior. When
#' \eqn{\boldsymbol{\Lambda}} is rank-deficient, the parameter space separates:
#' unpenalized coefficients are projected out in the REML sense while penalized
#' coefficients are integrated out through their prior, and additional linear
#' constraints in \eqn{\mathbf{A}} further reduce the effective parameter
#' dimension accordingly.
#'
#' Setting \eqn{\mathbf{U} = \mathbf{I}}, \eqn{\boldsymbol{\Lambda} = \mathbf{0}},
#' and \eqn{\tilde{\mathbf{W}} = \mathbf{I}} recovers the classical Gaussian
#' linear mixed model REML criterion.
#'
#' During penalty tuning, the block-diagonal approximation is retained for
#' computing GCV criteria and gradients, for efficiency. Since GCV is
#' rotation-invariant, the practical effect on automatic penalty selection is
#' expected to be negligible, though this has not been formally confirmed. The
#' tuned penalties can always be overridden by the user.
#' }
#'
#' \subsection{Built-In Correlation Structures}{
#' The package provides several built-in correlation structures for modeling
#' spatial and temporal dependence. These are specified via
#' \code{correlation_structure} with group membership in \code{correlation_id}
#' and spatial or temporal coordinates in \code{spacetime} (an \eqn{N}-row
#' matrix). Exchangeable correlation does not require \code{spacetime}.
#'
#' All positive scale parameters are estimated on the log scale,
#' with back-transform \eqn{\exp(\cdot)}. Parameters constrained to
#' \eqn{(0, 1)} use a double-exponential back-transform of the form
#' \eqn{\exp(-\exp(\eta))}, so optimization still occurs on the
#' unconstrained real line while the correlation remains bounded.
#'
#' \describe{
#'   \item{Exchangeable}{
#'     Aliases: \code{'exchangeable'}, \code{'cs'}, \code{'CS'},
#'     \code{'compoundsymmetric'}, \code{'compound-symmetric'}.
#'
#'     A constant correlation \eqn{\nu} between any two observations within
#'     the same cluster. Parameterization: \eqn{\nu = \exp(-\exp(\rho))},
#'     so \eqn{\nu \in (0, 1)}. Only positive within-cluster correlation is
#'     supported under this parameterization.
#'   }
#'   \item{Spatial Exponential}{
#'     Aliases: \code{'spatial-exponential'}, \code{'spatialexponential'},
#'     \code{'exp'}, \code{'exponential'}.
#'
#'     Correlation decays exponentially with distance:
#'     \eqn{\exp(-\omega d)} where \eqn{d} is Euclidean distance and
#'     \eqn{\omega > 0}. Parameterization: \eqn{\omega = \exp(\rho)}.
#'     Mathematically equivalent to the power correlation \eqn{\theta^{d}}
#'     with \eqn{\theta = e^{-\omega}}, but with better numerical properties
#'     during optimization.
#'   }
#'   \item{AR(1)}{
#'     Aliases: \code{'ar1'}, \code{'ar(1)'}, \code{'AR(1)'}, \code{'AR1'}.
#'
#'     Correlation depends on rank difference between observations:
#'     \eqn{\nu^{r}} where \eqn{r} is the rank difference within cluster.
#'     Parameterization: \eqn{\nu = \exp(-\exp(\rho))},
#'     so \eqn{\nu \in (0, 1)}. Only positive autocorrelation is supported.
#'   }
#'   \item{Gaussian / Squared Exponential}{
#'     Aliases: \code{'gaussian'}, \code{'rbf'}, \code{'squared-exponential'}.
#'
#'     Smooth decay with squared distance:
#'     \eqn{\exp(-d^{2}/(2\ell^{2}))} where \eqn{\ell} is the length scale.
#'     Parameterization: \eqn{\ell = \exp(\rho)}.
#'   }
#'   \item{Spherical}{
#'     Aliases: \code{'spherical'}, \code{'Spherical'}, \code{'cubic'},
#'     \code{'sphere'}.
#'
#'     Polynomial decay with a hard cutoff at range \eqn{r}:
#'     \eqn{1 - 1.5(d/r) + 0.5(d/r)^{3}} for \eqn{d \le r}, and \eqn{0}
#'     otherwise. Parameterization: \eqn{r = \exp(\rho)}.
#'   }
#'   \item{Matern}{
#'     Aliases: \code{'matern'}, \code{'Matern'}.
#'
#'     Flexible correlation with adjustable smoothness:
#'     \eqn{(2^{1-\nu}/\Gamma(\nu))(\sqrt{2\nu}\,d/\ell)^{\nu}K_{\nu}(\sqrt{2\nu}\,d/\ell)}.
#'     Two parameters: length scale \eqn{\ell = \exp(\rho_1)} and smoothness
#'     \eqn{\nu = \exp(\rho_2)}. No analytical gradient is available for
#'     \eqn{\nu} due to the difficulty of differentiating the modified Bessel
#'     function \eqn{K_{\nu}} with respect to its order, so finite differences
#'     are used; this makes Matern slower and potentially less stable than
#'     other structures.
#'   }
#'   \item{Gamma-Cosine}{
#'     Aliases: \code{'gamma-cosine'}, \code{'gammacosine'},
#'     \code{'GammaCosine'}.
#'
#'     Oscillatory dependence:
#'     \eqn{(d^{\alpha-1}e^{-\gamma d})/(\Gamma(\alpha)/\gamma^{\alpha})\cdot\cos(\omega d)}.
#'     Three parameters: shape \eqn{\alpha = \exp(\rho_1)}, rate
#'     \eqn{\gamma = \exp(\rho_2)}, frequency \eqn{\omega = \exp(\rho_3)}.
#'     Reduces to exponential when \eqn{\alpha = 1} and \eqn{\omega \approx 0}.
#'   }
#'   \item{Gaussian-Cosine}{
#'     Aliases: \code{'gaussian-cosine'}, \code{'gaussiancosine'},
#'     \code{'GaussianCosine'}.
#'
#'     Smooth oscillatory correlation:
#'     \eqn{\exp(-d^{2}/(2\ell^{2}))\cdot\cos(\omega d)}.
#'     Two parameters: length scale \eqn{\ell = \exp(\rho_1)} and frequency
#'     \eqn{\omega = \exp(\rho_2)}. Reduces to Gaussian when
#'     \eqn{\omega \approx 0}.
#'   }
#' }
#' }
#'
#' \subsection{Interpreting Estimated Correlation Parameters}{
#' Correlation parameters are estimated on transformed scales; they must be
#' back-transformed for interpretation. When \code{\link{confint.lgspline}} is called and the
#' inverse Hessian from BFGS is available, confidence intervals are returned
#' on the untransformed (working) scale and should be back-transformed as
#' described in the examples for \code{\link{lgspline}}.
#' }
#'
#' \subsection{Custom Correlation Structures}{
#' Custom correlation structures can be specified through:
#' \itemize{
#'   \item \code{VhalfInv_fxn}: Creates \eqn{\mathbf{V}^{-1/2}}.
#'   \item \code{Vhalf_fxn}: Creates \eqn{\mathbf{V}^{1/2}}. When omitted,
#'     the code computes it by explicit inversion of \code{VhalfInv}.
#'   \item \code{REML_grad}: Provides the analytical gradient of the REML
#'     objective.
#'   \item \code{VhalfInv_logdet}: Efficient log-determinant computation.
#'   \item \code{custom_VhalfInv_loss}: Replaces the REML objective entirely.
#' }
#' These functions enter \code{\link{lgspline}} through
#' \code{correlation_structure}, \code{VhalfInv_fxn}, \code{Vhalf_fxn},
#' \code{REML_grad}, and \code{custom_VhalfInv_loss}, and the fitted object
#' retains the resulting correlation machinery in components such as
#' \code{VhalfInv_fxn}, \code{Vhalf_fxn}, and
#' \code{VhalfInv_params_estimates}. When \code{VhalfInv} is supplied but
#' \code{Vhalf} is not, \code{Vhalf} is computed unconditionally as the inverse
#' of \code{VhalfInv} for all family/link combinations, since both
#' \code{\link{get_B}} and \code{\link{blockfit_solve}} require it for GEE
#' estimation.
#' }
#'
#' @section Variance Estimation and Inference:
#'
#' Once the constrained estimate has been obtained, the next questions are how
#' much flexibility the fitted model effectively used and how uncertainty should
#' be propagated through the same constrained geometry. The quantities in this
#' section are therefore all built on the projected information matrices from
#' the previous sections.
#'
#' \subsection{Effective Degrees of Freedom and Dispersion}{
#' The effective degrees of freedom is the trace of the hat matrix. In the
#' Gaussian identity case with observation weights and correlation, the fitted
#' linear operator is built from the dense GLS analogue
#' \eqn{\mathbf{G}_{\mathrm{correct}}} and can be written schematically as
#' \deqn{\mathbf{H} =
#'   \mathbf{V}^{-1/2}(\mathbf{W}\mathbf{D})^{1/2}
#'   \mathbf{X}\mathbf{U}\mathbf{G}_{\mathrm{correct}}\mathbf{X}^{\top}
#'   (\mathbf{W}\mathbf{D})^{1/2}\mathbf{V}^{-1/2},}
#' where for Gaussian identity \eqn{\mathbf{W} = \mathbf{I}}. In the
#' no-correlation Gaussian case this reduces to the familiar
#' \eqn{\mathbf{H} = \mathbf{X}\mathbf{U}\mathbf{G}\mathbf{X}^{\top}\mathbf{D}}.
#'
#' For Gaussian identity fits, the dispersion estimate is computed as a
#' weighted mean squared residual, optionally scaled by
#' \eqn{N/(N - \mathrm{tr}(\mathbf{H}))} when
#' \code{unbias_dispersion = TRUE}:
#' \deqn{\tilde{\sigma}^{2} =
#'   \frac{1}{N - \mathrm{tr}(\mathbf{H})}\|\mathbf{y} - \mathbf{\tilde{y}}_i \|^{2}.}
#'
#' More generally with weights, a correlation structure and non-linear link function:
#' \deqn{\tilde{\sigma}^{2} =
#'   \frac{1}{N - \mathrm{tr}(\mathbf{H})}\| \mathbf{V}^{-1/2}\mathbf{W}^{-1/2}\mathbf{D}^{1/2}(\mathbf{y} - \tilde{\mathbf{y}})\|^{2}.}
#'
#' This estimated dispersion is returned as \code{sigmasq_tilde}, and the
#' corresponding effective degrees of freedom trace is returned as
#' \code{trace_XUGX}. For non-Gaussian families, the fitting code delegates
#' dispersion estimation to \code{dispersion_function}; thus the package does
#' not assume a single closed-form Pearson-style formula outside the Gaussian
#' identity setting. The hat-matrix trace itself is assembled by
#' \code{\link{compute_trace_H}} in the dense correlation-aware case and by the
#' same blockwise products summarized by \code{trace_XUGX} in the simpler
#' no-correlation paths.
#' A concrete built-in non-Gaussian example is the Weibull AFT path, which pairs
#' \code{\link{weibull_family}} with \code{\link{weibull_dispersion_function}},
#' \code{\link{weibull_glm_weight_function}}, and
#' \code{\link{weibull_schur_correction}}.
#' Users who want these quantities available for downstream inference should
#' keep \code{estimate_dispersion = TRUE} and \code{return_varcovmat = TRUE}
#' (the defaults), since \code{\link{wald_univariate}},
#' \code{\link{confint.lgspline}}, and the prediction-standard-error path in
#' \code{\link{predict.lgspline}} all rely on the post-fit dispersion and
#' covariance components documented here.
#' }
#'
#' \subsection{Variance--Covariance Matrix}{
#' The variance--covariance matrix of \eqn{\tilde{\boldsymbol{\beta}}} is
#' estimated as:
#' \deqn{\mathrm{Var}(\tilde{\boldsymbol{\beta}}) = \tilde{\sigma}^{2}(\mathbf{U}\mathbf{G}^{1/2})(\mathbf{U}\mathbf{G}^{1/2})^{\top}}
#' using the outer-product form for numerical stability. The result is returned
#' as \code{varcovmat} when \code{return_varcovmat = TRUE}. The algebraically
#' equivalent expression \eqn{\tilde{\sigma}^{2}\mathbf{U}\mathbf{G}} is not
#' used because \eqn{\mathbf{G}} is only positive semi-definite when the
#' penalty matrix \eqn{\boldsymbol{\Lambda}} has zero eigenvalues (e.g., the
#' intercept and linear terms under the smoothing spline penalty when
#' \code{flat_ridge_penalty} = 0), which can introduce negative diagonal
#' entries in finite precision arithmetic. The outer-product form also
#' guarantees symmetry.
#'
#' This is the Bayesian posterior covariance, treating the penalty as a
#' Gaussian prior on the coefficients. When \code{exact_varcovmat = TRUE},
#' a frequentist correction is additionally computed:
#' \deqn{\mathrm{Var}_{\mathrm{exact}}(\tilde{\boldsymbol{\beta}}) =
#'   \tilde{\sigma^2}\mathbf{U}\mathbf{G}^{1/2}(\mathbf{X}^{\top}\mathbf{W}\mathbf{D}\mathbf{V}^{-1}\mathbf{X})\mathbf{G}^{1/2}\mathbf{U}^{\top} =
#'   \tilde{\sigma}^{2}\mathbf{U}\mathbf{G}\mathbf{U}^{\top}
#'   - \tilde{\sigma}^{2}\mathbf{U}\mathbf{G}\boldsymbol{\Lambda}\mathbf{G}\mathbf{U}^{\top}.}
#' The first term is the Bayesian posterior covariance; the second is a
#' bias correction that accounts for the penalty-induced shrinkage.
#' For Gaussian identity link (with or without correlation), this is the
#' exact variance--covariance matrix of the constrained estimator.
#'
#' When a correlation structure is present (\code{VhalfInv} non-\code{NULL}),
#' the block-diagonal \eqn{\mathbf{G}} is replaced by the full weighted GLS
#' analogue
#' \deqn{\mathbf{G}_{\mathrm{correct}} =
#'   \left(\mathbf{X}^{\top}\mathbf{W}\mathbf{D}\mathbf{V}^{-1}\mathbf{X}
#'   + \boldsymbol{\Lambda}\right)^{-1},}
#' where \eqn{\mathbf{W} = \mathbf{I}} in the Gaussian identity case.
#' This dense matrix is what enters the correlation-aware \eqn{\mathbf{U}},
#' \eqn{\mathrm{Var}(\tilde{\boldsymbol{\beta}})}, and
#' \eqn{\mathrm{Var}_{\mathrm{exact}}(\tilde{\boldsymbol{\beta}})}
#' computations.
#'
#' In user-facing terms, \code{return_varcovmat} controls whether this matrix is
#' stored at all, while \code{exact_varcovmat} switches between the default
#' posterior/Laplace approximation and the exact frequentist correction in the
#' Gaussian-identity setting. The stored covariance is what powers
#' \code{\link{wald_univariate}}, \code{\link{confint.lgspline}}, and
#' \code{se.fit = TRUE} in \code{\link{predict.lgspline}}; the
#' \code{critical_value} argument supplied at fit time is carried forward as the
#' default cutoff for those interval-producing helpers.
#' }
#'
#' \subsection{Recomputation of G at Convergence}{
#' At the final iterate, \eqn{\mathbf{G}} is recomputed to reflect the
#' converged working weights and Schur corrections. The implementation
#' computes the weighted design
#' \eqn{\mathbf{X}_{w}^{(k)} = \mathbf{X}_k \cdot \mathrm{diag}(\sqrt{\mathbf{w}_k})},
#' forms the weighted Gram matrix
#' \eqn{\mathbf{X}_{w}^{(k)\top}\mathbf{X}_{w}^{(k)}}, adds the Schur
#' correction, and performs eigendecomposition via
#' \code{\link{compute_G_eigen}} to obtain \eqn{\mathbf{G}_k},
#' \eqn{\mathbf{G}_k^{1/2}}, and \eqn{\mathbf{G}_k^{-1/2}}. The relationship
#' \eqn{\mathbf{G}_k = \mathbf{G}_k^{1/2}(\mathbf{G}_k^{1/2})^{\top}} is
#' enforced exactly by construction, and the fitted object can retain these as
#' \code{G} and \code{Ghalf} when \code{return_G = TRUE} and
#' \code{return_Ghalf = TRUE}. The square-root factors are numerically
#' stabilized through the helper routines \code{\link{matsqrt}} and
#' \code{\link{matinvsqrt}}, which are also used elsewhere in the package when
#' dense GLS analogues of \eqn{\mathbf{G}^{1/2}} or \eqn{\mathbf{G}^{-1/2}} are
#' required.
#' }
#'
#' @section Bayesian Interpretation:
#'
#' The penalty has a natural Gaussian-prior interpretation, so once the
#' constrained estimator and its covariance are available, Bayesian-style
#' posterior simulation follows almost immediately. This section records the
#' interpretation that is already implicit in the fitted object and in the
#' package's posterior simulation helpers.
#'
#' A Bayesian interpretation follows from viewing the penalty as a Gaussian
#' prior on the coefficients. Conditional on the fitted smoothing parameters,
#' the code samples on the coefficient scale from
#' \deqn{\boldsymbol{\beta}^{(m)} =
#'   \tilde{\boldsymbol{\beta}} +
#'   \sqrt{\tilde{\sigma}^{2}}\,\mathbf{F}_{\mathrm{post}}\mathbf{z}^{(m)},
#'   \qquad \mathbf{z}^{(m)} \sim \mathcal{N}(\mathbf{0}, \mathbf{I}),}
#' where \eqn{\mathbf{F}_{\mathrm{post}}} is a square-root factor of the
#' constrained covariance. In the no-correlation paths one can take
#' \eqn{\mathbf{F}_{\mathrm{post}} = \mathbf{U}\mathbf{G}^{1/2}}; when a
#' correlation structure is present the same role is played by the dense GLS
#' analogue built from \eqn{\mathbf{U}\mathbf{G}_{\mathrm{correct}}^{1/2}}.
#' The coefficients are then back-transformed to the original response and
#' predictor scales.
#'
#' When inequality constraints are absent, these draws are i.i.d. Gaussian
#' posterior draws around the fitted mode.
#' The underlying coefficient-draw closure also contains an elliptical slice
#' sampling route for active inequality constraints, using the same covariance
#' factor to keep retained draws in the feasible region.
#'
#' At the implementation level, standard Gaussian posterior draws may place
#' positive mass on the infeasible region, so the constrained-draw closure
#' instead targets the truncated posterior
#' \deqn{\pi(\boldsymbol{\beta} \mid \mathbf{y}) \propto
#'   \exp\!\left(-\frac{1}{2}(\boldsymbol{\beta} - \tilde{\boldsymbol{\beta}})^{\top}
#'   \mathbf{G}^{-1}(\boldsymbol{\beta} - \tilde{\boldsymbol{\beta}})\right)
#'   \mathbf{1}(\mathbf{C}^{\top}\boldsymbol{\beta} \succeq \mathbf{c}),}
#' yielding credible intervals that respect the constraint boundaries.
#' The public \code{\link{generate_posterior}} wrapper forwards
#' \code{enforce_qp_constraints} to the stored constrained-draw closure, so
#' constrained draws can be requested directly from the user-facing interface.
#' When a working correlation structure is present, the companion helper
#' \code{\link{generate_posterior_correlation}} extends this idea by propagating
#' uncertainty in the fitted correlation parameters through the same
#' \code{VhalfInv_fxn}/\code{Vhalf_fxn} machinery described in the correlation
#' section, rather than conditioning only on fixed covariance parameters. Correlation
#' parameters are drawn from a multivariate normal distribution centered about their estimates
#' with the inverse approxiamte BFGS Hessian of the REML optimization problem
#' \code{VhalfInv_params_vcov} used by default (or a custom alternative as
#' supplied to the argument \code{correlation_param_vcov_sc}).
#'
#'
#' @section Inequality Constraints via Sequential Quadratic Programming:
#'
#' The equality constraints that define smoothness are only part of the story.
#' In many applications, one also wants the fitted surface to satisfy shape
#' restrictions such as monotonicity, convexity, or boundedness. Because the
#' package works directly in a polynomial coefficient basis, many of these
#' restrictions remain linear in the coefficients and can therefore be handled
#' within the same constrained optimization framework.
#'
#' \subsection{Overview}{
#' Sometimes we want to impose not just smoothness (equality constraints) but
#' also shape restrictions: the fitted curve should be monotone, or stay
#' within certain bounds, or be concave-decreasing. These are expressible as
#' inequality constraints on the coefficients.
#'
#' Inequality constraints of the form
#' \eqn{\mathbf{C}^{\top}\boldsymbol{\beta} \succeq \mathbf{c}} are handled
#' via damped sequential quadratic programming (SQP). The idea is to
#' repeatedly solve a quadratic approximation of the constrained
#' log-likelihood, where each subproblem is a standard QP that respects both
#' the smoothness equalities and the shape inequalities simultaneously. SQP
#' is the natural generalization of damped Newton--Raphson to problems with
#' inequality constraints. In the implementation, the built-in inequality
#' pieces are first assembled by \code{\link{process_qp}}, which returns the
#' \code{qp_Amat}, \code{qp_bvec}, and \code{qp_meq} objects ultimately passed
#' into \code{\link[quadprog]{solve.QP}}, together with a \code{quadprog} flag
#' indicating whether any inequality constraints are active at all.
#'
#' In the monomial basis, shape constraints are linear in
#' \eqn{\boldsymbol{\beta}}. For instance, monotonicity at a grid of points
#' \eqn{\{x_1^{*}, \ldots, x_M^{*}\}} requires
#' \eqn{\mathbf{C}^{\top}\boldsymbol{\beta} \succeq \mathbf{0}}, where column
#' \eqn{m} of \eqn{\mathbf{C}} evaluates the derivative polynomial at
#' \eqn{x_m^{*}}. Range constraints and sign constraints on second
#' derivatives (convexity: \eqn{f''(x) \geq 0}) are similarly linear in
#' \eqn{\boldsymbol{\beta}}. For the built-in derivative-sign constraints,
#' \code{\link{process_qp}} delegates to \code{.build_deriv_qp}, which calls
#' \code{\link{make_derivative_matrix}} on the expansion-standardized design and then
#' maps those derivative rows back into the full \eqn{P}-dimensional
#' coefficient space one partition block at a time.
#'
#' Writing the penalized objective as
#' \eqn{\ell(\boldsymbol{\beta}) - \frac{1}{2\sigma^{2}}\boldsymbol{\beta}^{\top}\boldsymbol{\Lambda}\boldsymbol{\beta}},
#' take a second-order Taylor expansion of \eqn{\ell} around a current iterate
#' \eqn{\boldsymbol{\beta}^{*}}. Since the Hessian of the penalized
#' log-likelihood is \eqn{-\sigma^{-2}\mathbf{G}^{*-1}}, the expansion yields
#' a quadratic objective in \eqn{\boldsymbol{\beta}}. Collecting terms, the
#' subproblem at each iteration is:
#' \deqn{\tilde{\boldsymbol{\beta}} = \arg\min_{\boldsymbol{\beta}}
#'   \left\{-\mathbf{d}^{\top}\boldsymbol{\beta} +
#'   \frac{1}{2\sigma^{2}}\boldsymbol{\beta}^{\top}\mathbf{G}^{-1}\boldsymbol{\beta}\right\}
#'   \quad \text{s.t.} \quad \mathbf{A}^{\top}\boldsymbol{\beta} = \mathbf{0}, \quad
#'   \mathbf{C}^{\top}\boldsymbol{\beta} \succeq \mathbf{c}}
#' where
#' \eqn{\mathbf{d} = \nabla_{\boldsymbol{\beta}}\ell(\boldsymbol{\beta}^{*}) + \sigma^{-2}\mathbf{G}^{*-1}\boldsymbol{\beta}^{*}}.
#' In implementation, the linear term is built from \code{qp_score_function},
#' whose default is the canonical GLM score \eqn{\mathbf{X}^{\top}(\mathbf{y}-\boldsymbol{\mu})}
#' documented in \code{\link{lgspline}}. For non-canonical or custom models a
#' different score can be supplied, for example
#' \code{\link{weibull_qp_score_function}} together with the matching
#' \code{glm_weight_function}, \code{dispersion_function},
#' \code{schur_correction_function}, and, when needed,
#' \code{unconstrained_fit_fxn}. The inequality side of the subproblem is built
#' from the objects returned by \code{\link{process_qp}}. In the current
#' implementation, custom assembled constraints should be supplied through
#' \code{qp_Amat_fxn}, \code{qp_bvec_fxn}, and \code{qp_meq_fxn}; the low-level
#' \code{qp_Amat}, \code{qp_bvec}, and \code{qp_meq} arguments currently serve
#' only as advanced placeholders / activation markers rather than being merged
#' into the built-in constructor.
#' This is a convex quadratic program, solvable by active-set methods
#' (Goldfarb and Idnani, 1983). The equality constraints
#' \eqn{\mathbf{A}^{\top}\boldsymbol{\beta} = \mathbf{0}} (smoothness)
#' are incorporated alongside the inequality constraints
#' \eqn{\mathbf{C}^{\top}\boldsymbol{\beta} \succeq \mathbf{c}} (shape) in
#' a single call to \code{\link[quadprog]{solve.QP}}. The combined constraint
#' matrix is \eqn{[\mathbf{A} \mid \mathbf{C}^{\top}]} with the first
#' \eqn{r} columns declared as equalities. When the inequality system is
#' partition-local, \code{\link{get_B}} and \code{\link{blockfit_solve}} first
#' try a partition-wise active-set refinement that treats active inequalities
#' as temporary equalities and reuses the Lagrangian projection machinery. If
#' that active-set loop does not converge within its fixed iteration limit, the
#' code falls back automatically to dense SQP.
#' }
#'
#' \subsection{Damped SQP Iteration}{
#' When the link is non-identity, the QP subproblem is embedded in a damped
#' SQP outer loop. At each SQP iteration \eqn{s}:
#' \enumerate{
#'   \item Compute the current information matrix
#'     \eqn{\mathbf{M}^{(s)} = \mathbf{X}^{\top}\mathbf{W}^{(s)}\mathbf{X} + \boldsymbol{\Lambda}_{\mathrm{block}} + \mathbf{S}^{(s)}},
#'     where \eqn{\mathbf{S}^{(s)}} is the Schur complement correction.
#'   \item Compute the score vector
#'     \eqn{\mathbf{s}^{(s)} = \mathbf{X}^{\top}\mathbf{W}^{(s)}(\mathbf{y} - \boldsymbol{\mu}^{(s)}) / g'(\boldsymbol{\mu}^{(s)})}.
#'   \item Solve the QP with Hessian \eqn{\mathbf{M}^{(s)}} and linear term
#'     \eqn{\mathbf{s}^{(s)} - \boldsymbol{\Lambda}_{\mathrm{block}}\boldsymbol{\beta}^{(s)} + \mathbf{M}^{(s)}\boldsymbol{\beta}^{(s)}}.
#'   \item Apply damped update:
#'     \eqn{\boldsymbol{\beta}^{(s+1)} = (1-\alpha)\boldsymbol{\beta}^{(s)} + \alpha\boldsymbol{\beta}_{\mathrm{QP}}},
#'     where \eqn{\alpha = 2^{-d}} and \eqn{d} is incremented upon deviance
#'     increase.
#' }
#' A rescaling factor \eqn{\mathrm{sc} = \sqrt{\mathrm{mean}(|\mathbf{M}^{(s)}|)}}
#' is applied to the Hessian and linear term before calling the QP solver,
#' improving numerical stability.
#'
#' The equality-constrained estimate (from the Lagrangian projection) serves
#' as a warm start for the first SQP iteration. The vector \eqn{\mathbf{d}}
#' is updated iteratively using the current estimate of
#' \eqn{\boldsymbol{\beta}} until convergence.
#'
#' When blockfitting is active, the SQP refinement loop runs after backfitting
#' convergence, using the backfitting solution as a warm start. In all of these
#' cases, the QP data themselves can be thinned to a user-specified subset of
#' rows through \code{qp_observations}, which \code{\link{process_qp}} applies
#' before derivative, range, and monotonicity constraints are assembled.
#' }
#'
#' \subsection{Active Set and Lagrange Multipliers}{
#' The active set at the solution identifies the binding inequality
#' constraints, and the associated Lagrange multipliers quantify the ``cost''
#' of each constraint. A multiplier of zero indicates that the constraint is
#' not binding at the solution. The implementation stores the active constraint
#' indices, the corresponding submatrix of the constraint matrix, and the
#' Lagrange multiplier vector in a \code{qp_info} list returned alongside the
#' coefficient estimates, with components such as \code{lagrangian},
#' \code{iact}, and \code{Amat_active}. When
#' \code{return_lagrange_multipliers = TRUE}, the fitted object also returns the
#' final multiplier vector directly as \code{lagrange_multipliers}.
#'
#' The original assembled inequality data are also retained in the fitted
#' object's \code{quadprog_list} component, so the final active set can be
#' interpreted relative to the full \code{qp_Amat}, \code{qp_bvec}, and
#' \code{qp_meq} specification returned by \code{\link{process_qp}}.
#'
#' The final \eqn{\mathbf{U}} returned and used in the construction of the
#' posterior variance--covariance matrix is constructed from equality
#' constraints and active inequality constraints, excluding inactive
#' constraints.
#' }
#'
#' \subsection{Built-In Constraints}{
#' Built-in inequality constraints include:
#' \itemize{
#'   \item Monotonicity: \code{qp_monotonic_increase},
#'     \code{qp_monotonic_decrease}. Enforced by requiring consecutive
#'     fitted values to be non-decreasing (or non-increasing):
#'     \eqn{(\mathbf{x}_i - \mathbf{x}_{i-1})^{\top}\boldsymbol{\beta} \geq 0}
#'     for all \eqn{i}. These are constructed in \code{\link{process_qp}} from
#'     the partition-stacked block design reordered back to observation order.
#'   \item Derivative sign: \code{qp_positive_derivative},
#'     \code{qp_negative_derivative}. Enforced through the first-derivative
#'     design matrix from \code{\link{make_derivative_matrix}}. These arguments
#'     may be \code{TRUE}/\code{FALSE} or a character/integer vector selecting
#'     specific predictors, with character names resolved against the original
#'     columns via \code{og_cols} inside \code{\link{process_qp}}.
#'   \item Second-derivative sign: \code{qp_positive_2ndderivative},
#'     \code{qp_negative_2ndderivative}. Enforced through the
#'     second-derivative design matrix, again with optional per-variable
#'     targeting handled by \code{\link{process_qp}} and \code{.build_deriv_qp}.
#'   \item Response range: \code{qp_range_lower}, \code{qp_range_upper}.
#'     Constrains \eqn{\mathbf{x}^{\top}\boldsymbol{\beta}} to lie within
#'     bounds. For non-identity links, the bounds are transformed to the
#'     link scale automatically inside \code{\link{process_qp}}, and then
#'     rescaled to match the internally standardized response.
#'   \item Custom constraints via \code{qp_Amat_fxn}, \code{qp_bvec_fxn},
#'     and \code{qp_meq_fxn}, which receive the design matrix structure and
#'     return the constraint matrix, bound vector, and number of equalities.
#'     These functions are commonly paired with a custom \code{qp_score_function}
#'     when the quadratic approximation is built from a non-default likelihood
#'     or score. The pre-built low-level objects \code{qp_Amat},
#'     \code{qp_bvec}, and \code{qp_meq} remain documented in
#'     \code{\link{lgspline}} and \code{\link{lgspline.fit}}, but in the current
#'     implementation they are not merged into the constraint set assembled by
#'     \code{\link{process_qp}}.
#' }
#' }
#'
#' @section Blockfit Backfitting for Linear Non-Interactive Effects:
#'
#' One of the practical advantages of keeping the polynomial basis explicit is
#' that special coefficient structures can be exploited directly. The blockfit
#' path is the clearest example: rather than carrying redundant copies of flat
#' linear coefficients and forcing them to agree through generic equality
#' constraints, the solver can pool them structurally from the start.
#'
#' \subsection{Motivation}{
#' When a model contains both spline terms (which receive \eqn{K+1}
#' partition-specific coefficient vectors constrained to smoothness at knots)
#' and non-interactive linear terms (``flat'' terms, specified via
#' \code{just_linear_without_interactions}, which receive a single shared
#' coefficient vector \eqn{\mathbf{v}} across all partitions), the standard
#' estimation procedure treats \eqn{\mathbf{v}} as \eqn{K+1} copies linked by
#' equality constraints. This inflates the effective parameter count and forces
#' the generic solver to carry redundant copies of coefficients that are
#' conceptually shared.
#'
#' Backfitting avoids this inflation by alternating between a spline step
#' and a flat step, each solving a lower-dimensional problem. Write the
#' partition-\eqn{k} design as
#' \eqn{\mathbf{X}_k = [\mathbf{Z}_k \mid \mathbf{X}_{\mathrm{flat}}^{(k)}]},
#' where \eqn{\mathbf{Z}_k} contains the spline columns and
#' \eqn{\mathbf{X}_{\mathrm{flat}}^{(k)}} contains the flat columns. Let
#' \eqn{n_c^{(s)}} and \eqn{n_c^{(f)}} denote the number of spline and flat
#' columns, respectively. This is invoked when \code{blockfit = TRUE}, flat
#' columns are non-empty, and \eqn{K > 0}.
#' }
#'
#' \subsection{Block-Coordinate Descent}{
#' The design matrix columns are split into a \emph{spline block} (receiving
#' full spline treatment with partition-specific coefficients) and a
#' \emph{flat block} (columns whose coefficients are pooled identically across
#' all \eqn{K+1} partitions). Penalty matrices and the constraint matrix
#' \eqn{\mathbf{A}} are partitioned accordingly.
#'
#' \strong{Spline step.} Holding \eqn{\mathbf{v}} fixed, the code forms the
#' response adjusted for the current flat contribution and applies the same
#' spline-only Lagrangian projection used in the main solver:
#' \deqn{\tilde{\boldsymbol{\beta}}_{\mathrm{spline}}^{(k)} =
#'   \mathbf{U}_{\mathrm{spline}}\mathbf{G}_{\mathrm{spline}}
#'   \mathbf{Z}_k^{\top}\mathbf{W}_k\mathbf{D}_k
#'   (\mathbf{y}_k - \mathbf{X}_{\mathrm{flat}}^{(k)}\mathbf{v}),}
#' for \eqn{k = 0, \ldots, K}. The spline-only constraint matrix is obtained
#' by extracting the spline rows of \eqn{\mathbf{A}}, dropping null columns,
#' and rank-reducing by QR.
#'
#' \strong{Flat step.} Holding the spline coefficients fixed, the shared flat
#' coefficients are updated by pooled penalized regression:
#' \deqn{\mathbf{v} = \left(\sum_{k=0}^{K}
#'   \mathbf{X}_{\mathrm{flat}}^{(k)\top}\mathbf{W}_k\mathbf{D}_k\mathbf{X}_{\mathrm{flat}}^{(k)}
#'   + \boldsymbol{\Lambda}_{\mathrm{flat}}\right)^{-1}
#'   \sum_{k=0}^{K}\mathbf{X}_{\mathrm{flat}}^{(k)\top}\mathbf{W}_k\mathbf{D}_k
#'   (\mathbf{y}_k - \mathbf{Z}_k\tilde{\boldsymbol{\beta}}_{\mathrm{spline}}^{(k)}).}
#'
#' When the full constraint matrix \eqn{\mathbf{A}} has columns with nonzero
#' entries on both spline and flat rows, the flat update is instead obtained by
#' solving a KKT system enforcing the residual equality constraint
#' \eqn{\mathbf{A}_{\mathrm{flat}}^{\top}\mathbf{v} =
#'   \mathbf{c} - \mathbf{A}_{\mathrm{spline}}^{\top}\tilde{\boldsymbol{\beta}}_{\mathrm{spline}}},
#' rather than using the unconstrained pooled solve shown above.
#'
#' \strong{Convergence.} In the Gaussian no-correlation case, convergence is
#' checked using the maximum absolute change across spline and flat coefficients.
#' In the weighted inner loop used by the GLM blockfit solvers, the code stops
#' when the flat-block update is smaller than \code{tol}; the spline block is
#' re-solved at each step from the current flat coefficients.
#' }
#'
#' \subsection{Four Estimation Cases}{
#' The backfitting solver dispatches to one of four paths based on the model
#' configuration.
#'
#' \strong{Case (a): Gaussian identity + GEE.}
#' Whitening destroys the block-diagonal structure needed for backfitting, so
#' the code skips the block iterations and performs the same full-system
#' Gaussian GEE projection used by \code{\link{get_B}}.
#'
#' \strong{Case (b): Gaussian identity, no correlation.}
#' Standard block-coordinate descent as described above. The spline-only
#' \eqn{\mathbf{G}_{\mathrm{spline}}} factors are precomputed once, and the pooled
#' flat penalized inverse is reused across iterations.
#'
#' \strong{Case (c): GLM + GEE (two-stage).}
#' This is the \code{.bf_case_glm_gee} path inside
#' \code{\link{blockfit_solve}}. Stage 1 forms a warm start by repeatedly
#' computing working responses and weights on the original scale and running the
#' weighted backfitting inner loop. Stage 2 then refines that warm start on the
#' full whitened system with the same dense SQP loop used by the main
#' non-Gaussian GEE solver.
#'
#' \strong{Case (d): GLM without GEE.}
#' This is the \code{.bf_case_glm_no_corr} path inside
#' \code{\link{blockfit_solve}}. An outer damped Newton--Raphson iterate updates
#' working responses and weights, while the inner loop alternates between
#' weighted spline and weighted flat updates. Deviance is monitored across outer
#' iterations, but there is no separate line-search damping coefficient applied
#' to the blockfit inner update itself.
#' }
#'
#' \subsection{Constraint Preservation and Coefficient Reassembly}{
#' Because the flat coefficients are shared by construction, the associated
#' flat-equality constraints are satisfied exactly, not approximately. The
#' smoothness constraints on the spline block are likewise enforced by the
#' spline-only Lagrangian projection. After convergence, the shared flat vector
#' \eqn{\mathbf{v}} is copied into each partition's coefficient vector,
#' yielding
#' \eqn{\boldsymbol{\beta}_k = [\tilde{\boldsymbol{\beta}}_{\mathrm{spline}}^{(k)\top}, \mathbf{v}^{\top}]^{\top}}
#' for downstream inference.
#'
#' When inequality constraints are present, they are enforced after
#' backfitting convergence through the same partition-wise active-set or dense
#' SQP refinement used in the main solver. For GEE blockfit (Case (c)), that
#' inequality handling occurs inside Stage 2 on the whitened system.
#'
#' If \code{\link{blockfit_solve}} throws an error, a warning is issued and the
#' code falls back to \code{\link{get_B}}.
#' }
#'
#' @section Knot Selection and Partitioning:
#'
#' The topic of knot selection is not the main focus of the package, but the
#' partition structure is central because every later design matrix, penalty,
#' and smoothness constraint depends on it. The defaults in \pkg{lgspline} are
#' therefore meant to be practical and transparent rather than theoretically
#' final.
#'
#' \subsection{Univariate Case}{
#' For a single predictor, the default partitioning is now handled by
#' \code{\link{make_partitions}} in the same \eqn{k}-means framework used more
#' generally: \eqn{K+1} centers are fit on an internally standardized copy
#' of the predictor, controlled by \code{standardize_predictors_for_knots}, and
#' then returned on the raw scale. Custom knots can still be supplied via
#' \code{custom_knots}, in which case partition assignment is built directly
#' from those raw-scale breakpoints. The default number of knots \eqn{K} is
#' chosen adaptively based on \eqn{N}, \eqn{p}, \eqn{q}, and the GLM family.
#' For multivariate fits, the resulting partition metadata are returned as
#' \code{make_partition_list} and can be re-used in later calls to
#' \code{\link{lgspline}}. This is particularly useful when one wants to hold
#' the partition geometry fixed across repeated fits, for example while varying
#' penalties, families, or correlation structures.
#' }
#' \subsection{Multivariate Case}{
#' For multiple predictors, \eqn{K+1} cluster centers are identified by
#' \eqn{k}-means on an internally standardized predictor matrix via
#' \code{\link{make_partitions}}. This is the partitioning mechanism used to
#' determine the multivariate spline regions; see MacQueen (1967) for the
#' classical clustering formulation and Kisi et al. (2025) for a recent applied
#' example of \eqn{k}-means-driven partitioning in a nonlinear prediction
#' setting. Midpoints between neighboring centers
#' (those whose midpoint  does not fall into a third cluster) serve as knot
#' locations. Observations are assigned to the nearest cluster center using
#' \code{\link[FNN]{get.knnx}}, and the returned centers and knots are on
#' the original predictor scale. The resulting partition structure is a type
#' of Voronoi diagram and is stored in the fitted object as
#' \code{make_partition_list}. The \code{do_not_cluster_on_these} argument can
#' exclude certain predictors from clustering (e.g., a treatment indicator that
#' should not drive partitioning). The lower-level clustering behavior can be
#' further controlled by \code{cluster_args} and \code{neighbor_tolerance},
#' while \code{cluster_on_indicators} determines whether binary predictors are
#' allowed to influence the partition geometry at all.
#' }
#'
#' \subsection{Standardizing Predictors}{
#' Higher-order polynomial terms can dramatically inflate or deflate the
#' magnitude of basis expansions, introducing numerical instability. All
#' polynomial basis expansions are scaled by
#' \eqn{q_{0.69} - q_{0.31}}, where \eqn{q_{\zeta}} is the \eqn{\zeta}-th
#' quantile of the expansion. For a standard normal distribution this quantity
#' is approximately 1, so the scaling is close to one standard deviation for
#' symmetric distributions. This fitting-stage rescaling is controlled by
#' \code{standardize_expansions_for_fitting}, while knot construction is
#' controlled separately by \code{standardize_predictors_for_knots}. The same
#' scaling is applied to the constraint matrix to maintain smoothness, and
#' coefficients are back-transformed to the original scale after fitting.
#' }
#'
#' @section Smoothing Spline Penalty:
#'
#' \subsection{Penalty Construction}{
#' The penalty matrix \eqn{\boldsymbol{\Lambda}_s} penalizes the integrated
#' squared total curvature of the fitted function over the observed predictor
#' ranges. This is the step that makes the piecewise polynomial fit genuinely
#' behave like a smoothing spline rather than merely a constrained regression
#' spline. The package computes this penalty directly from the monomial
#' structure of the basis rather than by appealing to a pre-tabulated spline
#' basis. For a single partition \eqn{k} with basis expansion
#' \eqn{\mathbf{x}_k = (\phi_1(\mathbf{t}), \ldots, \phi_p(\mathbf{t}))^{\top}}
#' where each \eqn{\phi_i(\mathbf{t}) = \prod_{j=1}^{q} t_j^{\alpha_{ij}}}
#' is a multivariate monomial:
#' \deqn{\boldsymbol{\beta}_k^{\top}\boldsymbol{\Lambda}_s\boldsymbol{\beta}_k
#'   = \int_{\mathbf{a}}^{\mathbf{b}}
#'     \|\tilde{f}_k''(\mathbf{t})\|^{2}\,d\mathbf{t},}
#' where \eqn{\mathbf{a}} and \eqn{\mathbf{b}} are the observed predictor
#' minimums and maximums (computed globally from the data, not
#' partition-specific), and
#' \eqn{\tilde{f}_k(\mathbf{t}) = \mathbf{x}_k^{\top}\boldsymbol{\beta}_k}
#' is the fitted function for partition \eqn{\mathcal{P}_k}.
#'
#' \strong{Total curvature operator.}
#' The integrated squared second derivative decomposes into \eqn{q}
#' curvature operators, one per predictor. For predictor \eqn{v}, the
#' curvature operator \eqn{D_v} is defined as
#' \deqn{D_v = \frac{\partial^{2}}{\partial t_v^{2}}
#'   + \sum_{s \neq v}\frac{\partial^{2}}{\partial t_v\,\partial t_s}.}
#' That is, \eqn{D_v} captures both the pure second derivative with respect
#' to \eqn{t_v} and all mixed second partial derivatives involving \eqn{t_v}.
#' The penalty matrix entries are then
#' \deqn{[\boldsymbol{\Lambda}_s]_{ij}
#'   = \sum_{v=1}^{q}\int_{\mathbf{a}}^{\mathbf{b}}
#'     D_v(\phi_i)\,D_v(\phi_j)\,d\mathbf{t}.}
#'
#' \strong{Monomial derivative rule.}
#' For a monomial \eqn{\phi(\mathbf{t}) = \prod_j t_j^{\alpha_j}}, the
#' derivatives entering \eqn{D_v} have closed forms. The pure second
#' derivative is
#' \deqn{\frac{\partial^{2}}{\partial t_v^{2}}\prod_j t_j^{\alpha_j}
#'   = \alpha_v(\alpha_v - 1)\,t_v^{\alpha_v - 2}\prod_{j \neq v} t_j^{\alpha_j},}
#' which is zero when \eqn{\alpha_v < 2}. The mixed second derivative is
#' \deqn{\frac{\partial^{2}}{\partial t_v\,\partial t_s}\prod_j t_j^{\alpha_j}
#'   = \alpha_v\alpha_s\,t_v^{\alpha_v - 1}t_s^{\alpha_s - 1}
#'   \prod_{j \neq v,s} t_j^{\alpha_j},}
#' which is zero when \eqn{\alpha_v < 1} or \eqn{\alpha_s < 1}. Applying
#' \eqn{D_v} to a monomial \eqn{\phi_i} produces a sum of monomials with
#' known coefficients and exponent vectors.
#'
#' \strong{Factorized integration.}
#' Because every \eqn{D_v(\phi_i)} is polynomial, the product
#' \eqn{D_v(\phi_i)\,D_v(\phi_j)} is also polynomial and the multivariate
#' integral factorizes over predictors:
#' \deqn{\int_{\mathbf{a}}^{\mathbf{b}}\prod_{j=1}^{q} t_j^{e_j}\,d\mathbf{t}
#'   = \prod_{j=1}^{q}\frac{b_j^{e_j+1} - a_j^{e_j+1}}{e_j + 1}.}
#' Crucially, this integral runs over \emph{all} \eqn{q} predictor ranges,
#' including predictors that do not appear in the integrand (for which
#' \eqn{e_j = 0} and the factor reduces to \eqn{b_j - a_j}). This ensures
#' that the penalty is properly scaled relative to the volume of the
#' predictor space.
#'
#' \strong{Single-predictor verification.}
#' For \eqn{q = 1} with expansion
#' \eqn{\mathbf{x} = (1, t, t^{2}, t^{3})^{\top}} on \eqn{[a, b]}, the
#' curvature operator reduces to \eqn{D_1 = \partial^{2}/\partial t^{2}} (no
#' mixed partials exist), and the penalty matrix reduces to
#' \deqn{\boldsymbol{\Lambda}_s
#'   = \int_a^b \mathbf{x}''\mathbf{x}''^{\top}\,dt
#'   = \begin{pmatrix}
#'       0 & 0 & 0 & 0 \\
#'       0 & 0 & 0 & 0 \\
#'       0 & 0 & 4(b - a) & 6(b^{2} - a^{2}) \\
#'       0 & 0 & 6(b^{2} - a^{2}) & 12(b^{3} - a^{3})
#'     \end{pmatrix},}
#' recovering the classical cubic smoothing spline penalty.
#'
#' \strong{Handling of non-spline predictors.}
#' Predictors specified via \code{just_linear_without_interactions} or
#' \code{just_linear_with_interactions} do not receive higher-order
#' polynomial expansions in the design matrix. To ensure their curvature
#' contributions are still correctly computed (particularly through
#' interaction terms), the implementation temporarily appends phantom
#' higher-order columns (with zero data) for these predictors, computes
#' the full curvature penalty on the augmented basis, and then subsets the
#' result back to the original \eqn{p \times p} dimensions. This ensures
#' that interaction terms involving non-spline predictors receive appropriate
#' penalty contributions without affecting the rest of the estimation
#' pipeline.
#'
#' \strong{Parallel computation.}
#' Because the total penalty is an additive sum over predictors
#' (\eqn{\boldsymbol{\Lambda}_s = \sum_{v=1}^{q}\boldsymbol{\Lambda}_{s,v}}),
#' the computation can be parallelized by distributing the per-predictor
#' curvature matrices across workers via \code{parallel::parLapply} and
#' summing the results. This is controlled by the \code{parallel_penalty}
#' argument and is beneficial when \eqn{q} is large.
#'
#' The penalty is computed by \code{\link{get_2ndDerivPenalty}} (single
#' predictor or subset) and \code{\link{get_2ndDerivPenalty_wrapper}}
#' (full assembly with optional parallelism and non-spline handling).
#'
#' Because the smoothing penalty has zero eigenvalues for the intercept and
#' linear terms (whose second derivatives vanish), an optional ridge penalty on
#' lower-order terms is added for computational stability. The full penalty
#' block for partition \eqn{k} is:
#' \deqn{\boldsymbol{\Lambda}_k = \lambda_w\bigl(\boldsymbol{\Lambda}_s + \lambda_r\boldsymbol{\Lambda}_r + \sum_{m=1}^{M}\xi_{mk}\mathbf{P}_{mk}\bigr)}
#' where \eqn{\lambda_w} is the global wiggle penalty (\code{wiggle_penalty}),
#' \eqn{\lambda_r} is ridge penalty on linear and intercept terms
#' (\code{flat_ridge_penalty}) multiplied by the wiggle penalty, and
#' \eqn{\xi_{mk}} and \eqn{\mathbf{P}_{mk}} denote optional additional
#' penalty multipliers and matrices, including the predictor- and
#' partition-specific components activated through
#' \code{unique_penalty_per_predictor}, \code{unique_penalty_per_partition},
#' \code{predictor_penalties}, and \code{partition_penalties}. This assembly is handled by
#' \code{\link{compute_Lambda}}.
#'
#' The penalty matrix \eqn{\boldsymbol{\Lambda}} is stored as a list
#' of \eqn{K+1} \eqn{p \times p} square, symmetric, positive semi-definite
#' matrices.
#' }
#'
#' \subsection{Penalty Optimization via Generalized Cross-Validation}{
#' After the structural pieces of the model are fixed, the main remaining
#' question is how much smoothing to apply. In \pkg{lgspline}, that tuning is
#' performed with generalized cross-validation, but it is carried out using the
#' same constrained estimator that will be used in the final model fit.
#' Penalty parameters are estimated on the log scale via exponential
#' parameterization (\eqn{\lambda = \exp(\theta)},
#' \eqn{\theta \in \mathbb{R}}), ensuring positivity. The chain rule factor
#' \eqn{\partial\exp(\theta)/\partial\theta = \exp(\theta) = \lambda} is
#' applied throughout. User-facing arguments (\code{initial_wiggle},
#' \code{initial_flat}, \code{predictor_penalties},
#' \code{partition_penalties}) accept values on the raw, natural scale;
#' conversion to log scale is handled internally. The final tuned values and
#' assembled penalty pieces are returned in the fitted object's
#' \code{penalties} component.
#'
#' The total penalty matrix \eqn{\boldsymbol{\Lambda}} is constructed as:
#' \deqn{\boldsymbol{\Lambda} = \lambda_w\mathbf{P}_{w} + \lambda_r\mathbf{P}_{r}
#'   + \sum_{j}\nu_j\mathbf{P}_j^{(\mathrm{pred})}
#'   + \sum_{k}\tau_k\mathbf{P}_k^{(\mathrm{part})},}
#' where \eqn{\mathbf{P}_{w}} is the integrated squared second-derivative
#' penalty (i.e., \eqn{\boldsymbol{\Lambda}_s} above), \eqn{\mathbf{P}_{r}} is
#' a ridge penalty on intercept and linear coefficients,
#' \eqn{\mathbf{P}_j^{(\mathrm{pred})}} are predictor-specific penalties, and
#' \eqn{\mathbf{P}_k^{(\mathrm{part})}} are partition-specific penalties. The
#' scalars \eqn{\lambda_w} (\code{wiggle_penalty}), \eqn{\lambda_r}
#' (\code{flat_ridge_penalty}), \eqn{\{\nu_j\}}
#' (\code{predictor_penalties}), and \eqn{\{\tau_k\}}
#' (\code{partition_penalties}) are tuned.
#' The unbiased generalized cross-validation criterion is
#' \deqn{\mathrm{GCV}_{u} = \frac{\sum_{i=1}^{N}D_{ii}\,r_i^{2}}{N(1 - \bar{W})^{2}},}
#' where \eqn{r_i} are residuals on the link scale and
#' \eqn{\bar{W} = \mathrm{tr}(\mathbf{H})/N} is the mean of the hat-matrix
#' diagonal. For identity link, \eqn{r_i = y_i - \hat{\eta}_i}. For
#' non-identity links, the residuals are
#' \eqn{r_i = g((y_i + \delta)/(1+2\delta)) - (\hat{\eta}_i + \delta)/(1+2\delta)},
#' where \eqn{\delta \geq 0} is a pseudocount that stabilizes the link
#' transformation, automatically tuned within \code{\link{tune_Lambda}} if not supplied
#' to \code{delta}. Non-Gaussian families with observation weights
#' \eqn{\omega_i} have their residuals scaled by \eqn{\omega_i}. When the
#' family provides a custom deviance residual function, that function is used
#' in place of the link-scale residuals.
#'
#' \strong{Pseudocount selection.} The pseudocount \eqn{\delta} is chosen to
#' make the transformed response distribution most closely approximate a
#' \eqn{t}-distribution with \eqn{N-1} degrees of freedom, in the sense of
#' minimizing the (optionally weighted) mean absolute deviation between the
#' sorted standardized transformed responses and the corresponding
#' \eqn{t}-quantiles. This is solved via Brent's method over
#' \eqn{[10^{-64}, 1]}. When the link is identity, or when the response
#' naturally lies in the domain of the link function, \eqn{\delta = 0}.
#' This behavior is exposed through the \code{delta} argument in
#' \code{\link{lgspline}}: supplying a fixed numeric value bypasses the internal
#' search, while leaving it \code{NULL} allows the tuning code to choose the
#' stabilizing pseudocount automatically when needed.
#'
#' \strong{Meta-penalty regularization.} A regularization term pulls the
#' predictor- and partition-specific penalty parameters toward 1 on the raw
#' scale:
#' \deqn{P_{\mathrm{meta}}(\lambda_w, \{\nu_j\}, \{\tau_k\})
#'   = \frac{1}{2}c_{\mathrm{meta}}\sum_j(\nu_j - 1)^{2}
#'   + \frac{1}{2}\cdot 10^{-32}(\lambda_w - 1)^{2},}
#' where \eqn{c_{\mathrm{meta}}} is a user-specified coefficient
#' (\code{meta_penalty}). The gradient of \eqn{P_{\mathrm{meta}}} on the log
#' scale, incorporating the exp chain rule, is
#' \eqn{\partial P_{\mathrm{meta}}/\partial\theta_j = c_{\mathrm{meta}}(\nu_j - 1)\nu_j}
#' and
#' \eqn{\partial P_{\mathrm{meta}}/\partial\theta_1 = 10^{-32}(\lambda_w - 1)\lambda_w}.
#' The total objective is \eqn{\mathrm{GCV}_{u} + P_{\mathrm{meta}}}.
#' }
#'
#' \subsection{Closed-Form Gradient of GCV}{
#' The gradient of \eqn{\mathrm{GCV}_{u}} with respect to
#' \eqn{\theta_1 = \log\lambda_w} is computed analytically via the quotient
#' rule:
#' \deqn{\frac{\partial\mathrm{GCV}_{u}}{\partial\theta_1}
#'   = \frac{1}{D^{2}}\left(\frac{\partial\mathcal{N}}{\partial\theta_1}D
#'   - \mathcal{N}\frac{\partial D}{\partial\theta_1}\right),}
#' where \eqn{\mathcal{N} = \sum r_i^{2}} (numerator) and
#' \eqn{D = N(1 - \bar{W})^{2}} (denominator). The key intermediates are:
#' \itemize{
#'   \item \eqn{\partial\mathbf{G}/\partial\lambda_w}, computed from the
#'     matrix identity
#'     \eqn{\partial(\mathbf{X}^{\top}\mathbf{X} + \boldsymbol{\Lambda})^{-1}/\partial\lambda = -\mathbf{G}(\partial\boldsymbol{\Lambda}/\partial\lambda)\mathbf{G}}.
#'   \item \eqn{\partial\mathbf{G}^{1/2}/\partial\lambda_w}, derived from
#'     \eqn{\partial\mathbf{G}/\partial\lambda_w} via the eigendecomposition
#'     chain rule.
#'   \item \eqn{\partial\bar{W}/\partial\lambda_w}, the derivative of the
#'     trace of the hat matrix
#'     \eqn{\mathbf{H} = \mathbf{X}\mathbf{U}\mathbf{G}\mathbf{X}^{\top}},
#'     which depends on both
#'     \eqn{\partial\mathbf{G}/\partial\lambda_w} and
#'     \eqn{\partial\mathbf{G}^{1/2}/\partial\lambda_w}.
#'   \item \eqn{\partial\mathcal{N}/\partial\theta_1 = -2\mathbf{r}^{\top}\mathbf{X}(\partial(\mathbf{U}\mathbf{G})/\partial\lambda_w)\mathbf{X}^{\top}\mathbf{y}\cdot\lambda_w},
#'     via the chain rule applied to the residual vector.
#'   \item \eqn{\partial D/\partial\theta_1 = 2(1 - \bar{W})(-\partial\bar{W}/\partial\lambda_w)\cdot\lambda_w}.
#' }
#' In the implementation, these quantities are assembled by a small set of
#' helper routines: \code{\link{compute_dG_dlambda}} for
#' \eqn{\partial\mathbf{G}/\partial\lambda}, \code{\link{compute_dGhalf}}
#' for \eqn{\partial\mathbf{G}^{1/2}/\partial\lambda},
#' \code{\link{compute_dW_dlambda_wrapper}} for derivatives of the effective
#' degrees-of-freedom term, \code{\link{compute_trace_UGXX_wrapper}} for the
#' trace pieces entering GCV, and \code{\link{compute_dG_u_dlambda_xy}} for
#' the derivative of the fitted-value quadratic form.
#'
#' The full gradient is scaled by \eqn{N} before adding the meta-penalty
#' gradient.
#'
#' For the ridge penalty and predictor-/partition-specific penalties, a
#' trace-ratio heuristic is used:
#' \deqn{\frac{\partial\mathrm{GCV}_{u}}{\partial\lambda_l} \approx
#'   \frac{\mathrm{mean}(\mathrm{diag}(\mathbf{P}_l))}{\mathrm{mean}(\mathrm{diag}(\boldsymbol{\Lambda}))}
#'   \frac{\partial\mathrm{GCV}_{u}}{\partial\lambda_w},}
#' and analogously for predictor- and partition-specific penalties, where
#' \eqn{\mathbf{P}_j^{(\mathrm{pred})}} or \eqn{\mathbf{P}_k^{(\mathrm{part})}}
#' replaces \eqn{\mathbf{P}_l} in the numerator. This follows from a
#' chain-rule argument: by the Leibniz rule and the inverse derivative,
#' \eqn{\partial\lambda_w/\partial\lambda_l = (\partial\boldsymbol{\Lambda}/\partial\lambda_l)(\partial\boldsymbol{\Lambda}/\partial\lambda_w)^{-1}}.
#' Since the derivative appears as a matrix rather than a scalar, the
#' mean-diagonal ratio provides a scalar summary. Once the derivative for
#' \eqn{\lambda_w} is in hand, the derivatives of other penalties are cheap
#' to compute. The exp chain rule is then applied:
#' \eqn{\partial/\partial\theta = (\partial/\partial\lambda)\cdot\lambda}.
#' }
#'
#' \subsection{Optimization Procedure}{
#' \strong{Grid search initialization.} The \eqn{\mathrm{GCV}_{u}} criterion
#' is evaluated over a grid of candidate values for
#' \eqn{(\lambda_w, \lambda_r)} on the log scale. All combinations of
#' user-supplied candidate vectors (\code{initial_wiggle} and
#' \code{initial_flat}) are formed, and the combination yielding
#' the smallest finite \eqn{\mathrm{GCV}_{u}} is selected as the starting
#' point for BFGS optimization. Grid points producing non-finite
#' \eqn{\mathrm{GCV}_{u}} are discarded. If all grid points fail, an error
#' is raised advising the user to check the data or adjust the grid.
#'
#' \strong{Damped BFGS optimizer.} A custom damped BFGS quasi-Newton optimizer,
#' implemented in \code{\link{efficient_bfgs}}, minimizes
#' \eqn{\mathrm{GCV}_{u} + P_{\mathrm{meta}}}. When analytic gradients are not
#' usable, the fallback finite-difference helper is \code{\link{approx_grad}}.
#'
#' \emph{Iterations 1--2: steepest descent.} The first two iterations use
#' steepest descent with a damping factor \eqn{\alpha}:
#' \eqn{\boldsymbol{\theta}^{(t+1)} = \boldsymbol{\theta}^{(t)} - \alpha\nabla_{\boldsymbol{\theta}}}.
#'
#' \emph{Iterations 3+: BFGS.} From iteration 3, an inverse Hessian
#' approximation \eqn{\mathbf{Q}^{(t)}} is maintained via the standard secant
#' update. Let
#' \eqn{\mathbf{s}^{(t)} = \boldsymbol{\theta}^{(t)} - \boldsymbol{\theta}^{(t-1)}}
#' and
#' \eqn{\mathbf{v}^{(t)} = \nabla^{(t)} - \nabla^{(t-1)}}. The BFGS update
#' is:
#' \deqn{\mathbf{Q}^{(t+1)}
#'   = (\mathbf{I} - \rho\mathbf{s}\mathbf{v}^{\top})
#'   \mathbf{Q}^{(t)}
#'   (\mathbf{I} - \rho\mathbf{v}\mathbf{s}^{\top})
#'   + \rho\mathbf{s}\mathbf{s}^{\top},
#'   \qquad \rho = (\mathbf{v}^{\top}\mathbf{s})^{-1}.}
#' When \eqn{|\mathbf{v}^{\top}\mathbf{s}| < 10^{-64}}, the approximation is
#' reset to \eqn{\mathbf{I}} and the iteration is flagged for restart. The
#' search direction is
#' \eqn{\mathbf{d}^{(t)} = -\mathbf{Q}^{(t)}\nabla^{(t)}}.
#'
#' \emph{Step acceptance.} A step is accepted if
#' \eqn{\mathrm{GCV}_{u}^{(\mathrm{new})} \leq \mathrm{GCV}_{u}^{(\mathrm{old})}}.
#' On rejection, \eqn{\alpha} is halved. If \eqn{\alpha < 2^{-10}} (early
#' iterations) or \eqn{\alpha < 2^{-12}} (later iterations), the optimizer
#' terminates with the best solution found.
#'
#' \emph{Convergence.} The optimizer terminates when
#' \eqn{|\mathrm{GCV}_{u}^{(t)} - \mathrm{GCV}_{u}^{(t-1)}| < \epsilon} or
#' \eqn{\|\boldsymbol{\theta}^{(t)} - \boldsymbol{\theta}^{(t-1)}\|_{\infty} < \epsilon},
#' provided at least 10 iterations have elapsed.
#'
#' \emph{Fallback.} If the custom BFGS fails, a base-R
#' \code{\link[stats]{optim}} call with method \code{"BFGS"} and
#' finite-difference gradients is used instead. If both fail, the best
#' grid-search point is used.
#'
#' \strong{Post-optimization inflation.} After optimization, the penalty
#' parameters are inflated by a factor \eqn{((N+2)/(N-2))^{2}} to counteract
#' the in-sample bias toward underpenalization inherent in GCV-type criteria.
#'
#' The tuning loop is implemented in \code{\link{tune_Lambda}}.
#' }
#'
#' @section Incorporating Non-Spline Effects:
#'
#' Multiple fixed effects are accommodated naturally in the LMSS framework
#' because spline effects, linear effects, and many interaction terms all live
#' in the same partition-wise polynomial expansion. The distinction is therefore
#' not whether a term is "allowed" by the solver, but whether it receives full
#' spline treatment or remains structurally linear across partitions.
#'
#' The constrained framework naturally accommodates non-spline terms. If only
#' linear terms are included for a predictor (via
#' \code{just_linear_without_interactions} or
#' \code{just_linear_with_interactions}), the first-derivative smoothing
#' constraint forces the linear coefficient to be identical across all
#' partitions, since the derivative of a linear function is its slope. This
#' is not an algorithmic modification but a natural consequence of the
#' constraint structure.
#'
#' For example, a model with one spline effect and a linear treatment indicator
#' interaction will naturally keep the treatment--time interaction coefficient
#' constant across partitions while allowing the time effect to vary
#' nonlinearly. This conveniently extends to arbitrary combinations of spline
#' and linear terms without requiring special handling.
#'
#' When \code{blockfit = TRUE} is specified alongside
#' \code{just_linear_without_interactions}, the flat-block path provides an
#' alternative enforcement mechanism. Rather than relying on constraint
#' projection, flat coefficients are pooled structurally across partitions
#' during backfitting. The two approaches agree at the point estimate but
#' differ in their uncertainty quantification; see the Blockfit section above.
#'
#' @section Integration:
#'
#' Because the fitted object retains an explicit polynomial representation in
#' each partition, numerical integration can be carried out in a fairly direct
#' way. The package wraps that calculation in a user-facing S3 method so the
#' user does not need to manage knot boundaries or partition membership by hand.
#'
#' In the user-facing interface, numerical integration is applied through
#' \code{\link{integrate.lgspline}}, which applies Gauss--Legendre quadrature to
#' predictions from the fitted model produced by \code{\link{predict.lgspline}}.
#'
#' \subsection{Implementation}{
#' For a user-supplied rectangular domain, \code{\link{integrate.lgspline}} constructs a
#' tensor-product grid of Gauss--Legendre nodes, evaluates the fitted model at
#' those points, and forms the weighted sum. This works for both univariate and
#' multivariate models, respects the fitted partition structure automatically,
#' and avoids requiring the user to keep track of knot boundaries by hand.
#'
#' The \code{vars} argument selects which predictors are integrated over.
#' Predictors not listed in \code{vars} are held fixed at
#' \code{initial_values} when supplied, or otherwise at the midpoint of their
#' observed training range. The optional \code{B_predict} argument makes it
#' possible to integrate posterior draws or other alternate coefficient sets,
#' and \code{n_quad} controls the number of Gauss--Legendre nodes used per
#' integrated dimension.
#'
#' Integration is performed on the response scale by default. Setting
#' \code{link_scale = TRUE} instead integrates the linear predictor
#' \eqn{\eta = f(\mathbf{t})}. For identity-link Gaussian models the two scales
#' coincide.
#' }
#'
#'
#' @section Lagrange Multipliers:
#'
#' When \code{return_lagrange_multipliers = TRUE}, the multiplier vector
#' \deqn{\boldsymbol{\lambda} = (\mathbf{A}^{\top}\mathbf{G}\mathbf{A})^{-1}\mathbf{A}^{\top}\hat{\boldsymbol{\beta}}}
#' is returned. These quantify the sensitivity of the penalized objective to
#' relaxing each smoothness or user-supplied equality constraint. When
#' constraint target values are nonzero
#' (\eqn{\mathbf{A}^{\top}\boldsymbol{\beta}_0 \neq \mathbf{0}}), the modified
#' formulation is used:
#' \deqn{\boldsymbol{\lambda} = (\mathbf{A}^{\top}\mathbf{G}\mathbf{A})^{-1}\mathbf{A}^{\top}(\hat{\boldsymbol{\beta}} - \boldsymbol{\beta}_0)}
#' where \eqn{\mathbf{A}^{\top}\boldsymbol{\beta}_0} is the vector of
#' constraint target values. Multipliers are \code{NULL} when no constraints
#' are active (\eqn{\mathbf{A}} is \code{NULL} or \eqn{K = 0}).
#'
#' For inequality constraints, multipliers are returned as computed by
#' \code{\link[quadprog]{solve.QP}}. The Lagrange multipliers for active
#' inequality constraints can be used diagnostically to identify which shape
#' constraints are most costly in terms of goodness of fit.
#'
#' @section S3 Methods:
#'
#' Standard S3 methods are provided for objects of class \code{lgspline}:
#' \itemize{
#'   \item \code{\link{print.lgspline}} and \code{\link{summary.lgspline}}:
#'     Provide concise model summaries, with
#'     \code{\link{print.summary.lgspline}} formatting coefficient tables in a
#'     familiar regression-style layout.
#'
#'   \item \code{\link{logLik.lgspline}}: Returns a standard \code{logLik}
#'     object. For Gaussian responses with identity link, the exact
#'     log-likelihood is computed. When a correlation structure is present via
#'     \code{VhalfInv}, the log-likelihood includes the
#'     \eqn{\log|\mathbf{V}^{-1/2}|} adjustment and the corresponding
#'     whitened quadratic form. For other families, the method falls back to
#'     \code{family$aic} or a deviance-based approximation. An
#'     \code{include_prior} argument (default \code{TRUE}) optionally adds
#'     the Gaussian prior penalty interpretation of the smoothing spline
#'     penalty
#'     \eqn{-\frac{1}{2\sigma^{2}}\tilde{\boldsymbol{\beta}}^{\top}\boldsymbol{\Lambda}\tilde{\boldsymbol{\beta}}}
#'     to obtain a penalized MAP log-likelihood.
#'
#'   \item \code{\link{predict.lgspline}}: Produces fitted values and related
#'     quantities (e.g., derivatives and standard errors through
#'     \code{se.fit = TRUE}), lets \code{new_predictors} override
#'     \code{newdata}, accepts alternate coefficient lists through
#'     \code{B_predict}, and supports prediction on new predictor matrices
#'     consistent with the original spline expansions.
#'
#'   \item \code{\link{coef.lgspline}}: Extracts partition-specific
#'     coefficient vectors.
#'
#'   \item \code{\link{confint.lgspline}}: Extracts confidence intervals.
#'     When the inverse Hessian from BFGS optimization is available for
#'     correlation parameters, intervals for those correlation parameters are
#'     returned on the working (transformed) scale and should be
#'     back-transformed as described in the correlation section.
#'   \item \code{\link{plot.lgspline}}: For one-dimensional fits, produces
#'     base R graphics showing the fitted function (with optional partition-wise
#'     formulas) and supports overlay via \code{add = TRUE}. For two or more
#'     predictors, an interactive \pkg{plotly}-based visualization is
#'     returned. Specific predictors may be selected via \code{vars}.
#'
#'   \item \code{\link{integrate.lgspline}}: Computes definite integrals of
#'     the fitted surface over rectangular domains by Gauss--Legendre
#'     quadrature.
#' }
#'
#' Additional user-facing helpers include \code{\link{wald_univariate}} for
#' coefficient-wise Wald inference, \code{\link{generate_posterior}} for
#' posterior and posterior-predictive sampling,
#' \code{\link{generate_posterior_correlation}} for correlation-aware posterior
#' simulation, \code{\link{equation}} for closed-form display of the fitted
#' partition formulas, and
#' \code{\link{find_extremum}} for optimizing the fitted surface or a custom
#' acquisition function built from it.
#'
#' @references
#'
#' Buse, A. and Lim, L. (1977). Cubic Splines as a Special Case of
#' Restricted Least Squares. \emph{Journal of the American Statistical
#' Association}, 72, 64--68.
#'
#' Eilers, P. H. and Marx, B. D. (1996). Flexible Smoothing with B-splines
#' and Penalties. \emph{Statistical Science}, 11(2), 89--121.
#'
#' Ezhov, N., Neitzel, F. and Petrovic, S. (2018). Spline Approximation,
#' Part 1: Basic Methodology. \emph{Journal of Applied Geodesy}, 12(2),
#' 139--155.
#'
#' Goldfarb, D. and Idnani, A. (1983). A Numerically Stable Dual Method for
#' Solving Strictly Convex Quadratic Programs. \emph{Mathematical
#' Programming}, 27(1), 1--33.
#'
#' Harville, D. A. (1977). Maximum Likelihood Approaches to Variance
#' Component Estimation and to Related Problems. \emph{Journal of the
#' American Statistical Association}, 72(358), 320--338.
#'
#' Hastie, T. J. and Tibshirani, R. J. (1990). \emph{Generalized Additive
#' Models}. Chapman & Hall/CRC.
#'
#' Kisi, O., Heddam, S., Parmar, K. S., Petroselli, A., K\"ulls, C. and
#' Zounemat-Kermani, M. (2025). Integration of Gaussian Process Regression
#' and K Means Clustering for Enhanced Short Term Rainfall Runoff Modeling.
#' \emph{Scientific Reports}, 15, 7444.
#'
#' MacQueen, J. B. (1967). Some Methods for Classification and Analysis of
#' Multivariate Observations. In \emph{Proceedings of the Fifth Berkeley
#' Symposium on Mathematical Statistics and Probability}, Volume 1,
#' 281--297. University of California Press.
#'
#' McCullagh, P. and Nelder, J. A. (1989). \emph{Generalized Linear Models}.
#' Chapman & Hall, 2nd edition.
#'
#' Murray, I., Adams, R. P. and MacKay, D. J. C. (2010). Elliptical Slice
#' Sampling. \emph{Proceedings of the 13th International Conference on
#' Artificial Intelligence and Statistics (AISTATS)}, 9, 541--548.
#'
#' Nocedal, J. and Wright, S. J. (2006). \emph{Numerical Optimization}
#' (2nd ed.). Springer.
#'
#' Patterson, H. D. and Thompson, R. (1971). Recovery of Inter-Block
#' Information When Block Sizes Are Unequal. \emph{Biometrika}, 58,
#' 545--554.
#'
#' Pya, N. and Wood, S. N. (2015). Shape Constrained Additive Models.
#' \emph{Statistics and Computing}, 25(3), 543--559.
#'
#' Reinsch, C. H. (1967). Smoothing by Spline Functions. \emph{Numerische
#' Mathematik}, 10, 177--183.
#'
#' Ruppert, D., Wand, M. P. and Carroll, R. J. (2003).
#' \emph{Semiparametric Regression}. Cambridge University Press.
#'
#' Searle, S. R., Casella, G. and McCulloch, C. E. (2006). \emph{Variance
#' Components}. Wiley.
#'
#' Wahba, G. (1990). \emph{Spline Models for Observational Data}. SIAM.
#'
#' Wood, S. N. (2006). On Confidence Intervals for Generalized Additive
#' Models Based on Penalized Regression Splines. \emph{Australian & New
#' Zealand Journal of Statistics}, 48(4), 445--464.
#'
#' Wood, S. N. (2011). Fast Stable Restricted Maximum Likelihood and Marginal
#' Likelihood Estimation of Semiparametric Generalized Linear Models.
#' \emph{Journal of the Royal Statistical Society: Series B}, 73(1), 3--36.
#'
#' Wood, S. N. (2017). \emph{Generalized Additive Models: An Introduction
#' with R}. CRC Press, 2nd edition.
#'
#' @name Details
NULL











