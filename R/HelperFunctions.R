#' @useDynLib lgspline
#' @importFrom Rcpp sourceCpp
NULL

#' Create One-Hot Encoded Matrix
#'
#' @description
#' Converts a categorical vector into a one-hot encoded matrix where each unique value
#' becomes a binary column.
#'
#' @param x A vector containing categorical values (factors, character, etc.)
#'
#' @return A data frame containing the one-hot encoded binary columns with cleaned column names
#'
#' @details
#' The function creates dummy variables for each unique value in the input vector using
#' \code{model.matrix()} with dummy-intercept coding. Column names are cleaned by removing the
#' 'x' prefix added by \code{model.matrix()}.
#'
#' @examples
#'
#' ## lgspline will not accept this format of "catvar", because inputting data
#' # this way can cause difficult-to-diagnose issues in formula parsing
#' # all variables must be numeric
#' df <- data.frame(numvar = rnorm(100),
#'                  catvar = rep(LETTERS[1:4],
#'                               25))
#' print(head(df))
#'
#' ## Instead, replace with dummy-intercept coding by
#' # 1) applying one-hot encoding
#' # 2) dropping the first column
#' # 3) appending to our data
#'
#' dummy_intercept_coding <- create_onehot(df$catvar)[,-1]
#' df$catvar <- NULL
#' df <- cbind(df, dummy_intercept_coding)
#' print(head(df))
#'
#'
#' @export
create_onehot <- function(x, drop_first = FALSE) {
  ## [Change 2026-02-12] C5: Added drop_first parameter for convenience;
  # coerce to factor to ensure stable level ordering
  if(!is.factor(x)) x <- as.factor(x)
  mat <- model.matrix(~ x - 1)
  colnames(mat) <- sub("^x", "", colnames(mat))
  result <- as.data.frame(mat)
  if(drop_first && ncol(result) > 1){
    result <- result[, -1, drop = FALSE]
  }
  return(result)
}

#' Standardize Vector to Z-Scores
#'
#' @description
#' Centers a vector by its sample mean, then scales it by its sample standard deviation
#' \eqn{(\text{x}-\text{mean}(\text{x}))/\text{sd}(\text{x})}.
#'
#'
#' @param x Numeric vector to standardize
#'
#' @return Standardized vector with sample mean 0 and standard deviation 1
#'
#' @examples
#' x <- c(1, 2, 3, 4, 5)
#' std(x)
#' print(mean(x))
#' print(sd(x))
#'
#' @keywords internal
#' @export
std <- function(x){
  (x-mean(x))/sd(x)
}

#' Compute softplus transform
#'
#' @description
#' Computes the softplus transform, equivalent to the cumulant generating function
#' of a logistic regression model: \eqn{\log(1+e^x)}.
#'
#'
#' @param x Numeric vector to apply softplus to
#'
#' @return Softplus transformed vector
#'
#' @examples
#' x <- runif(5)
#' softplus(x)
#'
#' @keywords internal
#' @export
softplus <- function(x){
  log(1+exp(x))
}

#' Efficient Matrix Multiplication Operator
#'
#' @description
#' Operator wrapper around C++ \code{efficient_matrix_mult()} for matrix multiplication syntax.
#'
#' This is an internal function meant to provide improvement over base R's operator for
#' certain large matrix operations, at a cost of potential slight slowdown for
#' smaller problems.
#'
#' @param x Left matrix
#' @param y Right matrix
#'
#' @return Matrix product of x and y
#'
#' @examples
#' M1 <- matrix(1:4, 2, 2)
#' M2 <- matrix(5:8, 2, 2)
#' M1 %**% M2
#'
#' @keywords internal
`%**%` <- function(x, y) {
  efficient_matrix_mult(x, y)
}

#' Matrix Inversion with Fallback Methods
#'
#' @description
#' Attempts matrix inversion using multiple methods, falling back to more robust
#' approaches if standard inversion fails.
#'
#' @param mat Square matrix to invert
#' @param include_warnings Logical; default FALSE for current implementation.
#'
#' @return Inverted matrix or identity matrix if all methods fail
#'
#' @details
#' Tries methods in order:
#'
#' 1. Cholesky decomposition via \code{chol2inv(chol(...))} for symmetric
#'    positive-definite matrices (fastest for SPD)
#'
#' 2. Direct inversion using \code{armaInv()} as first fallback
#'
#' 3. Generalized inverse using eigendecomposition with small ridge
#'
#' 4. Returns identity matrix with warning if all methods fail
#'
#' For eigendecomposition, uses a small ridge penalty (\code{1e-16}) for stability and
#' zeroes eigenvalues below machine precision.
#'
#' @examples
#' ## Well-conditioned matrix
#' A <- matrix(c(4,2,2,4), 2, 2)
#' invert(A) %*% A
#'
#' ## Singular matrix falls back to M.P. generalized inverse
#' B <- matrix(c(1,1,1,1), 2, 2)
#' invert(B) %*% B
#'
#' @keywords internal
#' @export
invert <- function(mat, include_warnings = FALSE){

  ## [Change 2026-02-12] C2 (B2): Try Cholesky first for SPD matrices
  ## chol2inv(chol(M)) is faster than solve(M) for SPD M
  t <- try({
    chol2inv(chol(mat))
  }, silent = TRUE)

  if(!any(inherits(t, 'try-error'))){
    return(t)
  }

  ## Fallback 1: direct inversion via Armadillo
  t <- try({
    armaInv(mat)
  }, silent = TRUE)

  if(!any(inherits(t, 'try-error'))){
    return(t)
  }

  ## Fallback 2: generalized inverse with small ridge penalty
  # [Change 2026-02-12] C1: replaced %**% with crossprod/tcrossprod
  t <- try({
    eig <- eigen(mat + 1e-16*diag(nrow(mat)), symmetric = TRUE)
    d_inv <- ifelse(eig$values <= sqrt(.Machine$double.eps),
                    0,
                    1/eig$values)
    tcrossprod(eig$vectors * rep(d_inv, each = nrow(mat)), eig$vectors)
  }, silent = TRUE)

  if(!any(inherits(t, 'try-error'))){
    return(t)
  }

  ## Fallback 3: return identity matrix with warning
  # Justification for identity is that for Newton-Raphson, this reduces
  # to approx. gradient-descent. But this is a terrible option for inference!
  if(include_warnings) warning('Matrix not inverted, returning identity: ',
                               print(t))
  return(diag(nrow(mat)))
}

#' Multiply Block Diagonal Matrices in Parallel
#'
#' @description
#' Multiplies two lists of matrices that form block diagonal structures, with optional
#' parallel processing.
#'
#' @param A List of matrices forming first block diagonal matrix
#' @param B List of matrices forming second block diagonal matrix
#' @param K Number of blocks minus 1 (\eqn{K})
#' @param parallel Logical; whether to use parallel processing
#' @param cl Cluster object for parallel processing
#' @param chunk_size Number of blocks per chunk for parallel processing
#' @param num_chunks Number of chunks for parallel processing
#' @param rem_chunks Remaining blocks after chunking
#'
#' @return List containing products of corresponding blocks
#'
#' @details
#' When \code{parallel=TRUE}, splits computation into chunks processed in parallel.
#' Handles remainder chunks separately. Uses \code{matmult_block_diagonal_cpp()} for
#' actual multiplication.
#'
#' The function expects A and B to contain corresponding blocks that can be
#' matrix multiplied.
#'
#' @examples
#' A <- list(matrix(1:4,2,2), matrix(5:8,2,2))
#' B <- list(matrix(1:4,2,2), matrix(5:8,2,2))
#' matmult_block_diagonal(A, B, K=1, parallel=FALSE, cl=NULL,
#'                        chunk_size=1, num_chunks=1, rem_chunks=0)
#'
#' @keywords internal
#' @export
matmult_block_diagonal <- function(A,
                                   B,
                                   K,
                                   parallel,
                                   cl,
                                   chunk_size,
                                   num_chunks,
                                   rem_chunks){
  if(parallel & !is.null(cl)){
    ## Start with remainder chunks
    if(rem_chunks > 0){
      rem <- matmult_block_diagonal_cpp(A[num_chunks*chunk_size + 1:rem_chunks],
                                        B[num_chunks*chunk_size + 1:rem_chunks],
                                        rem_chunks-1)
    } else {
      rem <- list()
    }
    ## Handle the rest in parallel
    c(
      Reduce("c",parallel::parLapply(cl, 1:num_chunks, function(k){
        inds <- 1:chunk_size + (k-1)*chunk_size
        matmult_block_diagonal_cpp(A[inds],
                                   B[inds],
                                   chunk_size-1)
      })),
      rem
    )
  } else {
    matmult_block_diagonal_cpp(A, B, K)
  }
}

#' Generate Design Matrix with Polynomial and Interaction Terms
#'
#' @description
#' Internal function for creating a design matrix containing polynomial
#' expansions and interaction terms for predictor variables. Supports
#' customizable term generation including polynomial degrees up to quartic
#' terms, interaction types, and selective term exclusion.
#'
#' Column names take on the form "_v_" for linear terms, "_v_^d" for polynomial
#' powers up to d = 4, and "_v_x_w_" for interactions between variables v and w,
#' where v and w are column indices of the input predictor matrix.
#'
#' The \code{custom_basis_fxn} argument, if supplied, requires the same arguments
#' as this function, in the same order, minus the eponymous argument,
#' "custom_basis_fxn".
#'
#' @param predictors Numeric matrix of predictor variables
#' @param numerics Integer vector; column indices for variables to expand as polynomials
#' @param just_linear_with_interactions Integer vector; column indices for variables to keep linear but allow interactions
#' @param just_linear_without_interactions Integer vector; column indices for variables to keep linear without interactions
#' @param exclude_interactions_for Integer vector; column indices to exclude from all interactions
#' @param include_quadratic_terms Logical; whether to include squared terms (default TRUE)
#' @param include_cubic_terms Logical; whether to include cubic terms (default TRUE)
#' @param include_quartic_terms Logical; whether to include 4th degree terms (default FALSE)
#' @param include_2way_interactions Logical; whether to include two-way interactions (default TRUE)
#' @param include_3way_interactions Logical; whether to include three-way interactions (default TRUE)
#' @param include_quadratic_interactions Logical; whether to include interactions with squared terms (default TRUE)
#' @param exclude_these_expansions Character vector; names of specific terms to exclude from final matrix
#' @param custom_basis_fxn Function; optional custom basis expansion function that accepts all arguments listed here except itself
#' @param ... Additional arguments passed to \code{custom_basis_fxn}
#'
#' @return Matrix with columns for intercept, polynomial terms, and specified interactions
#'
#' @keywords internal
#' @export
get_polynomial_expansions <- function(predictors,
                                      numerics,
                                      just_linear_with_interactions,
                                      just_linear_without_interactions,
                                      exclude_interactions_for = NULL,
                                      include_quadratic_terms = TRUE,
                                      include_cubic_terms = TRUE,
                                      include_quartic_terms = FALSE,
                                      include_2way_interactions = TRUE,
                                      include_3way_interactions = TRUE,
                                      include_quadratic_interactions = FALSE,
                                      exclude_these_expansions = NULL,
                                      custom_basis_fxn = NULL,
                                      ...) {
  if(any(!is.null(custom_basis_fxn))){
    ## Custom basis-expansion function, if desired
    return(custom_basis_fxn(predictors,
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
                            ...))
  } else {

    ## Dimension setup
    n_rows <- nrow(predictors)
    n_numerics <- length(numerics)
    n_linear_with <- length(just_linear_with_interactions)
    n_linear_without <- length(just_linear_without_interactions)
    n_vars <- ncol(predictors)

    ## Define which variables can participate in interactions
    # Start with potential interaction variables
    vars_that_interact <- unique(c(numerics, just_linear_with_interactions))

    ## Remove any variables specifically excluded from interactions
    if(!is.null(exclude_interactions_for)) {
      vars_that_interact <- setdiff(vars_that_interact,
                                    exclude_interactions_for)
    }

    n_interact_vars <- length(vars_that_interact)

    ## Multiplier for degree of polynomial expansions (3 or 4) we use
    if(include_quartic_terms){
      mult <- 4
    } else {
      mult <- 3
    }
    if(!include_quadratic_terms){
      mult <- mult - 1
    }
    if(!include_cubic_terms){
      mult <- mult - 1
    }

    ## Number of columns pre-allocated
    n_cols <- 1 +
      mult * n_numerics +  # polynomial terms for numerics
      n_linear_with +      # linear terms that can interact
      n_linear_without +   # linear terms that can't interact
      choose(n_interact_vars, 2)*
      include_2way_interactions +  # 2-way interactions
      choose(n_interact_vars, 3)*
      include_3way_interactions +  # 3-way interactions
      max(0, n_numerics * (length(vars_that_interact) - 1))*
      include_quadratic_interactions
    # Note: each numeric can interact quadratically with all other vars


    ## Allocate matrix
    result <- matrix(0, nrow = n_rows, ncol = n_cols)
    col_names <- character(n_cols)
    col_names[1] <- "intercept"
    col_index <- 2

    ## Intercept
    result[, 1] <- 1

    ## Numeric variables and their powers
    if (n_numerics > 0) {
      numeric_data <- predictors[, numerics, drop = FALSE]

      ## Linear terms
      result[, col_index:(col_index + n_numerics - 1)] <-
        numeric_data
      col_names[col_index:(col_index + n_numerics - 1)] <-
        paste0("_", numerics, "_")
      col_index <- col_index + n_numerics

      ## Quadratic terms
      if(include_quadratic_terms){
        result[, col_index:(col_index + n_numerics - 1)] <-
          numeric_data^2
        col_names[col_index:(col_index + n_numerics - 1)] <-
          paste0("_", numerics, "_^2")
        col_index <- col_index + n_numerics
      }

      ## Cubic terms
      if(include_cubic_terms){
        result[, col_index:(col_index + n_numerics - 1)] <-
          numeric_data^3
        col_names[col_index:(col_index + n_numerics - 1)] <-
          paste0("_", numerics, "_^3")
        col_index <- col_index + n_numerics
      }

      ## Quartic terms
      if(include_quartic_terms){
        result[, col_index:(col_index + n_numerics - 1)] <-
          numeric_data^4
        col_names[col_index:(col_index + n_numerics - 1)] <-
          paste0("_", numerics, "_^4")
        col_index <- col_index + n_numerics
      }
    }

    ## Linear terms that can interact
    if (length(just_linear_with_interactions) > 0) {
      result[, col_index:(col_index + n_linear_with - 1)] <-
        predictors[, just_linear_with_interactions, drop = FALSE]
      col_names[col_index:(col_index + n_linear_with - 1)] <-
        paste0("_", just_linear_with_interactions, "_")
      col_index <- col_index + n_linear_with
    }

    ## Linear terms that cannot interact
    if (n_linear_without > 0) {
      result[, col_index:(col_index + n_linear_without - 1)] <-
        predictors[, just_linear_without_interactions, drop = FALSE]
      col_names[col_index:(col_index + n_linear_without - 1)] <-
        paste0("_", just_linear_without_interactions, "_")
      col_index <- col_index + n_linear_without
    }

    ## 2-way Interactions (among variables that can interact)
    if(n_interact_vars > 1 & include_2way_interactions){
      for (i in 1:(length(vars_that_interact) - 1)) {
        for (j in (i+1):length(vars_that_interact)) {
          var_i <- vars_that_interact[i]
          var_j <- vars_that_interact[j]
          result[, col_index] <- predictors[, var_i] * predictors[, var_j]
          col_names[col_index] <- paste0("_", var_i, "_x_", var_j, "_")
          col_index <- col_index + 1
        }
      }
    }

    ## 3-way Interactions (among variables that can interact)
    if(n_interact_vars > 2 & include_3way_interactions){
      for (i in 1:(length(vars_that_interact) - 2)) {
        for (j in (i+1):(length(vars_that_interact) - 1)) {
          for (k in (j+1):length(vars_that_interact)) {
            var_i <- vars_that_interact[i]
            var_j <- vars_that_interact[j]
            var_k <- vars_that_interact[k]
            result[, col_index] <-
              predictors[, var_i] *
              predictors[, var_j] *
              predictors[, var_k]
            col_names[col_index] <- paste0("_",
                                           var_i,
                                           "_x_",
                                           var_j,
                                           "_x_",
                                           var_k,
                                           "_")
            col_index <- col_index + 1
          }
        }
      }
    }

    ## Quadratic interactions
    # Note: only numerics (not in exclude list) can have quadratic terms here
    if(include_quadratic_interactions &
       (length(numerics) + length(just_linear_with_interactions)) > 0 &
       length(intersect(numerics, vars_that_interact)) > 0){
      numerics_that_interact <- intersect(numerics, vars_that_interact)
      for (i in numerics_that_interact) {
        for (j in vars_that_interact) {
          if (i != j) {
            result[, col_index] <- predictors[, i]^2 * predictors[, j]
            col_names[col_index] <- paste0("_", j, "_x_", i, "_^2")
            col_index <- col_index + 1
          }
        }
      }
    }

    ## Set column names and filter empty
    colnames(result) <- col_names
    result <- result[,colnames(result) != ""]

    ## remove custom expansions if desired
    if(any(!is.null(exclude_these_expansions))){
      if(length(exclude_these_expansions) > 0){
        excl <- which(unlist(colnames(result)) %in%
                        unlist(exclude_these_expansions))
        if(length(excl) > 0){
          result <- result[,-excl, drop=FALSE]
        }
      }
    }

    return(result)
  }
}

#' Calculate Derivatives of Polynomial Terms
#'
#' @description
#' Computes first or second derivatives of polynomial terms in a design matrix with
#' respect to a specified variable. Handles polynomial terms up to fourth degree.
#'
#' @param dat Numeric matrix; design matrix containing polynomial basis expansions
#' @param var Character; column name of variable to differentiate with respect to
#' @param second Logical; if TRUE compute second derivative, if FALSE compute first derivative (default FALSE)
#' @param scale Numeric; scaling factor for normalization
#'
#' @return Numeric matrix containing derivatives of polynomial terms, with same dimensions as input matrix
#'
#' @keywords internal
#' @export
take_derivative <- function(dat, var, second = FALSE, scale) {

  ## Initialize matrix for returning
  n_cols <- ncol(dat)
  dat_deriv <- matrix(0, nrow = nrow(dat), ncol = n_cols)
  colnames(dat_deriv) <- colnames(dat)
  variable_values <- dat[, var]

  ## For each colum nof the expansions
  for (i in 1:n_cols) {
    col_name <- colnames(dat)[i]
    match <- regexpr(var, col_name)[[1]]

    ## Check for quadratic terms of var
    qv <- paste0(var,'^2')
    prefix <- substr(col_name, 1, regexpr('_x', col_name)[[1]])
    col_name_sub <- substr(col_name, nchar(prefix)+2, nchar(col_name))
    match2 <- 1*(col_name_sub == qv) + 1*(col_name == qv)

    ## Derivative multiplier based on polynomial degree
    if(match2){
      mult <- 2
    } else if(substr(col_name, nchar(col_name)-1, nchar(col_name)) == '^3'){
      if(second){
        mult <- 6
      } else{
        mult <- 3
      }
    } else if(substr(col_name, nchar(col_name)-1, nchar(col_name)) == '^4'){
      if(second){
        mult <- 12
      } else{
        mult <- 4
      }
    } else {
      mult <- 1
    }

    ## Detect if basis expansion is a function of variable of interest
    # using column title
    if (match > 0) {

      ## Second derivative
      if(second){
        dat_deriv[, i] <- mult * dat[, i] /
          (variable_values^2 + 1*(variable_values == 0))

        ## First derivative
      } else {
        dat_deriv[, i] <- mult * dat[, i] /
          (variable_values + 1*(variable_values == 0))
      }
    }
  }

  ## For second derivatives, set linear terms to 0
  if(second){
    try(dat_deriv[,var] <- 0, silent = TRUE) # may be missing from manual excl.
  }

  return(dat_deriv)
}

#' Calculate Second Derivatives of Interaction Terms
#'
#' @description
#' Computes partial second derivatives for interaction terms including
#' two-way linear, quadratic, and three-way interactions. Handles special cases
#' for each type.
#'
#' This function is necessary to compute total second derivatives as the sum of
#' second partial "pure" derivatives (\eqn{d^2/dx^2}) plus second partial "mixed"
#' derivative (\eqn{d^2/dxdz}), for a predictor x and all other predictors z.
#'
#' @param dat Numeric matrix; design matrix containing basis expansions
#' @param var Character; variable name to differentiate with respect to
#' @param interaction_single_cols Integer vector; column indices for linear-linear interactions
#' @param interaction_quad_cols Integer vector; column indices for linear-quadratic interactions
#' @param triplet_cols Integer vector; column indices for three-way interactions
#' @param colnm_expansions Character vector; column names of expansions for each partition
#' @param power1_cols Integer vector; column indices of linear terms
#' @param power2_cols Integer vector; column indices of quadratic terms
#'
#' @return Numeric matrix of second derivatives, same dimensions as input
#'
#' @keywords internal
#' @export
take_interaction_2ndderivative <-
  function(dat,
           var,
           interaction_single_cols,
           interaction_quad_cols,
           triplet_cols,
           colnm_expansions,
           power1_cols,
           power2_cols) {

    ## Initialize output matrix
    n_cols <- ncol(dat)
    dat_deriv <- matrix(0, nrow = nrow(dat), ncol = n_cols)
    colnames(dat_deriv) <- colnames(dat)
    variable_values <- dat[,var]

    ## Index of linear term
    v <- which(colnm_expansions[power1_cols] == var)

    ## Detect interactions, if relevant
    if (length(interaction_single_cols) > 0) {
      interaction_singles <-
        interaction_single_cols[grep(paste0("_", var, "_"),
                                     colnm_expansions[interaction_single_cols])]
      if (length(interaction_singles) > 0) {
        ## 2nd derivative of two-way linear-linear interactions is always 1
        dat_deriv[, interaction_singles] <- 1
      }
    }

    ## For two-way linear-quadratic interactions, more work is needed
    if (length(interaction_quad_cols) > 0) {
      interaction_quads <-
        interaction_quad_cols[grep(paste0("_", var, "_"),
                                   colnm_expansions[interaction_quad_cols])]
      if (length(interaction_quads) > 0) {
        for (w in 1:length(power1_cols[-v])) {
          ## The other variable, with interactions affecting quadratic terms
          wvar <- c(power1_cols[-v])[w]

          ## Quadratic interaction indices
          interq <-
            interaction_quads[grep(colnm_expansions[wvar],
                                   colnm_expansions[interaction_quads])]
          if (length(interq) > 0) {
            ## This is the _var^2_x_w term
            if(length(power2_cols) > 0){
              nchv <- nchar(colnm_expansions[power2_cols[v]])
              interqv2 <-
                interq[substr(colnm_expansions[interq],
                              nchar(colnm_expansions[interq]) - nchv + 1,
                              nchar(colnm_expansions[interq])) ==
                         colnm_expansions[power2_cols[v]]]

              ## This is the _var_x_w_^2 term
              if(length(interqv2) > 0){
                interqv1 <- interq[-which(interq == interqv2)]
              }
            } else {
              ## This is the _var_x_w_^2 term
              interqv1 <- interq
            }

            ## 2nd derivative of each is 2*other, 2*other + 2*self respectively
            dat_deriv[, interqv1] <- 2 * dat[, colnm_expansions[wvar]]
            if(length(power2_cols) > 0){
              if(length(interqv2) > 0){
                dat_deriv[, interqv2] <-
                  dat_deriv[, interqv1] + 2 * dat[, var]
              }
            }
          }
        }
      }
    }

    ## 2nd deriv of triplets w.r.t.
    # one variable is always the sum of the other two variables
    if (length(triplet_cols) > 0) {
      triplets <-
        triplet_cols[grep(paste0("_", var, "_"),
                          colnm_expansions[triplet_cols])]
      if (length(triplets) > 0) {
        other2_vars <- lapply(triplets, function(tr) {
          vars <- unlist(strsplit(colnm_expansions[tr], 'x'))
          vars[vars != var]
        })
        for (tr in 1:length(other2_vars)) {
          dat_deriv[, triplets[tr]] <-
            dat[, other2_vars[[tr]][1]] +
            dat[, other2_vars[[tr]][2]]
        }
      }
    }

    return(dat_deriv)
  }


#' Expand Matrix into Partition Lists Based on Knot Boundaries
#'
#' @description
#' Takes an input \eqn{N \times p} matrix of polynomial expansions and outputs a list of
#' length \eqn{K+1}, isolating the rows of the input corresponding to assigned partition.
#'
#' @param partition_codes Numeric vector; values determining partition assignment for each row
#' @param partition_bounds Numeric vector; ordered knot locations defining partition boundaries
#' @param nr Integer; number of rows in input matrix
#' @param mat Numeric matrix; data to be partitioned
#' @param K Integer; number of interior knots (resulting in \eqn{K+1} partitions)
#'
#' @return List of length \eqn{K+1}, each element containing the submatrix for that partition
#'
#' @keywords internal
#' @export
knot_expand_list <- function(partition_codes,
                             partition_bounds,
                             nr,
                             mat,
                             K){
  ## No knots, return the input matrix
  if(K == 0){
    return(list(mat))
  } else {
    ## Lower and upper bounds are always infinite
    partition_bounds_app <- c(-Inf, partition_bounds, Inf)

    ## Isolate each partition
    expansions <- lapply(1:(K+1), function(partition){
      mat[partition_codes > partition_bounds_app[partition] &
            partition_codes <= partition_bounds_app[partition+1],,drop=FALSE]
    })
    return(expansions)
  }
}


#' Create Block Diagonal Matrix
#'
#' @param matrix_list List of matrices to arrange diagonally
#' @return Block diagonal matrix with input matrices on diagonal
#'
#' @details
#' Takes in a list of matrices, and returns a block-diagonal matrix with each
#' element of the list as one block. All off-diagonal elements are 0.
#' Matrices must have compatible dimensions.
#'
#' @keywords internal
#' @export
create_block_diagonal <- function(matrix_list) {

  ## Calculate the total dimensions of the resulting matrix
  total_dim <- sum(sapply(matrix_list, nrow))

  ## Create an empty matrix filled with zeros
  result <- matrix(0, nrow = total_dim, ncol = total_dim)

  ## Fill the diagonal blocks
  current_row <- 1
  current_col <- 1

  for (mat in matrix_list) {
    dim_mat <- nrow(mat)
    result[current_row:(current_row + dim_mat - 1),
           current_col:(current_col + dim_mat - 1)] <- mat
    current_row <- current_row + dim_mat
    current_col <- current_col + dim_mat
  }

  return(result)
}


#' Compute Gram Matrix for Block Diagonal Structure
#'
#' @param list_in List of matrices
#' @param parallel Use parallel processing
#' @param cl Cluster object
#' @param chunk_size Chunk size for parallel
#' @param num_chunks Number of chunks
#' @param rem_chunks Remaining chunks
#'
#' @return List of Gram matrices (\eqn{\textbf{X}^{T}\textbf{X}}) for each block
#'
#' @details
#' For a list of matrices, will compute the gram matrix of each element of the
#' list.
#'
#'
#' @keywords internal
#' @export
compute_gram_block_diagonal <- function(list_in,
                                        parallel,
                                        cl,
                                        chunk_size,
                                        num_chunks,
                                        rem_chunks) {
  if(parallel & !is.null(cl)) {

    # Handle remainder chunks first
    if(rem_chunks > 0) {
      rem_indices <- num_chunks*chunk_size + 1:rem_chunks
      rem <- lapply(rem_indices, function(k, func) {
        crossprod(list_in[[k]])
      })
    } else {
      rem <- list()
    }

    # Process main chunks in parallel
    result <- c(
      Reduce("c", parallel::parLapply(cl, 1:num_chunks, function(chunk) {
        inds <- (chunk - 1)*chunk_size + 1:chunk_size
        lapply(inds, function(k) {
          crossprod(list_in[[k]])
        })
      })),
      rem
    )
  } else {
    result <- lapply(list_in, crossprod)
  }
  return(result)
}

#' Compute First and Second Derivative Matrices
#'
#' @param nc Number of columns
#' @param Cpredictors Predictor matrix
#' @param power1_cols Indices of linear terms of spline effects
#' @param power2_cols Indices of quadratic terms of spline effects
#' @param nonspline_cols Indices of non-spline effects
#' @param interaction_single_cols Indices of first-order interactions
#' @param interaction_quad_cols Indices of quadratic interactions
#' @param triplet_cols Indices of three-way interactions
#' @param K Number of partitions minus 1 (\eqn{K})
#' @param include_2way_interactions Include 2-way interactions
#' @param include_3way_interactions Include 3-way interactions
#' @param include_quadratic_interactions Include quadratic interactions
#' @param colnm_expansions Column names
#' @param expansion_scales Scale factors
#' @param just_first_derivatives Only compute first derivatives
#' @param just_spline_effects Only compute derivatives for spline effects
#'
#' @return List containing first and second derivative matrices
#'
#' @keywords internal
#' @export
make_derivative_matrix  <-  function(
    nc,
    Cpredictors,
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
    just_first_derivatives = FALSE,
    just_spline_effects = TRUE){

  ## Include derivatives for non-spline effects, if desired
  if(!just_spline_effects){
    power1_cols <- c(power1_cols, nonspline_cols)
  }

  ## First derivative, for all numeric variables
  first_derivs <- lapply(colnm_expansions[power1_cols], function(v){
    take_derivative(dat = Cpredictors, var = v, scale = expansion_scales[v])
  })
  names(first_derivs) <- colnm_expansions[power1_cols]
  if(just_first_derivatives){
    return(list('first_derivatives' = first_derivs))
  }

  ## 2nd derivatives
  ## If interactions present
  if ((include_2way_interactions |
       include_3way_interactions |
       include_quadratic_interactions) &
      (length(power1_cols) > 1)){
    second_derivs <- lapply(
      colnm_expansions[power1_cols], function(v) {
        take_interaction_2ndderivative(
          dat = Cpredictors,
          var = v,
          interaction_single_cols,
          interaction_quad_cols,
          triplet_cols,
          colnm_expansions,
          power1_cols,
          power2_cols
        ) +
          take_derivative(dat = Cpredictors,
                          var = v,
                          second = TRUE,
                          scale = expansion_scales[v])
      }
    )
    names(second_derivs) <- colnm_expansions[power1_cols]
    ## No interactions present
  } else {
    second_derivs <- lapply(colnm_expansions[power1_cols], function(v){
      take_derivative(dat = Cpredictors,
                      var = v,
                      second = TRUE,
                      scale = expansion_scales[v])
    })
    names(second_derivs) <- colnm_expansions[power1_cols]
  }

  return(list('first_derivatives' = first_derivs,
              'second_derivatives' = second_derivs))
}

#' Create Smoothing Spline Constraint Matrix
#'
#' @description
#' Constructs constraint matrix \eqn{\textbf{A}} enforcing continuity and smoothness at knot boundaries
#' by constraining function values, derivatives, and interactions between partitions.
#'
#' @param nc Integer; number of columns in basis expansion
#' @param CKnots Matrix; basis expansions evaluated at knot points
#' @param power1_cols Integer vector; indices of linear terms
#' @param power2_cols Integer vector; indices of quadratic terms
#' @param nonspline_cols Integer vector; indices of non-spline terms
#' @param interaction_single_cols Integer vector; indices of linear interaction terms
#' @param interaction_quad_cols Integer vector; indices of quadratic interaction terms
#' @param triplet_cols Integer vector; indices of three-way interaction terms
#' @param K Integer; number of interior knots (\eqn{K+1} partitions)
#' @param include_constrain_fitted Logical; constrain function values at knots
#' @param include_constrain_first_deriv Logical; constrain first derivatives at knots
#' @param include_constrain_second_deriv Logical; constrain second derivatives at knots
#' @param include_constrain_interactions Logical; constrain interaction terms at knots
#' @param include_2way_interactions Logical; include two-way interactions
#' @param include_3way_interactions Logical; include three-way interactions
#' @param include_quadratic_interactions Logical; include quadratic interactions
#' @param colnm_expansions Character vector; column names for basis expansions
#' @param expansion_scales Numeric vector; scaling factors for standardization
#'
#' @return Matrix \eqn{\textbf{A}} of constraint coefficients. Columns correspond to
#' constraints, rows to coefficients across all \eqn{K+1} partitions.
#'
#' @keywords internal
#' @export
make_constraint_matrix <- function(nc,
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
                                   expansion_scales){

  ## First, create a checkered matrix-each column alternates 1s/-1s,
  # 0s on off-diagonals
  checkered = matrix(0, nrow = K, ncol = K+1)

  ## Handle matrix construction based on dimension and CKnots structure
  if(!(any(is.null(rownames(CKnots))))){
    rwnms <- rownames(CKnots)
    if(length(rwnms) == 1){
      if(rwnms == 'CKnots'){
        ## For univariate case: Alternate 1/-1 in adjacent columns
        for(j in 1:(K+1)){
          checkered[,j] = -(2*(j %% 2) - 1)
        }
        ## Zero out non-adjacent entries
        for(i in 1:K){
          checkered[i,-c(i,i+1)] = 0
        }
      }
    } else {
      ## For multivariate case: Set pairwise constraints from rownames
      # Rownames code the neighbors i.e. 2_4 implies partitions 2 and 4
      # are neighbors
      for(i in 1:nrow(CKnots)){
        rown <- unlist(strsplit(rwnms[i], "_"))
        col1 <- as.numeric(rown[1])
        col2 <- as.numeric(rown[2])
        checkered[i,col1] <- 1
        checkered[i,col2] <- -1
      }
    }
  } else {
    ## Default univariate construction when no rownames present
    for(j in 1:(K+1)){
      checkered[,j] = -(2*(j %% 2) - 1)
    }
    for(i in 1:K){
      checkered[i,-c(i,i+1)] = 0
    }
  }

  ## Special case handling for single knot with multiple variables
  # (Fix for K = 1, q > 1)
  if(length(checkered) == 2 & !any(unique(checkered) != 0)){
    checkered[1] <- 1
    checkered[2] <- -1
  }

  ## Expand out the checkered matrix to match dimensions of
  # basis expansion dimensions
  checkered_fitted_expand = checkered[, rep(1:(K+1), each = nc), drop = FALSE]

  ## Expand out the constraints, located at knots
  constrain_fitted <- CKnots[,rep(1:nc, K+1), drop = FALSE] *
    checkered_fitted_expand

  ## When non-spline and spline present, repeat fitted constraint for
  # spline-only
  if(length(nonspline_cols) > 0 & length(power1_cols) > 0){
    CKnots0 <- CKnots
    CKnots0[,nonspline_cols] <- 0
    constrain_fitted0 <- CKnots0[,rep(1:nc, K+1), drop = FALSE] *
      checkered_fitted_expand
    constrain_fitted <- rbind(constrain_fitted, constrain_fitted0)
  }

  ## Zero-out the fitted constraint if not desired
  if(!include_constrain_fitted){
    constrain_fitted <- 0 * constrain_fitted
  }

  ## First derivative, for all numeric variables
  if(include_constrain_first_deriv){
    first_derivs <- lapply(colnames(CKnots)[power1_cols], function(v){
      take_derivative(dat = CKnots, var = v, scale = expansion_scales[v])
    })
    first_deriv = Reduce('rbind',
                         first_derivs
    )
    constrain_first_deriv <- first_deriv[, rep(1:nc, K+1),
                                         drop = FALSE] *
      checkered_fitted_expand[

        rep(c(1:nrow(checkered_fitted_expand)),
            length(power1_cols)),
      ]
  }

  ## First derivatives for non-spline effects
  # = identical coefficients across partitions
  if(length(nonspline_cols) > 0 &
     include_constrain_first_deriv){
    first_derivs <- lapply(colnames(CKnots)[nonspline_cols], function(v){
      take_derivative(dat = CKnots,
                      var = v,
                      scale = expansion_scales[v])
    })
    first_derivs <- Reduce('rbind',
                           first_derivs
    )

    ## New constraint
    new_constr <- first_derivs[, rep(1:nc, K+1), drop = FALSE] *
      checkered_fitted_expand[
        rep(c(1:nrow(checkered_fitted_expand)),
            length(nonspline_cols)),
      ]
    new_constr <- unique(new_constr[rowSums(abs(new_constr)) > 0,,drop=FALSE])

    ## Constrain first derivative
    if(length(power1_cols) > 0){
      constrain_first_deriv <- rbind(
        constrain_first_deriv,
        new_constr
      )
    } else {
      constrain_first_deriv <- new_constr
    }
  }

  ## Second derivative, for all numeric variables
  if(include_constrain_second_deriv){
    second_deriv = Reduce('rbind',
                          lapply(colnames(CKnots)[c(power1_cols,
                                                    nonspline_cols)],
                                 function(v){
                                   take_derivative(dat = CKnots,
                                                   var = v,
                                                   second = TRUE,
                                                   scale = expansion_scales[v])
                                 })
    )
    constrain_second_deriv <- second_deriv[, rep(1:nc, K+1),
                                           drop = FALSE] *
      checkered_fitted_expand[
        rep(c(1:nrow(checkered_fitted_expand)),
            length(c(power1_cols, nonspline_cols))),
      ]
  }

  ## If interactions present
  if ((include_2way_interactions |
       include_3way_interactions |
       include_quadratic_interactions) &
      (length(c(power1_cols, nonspline_cols)) > 1) &
      include_constrain_interactions){
    second_deriv_interaction <-
      Reduce('rbind',
             lapply(colnames(CKnots)[
               c(power1_cols, nonspline_cols)],
               function(v) {
                 if(v %in% nonspline_cols){
                   take_derivative(dat = CKnots,
                                   var = v,
                                   scale = expansion_scales[v],
                                   second = TRUE)
                 } else {
                   take_interaction_2ndderivative(
                     dat = CKnots,
                     var = v,
                     interaction_single_cols,
                     interaction_quad_cols,
                     triplet_cols,
                     colnm_expansions,
                     power1_cols,
                     power2_cols
                   )
                 }
               }))
    constrain_second_deriv_interactions <-
      second_deriv_interaction[, rep(1:nc, K+1), drop = FALSE] *
      checkered_fitted_expand[
        rep(c(1:nrow(checkered_fitted_expand)),
            length(c(power1_cols,nonspline_cols))),
      ]
  }

  ## Combine the constraints into a single full-rank matrix A
  if ((include_2way_interactions |
       include_3way_interactions |
       include_quadratic_interactions) &
      (length(c(power1_cols, nonspline_cols)) > 1) &
      include_constrain_interactions){
    if (include_constrain_second_deriv) {
      if (include_constrain_first_deriv) {
        A = t(rbind(
          constrain_fitted,
          constrain_first_deriv,
          constrain_second_deriv +
            constrain_second_deriv_interactions
        ))
      } else {
        A = t(rbind(constrain_fitted,
                    constrain_second_deriv +
                      constrain_second_deriv_interactions))
      }
    } else {
      if (include_constrain_first_deriv) {
        A = t(rbind(constrain_fitted,
                    constrain_first_deriv,
                    constrain_second_deriv_interactions))
      } else {
        A = t(rbind(constrain_fitted,
                    constrain_second_deriv_interactions))
      }
    }
  } else {
    if (include_constrain_second_deriv) {
      if (include_constrain_first_deriv) {
        A = t(rbind(
          constrain_fitted,
          constrain_first_deriv,
          constrain_second_deriv
        ))
      } else {
        A = t(rbind(constrain_fitted,
                    constrain_second_deriv))
      }
    } else {
      if (include_constrain_first_deriv) {
        A = t(rbind(constrain_fitted,
                    constrain_first_deriv))
      } else {
        ## For compatibility
        A = t(rbind(constrain_fitted,
                    constrain_fitted))
      }
    }
  }

  ## If all 0s, return the 0s in conformable form
  if(!any(A != 0)){
    return(cbind(rep(0, (K+1)*nc)))
  } else {
    ## Otherwise, return A
    return(A)
  }
}

#' Test if Vector is Binary
#'
#' @param x Vector to test
#' @return Logical indicating if x has exactly 2 unique values
#'
#' @keywords internal
#' @export
is_binary <- function(x){
  if(length(unique(x)) > 2){
    return(FALSE)
  }
  TRUE
}


#' Compute Derivative of Penalty Matrix G with Respect to Lambda
#'
#' @description
#' Calculates the derivative of the penalty matrix \eqn{\textbf{G}} with respect to the
#' smoothing parameter lambda (\eqn{\lambda}), supporting both global and partition-specific penalties.
#' This is related to the derivative of the diagonal weight matrix \eqn{1/(1+\textbf{x}^{T}\textbf{U}\textbf{G}\textbf{x})} w.r.t. the penalty.
#'
#' @param G A list of penalty matrices \eqn{\textbf{G}} for each partition
#' @param L The base penalty matrix \eqn{\textbf{L}}
#' @param K Number of partitions minus 1 (\eqn{K})
#' @param lambda Smoothing parameter value \eqn{\lambda}
#' @param unique_penalty_per_partition Logical indicating partition-specific penalties
#' @param L_partition_list Optional list of partition-specific penalty matrices \eqn{\textbf{L}_k}
#' @param parallel Logical to enable parallel processing
#' @param cl Cluster object for parallel computation
#' @param chunk_size Size of chunks for parallel processing
#' @param num_chunks Number of chunks
#' @param rem_chunks Remainder chunks
#'
#' @return
#' A list of derivative matrices \eqn{d\textbf{G}/d\lambda} for each partition
#'
#' @keywords internal
#' @export
compute_dG_dlambda <- function(G,
                               L,
                               K,
                               lambda,
                               unique_penalty_per_partition,
                               L_partition_list,
                               parallel,
                               cl,
                               chunk_size,
                               num_chunks,
                               rem_chunks) {
  if(!unique_penalty_per_partition){
    ## Compute negL_lambda once
    negL_lambda <- -lambda * L
  } else {
    negL_lambda <- lapply(L_partition_list, function(L_l){
      -unlist(lambda) * (L + L_l)
    })
  }

  ## [2026-02-12] replaced %**% with crossprod/tcrossprod
  if(parallel & !is.null(cl)) {
    ## Handle remainder chunks first
    if(rem_chunks > 0) {
      rem_indices <- num_chunks*chunk_size + 1:rem_chunks
      rem <- lapply(rem_indices, function(k) {
        if(unique_penalty_per_partition){
          crossprod(t(G[[k]]), crossprod(t(-negL_lambda[[k]]), G[[k]]))
        } else {
          crossprod(t(G[[k]]), crossprod(t(negL_lambda), G[[k]]))
        }
      })
    } else {
      rem <- list()
    }

    ## Process main chunks in parallel
    result <- c(
      Reduce("c",parallel::parLapply(cl, 1:num_chunks, function(chunk) {
        inds <- (chunk - 1)*chunk_size + 1:chunk_size
        lapply(inds, function(k) {
          if(unique_penalty_per_partition){
            crossprod(t(G[[k]]), crossprod(t(-negL_lambda[[k]]), G[[k]]))
          } else {
            crossprod(t(G[[k]]), crossprod(t(negL_lambda), G[[k]]))
          }
        })
      })),
      rem
    )

  } else {
    ## Sequential computation
    result <- lapply(1:(K+1), function(k){
      if(unique_penalty_per_partition){
        crossprod(t(G[[k]]), crossprod(t(-negL_lambda[[k]]), G[[k]]))
      } else {
        crossprod(t(G[[k]]), crossprod(t(negL_lambda), G[[k]]))
      }
    })
  }

  return(result)
}

#' Compute Derivative of Penalty Matrix G with Respect to Lambda (Wrapper)
#'
#' @description
#' Wrapper function for computing the derivative of the weight matrix w.r.t lambda \eqn{\lambda}.
#' This involves computing terms related to the derivative of \eqn{1/(1+\textbf{x}^{T}\textbf{U}\textbf{G}\textbf{x})}.
#'
#' @param G A list of penalty matrices \eqn{\textbf{G}} for each partition
#' @param A Constraint matrix \eqn{\textbf{A}}
#' @param GXX List of \eqn{\textbf{G}\textbf{X}^{T}\textbf{X}} products
#' @param Ghalf List of \eqn{\textbf{G}^{1/2}} matrices
#' @param dG_dlambda List of \eqn{d\textbf{G}/d\lambda} matrices
#' @param dGhalf_dlambda List of \eqn{d\textbf{G}^{1/2}/d\lambda} matrices
#' @param AGAInv Inverse of \eqn{\textbf{A}^{T}\textbf{G}\textbf{A}}
#' @param nc Number of columns
#' @param K Number of partitions minus 1 (\eqn{K})
#' @param parallel Logical to enable parallel processing
#' @param cl Cluster object for parallel computation
#' @param chunk_size Size of chunks for parallel processing
#' @param num_chunks Number of chunks
#' @param rem_chunks Remainder chunks
#'
#' @return
#' Scalar value representing the trace derivative component.
#'
#' @keywords internal
#' @export
compute_dW_dlambda_wrapper <- function(G,
                                       A,
                                       GXX,
                                       Ghalf,
                                       dG_dlambda,
                                       dGhalf_dlambda,
                                       AGAInv,
                                       nc,
                                       K,
                                       parallel,
                                       cl,
                                       chunk_size,
                                       num_chunks,
                                       rem_chunks) {

  if(parallel & !is.null(cl)) {
    # Pre-calculate total size and indices
    total_chunks <- num_chunks + (rem_chunks > 0)
    chunk_indices <- vector("list", total_chunks)

    # Create all indices at once
    for(i in 1:num_chunks) {
      chunk_indices[[i]] <- (i-1)*chunk_size + 1:chunk_size
    }
    if(rem_chunks > 0) {
      chunk_indices[[total_chunks]] <- num_chunks*chunk_size + 1:rem_chunks
    }

    # Pre-compute dGXX for all partitions in one step
    dGXX <- matmult_block_diagonal_cpp(dG_dlambda, GXX, K)

    # Compute all traces and corrections in a single parallel operation
    results <- parallel::parLapply(cl, chunk_indices, function(inds) {
      start_idx <- min(inds) - 1

      # Compute initial trace
      trace_part <- sum(sapply(inds, function(k)
        sum(diag(crossprod(t(dG_dlambda[[k]]), GXX[[k]])))))

      # Get relevant chunk of A
      A_chunk <- A[(start_idx*nc + 1):min((max(inds))*nc, nrow(A)), ]
      len_inds <- length(inds) - 1

      # Compute both corrections at once
      correction1 <- compute_trace_correction(
        G[inds], A_chunk, GXX[inds],
        AGAInv, nc, len_inds)

      correction2 <- compute_trace_correction(
        G[inds], A_chunk, dGXX[inds],
        AGAInv, nc, len_inds)

      c(trace_part, correction1, correction2)
    })

    # Efficiently sum results
    all_results <- do.call(rbind, results)
    trace <- sum(all_results[,1])
    correction1 <- sum(all_results[,2])
    correction2 <- sum(all_results[,3])

    return(trace - correction1 - correction2)

  } else {
    return(compute_dW_dlambda(G,
                              A,
                              GXX,
                              Ghalf,
                              dG_dlambda,
                              dGhalf_dlambda,
                              AGAInv,
                              nc,
                              K))
  }
}


#' Calculate Trace of Matrix Product \eqn{\text{trace}(\textbf{X}\textbf{U}\textbf{G}\textbf{X}^{T})}
#'
#' @param G List of G matrices (\eqn{\textbf{G}})
#' @param A Constraint matrix (\eqn{\textbf{A}})
#' @param GXX List of \eqn{\textbf{G}\textbf{X}^{T}\textbf{X}} products
#' @param AGAInv Inverse of \eqn{\textbf{A}^{T}\textbf{G}\textbf{A}}
#' @param nc Number of columns
#' @param nca Number of constraint columns
#' @param K Number of partitions minus 1 (\eqn{K})
#' @param parallel Use parallel processing
#' @param cl Cluster object
#' @param chunk_size Size of parallel chunks
#' @param num_chunks Number of chunks
#' @param rem_chunks Remaining chunks
#'
#' @details
#' Computes \eqn{\text{trace}(\textbf{X}\textbf{U}\textbf{G}\textbf{X}^{T})} where \eqn{\textbf{U} = \textbf{I} - \textbf{G}\textbf{A}(\textbf{A}^{T}\textbf{G}\textbf{A})^{-1}\textbf{A}^{T}}.
#' Handles parallel computation by splitting into chunks.
#'
#' @return Trace value
#'
#' @keywords internal
#' @export
compute_trace_UGXX_wrapper <- function(G,
                                       A,
                                       GXX,
                                       AGAInv,
                                       nc,
                                       nca,
                                       K,
                                       parallel,
                                       cl,
                                       chunk_size,
                                       num_chunks,
                                       rem_chunks) {

  # P <- A %**% AGAInv %**% t(A)
  # U <- -t(matmult_U(P, G, nc, K)) + diag(nrow(P))
  # UGXX <- matmult_U(U, GXX, nc, K)
  # return(sum(diag(UGXX)))

  if(parallel & !is.null(cl)) {
    ## Pre-calculate total size and indices
    total_chunks <- num_chunks + (rem_chunks > 0)
    chunk_indices <- vector("list", total_chunks)

    ## Create all indices at once
    for(i in 1:num_chunks) {
      chunk_indices[[i]] <- (i-1)*chunk_size + 1:chunk_size
    }
    if(rem_chunks > 0) {
      chunk_indices[[total_chunks]] <- num_chunks*chunk_size + 1:rem_chunks
    }

    ## Compute traces in parallel with single parallel::parLapply call
    results <- parallel::parLapply(cl, chunk_indices, function(inds) {
      start_idx <- min(inds) - 1

      ## Compute first trace
      trace_part <- mean(sapply(GXX[inds], function(gxx)rowMeans(t(diag(gxx)))))

      ## Compute correction
      const <- length(inds) * nc
      A_chunk <- A[(start_idx*nc + 1):min((max(inds))*nc, nrow(A)), ]
      correction_part <- compute_trace_correction(G[inds],
                                                  A_chunk / sqrt(const),
                                                  GXX[inds],
                                                  AGAInv,
                                                  nc,
                                                  length(inds)-1)

      c(trace_part, correction_part) * const
    })

    ## Efficiently sum results
    all_results <- do.call(rbind, results)
    trace <- mean(all_results[,1])
    correction <- mean(all_results[,2])

    return(min(max((trace - correction)*nrow(all_results), 0), nc*(K+1)))

  } else {
    trace_part <- rowSums(
      rbind(sapply(GXX, function(gxx)rowSums(t(diag(gxx)))))
    )
    correction_part <-
      compute_trace_correction(G,
                               A,
                               GXX,
                               AGAInv,
                               nc,
                               K)
    return(min(max(trace_part - correction_part, 0), (K + 1)*nc))
  }
}

## [Change 2026-02-15] replace compute trace UGXX with more stable and efficient
# version below

#' Effective degrees of freedom via trace of the hat matrix
#'
#' Computes \eqn{\mathrm{tr}(\mathbf{X}\mathbf{U}\mathbf{G}\mathbf{X}^\top
#' \mathbf{W}^{1/2}\mathbf{V}^{-1}\mathbf{W}^{1/2})} without forming
#' the \eqn{N \times N} hat matrix, by reducing the problem to
#' \eqn{P \times P} and \eqn{r \times r} operations.
#'
#' @details
#' \subsection{Derivation}{
#' Let \eqn{\mathbf{G} = (\mathbf{X}^\top\mathbf{W}^{1/2}\mathbf{V}^{-1}\mathbf{W}^{1/2}\mathbf{X} + \boldsymbol{\Lambda})^{-1}}
#' and \eqn{\mathbf{U} = \mathbf{I} - \mathbf{G}\mathbf{A}(\mathbf{A}^\top\mathbf{G}\mathbf{A})^{-1}\mathbf{A}^\top}.
#' By the cyclic property of the trace:
#' \deqn{\mathrm{tr}(\mathbf{X}\mathbf{U}\mathbf{G}\mathbf{X}^\top\mathbf{W}^{1/2}\mathbf{V}^{-1}\mathbf{W}^{1/2}) = \mathrm{tr}(\mathbf{U}\mathbf{G}\,\mathbf{X}^\top\mathbf{W}^{1/2}\mathbf{V}^{-1}\mathbf{W}^{1/2}\mathbf{X})}
#'
#' Substituting \eqn{\mathbf{X}^\top\mathbf{W}^{1/2}\mathbf{V}^{-1}\mathbf{W}^{1/2}\mathbf{X} = \mathbf{G}^{-1} - \boldsymbol{\Lambda}}:
#' \deqn{= \mathrm{tr}(\mathbf{U}\mathbf{G}(\mathbf{G}^{-1} - \boldsymbol{\Lambda})) = \mathrm{tr}(\mathbf{U}) - \mathrm{tr}(\mathbf{U}\mathbf{G}\boldsymbol{\Lambda})}
#'
#' Since \eqn{\mathbf{U}} is idempotent,
#' \eqn{\mathrm{tr}(\mathbf{U}) = \mathrm{rank}(\mathbf{U}) = P - r}
#' where \eqn{r = \mathrm{rank}(\mathbf{A})}. For the second term,
#' expand \eqn{\mathbf{U} = \mathbf{I} - \mathbf{G}\mathbf{A}(\mathbf{A}^\top\mathbf{G}\mathbf{A})^{-1}\mathbf{A}^\top}:
#' \deqn{\mathrm{tr}(\mathbf{U}\mathbf{G}\boldsymbol{\Lambda}) = \mathrm{tr}(\mathbf{G}\boldsymbol{\Lambda}) - \mathrm{tr}\!\left((\mathbf{A}^\top\mathbf{G}\mathbf{A})^{-1}\mathbf{A}^\top\mathbf{G}\boldsymbol{\Lambda}\mathbf{G}\mathbf{A}\right)}
#'
#' using the cyclic property on the second term. Both \eqn{\mathbf{G}}
#' and \eqn{\boldsymbol{\Lambda}} are block-diagonal with \eqn{K+1}
#' blocks of dimension \eqn{p \times p}, so:
#' \deqn{\mathrm{tr}(\mathbf{G}\boldsymbol{\Lambda}) = \sum_{k=1}^{K+1}\mathrm{tr}(\mathbf{G}_k\boldsymbol{\Lambda}_k)}
#' and \eqn{\mathbf{G}\boldsymbol{\Lambda}\mathbf{G}} is also block-diagonal
#' with blocks \eqn{\mathbf{G}_k\boldsymbol{\Lambda}_k\mathbf{G}_k}.
#'
#' Combining:
#' \deqn{\mathrm{tr}(\mathbf{X}\mathbf{U}\mathbf{G}\mathbf{X}^\top\mathbf{W}^{1/2}\mathbf{V}^{-1}\mathbf{W}^{1/2}) = (P - r) - \sum_{k=1}^{K+1}\mathrm{tr}(\mathbf{G}_k\boldsymbol{\Lambda}_k) + \mathrm{tr}\!\left((\mathbf{A}^\top\mathbf{G}\mathbf{A})^{-1}\mathbf{A}^\top\mathbf{G}\boldsymbol{\Lambda}\mathbf{G}\mathbf{A}\right)}
#'
#' The correction term
#' \eqn{\mathrm{tr}((\mathbf{A}^\top\mathbf{G}\mathbf{A})^{-1}\mathbf{A}^\top\mathbf{G}\boldsymbol{\Lambda}\mathbf{G}\mathbf{A})}
#' is an \eqn{r \times r} trace and captures the degrees of freedom
#' recovered because the constraints pin certain linear combinations of
#' penalized coefficients, removing them from estimation. When
#' \eqn{\boldsymbol{\Lambda} = \mathbf{0}} (no penalty), both
#' \eqn{\mathrm{tr}(\mathbf{G}\boldsymbol{\Lambda})} and the correction
#' vanish, giving \eqn{\mathrm{edf} = P - r}. When \eqn{\mathbf{U} = \mathbf{I}}
#' (no constraints), \eqn{r = 0} and the correction vanishes, giving
#' \eqn{\mathrm{edf} = P - \mathrm{tr}(\mathbf{G}\boldsymbol{\Lambda})}.
#'
#' The result is invariant to the correlation structure \eqn{\mathbf{V}}
#' and GLM weights \eqn{\mathbf{W}}, which enter only through
#' \eqn{\mathbf{G}}.
#' }
#'
#' \subsection{Partition-specific penalties}{
#' When \code{unique_penalty_per_partition = TRUE}, each partition \eqn{k}
#' has an additional penalty matrix \eqn{\mathbf{L}_k} from
#' \code{L_partition_list}, so the effective per-partition penalty becomes
#' \eqn{\boldsymbol{\Lambda}_k = \boldsymbol{\Lambda} + \mathbf{L}_k}.
#' All block-wise traces and products use \eqn{\boldsymbol{\Lambda}_k}
#' in place of the shared \eqn{\boldsymbol{\Lambda}}.
#' }
#'
#' \subsection{Computational cost}{
#' The partition-wise traces \eqn{\mathrm{tr}(\mathbf{G}_k\boldsymbol{\Lambda}_k)}
#' cost \eqn{O(p^2)} each and are parallelizable. The correction term
#' requires forming \eqn{\mathbf{G}\boldsymbol{\Lambda}\mathbf{G}\mathbf{A}}
#' (block-diagonal times sparse, \eqn{O((K+1)p^2 r)}) and a single
#' \eqn{r \times r} trace (\eqn{O(r^2)}). The total cost is
#' \eqn{O((K+1)p^2 + r^2)}, compared to \eqn{O(N^2)} or \eqn{O(NP)}
#' for forming and tracing the full hat matrix.
#' }
#'
#' @param G List of \eqn{K+1} matrices \eqn{\mathbf{G}_k}, each \eqn{p \times p}.
#' @param Lambda Base penalty matrix \eqn{\boldsymbol{\Lambda}}, \eqn{p \times p},
#'   shared across all partitions. When \code{unique_penalty_per_partition = FALSE},
#'   this is the full per-block penalty. When \code{TRUE}, the effective penalty
#'   for partition \eqn{k} is \eqn{\boldsymbol{\Lambda} + \mathbf{L}_k}.
#' @param A Constraint matrix, \eqn{P \times r}.
#' @param AGAInv Precomputed \eqn{(\mathbf{A}^\top\mathbf{G}\mathbf{A})^{-1}},
#'   \eqn{r \times r}.
#' @param nc Number of columns per partition block (\eqn{p}).
#' @param nca Number of constraint columns (\eqn{r}).
#' @param K Number of interior knots (so there are \eqn{K+1} partitions).
#' @param parallel Logical; use parallel processing.
#' @param cl Cluster object for parallel computation.
#' @param chunk_size Size of parallel chunks.
#' @param num_chunks Number of full-size chunks.
#' @param rem_chunks Number of remaining partitions in the final chunk.
#' @param unique_penalty_per_partition Logical; if \code{TRUE}, each partition
#'   receives an additional penalty from \code{L_partition_list}.
#' @param L_partition_list List of partition-specific penalty matrices
#'   \eqn{\mathbf{L}_k}, each \eqn{p \times p}. Only used when
#'   \code{unique_penalty_per_partition = TRUE}; ignored otherwise.
#'
#' @return Scalar effective degrees of freedom, clamped to \eqn{[0, \, (K+1)p]}.
#'
#' @keywords internal
#' @export
compute_trace_H <- function(G,
                            Lambda,
                            A,
                            AGAInv,
                            nc,
                            nca,
                            K,
                            parallel,
                            cl,
                            chunk_size,
                            num_chunks,
                            rem_chunks,
                            unique_penalty_per_partition,
                            L_partition_list) {

  n_partitions <- K + 1
  P <- nc * n_partitions
  r <- nca

  ## Term 1: rank(U) = P - r
  rank_U <- P - r

  ## Helper: get the effective Lambda block for partition k
  #  When unique_penalty_per_partition, Lambda_k = Lambda + L_partition_list[[k]]
  #  Otherwise Lambda_k = Lambda for all k
  get_Lambda_k <- function(k) {
    if (unique_penalty_per_partition) {
      Lambda + L_partition_list[[k]]
    } else {
      Lambda
    }
  }

  ## Terms 2 and 3 computed jointly per partition
  # Term 2: sum_k tr(G_k Lambda_k)
  # Term 3: B = G Lambda G A (rows assembled block-wise)
  B <- matrix(0, P, r)

  if (parallel & !is.null(cl)) {
    ## Build chunk index lists
    total_chunks <- num_chunks + (rem_chunks > 0)
    chunk_indices <- vector("list", total_chunks)
    for (i in seq_len(num_chunks)) {
      chunk_indices[[i]] <- (i - 1) * chunk_size + seq_len(chunk_size)
    }
    if (rem_chunks > 0) {
      chunk_indices[[total_chunks]] <- num_chunks * chunk_size +
        seq_len(rem_chunks)
    }

    ## Parallel: each chunk computes its partition-wise traces and B rows
    results <- parallel::parLapply(cl, chunk_indices, function(inds) {
      tr_parts <- numeric(length(inds))
      B_chunk <- matrix(0, length(inds) * nc, r)
      for (j in seq_along(inds)) {
        k <- inds[j]
        Lk <- get_Lambda_k(k)
        tr_parts[j] <- sum(G[[k]] * Lk)
        rows_local <- (j - 1) * nc + seq_len(nc)
        rows_global <- (k - 1) * nc + seq_len(nc)
        A_k <- A[rows_global, , drop = FALSE]
        B_chunk[rows_local, ] <- G[[k]] %**% (Lk %**% (G[[k]] %**% A_k))
      }
      list(tr_parts = tr_parts, B_chunk = B_chunk, inds = inds)
    })

    ## Collect parallel results
    tr_GLambda <- 0
    for (res in results) {
      tr_GLambda <- tr_GLambda + sum(res$tr_parts)
      for (j in seq_along(res$inds)) {
        k <- res$inds[j]
        rows_global <- (k - 1) * nc + seq_len(nc)
        rows_local <- (j - 1) * nc + seq_len(nc)
        B[rows_global, ] <- res$B_chunk[rows_local, ]
      }
    }
  } else {
    ## Sequential
    tr_GLambda <- 0
    for (k in seq_len(n_partitions)) {
      Lk <- get_Lambda_k(k)
      tr_GLambda <- tr_GLambda + sum(G[[k]] * Lk)
      rows <- (k - 1) * nc + seq_len(nc)
      A_k <- A[rows, , drop = FALSE]
      B[rows, ] <- G[[k]] %**% (Lk %**% (G[[k]] %**% A_k))
    }
  }

  ## Term 3: tr((A'GA)^{-1} A' G Lambda G A)
  D <- crossprod(A, B)
  tr_correction <- sum(AGAInv * t(D))

  ## Combine
  edf <- rank_U - tr_GLambda + tr_correction

  return(min(max(edf, 0), P))
}

#' Derivative of Constrained Penalized Coefficients with Respect to Lambda
#'
#' Computes \eqn{d(\mathbf{U}\mathbf{G}\mathbf{X}^{T}\mathbf{y})/d\lambda},
#' the sensitivity of the constrained coefficient estimates
#' \eqn{\tilde{\boldsymbol{\beta}} = \mathbf{U}\mathbf{G}\mathbf{X}^T\mathbf{y}}
#' to changes in the smoothing parameter \eqn{\lambda}.
#'
#' @param AGAInv_AGXy Product of \eqn{(\mathbf{A}^{T}\mathbf{G}\mathbf{A})^{-1}} and \eqn{\mathbf{A}^{T}\mathbf{G}\mathbf{X}^{T}\mathbf{y}}
#' @param AGAInv Inverse of \eqn{\mathbf{A}^{T}\mathbf{G}\mathbf{A}}
#' @param G List of \eqn{\mathbf{G}} matrices
#' @param A Constraint matrix \eqn{\mathbf{A}}
#' @param dG_dlambda List of \eqn{d\mathbf{G}/d\lambda} matrices
#' @param nc Number of columns
#' @param nca Number of constraint columns
#' @param K Number of partitions minus 1 (\eqn{K})
#' @param Xy List of \eqn{\mathbf{X}^{T}\mathbf{y}} products
#' @param Ghalf List of \eqn{\mathbf{G}^{1/2}} matrices
#' @param dGhalf List of \eqn{d\mathbf{G}^{1/2}/d\lambda} matrices
#' @param GhalfXy_temp Temporary storage for \eqn{\mathbf{G}^{1/2}\mathbf{X}^{T}\mathbf{y}}
#' @param parallel Use parallel processing
#' @param cl Cluster object
#' @param chunk_size Size of parallel chunks
#' @param num_chunks Number of chunks
#' @param rem_chunks Remaining chunks
#'
#' @details
#' This function is called during GCV/penalty optimization and computes
#' how the constrained coefficient vector changes as the penalty weight
#' varies. Two implementations are provided depending on problem size.
#'
#' \subsection{Derivation}{
#' The constrained estimate is
#' \eqn{\tilde{\boldsymbol{\beta}} = \mathbf{U}\hat{\boldsymbol{\beta}}}
#' where
#' \eqn{\hat{\boldsymbol{\beta}} = \mathbf{G}\mathbf{X}^T\mathbf{y}}
#' and
#' \eqn{\mathbf{U} = \mathbf{I} - \mathbf{G}\mathbf{A}(\mathbf{A}^T\mathbf{G}\mathbf{A})^{-1}\mathbf{A}^T}.
#' Both \eqn{\mathbf{G}} and \eqn{\mathbf{U}} depend on \eqn{\lambda}
#' through \eqn{\mathbf{G} = (\mathbf{X}^T\mathbf{X} + \lambda\boldsymbol{\Lambda})^{-1}}.
#'
#' Differentiating the product \eqn{\mathbf{U}\mathbf{G}\mathbf{X}^T\mathbf{y}}
#' requires the chain rule applied to three \eqn{\lambda}-dependent
#' components:
#' \deqn{\frac{d}{d\lambda}(\mathbf{U}\mathbf{G}\mathbf{X}^T\mathbf{y})}
#' \deqn{= \frac{d\mathbf{U}}{d\lambda}\hat{\boldsymbol{\beta}} + \mathbf{U}\frac{d\hat{\boldsymbol{\beta}}}{d\lambda}}
#'
#' The unconstrained part is straightforward:
#' \eqn{d\hat{\boldsymbol{\beta}}/d\lambda = (d\mathbf{G}/d\lambda)\mathbf{X}^T\mathbf{y}},
#' computed partition-wise since \eqn{\mathbf{G}} is block-diagonal.
#'
#' The constraint projection derivative is more involved. Writing
#' \eqn{\mathbf{C} = (\mathbf{A}^T\mathbf{G}\mathbf{A})^{-1}} and
#' expanding \eqn{d\mathbf{U}/d\lambda} gives three terms (after
#' applying the product rule and the identity
#' \eqn{d\mathbf{C}/d\lambda = -\mathbf{C}(d(\mathbf{A}^T\mathbf{G}\mathbf{A})/d\lambda)\mathbf{C}}):
#' \describe{
#'   \item{term2a}{\eqn{(d\mathbf{G}/d\lambda)\mathbf{A}\mathbf{C}\mathbf{A}^T\hat{\boldsymbol{\beta}}}
#'     direct effect of \eqn{d\mathbf{G}/d\lambda} on the projection}
#'   \item{term2b}{\eqn{\mathbf{G}\mathbf{A}\mathbf{C}(d(\mathbf{A}^T\mathbf{G}\mathbf{A})/d\lambda)\mathbf{C}\mathbf{A}^T\hat{\boldsymbol{\beta}}}
#'     effect through the change in \eqn{(\mathbf{A}^T\mathbf{G}\mathbf{A})^{-1}}}
#'   \item{term2c}{\eqn{\mathbf{G}\mathbf{A}\mathbf{C}\mathbf{A}^T(d\hat{\boldsymbol{\beta}}/d\lambda)}
#'     projection of the unconstrained derivative back through the constraint}
#' }
#'
#' The full derivative is
#' \eqn{d\hat{\boldsymbol{\beta}}/d\lambda} minus the sum of these
#' three correction terms.
#' }
#'
#' \subsection{Large problem path (K >= 10, nc > 4)}{
#' Computes the three correction terms explicitly using the block
#' structure of \eqn{\mathbf{G}} and \eqn{d\mathbf{G}/d\lambda}.
#' The intermediate quantity
#' \eqn{d(\mathbf{A}^T\mathbf{G}\mathbf{A})/d\lambda = \mathbf{A}^T(d\mathbf{G}/d\lambda)\mathbf{A}}
#' is accumulated partition-wise. Shared vectors
#' \eqn{\mathbf{A}\mathbf{C}\mathbf{A}^T\hat{\boldsymbol{\beta}}}
#' and related products are precomputed once and reused across
#' partitions. Parallelism is over chunks of partitions.
#' }
#'
#' \subsection{Small problem path}{
#' Reformulates the derivative using the matrix square root
#' \eqn{\mathbf{G}^{1/2}} and its derivative
#' \eqn{d\mathbf{G}^{1/2}/d\lambda}. The constraint
#' \eqn{\mathbf{A}^T\boldsymbol{\beta} = 0} is imposed via least
#' squares projection: the residuals from regressing
#' \eqn{\mathbf{G}^{1/2}\mathbf{X}^T\mathbf{y}} onto
#' \eqn{\mathbf{G}^{1/2}\mathbf{A}} give the constrained component,
#' and differentiating this projection with respect to \eqn{\lambda}
#' yields the derivative. Uses \code{.lm.fit} for speed with a
#' stabilizing rescaling factor to prevent numerical issues when the
#' constraint matrix is poorly scaled.
#' }
#'
#' @return \eqn{P \times 1} vector of derivatives
#'   \eqn{d\tilde{\boldsymbol{\beta}}/d\lambda}
#'
#' @keywords internal
#' @export
compute_dG_u_dlambda_xy <- function(AGAInv_AGXy,
                                    AGAInv,
                                    G,
                                    A,
                                    dG_dlambda,
                                    nc,
                                    nca,
                                    K,
                                    Xy,
                                    Ghalf,
                                    dGhalf,
                                    GhalfXy_temp,
                                    parallel,
                                    cl,
                                    chunk_size,
                                    num_chunks,
                                    rem_chunks) {

  if((K >= 10) & (nc > 4)){

    ##  ## Large problem path: explicit three-term derivative

    ## [Change 2026-02-12] replaced %**% with crossprod/tcrossprod,
    # rigorously improved documentation

    ## Helper function for chunk processing: computes partition-wise
    # GXy = G_k X_k'y_k, dGXy = (dG_k/dl) X_k'y_k, and accumulates
    # A_k'(dG_k/dl)A_k into the r x r matrix AdGA for the chain rule
    # on (A'GA)^{-1}.
    process_chunk <- function(inds, G_chunk, dG_chunk, A_chunk, Xy_chunk) {
      ## Process GXy and dGXy for chunk
      GXy_chunk <- lapply(seq_along(inds), function(i) {
        crossprod(t(G_chunk[[i]]), Xy_chunk[[i]])
      })

      dGXy_chunk <- lapply(seq_along(inds), function(i) {
        crossprod(t(dG_chunk[[i]]), Xy_chunk[[i]])
      })

      ## Process AdGA for chunk: accumulates A_k'(dG_k/dl)A_k
      AdGA_chunk <- Reduce("+", lapply(seq_along(inds), function(i) {
        k <- inds[i]
        start_idx <- (k-1)*nc + 1
        end_idx <- k*nc
        crossprod(A[start_idx:end_idx,],
                  crossprod(t(dG_chunk[[i]]), A[start_idx:end_idx,]))
      }))

      list(GXy = GXy_chunk, dGXy = dGXy_chunk, AdGA = AdGA_chunk)
    }

    if(parallel & !is.null(cl)) {
      ## Handle remainder chunks
      if(rem_chunks > 0) {
        rem_indices <- num_chunks*chunk_size + 1:rem_chunks
        rem_result <- process_chunk(
          rem_indices,
          G[rem_indices],
          dG_dlambda[rem_indices],
          A[(num_chunks*chunk_size)*nc + 1:((rem_chunks)*nc), ],
          Xy[rem_indices]
        )
      } else {
        rem_result <- list(GXy = list(), dGXy = list(), AdGA = 0)
      }

      ## Process main chunks in parallel
      chunk_results <- parallel::parLapply(cl, 1:num_chunks, function(chunk) {
        start_idx <- (chunk - 1) * chunk_size
        inds <- (start_idx + 1):min(start_idx + chunk_size, K+1)

        process_chunk(
          inds,
          G[inds],
          dG_dlambda[inds],
          A[(start_idx*nc + 1):min((start_idx + chunk_size)*nc, nrow(A)), ],
          Xy[inds]
        )
      })

      ## Combine results
      GXy <- Reduce("c", c(lapply(chunk_results, `[[`, "GXy"),
                           rem_result$GXy))
      dGXy <- Reduce("c", c(lapply(chunk_results, `[[`, "dGXy"),
                            rem_result$dGXy))
      AdGA <- Reduce("+", c(lapply(chunk_results, `[[`, "AdGA"),
                            list(rem_result$AdGA)))

    } else {
      ## Sequential computation
      GXy <- Reduce("c", lapply(1:(K+1), function(k)
        crossprod(t(G[[k]]), Xy[[k]])))
      dGXy <- Reduce("c", lapply(1:(K+1), function(k)
        crossprod(t(dG_dlambda[[k]]), Xy[[k]])))
      AdGA <- Reduce("+", lapply(1:(K+1), function(k) {
        start_idx <- (k-1)*nc + 1
        end_idx <- k*nc
        crossprod(A[start_idx:end_idx,],
                  crossprod(t(dG_dlambda[[k]]), A[start_idx:end_idx,]))
      }))
    }

    ## Precompute shared vectors used by all three correction terms:
    #    A_AGAInv_AGXy  = A * C * A' * beta_hat
    #      (the constraint correction applied to the unconstrained estimate)
    #    A_term2b_mat   = A * (-C * AdGA * C) * A' * beta_hat
    #      (effect of dC/dl on the projection, pre-expanded by A)
    #    A_term2c_mat   = A * C * A' * d(beta_hat)/dl
    #      (the constraint correction applied to the unconstrained derivative)
    A_AGAInv_AGXy <- cbind(crossprod(t(A), AGAInv_AGXy))
    A_term2b_mat <- cbind(crossprod(t(A),
                                    crossprod(t(-crossprod(t(AGAInv), AdGA)), AGAInv_AGXy)))
    A_term2c_mat <- crossprod(t(A),
                              cbind(crossprod(t(AGAInv),
                                              crossprod(A, cbind(unlist(dGXy))))))

    ## Compute the three correction terms partition-wise:
    #    term2a_k = (dG_k/dl) * [A C A' beta_hat]_k
    #    term2b_k = G_k * [A (-C (AdGA) C) A' beta_hat]_k
    #    term2c_k = G_k * [A C A' (d beta_hat/dl)]_k
    if(parallel & !is.null(cl)) {
      # Handle remainder chunks
      if(rem_chunks > 0) {
        rem_indices <- num_chunks*chunk_size + 1:rem_chunks

        term2a_rem <- lapply(rem_indices, function(k) {
          crossprod(t(dG_dlambda[[k]]),
                    A_AGAInv_AGXy[(k-1)*nc + 1:nc,,drop=FALSE])
        })

        term2b_rem <- lapply(rem_indices, function(k) {
          idx <- (k-1)*nc + 1:nc
          crossprod(t(G[[k]]), A_term2b_mat[idx,, drop=FALSE])
        })

        term2c_rem <- lapply(rem_indices, function(k) {
          idx <- (k-1)*nc + 1:nc
          crossprod(t(G[[k]]), A_term2c_mat[idx,, drop=FALSE])
        })
      } else {
        term2a_rem <- term2b_rem <- term2c_rem <- list()
      }

      ## Process main chunks in parallel
      chunk_results <- parallel::parLapply(cl, 1:num_chunks, function(chunk) {
        start_idx <- (chunk - 1) * chunk_size
        inds <- (start_idx + 1):min(start_idx + chunk_size, K+1)

        terms <- lapply(inds, function(k){
          t2a <- crossprod(t(dG_dlambda[[k]]),
                           A_AGAInv_AGXy[(k-1)*nc + 1:nc,, drop=FALSE])
          idx <- (k-1)*nc + 1:nc
          t2b <- crossprod(t(G[[k]]), A_term2b_mat[idx,, drop=FALSE])
          t2c <- crossprod(t(G[[k]]), A_term2c_mat[idx,, drop=FALSE])
          list(t2a,t2b,t2c)
        })
        list(term2a = lapply(terms, `[[`, 1),
             term2b = lapply(terms, `[[`, 2),
             term2c = lapply(terms, `[[`, 3))
      })

      term2a <- Reduce("c", c(lapply(chunk_results, `[[`, "term2a"),
                              term2a_rem))
      term2b <- Reduce("c", c(lapply(chunk_results, `[[`, "term2b"),
                              term2b_rem))
      term2c <- Reduce("c", c(lapply(chunk_results, `[[`, "term2c"),
                              term2c_rem))

    } else {
      terms <- lapply(1:(K+1), function(k){
        t2a <- crossprod(t(dG_dlambda[[k]]),
                         A_AGAInv_AGXy[(k-1)*nc + 1:nc,,drop=FALSE])
        idx <- (k-1)*nc + 1:nc
        t2b <- crossprod(t(G[[k]]), A_term2b_mat[idx,, drop=FALSE])
        t2c <- crossprod(t(G[[k]]), A_term2c_mat[idx,, drop=FALSE])
        list(t2a,t2b,t2c)
      })
      term2a <- lapply(terms, `[[`, 1)
      term2b <- lapply(terms, `[[`, 2)
      term2c <- lapply(terms, `[[`, 3)
    }

    term2a <- Reduce("c", term2a)
    term2b <- Reduce("c", term2b)
    term2c <- Reduce("c", term2c)

    ## Final result: unconstrained derivative minus three correction terms
    return(cbind(unlist(dGXy) - (unlist(term2a) +
                                   unlist(term2b) +
                                   unlist(term2c))))
  } else {

    ## ## Small problem path: least-squares formulation via G^{1/2}

    ## The constrained estimate can be written as
    #    beta_tilde = G^{1/2} * residuals(G^{1/2}X'y ~ G^{1/2}A)
    #  because projecting out the column space of G^{1/2}A from
    #  G^{1/2}X'y and then multiplying by G^{1/2} is equivalent to
    #  applying the constraint projection U. Differentiating this
    #  representation with respect to lambda using the product rule on
    #  G^{1/2} and dG^{1/2}/dl yields two residual vectors (resids1
    #  from the derivative of the basis, resids2 from the derivative
    #  of the right-hand side) that combine into the final derivative.

    ## G^{1/2} A and its derivative (block-diagonal times A)
    GhalfA <- Reduce("rbind", GAmult_wrapper(Ghalf, A, K, nc, nca,
                                             parallel, cl, chunk_size,
                                             num_chunks, rem_chunks))

    dGhalfA <- Reduce("rbind", GAmult_wrapper(dGhalf, A, K, nc, nca,
                                              parallel, cl, chunk_size,
                                              num_chunks, rem_chunks))

    ## (dG^{1/2}/dl) X'y (block-diagonal multiply)
    dGhalfXy <- cbind(unlist(
      matmult_block_diagonal(dGhalf, Xy, K, parallel, cl,
                             chunk_size, num_chunks, rem_chunks)))

    ## Stabilize the least squares solve: rescale by the average
    # magnitude of G^{1/2}A to prevent numerical issues when the
    # constraint matrix is poorly conditioned
    comp_stab_sc <- sqrt(mean(abs(GhalfA)))
    comp_stab_sc <- 1*(comp_stab_sc == 0) + comp_stab_sc

    ## resids1: residuals from regressing (dG^{1/2}/dl)X'y onto
    #           (dG^{1/2}/dl)A -- captures derivative of the basis
    resids1 <- do.call('.lm.fit', list(x = GhalfA / comp_stab_sc,
                                       y = GhalfXy_temp / comp_stab_sc)
    )$residuals * comp_stab_sc

    ## resids2: residuals from regressing (dG^{1/2}/dl)X'y onto
    #           (dG^{1/2}/dl)A -- captures derivative of the RHS
    resids2 <- do.call('.lm.fit', list(x = dGhalfA / comp_stab_sc,
                                       y = dGhalfXy / comp_stab_sc)
    )$residuals * comp_stab_sc

    ## Combine: d(beta_tilde)/dl = (dG^{1/2}/dl)*resids2 + G^{1/2}*resids1
    # applied partition-wise since G^{1/2} and dG^{1/2}/dl are block-diagonal
    if(parallel & !is.null(cl)) {
      if(rem_chunks > 0) {
        rem_indices <- num_chunks * chunk_size + 1:rem_chunks
        b2_rem <- lapply(rem_indices, function(k) {
          crossprod(t(dGhalf[[k]]), cbind(resids2[(k-1)*nc + 1:nc])) +
            crossprod(t(Ghalf[[k]]), cbind(resids1[(k-1)*nc + 1:nc]))
        })
      } else {
        b2_rem <- list()
      }

      b2_result <- c(
        Reduce("c",parallel::parLapply(cl, 1:num_chunks, function(chunk) {
          inds <- (chunk - 1) * chunk_size + 1:chunk_size
          lapply(inds, function(k) {
            crossprod(t(dGhalf[[k]]), cbind(resids2[(k-1)*nc + 1:nc])) +
              crossprod(t(Ghalf[[k]]), cbind(resids1[(k-1)*nc + 1:nc]))
          })
        })),
        b2_rem
      )
      b2 <- unlist(b2_result)
    } else {
      b2 <- Reduce("c", lapply(1:(K+1), function(k) {
        crossprod(t(dGhalf[[k]]), cbind(resids2[(k-1)*nc + 1:nc])) +
          crossprod(t(Ghalf[[k]]), cbind(resids1[(k-1)*nc + 1:nc]))
      }))
    }

    return(cbind(b2))
  }
}


#' Compute Component \eqn{\textbf{G}^{1/2}\textbf{A}(\textbf{A}^{T}\textbf{G}\textbf{A})^{-1}\textbf{A}^{T}\textbf{G}\textbf{X}^{T}\textbf{y}}
#'
#' @param G List of \eqn{\textbf{G}} matrices
#' @param Ghalf List of \eqn{\textbf{G}^{1/2}} matrices
#' @param A Constraint matrix \eqn{\textbf{A}}
#' @param AGAInv Inverse of \eqn{\textbf{A}^{T}\textbf{G}\textbf{A}}
#' @param Xy List of \eqn{\textbf{X}^{T}\textbf{y}} products
#' @param nc Number of columns
#' @param K Number of partitions minus 1 (\eqn{K})
#' @param parallel Use parallel processing
#' @param cl Cluster object
#' @param chunk_size Size of parallel chunks
#' @param num_chunks Number of chunks
#' @param rem_chunks Remaining chunks
#'
#' @details
#' Computes \eqn{\textbf{G}^{1/2}\textbf{A}(\textbf{A}^{T}\textbf{G}\textbf{A})^{-1}\textbf{A}^{T}\textbf{G}\textbf{X}^{T}\textbf{y}} efficiently in parallel chunks.
#' Note: The description seems slightly off from the C++ helper functions called (e.g., `compute_AGXy`, `compute_result_blocks`). This computes a component needed for least-squares transformation involving \eqn{\textbf{G}^{1/2}}.
#' Returns both the result and intermediate \eqn{\textbf{A}^{T}\textbf{G}\textbf{A})^{-1}\textbf{A}^{T}\textbf{G}\textbf{X}^{T}\textbf{y}} product for reuse.
#'
#' @return List containing:
#' \itemize{
#'   \item Result vector
#'   \item AGAInvAGXy intermediate product \eqn{(\textbf{A}^{T}\textbf{G}\textbf{A})^{-1}\textbf{A}^{T}\textbf{G}\textbf{X}^{T}\textbf{y}}
#' }
#'
#' @keywords internal
#' @export
compute_G_eigen <- function(X_gram, Lambda, K, parallel, cl,
                            chunk_size, num_chunks, rem_chunks,
                            family, unique_penalty_per_partition,
                            L_partition_list, keep_G = TRUE,
                            schur_corrections) {
  ## [Changer 2026-02-12] replaced %**% with tcrossprod(V * d, V)
  ## helper to compute eigen-based matrix powers efficiently
  compute_from_eigen <- function(eig, keep_G, need_GhalfInv) {
    eigen_values <- eig$values
    eigen_values[eig$values <= 0] <- 1
    inv_eigen_values <- 1/eigen_values
    inv_eigen_values[eig$values <= 0] <- 0
    sqrt_inv <- sqrt(inv_eigen_values)

    Ghalf <- tcrossprod(eig$vectors * rep(sqrt_inv, each = nrow(eig$vectors)),
                        eig$vectors)
    G_out <- NULL
    if(keep_G){
      G_out <- tcrossprod(eig$vectors * rep(inv_eigen_values, each = nrow(eig$vectors)),
                          eig$vectors)
    }
    if(need_GhalfInv){
      GhalfInv <- tcrossprod(eig$vectors / rep(sqrt_inv, each = nrow(eig$vectors)),
                             eig$vectors)
      return(list(G = G_out, Ghalf = Ghalf, GhalfInv = GhalfInv))
    } else {
      return(list(G = G_out, Ghalf = Ghalf))
    }
  }

  need_GhalfInv <- (paste0(family)[2] != 'identity') |
    (paste0(family)[1] != 'gaussian')

  if(parallel & !is.null(cl)) {
    if(rem_chunks > 0) {
      rem_indices <- num_chunks * chunk_size + 1:rem_chunks
      rem <- lapply(rem_indices, function(k) {
        if(unique_penalty_per_partition){
          eig <- tryCatch(eigen(X_gram[[k]] + Lambda +
                                  L_partition_list[[k]] + schur_corrections[[k]],
                                symmetric = TRUE), error = function(e) NULL)
        } else {
          eig <- tryCatch(eigen(X_gram[[k]] + Lambda +
                                  schur_corrections[[k]], symmetric = TRUE),
                          error = function(e) NULL)
        }
        if(is.null(eig)) return(list(G = NULL, Ghalf = NULL))
        compute_from_eigen(eig, keep_G, need_GhalfInv)
      })
    } else {
      rem <- list()
    }

    result <- c(
      Reduce("c", parallel::parLapply(cl, 1:num_chunks, function(chunk) {
        inds <- (chunk - 1)*chunk_size + 1:chunk_size
        lapply(inds, function(k) {
          if(unique_penalty_per_partition){
            eig <- tryCatch(eigen(X_gram[[k]] + Lambda +
                                    L_partition_list[[k]] + schur_corrections[[k]],
                                  symmetric = TRUE), error = function(e) NULL)
          } else {
            eig <- tryCatch(eigen(X_gram[[k]] + Lambda +
                                    schur_corrections[[k]], symmetric = TRUE),
                            error = function(e) NULL)
          }
          if(is.null(eig)) return(list(G = NULL, Ghalf = NULL))
          compute_from_eigen(eig, keep_G, need_GhalfInv)
        })
      })),
      rem
    )
  } else {
    result <- lapply(1:(K+1), function(k) {
      if(unique_penalty_per_partition){
        eig <- tryCatch(eigen(X_gram[[k]] + Lambda +
                                L_partition_list[[k]] + schur_corrections[[k]],
                              symmetric = TRUE), error = function(e) NULL)
      } else {
        eig <- tryCatch(eigen(X_gram[[k]] + Lambda +
                                schur_corrections[[k]], symmetric = TRUE),
                        error = function(e) NULL)
      }
      if(is.null(eig)) return(list(G = NULL, Ghalf = NULL))
      compute_from_eigen(eig, keep_G, need_GhalfInv)
    })
  }

  ## Reorganize results by matrix type
  if(need_GhalfInv){
    result_processed <- list(
      G = lapply(result, `[[`, "G"),
      Ghalf = lapply(result, `[[`, "Ghalf"),
      GhalfInv = lapply(result, `[[`, "GhalfInv")
    )
  } else {
    result_processed <- list(
      G = lapply(result, `[[`, "G"),
      Ghalf = lapply(result, `[[`, "Ghalf")
    )
  }
  return(result_processed)
}


#' Compute Component \eqn{\textbf{G}^{1/2}\textbf{A}(\textbf{A}^{T}\textbf{G}\textbf{A})^{-1}\textbf{A}^{T}\textbf{G}\textbf{X}^{T}\textbf{y}}
#'
#' @param G List of \eqn{\textbf{G}} matrices
#' @param Ghalf List of \eqn{\textbf{G}^{1/2}} matrices
#' @param A Constraint matrix \eqn{\textbf{A}}
#' @param AGAInv Inverse of \eqn{\textbf{A}^{T}\textbf{G}\textbf{A}}
#' @param Xy List of \eqn{\textbf{X}^{T}\textbf{y}} products
#' @param nc Number of columns
#' @param K Number of partitions minus 1 (\eqn{K})
#' @param parallel Use parallel processing
#' @param cl Cluster object
#' @param chunk_size Size of parallel chunks
#' @param num_chunks Number of chunks
#' @param rem_chunks Remaining chunks
#'
#' @details
#' Computes \eqn{\textbf{G}^{1/2}\textbf{A}(\textbf{A}^{T}\textbf{G}\textbf{A})^{-1}\textbf{A}^{T}\textbf{G}\textbf{X}^{T}\textbf{y}} efficiently in parallel chunks.
#' Note: The description seems slightly off from the C++ helper functions called (e.g., `compute_AGXy`, `compute_result_blocks`). This computes a component needed for least-squares transformation involving \eqn{\textbf{G}^{1/2}}.
#' Returns both the result and intermediate \eqn{\textbf{A}^{T}\textbf{G}\textbf{A})^{-1}\textbf{A}^{T}\textbf{G}\textbf{X}^{T}\textbf{y}} product for reuse.
#'
#' @return List containing:
#' \itemize{
#'   \item Result vector
#'   \item AGAInvAGXy intermediate product \eqn{(\textbf{A}^{T}\textbf{G}\textbf{A})^{-1}\textbf{A}^{T}\textbf{G}\textbf{X}^{T}\textbf{y}}
#' }
#'
#' @keywords internal
#' @export
compute_GhalfXy_temp_wrapper <- function(G,
                                         Ghalf,
                                         A,
                                         AGAInv,
                                         Xy,
                                         nc,
                                         K,
                                         parallel,
                                         cl,
                                         chunk_size,
                                         num_chunks,
                                         rem_chunks) {
  if(parallel & !is.null(cl)) {
    ## First compute AGXy in parallel chunks
    if(rem_chunks > 0) {
      rem_start <- num_chunks * chunk_size + 1
      rem_end <- rem_start + rem_chunks - 1
      AGXy_rem <- compute_AGXy(G[rem_start:rem_end],
                               A,
                               Xy[rem_start:rem_end],
                               nc, K, rem_start-1, rem_end-1)
    } else {
      AGXy_rem <- matrix(0, ncol(A), 1)
    }

    ## Process main chunks for AGXy
    AGXy_chunks <- parallel::parLapply(cl, 1:num_chunks, function(chunk) {
      chunk_start <- (chunk - 1) * chunk_size + 1
      chunk_end <- min(chunk * chunk_size, K+1)
      compute_AGXy(G[chunk_start:chunk_end],
                   A,
                   Xy[chunk_start:chunk_end],
                   nc, K, chunk_start-1, chunk_end-1)
    })

    ## Sum up AGXy results
    AGXy <- Reduce('+', c(AGXy_chunks, list(AGXy_rem)))

    ## Compute AAGAInvAGXy once
    AGAInvAGXy <- AGAInv %**% AGXy
    AAGAInvAGXy <- A %**% AGAInvAGXy

    ## Now compute final blocks in parallel
    if(rem_chunks > 0) {
      rem_start <- num_chunks * chunk_size + 1
      rem_end <- rem_start + rem_chunks - 1
      result_rem <- compute_result_blocks(G[rem_start:rem_end],
                                          Ghalf[rem_start:rem_end],
                                          A,
                                          AAGAInvAGXy,
                                          nc,
                                          rem_start-1,
                                          rem_end-1)
    } else {
      result_rem <- numeric(0)
    }

    ## Process main chunks for final result
    result_chunks <- parallel::parLapply(cl, 1:num_chunks, function(chunk) {
      chunk_start <- (chunk - 1) * chunk_size + 1
      chunk_end <- min(chunk * chunk_size, K+1)
      compute_result_blocks(G[chunk_start:chunk_end],
                            Ghalf[chunk_start:chunk_end],
                            A,
                            AAGAInvAGXy,
                            nc,
                            chunk_start-1,
                            chunk_end-1)
    })

    ## Combine all results
    return(list(
      c(unlist(result_chunks), result_rem),
      AGAInvAGXy
    ))

  } else {

    ## Compute AAGAInvAGXy once
    AGAInvAGXy <- AGAInv %**% compute_AGXy(G,
                                           A,
                                           Xy,
                                           nc,
                                           K,
                                           0,
                                           K)

    ## Now compute final blocks in parallel
    result <- compute_result_blocks(G,
                                    Ghalf,
                                    A,
                                    A %**% AGAInvAGXy,
                                    nc,
                                    0,
                                    K)

    ## Combine all results
    return(list(
      c(unlist(result)),
      AGAInvAGXy
    ))
  }
}

## [Change 2026-02-15] Replace softplus with exponential transforms,
#  accept initial values on the raw rather than inverse-softplus scale

#' Tune Smoothing and Ridge Penalties via Generalized Cross Validation
#'
#' @description
#' Optimizes smoothing spline and ridge regression penalties by minimizing GCV criterion.
#' Uses BFGS optimization with analytical gradients or finite differences.
#'
#' @param y List; response vectors by partition
#' @param X List; design matrices by partition
#' @param X_gram List; Gram matrices by partition
#' @param smoothing_spline_penalty Matrix; integrated squared second derivative penalty
#' @param A Matrix; smoothness constraints at knots
#' @param K Integer; number of interior knots in 1-D, number of partitions - 1 in higher dimensions
#' @param nc Integer; columns per partition
#' @param nr Integer; total sample size
#' @param opt Logical; TRUE to optimize penalties, FALSE to use initial values
#' @param use_custom_bfgs Logical; TRUE for analytic gradient BFGS as natively implemented, FALSE for finite differences as implemented by \code{stats::optim()}.
#' @param family GLM family with optional custom tuning loss
#' @param wiggle_penalty,flat_ridge_penalty Fixed penalty values if provided
#' @param initial_wiggle,initial_flat Numeric vectors; candidate values for grid
#'   search initialization on the raw (non-negative) scale. Converted to log
#'   scale internally for optimization.
#' @param unique_penalty_per_predictor,unique_penalty_per_partition Logical; allow predictor/partition-specific penalties
#' @param penalty_vec Numeric vector; initial values for predictor/partition
#'   penalties on the raw (non-negative) scale. Converted to log scale
#'   internally for optimization. Use \code{c()} when no per-predictor or
#'   per-partition penalties are needed.
#' @param meta_penalty The "meta" ridge penalty, a regularization for predictor/partition penalties to pull them towards 1 on the raw scale
#' @param keep_weighted_Lambda,iterate Logical controlling GLM fitting
#' @param qp_score_function,quadprog,qp_Amat,qp_bvec,qp_meq Quadratic programming parameters (see arguments of \code{\link[lgspline]{lgspline}})
#' @param tol Numeric; convergence tolerance
#' @param sd_y,delta Response standardization parameters
#' @param constraint_value_vectors List; constraint values
#' @param parallel Logical; enable parallel computation
#' @param cl,chunk_size,num_chunks,rem_chunks Parallel computation parameters
#' @param custom_penalty_mat Optional custom penalty matrix
#' @param order_list List; observation ordering by partition
#' @param glm_weight_function,schur_correction_function Functions for GLM weights and corrections
#' @param need_dispersion_for_estimation,dispersion_function Control dispersion estimation
#' @param observation_weights Optional observation weights
#' @param homogenous_weights Logical; TRUE if all weights equal
#' @param blockfit Logical; when TRUE, block-fitting (not per-partition fitting) approach is used, analogous to quadratic programming.
#' @param just_linear_without_interactions Numeric; vector of columns of input predictor matrix that correspond to non-spline effects without interactions, used for block-fitting.
#' @param Vhalf,VhalfInv, Square root and inverse square root correlation structures for fitting GEEs.
#' @param verbose Logical; print progress
#' @param include_warnings Logical; print warnings/try-errors
#' @param ... Additional arguments passed to fitting functions
#'
#' @return List containing:
#' \itemize{
#'   \item Lambda - Final combined penalty matrix
#'   \item flat_ridge_penalty - Optimized ridge penalty
#'   \item wiggle_penalty - Optimized smoothing penalty
#'   \item other_penalties - Optimized predictor/partition penalties
#'   \item L_predictor_list - Predictor-specific penalty matrices
#'   \item L_partition_list - Partition-specific penalty matrices
#' }
#'
#' @details
#' Uses BFGS optimization to minimize GCV criterion for penalty selection.
#' Supports analytical gradients for efficiency with standard GLM families.
#' Can optimize unique penalties per predictor/partition.
#' Handles custom loss functions and GLM weights.
#' Parallel computation available for large problems.
#'
#' Parameterization: initial penalty values are accepted on the raw
#' (non-negative) scale and converted to natural log-scale internally,
#' i.e. raw_penalty = exp(theta), so that raw penalties are always positive.
#' The chain rule factor d(exp(theta))/d(theta) = exp(theta) = raw_penalty
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[stats]{optim}} for Hessian-free optimization
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
    nc,
    nr,
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
    qp_score_function, quadprog,  qp_Amat, qp_bvec, qp_meq,
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
){

  if(verbose){
    cat('    Starting tuning\n')
  }

  ## Convert raw-scale inputs to log scale for internal optimization
  log_initial_wiggle <- log(initial_wiggle)
  log_initial_flat <- log(initial_flat)
  if(length(penalty_vec) > 0){
    log_penalty_vec <- log(penalty_vec)
  } else {
    log_penalty_vec <- c()
  }

  ## No getting G estimates here
  return_G_getB <- FALSE

  if(verbose){
    cat('    Xy\n')
  }

  ## Computational efficiency
  snr <- sqrt(nr)
  Xy <- vectorproduct_block_diagonal(X, y, K)

  ## X'y sufficient statistic
  Xyr <- Reduce("rbind",Xy)
  nca <- ncol(A)
  snr <- sqrt(nr)
  Xt <- lapply(X, t)
  unl_y <- unlist(y)

  if(verbose){
    cat('    Getting puedocount delta\n')
  }
  ## pseudocount
  if(is.null(delta) &
     any(is.null(family$custom_dev.resids)) &
     opt){
    ## ignore pseudocount for y if not needed for link-fxn transform
    if((paste0(family)[2] == 'identity') |
       (paste0(family)[2] == 'log' & (min(unl_y) > 0 )) |
       ((any(paste0(family)[2] %in% c('inverse','1/mu^2'))) &
        (!(any(unl_y == 0)))) |
       (paste0(family)[2] == 'logit' & (!any(unl_y %in% c(0,1))))
    ){
      delta <- 0
      ## obtain pseudocount for y if needed for link-fxn transform
    } else {
      ## empirical, find delta such that t-quantiles are most closely matched
      t_quants <- qt((seq(0, 1, len = nr+2))[-c(1, nr+2)], df = nr - 1)
      delta <- stats::optim(
        1/16,
        ## minimize difference between normal quantile and observed quantiles
        fn = function(par){
          y_delta <- std(sort(family$linkfun((unl_y+par)/(1+2*par))))
          if(length(unlist(observation_weights)) == length(unl_y)){
            mean(abs(y_delta - t_quants)*unlist(observation_weights))
          } else {
            mean(abs(y_delta - t_quants))
          }
        },
        method = 'Brent',
        lower = 1e-64,
        upper = 1
      )$par
    }
  }

  if(verbose){
    cat('    GCV, gradient, and BFGS Function prep\n')
  }

  ## intialize gradient function
  gr_fxn <- NULL

  ## for glm/quadprog problems
  prevB <- lapply(1:(K+1),function(k)0)

  ## ensuring compatibility with no A
  if(any(is.null(A))){
    ## for compatibility, albeit inefficient
    A <- cbind(rep(0, (K+1)*nc))
    A <- cbind(A, A)
    nca <- 2
  }

  ## Compute GCV_u and matrix components
  gcvu_fxn <- function(par){
    if(verbose){
      cat('        gcvu_fxn start\n')
    }
    ## [Change: softplus -> exp]
    #  2 penalties of interest
    wiggle_penalty <- exp(par[1])       # smoothing spline f''(x)^2
    flat_ridge_penalty <- exp(par[2])   # flat ridge regression penalty
    if(unique_penalty_per_predictor | unique_penalty_per_partition){
      penalty_vec <- exp(c(par[-c(1:2)]))
    } else {
      penalty_vec <- c()
    }

    ## Reparameterize
    lambda_1 <- wiggle_penalty
    lambda_2 <- flat_ridge_penalty

    ## Compute penalty matrix Lambda and components
    if(verbose){
      cat('        compute_Lambda\n')
    }
    Lambda_list <- compute_Lambda(custom_penalty_mat,
                                  smoothing_spline_penalty,
                                  wiggle_penalty,
                                  flat_ridge_penalty,
                                  K,
                                  nc,
                                  unique_penalty_per_predictor,
                                  unique_penalty_per_partition,
                                  penalty_vec,
                                  colnm_expansions,
                                  just_Lambda = FALSE)
    Lambda <- Lambda_list[[1]]
    L1 <- Lambda_list[[2]]
    L2 <- Lambda_list[[3]]

    ## Compute
    if(verbose){
      cat('        compute_G_eigen\n')
    }

    ## Schur complements are not needed here, only in get_B
    # This is for initialization only
    schur_corrections <- lapply(1:(K+1), function(k)0)
    G_list <- compute_G_eigen(X_gram,
                              Lambda,
                              K,
                              parallel & parallel_eigen,
                              cl,
                              chunk_size,
                              num_chunks,
                              rem_chunks,
                              family,
                              unique_penalty_per_partition,
                              Lambda_list$L_partition_list,
                              keep_G = TRUE,
                              schur_corrections)

    if(verbose){
      cat('        gcvu_fxn get_B\n')
    }

    ## For getting updates of coefficient estimates and correlation matrix
    return_G_getB <- TRUE
    B_list <- get_B(
      X,
      X_gram,
      Lambda,
      keep_weighted_Lambda,
      unique_penalty_per_partition,
      Lambda_list$L_partition_list,
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
      iterate,
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
      constraint_value_vectors,
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
      ...)
    G_list <- B_list$G_list
    B <- B_list$B

    if(verbose){
      cat('        gcvu_fxn AGAmult_wrapper\n')
    }
    AGAInv <- invert(AGAmult_wrapper(G_list$G,
                                     A,
                                     K,
                                     nc,
                                     nca,
                                     parallel & parallel_aga,
                                     cl,
                                     chunk_size,
                                     num_chunks,
                                     rem_chunks) +
                       1e-16*diag(ncol(A)))

    if(verbose){
      cat('        gcvu_fxn matmult_block_diagonal for GXX\n')
    }
    GXX <- matmult_block_diagonal(G_list$G,
                                  X_gram,
                                  K,
                                  parallel & parallel_matmult,
                                  cl,
                                  chunk_size,
                                  num_chunks,
                                  rem_chunks)

    if(verbose){
      cat('        gcvu_fxn compute_trace_UGXX_wrapper\n')
    }
    ## [Change 2026-02-15] Replace with more stable compute trace H
    sum_W <-
      compute_trace_H(G_list$G,
                      Lambda,
                      A,
                      AGAInv,
                      nc,
                      nca,
                      K,
                      parallel,
                      cl,
                      chunk_size,
                      num_chunks,
                      rem_chunks,
                      unique_penalty_per_partition,
                      Lambda_list$L_partition_list)

    if(verbose)cat('        gcvu_fxn get predictions\n')
    preds <- matmult_block_diagonal(X,
                                    B,
                                    K,
                                    parallel & parallel_matmult,
                                    cl,
                                    chunk_size,
                                    num_chunks,
                                    rem_chunks)

    if(verbose){
      cat('        gcvu_fxn custom or default residuals\n')
    }
    ## Residuals for GCV
    if(paste0(family)[2] == 'identity' |
       is.null(family$custom_dev.resids)){

      ## If not canonical Gaussian and weights are present, use them
      # recall, for Gaussian, we transformed X and y to incorporate weights
      # prior to inclusion here
      if(any(!is.null(observation_weights[[1]])) &
         (paste0(family)[2] != 'identity' | paste0(family)[1] != 'gaussian')){
        residuals <- lapply(1:(K+1),function(k){
          (family$linkfun((y[[k]]+delta)/(1+2*delta)) -
             (preds[[k]]+delta)/(1+2*delta))*
            c(observation_weights[[k]])
        })
      } else{
        residuals <- lapply(1:(K+1),function(k){
          family$linkfun((y[[k]]+delta)/(1+2*delta)) -
            (preds[[k]]+delta)/(1+2*delta)
        })
      }
    } else {
      ## Can replace regular y-mu residuals with deviance/obj. function, etc.
      residuals <- lapply(1:(K+1),function(k){
        family$custom_dev.resids(y[[k]],
                                 family$linkinv(c(preds[[k]])),
                                 order_list[[k]],
                                 family,
                                 observation_weights[[k]],
                                 ...)
      })
    }

    ## Compute GCV_u components
    if(verbose){
      cat('        gcvu_fxn GCVu operations\n')
    }
    numerator <- sum(unlist(residuals)^2)
    mean_W <- sum_W / nr
    denominator <- nr * (1 - mean_W)^2
    denom_sq <- denominator^2
    GCV_u <- numerator / denominator

    ## Regularization penalty, pulls penalty terms towards 1
    if(verbose){
      cat('        gcvu_fxn penalization operations\n')
    }
    if(unique_penalty_per_partition | unique_penalty_per_predictor){
      meta_penalty <- 0.5*meta_penalty*sum((penalty_vec - 1)^2) +
        0.5*1e-32*((wiggle_penalty - 1))^2
    } else {
      meta_penalty <- 0.5*1e-32*((wiggle_penalty - 1))^2
    }

    if(verbose)cat('        done GCVu,', GCV_u, '\n')
    ## Output list, prevent from computing twice
    return(list(GCV_u = GCV_u + meta_penalty,
                B = B,
                GXX = GXX,
                G_list = G_list,
                mean_W = mean_W,
                sum_W = sum_W,
                Lambda = Lambda,
                L1 = L1,
                L2 = L2,
                L_predictor_list = Lambda_list$L_predictor_list,
                L_partition_list = Lambda_list$L_partition_list,
                numerator = numerator,
                denominator = denominator,
                residuals = residuals,
                denom_sq = denom_sq,
                AGAInv = AGAInv
    ))
  }

  ## Gradient function for BFGS
  gr_fxn <- function(par, outlist = NULL){
    if(verbose){
      cat('        gr_fxn start\n')
    }
    ## [Change: softplus -> exp]
    # 2 penalties of interest
    wiggle_penalty <- exp(par[1])       # smoothing spline f''(x)^2
    flat_ridge_penalty <- exp(par[2])   # flat ridge regression penalty
    if(unique_penalty_per_predictor | unique_penalty_per_partition){
      penalty_vec <- exp(c(par[-c(1:2)]))
    } else {
      penalty_vec <- c()
    }

    ## Reparameterize
    lambda_1 <- wiggle_penalty
    lambda_2 <- flat_ridge_penalty

    if(verbose){
      cat('        lambda_1, lambda_2: ', lambda_1, ', ', lambda_2, '\n')
    }
    if(any(is.null(outlist))){

      ## Compute necessary matrices, if not supplied
      if(verbose){
        cat('        Lambda list\n')
      }
      Lambda_list <- compute_Lambda(custom_penalty_mat,
                                    smoothing_spline_penalty,
                                    wiggle_penalty,
                                    flat_ridge_penalty,
                                    K,
                                    nc,
                                    unique_penalty_per_predictor,
                                    unique_penalty_per_partition,
                                    penalty_vec,
                                    colnm_expansions,
                                    just_Lambda = FALSE)
      Lambda <- Lambda_list[[1]]
      L1 <- Lambda_list[[2]]
      L2 <- Lambda_list[[3]]

      if(verbose){
        cat('        G list\n')
      }

      ## No need for corrections yet, this is just initialization
      schur_corrections <- lapply(1:(K+1), function(k)0)
      G_list <- compute_G_eigen(X_gram,
                                Lambda,
                                K,
                                parallel & parallel_eigen,
                                cl,
                                chunk_size,
                                num_chunks,
                                rem_chunks,
                                family,
                                unique_penalty_per_partition,
                                Lambda_list$L_partition_list,
                                keep_G = TRUE,
                                schur_corrections)

      if(verbose){
        cat('        gr fxn coefficients\n')
      }

      ## For getting updates of coefficient estimates and correlation matrix
      return_G_getB <- TRUE
      B_list <- get_B(
        X,
        X_gram,
        Lambda,
        keep_weighted_Lambda,
        unique_penalty_per_partition,
        Lambda_list$L_partition_list,
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
        iterate,
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
        constraint_value_vectors,
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
        ...)
      G_list <- B_list$G_list
      B <- B_list$B

      if(verbose){
        cat('        AGAmult_wrapper\n')
      }
      AGAInv <- invert(AGAmult_wrapper(G_list$G,
                                       A,
                                       K,
                                       nc,
                                       nca,
                                       parallel & parallel_aga,
                                       cl,
                                       chunk_size,
                                       num_chunks,
                                       rem_chunks) +
                         1e-16*diag(ncol(A)))


      if(verbose){
        cat('        GXX matmult_block_diagonal\n')
      }
      GXX <- matmult_block_diagonal(G_list$G,
                                    X_gram,
                                    K,
                                    parallel & parallel_matmult,
                                    cl,
                                    chunk_size,
                                    num_chunks,
                                    rem_chunks)

      if(verbose){
        cat('        sum_W compute_trace_UGXX_wrapper\n')
      }

      ## [Change 2026-02-15] Use more stable version of trace computation
      sum_W <-
        compute_trace_H(G_list$G,
                        Lambda,
                        A,
                        AGAInv,
                        nc,
                        nca,
                        K,
                        parallel,
                        cl,
                        chunk_size,
                        num_chunks,
                        rem_chunks,
                        unique_penalty_per_partition,
                        Lambda_list$L_partition_list)

      if(verbose){
        cat('        gr fxn preds\n')
      }
      preds <- matmult_block_diagonal(X,
                                      B,
                                      K,
                                      parallel & parallel_matmult,
                                      cl,
                                      chunk_size,
                                      num_chunks,
                                      rem_chunks)

      ## Residuals for GCV
      if(verbose){
        cat('        gr fxn residuals\n')
      }
      if(paste0(family)[2] == 'identity' |
         is.null(family$custom_dev.resids)){

        ## If not canonical Gaussian and weights are present, use them
        # recall, for Gaussian, we transformed X and y to incorporate weights
        # prior to inclusion here
        if(any(!is.null(observation_weights[[1]])) &
           (paste0(family)[2] != 'identity' | paste0(family)[1] != 'gaussian')){
          residuals <- lapply(1:(K+1),function(k){
            (family$linkfun((y[[k]]+delta)/(1+2*delta)) -
               (preds[[k]]+delta)/(1+2*delta))*
              c(observation_weights[[k]])
          })
        } else{
          residuals <- lapply(1:(K+1),function(k){
            family$linkfun((y[[k]]+delta)/(1+2*delta)) -
              (preds[[k]]+delta)/(1+2*delta)
          })
        }

      } else {
        residuals <- lapply(1:(K+1),function(k){
          family$custom_dev.resids(y[[k]],
                                   family$linkinv(c(preds[[k]])),
                                   order_list[[k]],
                                   family,
                                   observation_weights[[k]],
                                   ...)
        })
      }

      ## Compute GCV_u components
      if(verbose){
        cat('        gr fxn compute GCV_u\n')
      }
      numerator <- sum(unlist(residuals)^2)
      mean_W <- sum_W / nr
      denominator <- nr * (1 - mean_W)^2
      denom_sq <- denominator^2
      GCV_u <- numerator / denominator

      ## Regularization penalty - pull penalties to 1
      if(unique_penalty_per_partition | unique_penalty_per_predictor){
        meta_penalty <- 0.5*meta_penalty*sum((penalty_vec-1)^2) +
          0.5*1e-32*((wiggle_penalty - 1))^2
      } else {
        meta_penalty <- 0.5*1e-32*((wiggle_penalty - 1))^2
      }

      if(verbose){
        cat('        gr fxn outlist\n')
      }
      outlist <- list(GCV_u = GCV_u + meta_penalty,
                      B = B,
                      GXX = GXX,
                      G_list = G_list,
                      mean_W = mean_W,
                      sum_W = sum_W,
                      Lambda = Lambda,
                      L1 = L1,
                      L2 = L2,
                      L_predictor_list = Lambda_list$L_predictor_list,
                      L_partition_list = Lambda_list$L_partition_list,
                      numerator = numerator,
                      denominator = denominator,
                      residuals = residuals,
                      denom_sq = denom_sq,
                      AGAInv = AGAInv
      )
    }

    ## An important but computationally expensive-to-compute vector for
    # computing derivatives
    if(verbose){
      cat('        GhalfXy_temp_list \n')
    }
    GhalfXy_temp_list <- compute_GhalfXy_temp_wrapper(outlist$G_list$G,
                                                      outlist$G_list$Ghalf,
                                                      A,
                                                      outlist$AGAInv,
                                                      Xy,
                                                      nc,
                                                      K,
                                                      parallel &
                                                        parallel_aga,
                                                      cl,
                                                      chunk_size,
                                                      num_chunks,
                                                      rem_chunks)
    GhalfXy_temp <- GhalfXy_temp_list[[1]]
    AGAInvAGXy <- GhalfXy_temp_list[[2]]


    ## Compute derivatives for lambda
    if(verbose){
      cat('        compute_dG_dlambda \n')
    }
    dG_dlambda <- compute_dG_dlambda(outlist$G_list$G,
                                     outlist$Lambda,
                                     K,
                                     lambda_1,
                                     unique_penalty_per_partition,
                                     outlist$L_partition_list,
                                     parallel & parallel_matmult,
                                     cl,
                                     chunk_size,
                                     num_chunks,
                                     rem_chunks)

    if(verbose){
      cat('        Compute dGhalf \n')
    }
    dGhalf <- compute_dGhalf(dG_dlambda,
                             nc,
                             K,
                             parallel & parallel_eigen,
                             cl,
                             chunk_size,
                             num_chunks,
                             rem_chunks)

    if(verbose){
      cat('        compute_dG_u_dlambda_xy \n')
    }
    dG_u_dlambda1_Xyr <- compute_dG_u_dlambda_xy(
      AGAInvAGXy,
      outlist$AGAInv,
      outlist$G_list$G,
      A,
      dG_dlambda,
      nc,
      nca,
      K,
      Xy,
      outlist$G_list$Ghalf,
      dGhalf,
      GhalfXy_temp,
      parallel & parallel_matmult,
      cl,
      chunk_size,
      num_chunks,
      rem_chunks)

    if(verbose){
      cat('        compute_dW_dlambda_wrapper \n')
    }
    dW_dlambda <- compute_dW_dlambda_wrapper(outlist$G_list$G,
                                             A,
                                             outlist$GXX,
                                             outlist$G_list$Ghalf,
                                             dG_dlambda,
                                             dGhalf,
                                             outlist$AGAInv,
                                             nc,
                                             K,
                                             parallel & parallel_matmult,
                                             cl,
                                             chunk_size,
                                             num_chunks,
                                             rem_chunks)

    ## Compute other components)
    if(verbose){
      cat('        neg2tresidX\n')
    }
    neg2tresidX <- Reduce('cbind',
                          matmult_block_diagonal(
                            lapply(
                              outlist$residuals,
                              function(r)-2*t(r)),
                            X,
                            K,
                            parallel & parallel_matmult,
                            cl,
                            chunk_size,
                            num_chunks,
                            rem_chunks))

    dnumerator_dlambda1 <- c(neg2tresidX %**%
                               dG_u_dlambda1_Xyr)
    ddenominator_dlambda1 <- 2 * (1 - outlist$mean_W) * -dW_dlambda
    dGCV_u_dlambda1 <- (dnumerator_dlambda1 *
                          outlist$denominator -
                          outlist$numerator *
                          ddenominator_dlambda1) /
      outlist$denom_sq

    ## Trace-ratio approximation for computing gradients of other elements
    #  of the penalty. This heuristic operates on the raw penalty scale and
    #  is invariant to the parameterization (softplus vs exp). It maps
    #  dGCV/d(raw_lambda_w) to dGCV/d(raw_lambda_r) via the ratio of traces:
    #    dGCV/d(raw_lambda_r) ~ trace(L2)/trace(Lambda) * dGCV/d(raw_lambda_w)
    #  The parameterization-specific chain rule is then applied below.
    dGCV_u_dlambda2 <- mean(diag(outlist$L2))/
      mean(diag(outlist$Lambda)) *
      dGCV_u_dlambda1

    ## [Change 2026-02-15] softplus -> exp chain rule
    #
    #  Under exp parameterization: raw_penalty = exp(theta)
    #  Chain rule: d(GCV)/d(theta) = d(GCV)/d(raw_penalty) * d(exp(theta))/d(theta)
    #            = d(GCV)/d(raw_penalty) * exp(theta)
    #            = d(GCV)/d(raw_penalty) * raw_penalty
    if(verbose){
      cat('        Gradient start \n')
    }
    gradient <- cbind(c(dGCV_u_dlambda1 * lambda_1,
                        dGCV_u_dlambda2 * lambda_2))

    ## Gradients for unique penalties per predictor
    if(unique_penalty_per_predictor){
      predictor_penalties <-
        penalty_vec[grep('predictor',names(penalty_vec))]

      predictor_penalty_gradient <-
        sapply(1:length(predictor_penalties),function(j){
          mean(diag(outlist$L_predictor_list[[j]]))/
            mean(diag(outlist$Lambda)) *     # trace-ratio (raw scale, unchanged)
            dGCV_u_dlambda1 *
            predictor_penalties[j]           # [Change: exp chain rule, was /(1 + .)]
        })
      gradient <- cbind(c(c(gradient),
                          predictor_penalty_gradient))
    }

    ## Gradients for unique penalties per partition
    if(unique_penalty_per_partition){
      partition_penalties <-
        penalty_vec[grep('partition',names(penalty_vec))]

      partition_penalty_gradient <-
        sapply(1:length(partition_penalties),function(j){
          mean(diag(outlist$L_partition_list[[j]]))/
            mean(diag(outlist$Lambda +                  # ratio of trace requires
                        outlist$L_partition_list[[j]])) * # addition because unlike
            dGCV_u_dlambda1 *                           # predictor penalties,
            partition_penalties[j]          # [Change: exp chain rule, was /(1 + .)]
        })                                            # partition penalties
      # have not been added to
      ## Combined gradient                            # Lambda yet
      gradient <- cbind(c(gradient, partition_penalty_gradient))
    }

    ## [Change: exp chain rule for regularizer]
    #  Regularization penalty - pull penalty terms to 1
    #  Objective: 0.5 * c * (exp(theta) - 1)^2
    #  Gradient:  c * (exp(theta) - 1) * exp(theta) = c * (lambda - 1) * lambda
    #
    #  Note: the original softplus code had c*(lambda - 1) without the chain
    #  rule factor
    #  Under exp, the exact chain rule factor is simply lambda.
    if(unique_penalty_per_partition | unique_penalty_per_predictor){
      regulizer <- c(1e-32*(wiggle_penalty - 1)*wiggle_penalty,
                     0,
                     meta_penalty*(penalty_vec - 1)*penalty_vec)
      meta_penalty <- 0.5*meta_penalty*sum((penalty_vec - 1)^2) +
        0.5*1e-32*((wiggle_penalty - 1)^2)

    } else {
      regulizer <- c(1e-32*(wiggle_penalty - 1)*wiggle_penalty,
                     0)
      meta_penalty <- 0.5*1e-32*((wiggle_penalty - 1)^2)
    }
    if(verbose){
      cat('        Gradient end \n')
    }

    ## Return output
    return(list(GCV_u = outlist$GCV_u + meta_penalty,
                gradient = nr*gradient + regulizer,
                outlist = outlist))
  }

  ## Damped BFGS optimization for minimizing GCV criterion
  quasi_nr_fxn <- function(par){

    ## Initialize optimization parameters and storage
    lambda <- par # Current lambda values
    old_lambda <- lambda # Previous iteration lambda
    new_lambda <- lambda # Proposed next lambda
    n_params <- length(lambda)
    Id <- diag(n_params) # Identity matrix for BFGS updates
    initial_damp <- 0.5 # Initial damping factor for step size
    damp <- initial_damp
    prev_gradient <- lambda*0 # Previous gradient
    gradient <- prev_gradient # Current gradient
    best_gradient <- Inf
    old_gradient <- prev_gradient
    outlist <- NULL # Storage for GCV computations
    prev_outlist <- NULL
    gcv_u <- Inf # Current GCV value
    ridge <- NULL # Ridge term for numerical stability
    dont_skip_gr <- TRUE # Flag to compute gradient
    rho <- NULL # BFGS scaling parameter
    Inv <- diag(n_params) # Approximated inverse Hessian
    best_lambda <- lambda # Best lambda found so far
    best_gcv_u <- gcv_u # Best GCV value found
    prev_lambda <- old_lambda
    restart <- TRUE # Flag to restart BFGS approximation

    ## Main optimization loop with max 100 iterations
    for (iter in 1:100) {
      ## Compute gradient if needed
      if(dont_skip_gr){
        result <- gr_fxn(c(lambda[1],lambda[2], log_penalty_vec), outlist)
        gradient <- result$gradient
        outlist <- result$outlist
      }

      ## For first two iterations, use steepest descent
      if(iter <= 2){
        new_lambda <- lambda - damp * gradient

        ## [Change: softplus -> exp] Reset to best solution if numerical issues
        if(any(!is.finite(exp(new_lambda))) |
           any(is.nan(exp(new_lambda))) |
           any(is.na(exp(new_lambda)))){
          new_lambda <- best_lambda
        }

      } else {
        ## Add small ridge term for stability if needed
        if(any(is.null(ridge))) ridge <-
            1e-8*diag(length(lambda))

        ##s Update BFGS approximation
        if(dont_skip_gr){
          ## Initial BFGS approximation on 3rd iteration or after restart
          if(iter == 3 | restart){
            ## Sherman-Morrison-Woodbury update for initial inverse
            # Hessian approximation
            diff_grad <- gradient - prev_gradient
            diff_lam <- cbind(lambda - prev_lambda)
            denom <- as.numeric(t(diff_grad) %**% diff_lam)
            Inv <- Inv + (denom + as.numeric(t(diff_grad) %**%
                                               Inv %**%
                                               diff_grad)) * (diff_lam %**%
                                                                t(diff_lam)) /
              (denom^2) -
              (Inv %**% cbind(diff_grad) %**%
                 t(diff_lam) +
                 diff_lam %**%
                 t(diff_grad) %**%
                 Inv) /
              denom
            restart <- FALSE
          }

          ## Standard BFGS update for inverse Hessian approximation
          if(iter > 3){
            diff_grad <- gradient - prev_gradient
            diff_lam <- cbind(lambda - prev_lambda)
            denom <- as.numeric(t(diff_grad) %**% diff_lam)
            if(!is.na(denom) && abs(denom) > 1e-64){
              rho <- 1 / denom
              term1 <- Id - rho * (diff_lam %**% t(diff_grad))
              term2 <- Id - rho * (cbind(diff_grad) %**% t(diff_lam))
              Inv <- term1 %**% Inv %**% term2 + rho * (diff_lam %**%
                                                          t(diff_lam))
            } else {
              ## Reset if update is numerically unstable
              Inv <- diag(length(gradient))
              restart <- TRUE
            }
          }
        }

        ## Compute BFGS step direction
        new_lambda <- lambda - damp * Inv %**% cbind(gradient)

        ## [Change: softplus -> exp] Reset to best solution if numerical issues
        if(any(!is.finite(exp(new_lambda))) |
           any(is.nan(exp(new_lambda))) |
           any(is.na(exp(new_lambda)))){
          new_lambda <- best_lambda
        }
      }

      ## Evaluate GCV at new point
      if(any(is.na(new_lambda))){
        ## Backtrack if invalid step
        lambda <- old_lambda
        gradient <- old_gradient
        dont_skip_gr <- FALSE
        damp <- damp / 2
        if(damp < 2^-12 & iter > 9) return(list(par = best_lambda,
                                                gcv_u = best_gcv_u,
                                                iterations = iter))
        next
      } else {
        prev_outlist <- outlist
        outlist <- gcvu_fxn(c(new_lambda[1],
                              new_lambda[2],
                              log_penalty_vec))
        new_gcv_u <- outlist$GCV_u
        if(any(is.na(new_gcv_u))){
          new_gcv_u <- gcv_u
        }
      }

      ## Accept step if improvement or early iterations
      if (new_gcv_u <= gcv_u | iter <= 2) {
        ## Update solution history
        old_gradient <- prev_gradient
        prev_gradient <- gradient
        dont_skip_gr <- TRUE
        prev_outlist <- outlist
        old_lambda <- prev_lambda
        prev_lambda <- lambda
        lambda <- new_lambda
        damp <- 1

        ## Track best solution
        if(new_gcv_u <= gcv_u | iter == 1){
          best_gcv_u <- new_gcv_u
          best_lambda <- lambda
        }

        ## Check convergence criteria
        if(((abs(new_gcv_u - gcv_u) < tol) |
            (max(abs(lambda - prev_lambda)) < tol))  &
           (iter > 9)) return(list(par = best_lambda,
                                   gcv_u = best_gcv_u,
                                   iterations = iter))
        gcv_u <- new_gcv_u
      } else {
        ## Reject step and backtrack
        dont_skip_gr <- FALSE
        outlist <- prev_outlist
        damp <- damp / 2
        if(damp < 2^(-10) & iter > 0) return(
          list(par = best_lambda,
               gcv_u = best_gcv_u,
               iterations = iter))
      }
    }
    return(list(par = best_lambda, gcv_u = best_gcv_u, iterations = iter))
  }
  if(verbose){
    cat('    Starting grid search for initialization\n')
  }

  ## If optimization is desired for penalties
  if(opt){

    ## Create all combinations of the grid values
    initial_grid <- expand.grid(wiggle = log_initial_wiggle,
                                flat = log_initial_flat)

    ## Function to safely evaluate gcv_u
    safe_gcvu <- function(par) {
      tryCatch({
        result <-
          gcvu_fxn(c(unlist(par), log_penalty_vec))$GCV_u
        if(is.na(result) | is.nan(result)) {
          return(Inf)
        }
        return(result)
      }, error = function(e) {
        if(include_warnings)
          return(Inf)
      })
    }

    ## Evaluate GCV_u for each grid point
    gcv_values <- apply(initial_grid, 1, safe_gcvu)
    bads <- which(is.na(gcv_values) |
                    is.nan(gcv_values) |
                    !is.finite(gcv_values))
    if(verbose){
      cat('    Finished grid evaluations\n')
    }
    if(length(bads) == length(gcv_values)){
      stop('All GCV criteria for the initial tuning grid were computed as NA,',
           ' NaN, or non-finite: check your data for corrupt or missing values,',
           ' try changing initial tuning grid, or try manual tuning instead.',
           ' If you are setting no_intercept = TRUE, try experimenting with',
           ' standardize_response = FALSE vs. TRUE.')
    } else if(length(bads) > 0){
      gcv_values <- gcv_values[-c(bads)]
    }

    ## Find the best starting point
    best_index <- which.min(gcv_values)[1]
    best_start <- as.numeric(initial_grid[best_index, ])

    ## Run quasi-Newton method with the best starting point
    if(verbose){
      cat('    Best from grid search: ', cbind(c(best_start)), '\n')
    }
    ## lgspline BFGS custom implementation (closed-form gradient)
    if(use_custom_bfgs){
      res <- withCallingHandlers(
        try(quasi_nr_fxn(c(best_start, log_penalty_vec)), silent = TRUE),
        warning = function(w) if (include_warnings) warning(w) else invokeRestart("muffleWarning"),
        message = function(m) if (include_warnings) message(m) else invokeRestart("muffleMessage")
      )
      if(any(inherits(res, 'try-error'))){
        if(include_warnings) print(res)
        if(include_warnings) warning('Custom BFGS implementation failed. Try use_custom_bfgs =',
                                     ' FALSE, or manual tuning. Resorting to best as selected from',
                                     ' grid search.')
        par <- c(best_start, log_penalty_vec)
      } else {
        par <- res$par
      }
    } else {
      ## vs. base R (finite-difference approx.)
      res <- withCallingHandlers(
        try({optim(c(best_start, log_penalty_vec),
                   fn = function(par){
                     gcvu_fxn(par)$GCV_u
                   },
                   method = 'BFGS'
        )}, silent = TRUE),
        warning = function(w) if (include_warnings) warning(w) else invokeRestart("muffleWarning"),
        message = function(m) if (include_warnings) message(m) else invokeRestart("muffleMessage")
      )
      if(any(inherits(res, 'try-error'))){
        if(include_warnings) print(res)
        if(include_warnings) warning('Base R BFGS failed. Try use_custom_bfgs = TRUE, or manual',
                                     ' tuning. Resorting to best as selected from grid search.')
        par <- c(best_start, log_penalty_vec)
      } else {
        par <- res$par
      }
    }
    if(verbose){
      cat('    Finished tuning penalties\n')
    }

    ## Inflate the penalties, since there's an inherent bias towards
    # no penalization in-sample
    ## [Change: softplus -> exp]
    infl <- ((nr+2)/((nr-2)))^2
    wiggle_penalty <- exp(par[1])*infl
    flat_ridge_penalty <- exp(par[2])*infl

    ## Update penalty vec for predictor-and-partition specific penalties
    if(length(log_penalty_vec) > 0){
      penalty_vec[1:length(penalty_vec)] <-
        exp(c(par[-c(1,2)]))*infl
    }
  } else if(length(penalty_vec) == 0){
    ## No per-predictor or per-partition penalties; nothing to do
    # penalty_vec remains c() from input
  }
  ## When !opt and length(penalty_vec) > 0, penalty_vec is already on raw
  #  scale from input and used as-is


  if(verbose){
    cat('    Final update\n')
  }

  ## Update prior precision
  Lambda_list <- compute_Lambda(custom_penalty_mat,
                                smoothing_spline_penalty,
                                wiggle_penalty,
                                flat_ridge_penalty,
                                K,
                                nc,
                                unique_penalty_per_predictor,
                                unique_penalty_per_partition,
                                penalty_vec,
                                colnm_expansions,
                                just_Lambda = FALSE)

  return(list("Lambda" = Lambda_list$Lambda,
              "flat_ridge_penalty" = flat_ridge_penalty,
              "wiggle_penalty" = wiggle_penalty,
              "other_penalties" = penalty_vec,
              "L_predictor_list" = Lambda_list$L_predictor_list,
              "L_partition_list" = Lambda_list$L_partition_list))
}


#' Efficient Matrix Multiplication of G and A Matrices
#'
#' @param G List of G matrices
#' @param A Constraint matrix
#' @param K Number of partitions minus 1
#' @param nc Number of columns per partition
#' @param nca Number of constraint columns
#' @param parallel Use parallel processing
#' @param cl Cluster object
#' @param chunk_size Size of parallel chunks
#' @param num_chunks Number of chunks
#' @param rem_chunks Remaining chunks
#'
#' @details
#' Computes G %**% A when G has block diagonal structure and A is a matrix
#' Processes in parallel chunks if enabled. Avoids unnecessary operations.
#'
#' @return Matrix product
#'
#' @keywords internal
#' @export
GAmult_wrapper <- function(G,
                           A,
                           K,
                           nc,
                           nca,
                           parallel,
                           cl,
                           chunk_size,
                           num_chunks,
                           rem_chunks) {

  if(parallel & !is.null(cl)) {

    ## Handle remainder chunks first
    if(rem_chunks > 0) {
      rem_indices <- num_chunks * chunk_size + 1:rem_chunks
      G_rem <- G[rem_indices]
      rem <- GAmult(G_rem, A[(num_chunks*chunk_size)*nc + 1:((rem_chunks)*nc),],
                    rem_chunks-1, nc, nca)
    } else {
      rem <- list()
    }

    ## Process main chunks in parallel
    result <- c(
      Reduce("c",parallel::parLapply(cl, 1:num_chunks, function(chunk) {
        start_idx <- (chunk - 1) * chunk_size
        G_chunk <- G[(start_idx + 1):min(start_idx + chunk_size, length(G))]
        A_chunk <- A[(start_idx*nc + 1):min((start_idx + chunk_size)*nc,
                                            nrow(A)), ]
        GAmult(G_chunk,
               A_chunk,
               chunk_size-1,
               nc,
               nca)
      })),
      rem
    )

  } else {
    ## Sequential computation using original C++ function
    result <- GAmult(G, A, K, nc, nca)
  }

  return(result)
}

## [Change 2026-02-16] Removed blockfit path (old Path 2) entirely;
#  blockfit is now handled by a separate blockfit_solve function called
#  by lgspline.fit before dispatching to get_B. Parameters blockfit,
#  just_linear_without_interactions, prevB, prevUnconB, iter_count,
#  and prev_diff are retained in the signature for call-site
#  compatibility but are no longer used internally. IRWLS recursion
#  replaced with a stable for-loop. Documentation rewritten.

#' Compute Constrained GLM Coefficient Estimates via Lagrangian Multipliers
#'
#' @description
#' Core estimation function for Lagrangian multiplier smoothing splines.
#' Computes penalized coefficient estimates subject to smoothness constraints
#' (continuity, differentiability at knots) and optional user-supplied linear
#' equality or inequality constraints. Dispatches to one of three computational
#' paths depending on the model structure:
#'
#' \bold{Path 1. GEE without blockfitting:}
#' When \code{Vhalf} and \code{VhalfInv} are provided, the function uses
#' generalized estimating equations. X and y arrive UNWHITENED; whitening
#' is applied internally to the full N x P block-diagonal design matrix
#' to preserve cross-partition correlation contributions.
#' Two sub-paths:
#' \itemize{
#'   \item \bold{Path 1a: Gaussian identity + GEE}: closed-form full-system
#'     solve in whitened space.
#'   \item \bold{Path 1b: Non-Gaussian GEE}: damped SQP with full whitened
#'     design X_tilde = V^{-1/2} X_block.
#' }
#'
#' \bold{Path 2. Gaussian identity link, no correlation:}
#' For the canonical Gaussian case, the constrained estimate has a closed-form
#' solution via projection. No unconstrained fitting step or iteration is
#' needed; \eqn{\hat{\beta} = G X^T y} projected onto the constraint null
#' space via an OLS-style residual computation.
#'
#' \bold{Path 3. Non-Gaussian GLMs, no correlation:}
#' For GLMs with non-identity links or non-Gaussian families, unconstrained
#' estimates are first obtained per partition, then projected onto the
#' constraint space. For non-canonical links, iterative (IRWLS-like) fitting
#' is performed via a for-loop until convergence.
#'
#' In Paths 2 and 3, the Lagrangian projection is computed efficiently via
#' an OLS reformulation: letting \eqn{y^* = G^{1/2} \hat{\beta}_{unconstrained}}
#' and \eqn{X^* = G^{1/2} A}, the constrained estimate equals
#' \eqn{G^{1/2}} times the residuals from regressing \eqn{y^*} on \eqn{X^*}.
#' This avoids explicitly forming and inverting \eqn{A^T G A}.
#'
#' @param X List of length \eqn{K+1}; partition-specific design matrices.
#' @param X_gram List of Gram matrices \eqn{X_k^T W_k X_k} by partition.
#' @param Lambda Combined penalty matrix (smoothing spline + ridge).
#' @param keep_weighted_Lambda Logical; if \code{TRUE}, retain GLM weights
#'   in the smoothing spline penalty during IRWLS updates.
#' @param unique_penalty_per_partition Logical; if \code{TRUE}, add
#'   partition-specific penalties from \code{L_partition_list} to
#'   \code{Lambda}.
#' @param L_partition_list List of partition-specific penalty matrices.
#' @param A Constraint matrix. Columns encode smoothness constraints
#'   (continuity and derivative matching at knots) and any user-supplied
#'   equality constraints. The system \eqn{A^T \beta = c} is enforced,
#'   where \eqn{c} comes from \code{constraint_value_vectors}.
#' @param Xy List of cross-products \eqn{X_k^T y_k} by partition.
#' @param y List of response vectors by partition (UNWHITENED even when GEE).
#' @param K Integer; number of interior knots. The predictor domain is
#'   partitioned into \eqn{K+1} segments.
#' @param nc Integer; number of coefficients per partition.
#' @param nca Integer; number of columns in \code{A}.
#' @param Ghalf List of \eqn{G^{1/2}} matrices by partition, where
#'   \eqn{G = (X^T W X + \Lambda)^{-1}}.
#' @param GhalfInv List of \eqn{G^{-1/2}} matrices by partition.
#' @param parallel_eigen,parallel_aga,parallel_matmult,parallel_unconstrained
#'   Logical flags for parallelizing eigendecompositions, \eqn{A^T G A}
#'   computation, matrix multiplications, and unconstrained partition fits.
#' @param cl Cluster object from \code{parallel::makeCluster}.
#' @param chunk_size,num_chunks,rem_chunks Partition distribution
#'   parameters for parallel workers.
#' @param family GLM family object (link, variance, and optionally
#'   \code{custom_dev.resids}).
#' @param unconstrained_fit_fxn Function for partition-wise unconstrained
#'   coefficient estimation (non-Gaussian GLMs only).
#' @param iterate Logical; if \code{TRUE}, iterate (IRWLS) for
#'   non-canonical links.
#' @param qp_score_function Score function for quadratic programming.
#'   Signature: \code{function(X, y, mu, order_list, dispersion,
#'   VhalfInv, ...)}.
#' @param quadprog Logical; if \code{TRUE}, use \code{quadprog::solve.QP}
#'   to handle inequality constraints.
#' @param qp_Amat,qp_bvec,qp_meq Inequality constraint specification
#'   for \code{quadprog::solve.QP}.
#' @param prevB Retired; kept for call-site compatibility. Previously
#'   held the constrained estimate from the prior IRWLS iteration.
#'   Ignored internally; IRWLS state is now managed by a local for-loop.
#' @param prevUnconB Retired; kept for call-site compatibility. Previously
#'   held the unconstrained estimates, reused across IRWLS iterations.
#'   Ignored internally.
#' @param iter_count Retired; kept for call-site compatibility. Previously
#'   the IRWLS iteration counter. Ignored internally.
#' @param prev_diff Retired; kept for call-site compatibility. Previously
#'   the mean absolute coefficient change from the prior IRWLS step.
#'   Ignored internally.
#' @param tol Convergence tolerance for coefficient change and deviance.
#' @param constraint_value_vectors List encoding the right-hand side
#'   \eqn{c} in \eqn{A^T \beta = c}. When empty or all zero, the
#'   homogeneous system is solved.
#' @param order_list List of index vectors mapping partition rows to
#'   original data ordering.
#' @param glm_weight_function Function computing GLM working weights.
#' @param schur_correction_function Function computing Schur complement
#'   corrections to the information matrix for nuisance parameter
#'   uncertainty.
#' @param need_dispersion_for_estimation Logical; if \code{TRUE},
#'   dispersion enters weight and score functions.
#' @param dispersion_function Dispersion estimation function.
#' @param observation_weights List of observation weights by partition.
#' @param homogenous_weights Logical; constant weights across observations.
#' @param return_G_getB Logical; if \code{TRUE}, return covariance
#'   components alongside coefficients.
#' @param blockfit Retired; kept for call-site compatibility. Blockfit
#'   estimation is now handled by a separate \code{blockfit_solve}
#'   function before dispatching to \code{get_B}. Ignored internally.
#' @param just_linear_without_interactions Retired; kept for call-site
#'   compatibility. Was used to identify covariates pooled across
#'   partitions in blockfit mode. Ignored internally.
#' @param Vhalf,VhalfInv Square root and inverse square root of the
#'   working correlation matrix for GEE estimation. X and y are
#'   UNWHITENED; whitening is applied internally.
#' @param ... Passed to fitting, weight, correction, and dispersion
#'   functions.
#'
#' @return
#' If \code{return_G_getB = FALSE}: a list of coefficient column vectors,
#' one per partition.
#'
#' If \code{return_G_getB = TRUE}: a list with elements:
#' \describe{
#'   \item{B}{List of coefficient column vectors by partition.}
#'   \item{G_list}{List with components \code{G}, \code{Ghalf}, and
#'     \code{GhalfInv} (partition-wise covariance matrices and their
#'     matrix square roots). For GEE paths, G contains the diagonal
#'     blocks of the full information matrix inverse; the full-system
#'     G^{1/2} is used for the solve but per-partition diagonal blocks
#'     are returned for downstream inference compatibility.}
#' }
#'
#' @keywords internal
#' @export
get_B <- function(X,
                  X_gram,
                  Lambda,
                  keep_weighted_Lambda,
                  unique_penalty_per_partition,
                  L_partition_list,
                  A,
                  Xy,
                  y,
                  K,
                  nc,
                  nca,
                  Ghalf,
                  GhalfInv,
                  parallel_eigen,
                  parallel_aga,
                  parallel_matmult,
                  parallel_unconstrained,
                  cl,
                  chunk_size,
                  num_chunks,
                  rem_chunks,
                  family,
                  unconstrained_fit_fxn,
                  iterate,
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
                  constraint_value_vectors,
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
                  ...){

  ## Path 1: GEE (correlation structures present)
  #  [Change 2026-02-16] X and y are now UNWHITENED (the buggy per-
  #  partition whitening in lgspline.fit has been removed). This path
  #  forms the full whitened design X_tilde = V^{-1/2} X_block
  #  internally, preserving cross-partition correlation contributions.
  #  Two sub-paths:
  #    1a. Gaussian identity + GEE: closed-form solve in full P-space.
  #    1b. Non-Gaussian GEE: damped SQP with full whitened design.
  ## [Change 2026-02-16] Simplified condition: blockfit guard removed.
  if(!is.null(Vhalf) & !is.null(VhalfInv)){

    ## Permute correlation matrices to partition ordering
    perm <- unlist(order_list)
    Vhalf_perm    <- Vhalf[perm, perm]
    VhalfInv_perm <- VhalfInv[perm, perm]

    if(is.null(observation_weights) |
       length(unlist(observation_weights)) == 0){
      observation_weights <- 1
    }

    ## Assemble the unwhitened block-diagonal design
    X_block <- collapse_block_diagonal(X)
    y_block <- cbind(unlist(y))

    ## Full whitened design (NOT block-diagonal) and whitened response
    X_tilde <- VhalfInv_perm %**% X_block
    y_tilde <- VhalfInv_perm %**% y_block

    ## Block-diagonal penalty matrix
    if(unique_penalty_per_partition){
      Lambda_block <- Reduce("rbind", lapply(1:(K+1), function(k){
        Reduce("cbind", lapply(1:(K+1), function(j){
          if(j == k) Lambda + L_partition_list[[k]] else 0 * Lambda
        }))
      }))
    } else {
      Lambda_block <- Reduce("rbind", lapply(1:(K+1), function(k){
        Reduce("cbind", lapply(1:(K+1), function(j){
          if(j == k) Lambda else 0 * Lambda
        }))
      }))
    }

    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    ## Path 1a: Gaussian identity + GEE (closed-form full-system solve)
    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    is_gauss_id <- (paste0(family)[1] == 'gaussian' &
                      paste0(family)[2] == 'identity')

    if(is_gauss_id){

      ## Full P x P Gram matrix (includes cross-partition blocks)
      Gram_full <- crossprod(X_tilde)

      ## Full G = (X_tilde^T X_tilde + Lambda_block)^{-1}
      G_full_inv <- Gram_full + Lambda_block
      G_full <- invert(G_full_inv)

      ## G^{1/2} via eigendecomposition of the full P x P G
      eig_G <- eigen(G_full, symmetric = TRUE)
      vals_G <- eig_G$values
      vals_G[vals_G <= 0] <- 0
      G_full_half <- eig_G$vectors %**%
        (t(eig_G$vectors) * sqrt(vals_G))

      ## G^{-1/2} for per-partition extraction
      inv_sqrt_vals_G <- ifelse(sqrt(pmax(eig_G$values, 0)) > 0,
                                1 / sqrt(pmax(eig_G$values, 0)), 0)
      G_full_half_inv <- eig_G$vectors %**%
        (t(eig_G$vectors) * inv_sqrt_vals_G)

      ## Cross-product X_tilde^T y_tilde (full P x 1)
      Xy_tilde <- crossprod(X_tilde, y_tilde)

      ## Lagrangian projection in full P-dimensional space
      if(K == 0 | any(all.equal(unique(A), 0) == TRUE)){
        A_proj <- cbind(rep(0, nc * (K + 1)))
      } else {
        A_proj <- A
      }

      ## Use the G^{1/2}r* trick (see details)
      y_star <- G_full_half %**% Xy_tilde
      X_star <- G_full_half %**% A_proj

      comp_stab_sc <- 1 / sqrt(K + 1)
      resids_star <- .lm.fit(
        X_star * comp_stab_sc,
        y_star * comp_stab_sc
      )$residuals / comp_stab_sc

      ## Nonzero constraint values: add particular solution
      if(length(constraint_value_vectors) > 0){
        if(any(unlist(constraint_value_vectors) != 0)){
          preds_star <- X_star %**%
            (invert(crossprod(X_star) * comp_stab_sc) %**%
               crossprod(A_proj,
                         Reduce("rbind",
                                constraint_value_vectors) *
                           comp_stab_sc))
          resids_star <- resids_star + c(preds_star)
        }
      }

      beta_block <- G_full_half %**% cbind(resids_star)

      ## Unpack into per-partition result list
      result <- lapply(1:(K+1), function(k){
        cbind(beta_block[(k-1)*nc + 1:nc])
      })

      ## Optional QP refinement for inequality constraints
      if(quadprog){
        ## Merge equality + inequality constraints
        if(K == 0 | any(all.equal(unique(A), 0) == TRUE)){
          A_qp <- cbind(rep(0, nc * (K + 1)))
        } else {
          A_qp <- A
        }
        if(!is.null(qp_Amat)){
          qp_Amat_c <- cbind(A_qp, qp_Amat)
        } else {
          qp_Amat_c <- A_qp
        }
        if(length(constraint_value_vectors) > 0){
          constr_rhs <- Reduce('rbind', constraint_value_vectors)
          if(nrow(constr_rhs) < ncol(A_qp)){
            constr_rhs <- c(rep(0, ncol(A_qp) - nrow(constr_rhs)),
                            constr_rhs)
          }
        } else {
          constr_rhs <- rep(0, ncol(A_qp))
        }
        if(!is.null(qp_bvec)){
          qp_bvec_c <- c(constr_rhs, qp_bvec)
        } else {
          qp_bvec_c <- constr_rhs
        }
        if(!is.null(qp_meq)){
          qp_meq_c <- ncol(A_qp) + qp_meq
        } else {
          qp_meq_c <- ncol(A_qp)
        }

        ## Use full whitened X_tilde for info matrix
        info <- Gram_full + Lambda_block
        sc <- sqrt(mean(abs(info)))

        qp_score <- qp_score_function(
          X_tilde, y_tilde,
          VhalfInv_perm %**% cbind(family$linkinv(
            c(X_block %**% beta_block))),
          unlist(order_list), 1, VhalfInv_perm,
          unlist(observation_weights), ...
        )
        beta_new <- try({quadprog::solve.QP(
          Dmat = info / sc,
          dvec = (qp_score -
                    Lambda_block %**% beta_block +
                    info %**% beta_block) / sc,
          Amat = qp_Amat_c,
          bvec = qp_bvec_c,
          meq  = qp_meq_c
        )$solution}, silent = TRUE)

        if(!any(inherits(beta_new, 'try-error'))){
          beta_block <- cbind(beta_new)
          result <- lapply(1:(K+1), function(k){
            cbind(beta_block[(k-1)*nc + 1:nc])
          })
        }
      }

      ## Return with G_list
      if(return_G_getB){
        ## Extract per-partition diagonal blocks of G_full for
        #  downstream inference. Off-diagonal blocks are needed for the
        #  correct projection but get_U, compute_trace_H, etc. only
        #  use diagonal blocks -- a pre-existing approximation in the
        #  GEE pipeline.
        G_diag <- lapply(1:(K+1), function(k){
          G_full[(k-1)*nc + 1:nc, (k-1)*nc + 1:nc]
        })
        Ghalf_diag <- lapply(1:(K+1), function(k){
          G_full_half[(k-1)*nc + 1:nc, (k-1)*nc + 1:nc]
        })
        GhalfInv_diag <- lapply(1:(K+1), function(k){
          G_full_half_inv[(k-1)*nc + 1:nc, (k-1)*nc + 1:nc]
        })
        ## [Change 2026-02-16] Recompute G_diag as tcrossprod(Ghalf_diag)
        #  so that downstream varcov (sigma^2 * UG) is computed from the
        #  same Ghalf that will be used in matmult_U. This ensures
        #  G = Ghalf %*% t(Ghalf) exactly, avoiding floating-point
        #  asymmetry between G_full diagonal blocks and their square roots.
        G_diag_from_half <- lapply(1:(K+1), function(k){
          tcrossprod(Ghalf_diag[[k]])
        })
        return(list(
          B = result,
          G_list = list(G = G_diag_from_half,
                        Ghalf = Ghalf_diag,
                        GhalfInv = GhalfInv_diag)
        ))
      } else {
        return(result)
      }
    }

    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    ## Path 1b: Non-Gaussian GEE (damped SQP with full whitened design)
    #  [Change 2026-02-16] Uses X_tilde (full whitened design) for the
    #  information matrix. Linear predictor is X_block %*% beta
    #  (original scale) since X_block is unwhitened.
    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    beta_block <- cbind(rep(0, ncol(X_block)))

    ## SQP iteration control
    damp_cnt <- 0
    master_cnt <- 0
    err <- Inf

    ## Original-scale linear predictor (X_block is unwhitened)
    XB <- X_block %**% beta_block

    ## Build combined constraint matrix for solve.QP
    if(K == 0 | any(all.equal(unique(A), 0) == TRUE)){
      A <- cbind(rep(0, nc*(K+1)))
    }
    if(!is.null(qp_Amat)){
      qp_Amat <- cbind(A, qp_Amat)
    } else {
      qp_Amat <- A
    }
    if(length(constraint_value_vectors) > 0){
      constr_A <- Reduce('rbind', constraint_value_vectors)
      if(nrow(constr_A) < ncol(A)){
        constr_A <- c(rep(0, ncol(A) - nrow(constr_A)), constr_A)
      }
    } else {
      constr_A <- rep(0, ncol(A))
    }
    if(!is.null(qp_bvec)){
      qp_bvec <- c(constr_A, qp_bvec)
    } else {
      qp_bvec <- constr_A
    }
    if(!is.null(qp_meq)){
      qp_meq <- ncol(A) + qp_meq
    } else {
      qp_meq <- ncol(A)
    }

    ## Damped SQP loop
    while(err > tol & damp_cnt < 10 & master_cnt < 100){
      master_cnt <- master_cnt + 1
      damp <- 2^(-(damp_cnt))

      ## Dispersion on original scale
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

      ## GLM working weights on original scale
      W <- c(glm_weight_function(family$linkinv(XB),
                                 y_block,
                                 unlist(order_list),
                                 family, dispersion_temp,
                                 unlist(observation_weights), ...))

      ## Schur correction (X_block is original-scale design)
      result <- lapply(1:(K+1), function(k){
        cbind(beta_block[(k-1)*nc + 1:nc])
      })
      schur_correction <-
        schur_correction_function(
          list(X_block), list(y_block),
          list(cbind(unlist(result))),
          dispersion_temp, list(unlist(order_list)),
          0, family, unlist(observation_weights), ...
        )
      if(!(any(unlist(schur_correction) != 0))){
        schur_correction_collapsed <- 0
      } else {
        schur_correction_collapsed <-
          collapse_block_diagonal(schur_correction)
      }

      ## Info matrix using full whitened design X_tilde
      info <- crossprod(X_tilde, W * X_tilde) +
        Lambda_block +
        schur_correction_collapsed

      sc <- sqrt(mean(abs(info)))

      if(master_cnt == 1){
        ## First iteration: constrained Newton step
        infoinv_block <- sc * invert(sc * info)
        U <- (diag(1, nrow(info)) -
                tcrossprod(infoinv_block %**%
                             A %**%
                             invert(crossprod(A,
                                              infoinv_block %**% A)),
                           A))
        ## Score: X_tilde^T (y_tilde - V^{-1/2} mu)
        mu_vec <- cbind(family$linkinv(c(XB)))
        beta_new <- c(
          U %**% infoinv_block %**%
            crossprod(X_tilde,
                      y_tilde - VhalfInv_perm %**% mu_vec)
        )
        infoinv_block <- NULL
      } else {
        ## Subsequent iterations: solve.QP
        qp_score <- qp_score_function(
          X_tilde, y_tilde,
          VhalfInv_perm %**% cbind(family$linkinv(XB)),
          unlist(order_list), dispersion_temp,
          VhalfInv_perm,
          unlist(observation_weights), ...
        )
        beta_new <- try({quadprog::solve.QP(
          Dmat = info / sc,
          dvec = (qp_score -
                    Lambda_block %**% beta_block +
                    info %**% beta_block) / sc,
          Amat = qp_Amat,
          bvec = qp_bvec,
          meq = qp_meq
        )$solution}, silent = TRUE)
      }

      if(any(inherits(beta_new, 'try-error'))){
        beta_new <- 0 * beta_block
      }

      ## Non-iterative: accept after Newton step
      if(!iterate & master_cnt > 2){
        beta_block <- beta_new
        damp_cnt <- 11
        master_cnt <- 101
        err <- tol - 1
      } else {
        ## Damped update; deviance on original scale
        beta_new <- (1 - damp) * beta_block + damp * beta_new
        XB <- X_block %**% cbind(beta_new)

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
        W <- c(glm_weight_function(family$linkinv(XB),
                                   y_block,
                                   unlist(order_list),
                                   family, dispersion_temp,
                                   unlist(observation_weights), ...))

        ## Deviance on original scale
        mu_new <- family$linkinv(c(XB))
        if(!is.null(family$custom_dev.resids)){
          raw <- family$custom_dev.resids(
            y_block, mu_new, unlist(order_list),
            family, unlist(observation_weights), ...
          )
          W_safe <- pmax(W, sqrt(.Machine$double.eps))
          err_new <- mean((
            VhalfInv_perm %**%
              cbind(sign(raw) * sqrt(abs(raw)) / sqrt(c(W_safe)))
          )^2)
        } else if(is.null(family$dev.resids)){
          err_new <- mean((unlist(observation_weights) *
                             (y_block - cbind(mu_new)))^2)
        } else {
          err_new <-
            mean(unlist(observation_weights) *
                   family$dev.resids(y_block, cbind(mu_new),
                                     wt = 1 / W))
        }

        ## Step acceptance
        if(is.null(err_new) | is.na(err_new) | !is.finite(err_new)){
          damp_cnt <- damp_cnt + 1
        } else if(err_new <= err){
          prev_err <- err
          err <- err_new
          abs_diff <- max(abs(beta_new - beta_block))
          beta_block <- beta_new
          damp_cnt <- 0
          if((abs_diff < tol) & (prev_err - err < tol) &
             (master_cnt > 10)){
            damp_cnt <- 11
            master_cnt <- 101
            err <- tol - 1
          }
        } else {
          damp_cnt <- damp_cnt + 1
        }
      }

      ## Unpack per-partition
      result <- lapply(1:(K+1), function(k){
        cbind(beta_block[(k-1)*nc + 1:nc])
      })
    }

    ## Return with covariance components
    if(return_G_getB){
      ## Recompute G at final estimates using full whitened info
      XB <- X_block %**% cbind(unlist(result))
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
      W <- c(glm_weight_function(family$linkinv(XB),
                                 y_block,
                                 unlist(order_list),
                                 family, dispersion_temp,
                                 unlist(observation_weights), ...))

      schur_correction <-
        schur_correction_function(
          list(X_block), list(y_block),
          list(cbind(unlist(result))),
          dispersion_temp, list(unlist(order_list)),
          0, family, unlist(observation_weights), ...
        )
      if(!(any(unlist(schur_correction) != 0))){
        schur_correction_collapsed <- 0
      } else {
        schur_correction_collapsed <-
          collapse_block_diagonal(schur_correction)
      }

      info <- crossprod(X_tilde, W * X_tilde) +
        Lambda_block +
        schur_correction_collapsed

      ## Extract per-partition diagonal blocks of the full info
      #  matrix and compute G, Ghalf, GhalfInv via eigendecomposition.
      #  [Change 2026-02-16] G_diag is computed as tcrossprod(Ghalf_diag)
      #  to ensure G = Ghalf %*% t(Ghalf) exactly for downstream varcov.
      G <- lapply(1:(K+1), function(k) NA)
      Ghalf <- lapply(1:(K+1), function(k) NA)
      GhalfInv <- lapply(1:(K+1), function(k) NA)
      for(k in 1:(K+1)){
        info_kk <- info[(k-1)*nc + 1:nc, (k-1)*nc + 1:nc]
        eig <- eigen(info_kk, symmetric = TRUE)
        vals <- eig$values
        vals_safe <- pmax(vals, 0)
        vals_safe[vals <= sqrt(.Machine$double.eps)] <- 0
        sqrt_vals <- sqrt(vals_safe)
        Ghalf[[k]] <- eig$vectors %**%
          (t(eig$vectors) * sqrt_vals)
        inv_sqvals <- ifelse(sqrt_vals > 0, 1/sqrt_vals, 0)
        GhalfInv[[k]] <- eig$vectors %**%
          (t(eig$vectors) * inv_sqvals)
        ## G = Ghalf %*% t(Ghalf) exactly
        G[[k]] <- tcrossprod(Ghalf[[k]])
      }

      return(list(
        B = result,
        G_list = list(Ghalf = Ghalf,
                      GhalfInv = GhalfInv,
                      G = G)
      ))
    } else {
      return(result)
    }
  }

  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  ## Paths 2 and 3: No correlation structure
  #  The code below handles both the Gaussian identity case (Path 2, which
  #  has a closed-form solution) and the general GLM case (Path 3, which
  #  requires unconstrained fitting followed by Lagrangian projection).
  #  The family check determines which sub-path to take.
  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

  ## Path 3: Non-Gaussian or non-identity link (requires iterative fitting)
  if(!(paste0(family)[2] == 'identity') |
     !(paste0(family)[1] == 'gaussian')){
    use_lm <- FALSE

    ## Compute Lambda^{1/2} for augmented regression in unconstrained fits
    eig <- eigen(Lambda, symmetric = TRUE)
    LambdaHalf <- eig$vectors %**% (t(eig$vectors) * sqrt(ifelse(
      eig$values > 0,
      eig$values,
      0)))

    ## Obtain unconstrained (partition-wise) estimates, either in parallel
    #  or sequentially.
    if(parallel_unconstrained){
      ## Parallel unconstrained fitting across partitions
      unconstrained_estimate <-
        parallel::parLapply(cl,
                            1:(K+1),
                            function(k,
                                     unique_penalty_per_partition,
                                     unconstrained_fit_fxn,
                                     keep_weighted_Lambda,
                                     family,
                                     tol,
                                     K,
                                     parallel_unconstrained,
                                     cl,
                                     chunk_size,
                                     num_chunks,
                                     rem_chunks,
                                     observation_weights){
                              ## Add partition-specific penalty if requested
                              if(unique_penalty_per_partition){
                                Lambda_temp <- Lambda + L_partition_list[[k]]
                                eig <- eigen(Lambda_temp, symmetric = TRUE)
                                LambdaHalf_temp <- eig$vectors %**%
                                  (t(eig$vectors) * sqrt(ifelse(eig$values <= 0,
                                                                0,
                                                                eig$values)))
                              } else {
                                LambdaHalf_temp <- LambdaHalf
                                Lambda_temp <- Lambda
                              }
                              ## Fit unconstrained model to this partition
                              cbind(
                                c(unconstrained_fit_fxn(X[[k]],
                                                        y[[k]],
                                                        LambdaHalf_temp,
                                                        Lambda_temp,
                                                        keep_weighted_Lambda,
                                                        family,
                                                        tol,
                                                        K,
                                                        parallel_unconstrained,
                                                        cl,
                                                        chunk_size,
                                                        num_chunks,
                                                        rem_chunks,
                                                        order_list[[k]],
                                                        observation_weights[[k]],
                                                        ...)))
                            },
                            unique_penalty_per_partition,
                            unconstrained_fit_fxn,
                            keep_weighted_Lambda,
                            family,
                            tol,
                            K,
                            parallel_unconstrained,
                            cl,
                            chunk_size,
                            num_chunks,
                            rem_chunks,
                            observation_weights)
    } else {
      ## Sequential unconstrained fitting across partitions
      unconstrained_estimate <- lapply(1:(K+1), function(k){
        if(unique_penalty_per_partition){
          Lambda_temp <- Lambda + L_partition_list[[k]]
          eig <- eigen(Lambda_temp, symmetric = TRUE)
          LambdaHalf_temp <- eig$vectors %**% (t(eig$vectors) *
                                                 (sqrt(ifelse(eig$values <= 0,
                                                              0,
                                                              eig$values))))
        } else {
          LambdaHalf_temp <- LambdaHalf
          Lambda_temp <- Lambda
        }
        cbind(c(unconstrained_fit_fxn(X[[k]],
                                      y[[k]],
                                      LambdaHalf_temp,
                                      Lambda_temp,
                                      keep_weighted_Lambda,
                                      family,
                                      tol,
                                      K,
                                      parallel_unconstrained,
                                      cl,
                                      chunk_size,
                                      num_chunks,
                                      rem_chunks,
                                      order_list[[k]],
                                      observation_weights[[k]],
                                      ...)))
      })
    }

    ## Special case: no knots, no QP, no nonzero constraint values.
    #  The unconstrained estimate is already the solution.
    if(K == 0 & !quadprog & length(constraint_value_vectors) == 0){
      if(return_G_getB){
        ## Recompute G at the unconstrained estimate
        if(need_dispersion_for_estimation){
          mu <- family$linkinv(
            unlist(
              matmult_block_diagonal(X,
                                     unconstrained_estimate,
                                     K,
                                     parallel_matmult,
                                     cl,
                                     chunk_size,
                                     num_chunks,
                                     rem_chunks)))
          dispersion_temp <- dispersion_function(
            mu = mu,
            y = unlist(y),
            order_indices = unlist(order_list),
            family = family,
            observation_weights = unlist(observation_weights),
            VhalfInv = VhalfInv,
            ...
          )
        } else {
          dispersion_temp <- 1
        }
        ## Weighted design for G = (X^T W X + Lambda)^{-1}
        Xw <- lapply(1:(K+1),
                     function(k){
                       var <- glm_weight_function(family$linkinv(
                         X[[k]] %**%
                           cbind(c(
                             unconstrained_estimate[[k]]
                           ))
                       ),
                       y[[k]],
                       order_list[[k]],
                       family,
                       dispersion_temp,
                       observation_weights[[k]],
                       ...)
                       cbind(X[[k]] * c(sqrt(var)))
                     })
        X_gram <- compute_gram_block_diagonal(Xw,
                                              parallel_matmult,
                                              cl,
                                              chunk_size,
                                              num_chunks,
                                              rem_chunks)
        schur_corrections <- schur_correction_function(
          X,
          y,
          unconstrained_estimate,
          dispersion_temp,
          order_list,
          K,
          family,
          observation_weights,
          ...
        )
        G_list <- compute_G_eigen(X_gram,
                                  Lambda,
                                  K,
                                  parallel_eigen,
                                  cl,
                                  chunk_size,
                                  num_chunks,
                                  rem_chunks,
                                  family,
                                  unique_penalty_per_partition,
                                  L_partition_list,
                                  keep_G = TRUE,
                                  schur_corrections)
        ## [Change 2026-02-16] Ensure G = tcrossprod(Ghalf) exactly
        G_list$G <- lapply(G_list$Ghalf, function(mat) tcrossprod(mat))
        return(list(
          B = unconstrained_estimate,
          G_list = G_list
        ))
      } else {
        return(unconstrained_estimate)
      }
    }

    ## Transform unconstrained estimate to y* for the OLS projection:
    #  y* = G^{1/2} beta_unconstrained = G^{-1/2} G beta_unconstrained
    #  The identity G^{1/2} = G^{-1/2} G is used because
    #  G beta_unconstrained is already stored as unconstrained_estimate,
    #  and G^{-1/2} is available from the eigendecomposition.
    GhalfXy <- cbind(
      unlist(
        matmult_block_diagonal(
          GhalfInv,
          unconstrained_estimate,
          K,
          parallel_matmult,
          cl,
          chunk_size,
          num_chunks,
          rem_chunks)))

    ## Path 2: Gaussian identity link (closed-form solution)
  } else {
    use_lm <- TRUE

    ## No knots, no QP, no constraint values: direct solution beta = G X^T y
    if(K == 0 & !quadprog & length(constraint_value_vectors) == 0){
      ## [Change 2026-02-16] tcrossprod for symmetric Ghalf
      G <- list(tcrossprod(Ghalf[[1]]))
      result <- list(G[[1]] %**% Xy[[1]])
      if(return_G_getB){
        ## [Change 2026-02-16] Ensure G = tcrossprod(Ghalf) exactly
        Ghalf_exact <- lapply(G, function(Gk){
          eig_k <- eigen(Gk, symmetric = TRUE)
          vals_k <- pmax(eig_k$values, 0)
          eig_k$vectors %**% (t(eig_k$vectors) * sqrt(vals_k))
        })
        GhalfInv_exact <- lapply(Ghalf_exact, function(Gh){
          eig_k <- eigen(tcrossprod(Gh), symmetric = TRUE)
          vals_k <- pmax(eig_k$values, 0)
          sq_v <- sqrt(vals_k)
          inv_sq_v <- ifelse(sq_v > 0, 1/sq_v, 0)
          eig_k$vectors %**% (t(eig_k$vectors) * inv_sq_v)
        })
        return(list(
          B = result,
          G_list = list(G = G,
                        Ghalf = Ghalf_exact,
                        GhalfInv = GhalfInv_exact)
        ))
      } else {
        return(result)
      }
    } else {
      ## With constraints: transform to y* = G^{1/2} X^T y.
      #  For Gaussian identity, G^{1/2} (not G^{-1/2}) multiplies X^T y
      #  because the unconstrained estimate IS G X^T y = G^{1/2}(G^{1/2} X^T y),
      #  so y* = G^{1/2} X^T y is the appropriate transform.
      GhalfXy <- cbind(
        unlist(
          matmult_block_diagonal(Ghalf,
                                 Xy,
                                 K,
                                 parallel_matmult,
                                 cl,
                                 chunk_size,
                                 num_chunks,
                                 rem_chunks)))
    }
  }

  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  #  Common Lagrangian projection (Paths 2 and 3)
  #  The constrained estimate is obtained by projecting the unconstrained
  #  estimate onto the null space of A^T, reformulated as an OLS problem:
  #    y* = G^{1/2} beta_unconstrained  (computed above)
  #    X* = G^{1/2} A                   (computed below)
  #    resid = (I - X*(X*^T X*)^{-1} X*^T) y*
  #    beta_constrained = G^{1/2} resid
  #  This avoids explicitly forming and inverting A^T G A.
  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

  ## Compute X* = G^{1/2} A, stacking partition blocks vertically
  GhalfA <- Reduce("rbind",
                   GAmult_wrapper(Ghalf,
                                  A,
                                  K,
                                  nc,
                                  nca,
                                  parallel_aga,
                                  cl,
                                  chunk_size,
                                  num_chunks,
                                  rem_chunks))

  ## Replace any NA values in y* with the link-transformed zero
  #  (can arise from empty partitions)
  GhalfXy <- ifelse(is.na(GhalfXy), family$linkinv(0), GhalfXy)

  ## OLS residuals: project y* onto the orthogonal complement of col(X*).
  #  The scaling by 1/sqrt(K+1) improves numerical stability of .lm.fit
  #  when the constraint matrix has many rows.
  comp_stab_sc <- 1/sqrt(K + 1)
  resids_star <- do.call('.lm.fit',list(x = GhalfA * comp_stab_sc,
                                        y = GhalfXy * comp_stab_sc)
  )$residuals /  comp_stab_sc

  ## If constraint values are nonzero (A^T beta = c with c != 0),
  #  add the particular solution (I - U) b_0 where A^T b_0 = c.
  #  The full Lagrangian solution is: U beta_hat + (I - U) b_0
  if(length(constraint_value_vectors) > 0){
    if(any(unlist(constraint_value_vectors) != 0)){
      comp_stab_sc <- 1/sqrt(K + 1)
      ## [Change 2026-02-16] crossprod avoids explicit transpose of A
      preds_star <- GhalfA %**%
        (invert(crossprod(GhalfA) * comp_stab_sc) %**%
           crossprod(A,
                     Reduce("rbind", constraint_value_vectors) * comp_stab_sc)
        )
      resids_star <- resids_star + c(preds_star)
    }
  }

  ## Back-transform: beta_constrained = G^{1/2} resid (per partition)
  if(parallel_matmult & !is.null(cl)){
    ## Handle remainder partitions that don't fill a full chunk
    if(rem_chunks > 0) {
      rem_indices <- num_chunks * chunk_size + 1:rem_chunks
      rem <- lapply(rem_indices, function(k) {
        Ghalf[[k]] %**% cbind(resids_star[(k-1)*nc + 1:nc])
      })
    } else {
      rem <- list()
    }
    ## Parallel back-transformation across partition chunks
    result <- c(
      Reduce("c",parallel::parLapply(cl, 1:num_chunks, function(chunk) {
        start_idx <- (chunk - 1) * chunk_size + 1
        end_idx <- chunk * chunk_size
        lapply(start_idx:end_idx, function(k) {
          Ghalf[[k]] %**% cbind(resids_star[(k-1)*nc + 1:nc])
        })
      })),
      rem
    )
  } else {
    result <- lapply(1:(K+1), function(k) {
      Ghalf[[k]] %**% cbind(resids_star[(k-1)*nc + 1:nc])
    })
  }

  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  ## [Change 2026-02-16] IRWLS iteration for non-canonical links (Path 3 only).
  #  Replaced recursive calls to get_B with a for-loop that updates G and
  #  re-projects until convergence, coefficient divergence, or 100 iterations.
  #  The unconstrained estimates are computed once and reused across iterations.
  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

  if(iterate & !use_lm){
    prevB_IRWLS <- NULL
    prev_diff_IRWLS <- Inf

    for(IRWLS_iter in 1:100){
      ## Check convergence against previous constrained estimate
      if(!is.null(prevB_IRWLS)){
        diff_IRWLS <- mean(
          unlist(
            sapply(1:(K+1),
                   function(k)mean(
                     abs(result[[k]] - prevB_IRWLS[[k]])))))
        if(diff_IRWLS < tol) break
        ## If difference is growing, revert to previous and stop
        if(prev_diff_IRWLS <= diff_IRWLS){
          result <- prevB_IRWLS
          break
        }
        prev_diff_IRWLS <- diff_IRWLS
      }

      ## Save current constrained estimate for next convergence check
      prevB_IRWLS <- result

      ## Recompute G at the current constrained estimate
      if(need_dispersion_for_estimation){
        dispersion_temp <- dispersion_function(
          mu = family$linkinv(unlist(matmult_block_diagonal(X,
                                                       result,
                                                       K,
                                                       parallel_matmult,
                                                       cl,
                                                       chunk_size,
                                                       num_chunks,
                                                       rem_chunks))),
          y = unlist(y),
          order_indices = unlist(order_list),
          family = family,
          observation_weights = unlist(observation_weights),
          VhalfInv = VhalfInv,
          ...
        )
      } else {
        dispersion_temp <- 1
      }

      ## Weighted design matrix X W^{1/2} for computing X^T W X
      Xw <- lapply(1:(K+1),
                   function(k){
                     if(nrow(X[[k]]) == 0){
                       return(X[[k]])
                     }
                     var <- glm_weight_function(family$linkinv(
                       X[[k]] %**% cbind(c(result[[k]]))
                     ),
                     y[[k]],
                     order_list[[k]],
                     family,
                     dispersion_temp,
                     observation_weights[[k]],
                     ...)
                     ## [Change 2026-02-16] Direct row-scaling replaces
                     #  t(t(X) %**% D) pattern: avoids two transposes and
                     #  diag() construction.
                     if(length(var) == 1){
                       if(c(var) == 0){
                         return(X[[k]] * 0)
                       } else {
                         return(X[[k]] * c(sqrt(var)))
                       }
                     } else {
                       var <- c(sqrt(var))
                     }
                     X[[k]] * var
                   })
      X_gram_IRWLS <- compute_gram_block_diagonal(Xw,
                                                 parallel_matmult,
                                                 cl,
                                                 chunk_size,
                                                 num_chunks,
                                                 rem_chunks)

      ## Schur corrections at current estimate
      schur_corrections_IRWLS <- schur_correction_function(
        X,
        y,
        result,
        dispersion_temp,
        order_list,
        K,
        family,
        observation_weights,
        ...
      )
      G_list_IRWLS <- compute_G_eigen(X_gram_IRWLS,
                                     Lambda,
                                     K,
                                     parallel_eigen,
                                     cl,
                                     chunk_size,
                                     num_chunks,
                                     rem_chunks,
                                     family,
                                     unique_penalty_per_partition,
                                     L_partition_list,
                                     keep_G = FALSE,
                                     schur_corrections_IRWLS)

      ## Update Ghalf and GhalfInv for this iteration
      Ghalf <- G_list_IRWLS$Ghalf
      GhalfInv <- G_list_IRWLS$GhalfInv

      ## Recompute transformed unconstrained estimate with updated G^{-1/2}
      GhalfXy <- cbind(
        unlist(
          matmult_block_diagonal(
            GhalfInv,
            unconstrained_estimate,
            K,
            parallel_matmult,
            cl,
            chunk_size,
            num_chunks,
            rem_chunks)))

      ## Lagrangian projection with updated G^{1/2}
      GhalfA <- Reduce("rbind",
                       GAmult_wrapper(Ghalf,
                                      A,
                                      K,
                                      nc,
                                      nca,
                                      parallel_aga,
                                      cl,
                                      chunk_size,
                                      num_chunks,
                                      rem_chunks))

      GhalfXy <- ifelse(is.na(GhalfXy), family$linkinv(0), GhalfXy)

      comp_stab_sc <- 1/sqrt(K + 1)
      resids_star <- do.call('.lm.fit',list(x = GhalfA * comp_stab_sc,
                                            y = GhalfXy * comp_stab_sc)
      )$residuals / comp_stab_sc

      if(length(constraint_value_vectors) > 0){
        if(any(unlist(constraint_value_vectors) != 0)){
          comp_stab_sc <- 1/sqrt(K + 1)
          ## [Change 2026-02-16] crossprod avoids explicit transpose of A
          preds_star <- GhalfA %**%
            (invert(crossprod(GhalfA) * comp_stab_sc) %**%
               crossprod(A,
                         Reduce("rbind", constraint_value_vectors) * comp_stab_sc)
            )
          resids_star <- resids_star + c(preds_star)
        }
      }

      ## Back-transform: beta_constrained = G^{1/2} resid (per partition)
      if(parallel_matmult & !is.null(cl)){
        if(rem_chunks > 0) {
          rem_indices <- num_chunks * chunk_size + 1:rem_chunks
          rem <- lapply(rem_indices, function(k) {
            Ghalf[[k]] %**% cbind(resids_star[(k-1)*nc + 1:nc])
          })
        } else {
          rem <- list()
        }
        result <- c(
          Reduce("c",parallel::parLapply(cl, 1:num_chunks, function(chunk) {
            start_idx <- (chunk - 1) * chunk_size + 1
            end_idx <- chunk * chunk_size
            lapply(start_idx:end_idx, function(k) {
              Ghalf[[k]] %**% cbind(resids_star[(k-1)*nc + 1:nc])
            })
          })),
          rem
        )
      } else {
        result <- lapply(1:(K+1), function(k) {
          Ghalf[[k]] %**% cbind(resids_star[(k-1)*nc + 1:nc])
        })
      }
    }
  }

  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
  ## Optional QP refinement (Path 3 only, when quadprog = TRUE)
  #  After obtaining the Lagrangian projection estimate, further refine
  #  with inequality constraints via solve.QP in full block-diagonal form.
  #  This is not parallelized or memory-efficient.
  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

  if(quadprog){

    ## Construct full block-diagonal design and penalty matrices
    X_block <- Reduce("rbind", lapply(1:(K+1),function(k){
      Reduce("cbind",lapply(1:(K+1),function(j){
        if(nrow(X[[k]]) == 0){
          return(X[[k]])
        } else if(j == k) X[[k]] else 0*X[[k]]
      }))
    }))
    beta_block <- cbind(unlist(result))
    if(unique_penalty_per_partition){
      Lambda_block <- Reduce("rbind", lapply(1:(K+1),function(k){
        Reduce("cbind", lapply(1:(K+1), function(j){
          if(j == k) Lambda + L_partition_list[[k]] else 0*Lambda
        }))
      }))
    } else {
      Lambda_block <- Reduce("rbind", lapply(1:(K+1),function(k){
        Reduce("cbind", lapply(1:(K+1), function(j){
          if(j == k) Lambda else 0*Lambda
        }))
      }))
    }
    y_block <- cbind(unlist(y))

    ## Damped SQP loop for QP refinement
    damp_cnt <- 0
    master_cnt <- 0
    err <- Inf
    XB <- X_block %**% beta_block

    while(err > tol & damp_cnt < 10 & master_cnt < 100){
      master_cnt <- master_cnt + 1
      damp <- 2^(-(damp_cnt))

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
      W <- c(glm_weight_function(family$linkinv(XB),
                                 y_block,
                                 unlist(order_list),
                                 family,
                                 dispersion_temp,
                                 unlist(observation_weights),
                                 ...))

      ## Schur correction in partition form
      result <- lapply(1:(K+1),function(k){
        cbind(beta_block[(k-1)*nc + 1:nc])
      })
      schur_correction <-
        schur_correction_function(
          X,
          y,
          result,
          dispersion_temp,
          order_list,
          K,
          family,
          observation_weights,
          ...
        )
      if(!(any(unlist(schur_correction) != 0))){
        schur_correction <- 0
      } else {
        schur_correction <- collapse_block_diagonal(schur_correction)
      }
      ## [Change 2026-02-16] crossprod avoids explicit transpose
      info <- crossprod(X_block, W * X_block) + Lambda_block + schur_correction

      sc <- sqrt(mean(abs(info)))

      ## QP step with both equality (smoothness) and inequality constraints
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
      beta_new <- try({quadprog::solve.QP(
        Dmat = info/sc,
        dvec = (qp_score -
                  Lambda_block %**% beta_block +
                  info %**% beta_block)/sc,
        Amat = cbind(A, qp_Amat),
        bvec = c(rep(0, ncol(A)), qp_bvec),
        meq = ncol(A) + qp_meq
      )$solution}, silent = TRUE)

      if(any(inherits(beta_new, 'try-error'))){
        beta_new <- 0*beta_block
      }

      if(!iterate & master_cnt > 1){
        beta_block <- beta_new
        damp_cnt <- 11
        master_cnt <- 101
        err <- tol - 1
      } else {
        ## Damped update with deviance monitoring.
        #  NOTE: The custom_dev.resids check here requires BOTH
        #  !is.null(custom_dev.resids) AND is.null(dev.resids), which
        #  differs from the GEE path where custom_dev.resids takes
        #  unconditional priority.
        # [Change 2026-02-15] Important comment above for clarification
        beta_new <- (1-damp)*beta_block + damp*beta_new
        XB <- X_block %**% beta_new
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
          err_new <- mean(unlist(observation_weights)*
                            (y_block - family$linkinv(XB))^2)
        } else {
          err_new <-
            mean(unlist(observation_weights)*
                   family$dev.resids(y_block,
                                     family$linkinv(XB),
                                     wt = 1))
        }

        if(is.null(err_new) | is.na(err_new) | !is.finite(err_new)){
          damp_cnt <- damp_cnt + 1
        } else if(err_new <= err){

          prev_err <- err
          err <- err_new
          abs_diff <- max(abs(beta_new - beta_block))
          beta_block <- beta_new
          damp_cnt <- 0

          if((abs_diff < tol) &
             (prev_err - err < tol) &
             (master_cnt > 10)){
            damp_cnt <- 11
            master_cnt <- 101
            err <- tol - 1
          }

        } else {
          damp_cnt <- damp_cnt + 1
        }
      }

      ## Unpack into per-partition form
      result <- lapply(1:(K+1),function(k){
        cbind(beta_block[1:nc + (k-1)*nc])
      })
    }
  }

  ## Return coefficients, optionally with recomputed G
  if(return_G_getB){
    if(paste0(family)[1] == 'gaussian' &
       paste0(family)[2] == 'identity'){
      ## [Change 2026-02-16] Ensure G = tcrossprod(Ghalf) exactly
      #  so downstream varcov computation is consistent.
      G_exact <- lapply(Ghalf, function(mat) tcrossprod(mat))
      return(list(
        B = result,
        G_list = list(G = G_exact,
                      Ghalf = Ghalf,
                      GhalfInv = GhalfInv)
      ))
    }
    ## Non-Gaussian: recompute G at the final estimates
    if(need_dispersion_for_estimation){
      dispersion_temp <- dispersion_function(
        mu = family$linkinv(unlist(matmult_block_diagonal(X,
                                                     result,
                                                     K,
                                                     parallel_matmult,
                                                     cl,
                                                     chunk_size,
                                                     num_chunks,
                                                     rem_chunks))),
        y = unlist(y),
        order_indices = unlist(order_list),
        family = family,
        observation_weights = unlist(observation_weights),
        VhalfInv = VhalfInv,
        ...
      )
    } else {
      dispersion_temp <- 1
    }
    ## X W^{1/2} for computing X^T W X
    Xw <- lapply(1:(K+1),
                 function(k){
                   if(nrow(X[[k]]) == 0){
                     return(X[[k]])
                   }
                   var <- glm_weight_function(family$linkinv(
                     X[[k]] %**% cbind(c(result[[k]]))
                   ),
                   y[[k]],
                   order_list[[k]],
                   family,
                   dispersion_temp,
                   observation_weights[[k]],
                   ...)
                   ## [Change 2026-02-16] Direct row-scaling replaces
                   #  t(t(X) %**% D) pattern
                   if(length(var) == 1){
                     if(c(var) == 0){
                       return(X[[k]] * 0)
                     } else {
                       return(X[[k]] * c(sqrt(var)))
                     }
                   } else {
                     var <- c(sqrt(var))
                   }
                   X[[k]] * var
                 })
    X_gram <- compute_gram_block_diagonal(Xw,
                                          parallel_matmult,
                                          cl,
                                          chunk_size,
                                          num_chunks,
                                          rem_chunks)
    schur_corrections <- schur_correction_function(
      X,
      y,
      result,
      dispersion_temp,
      order_list,
      K,
      family,
      observation_weights,
      ...
    )
    G_list <- compute_G_eigen(X_gram,
                              Lambda,
                              K,
                              parallel_eigen,
                              cl,
                              chunk_size,
                              num_chunks,
                              rem_chunks,
                              family,
                              unique_penalty_per_partition,
                              L_partition_list,
                              keep_G = TRUE,
                              schur_corrections)
    ## [Change 2026-02-16] Ensure G = tcrossprod(Ghalf) exactly
    #  so downstream varcov computation is consistent regardless of
    #  what compute_G_eigen returns for G.
    G_list$G <- lapply(G_list$Ghalf, function(mat) tcrossprod(mat))
    return(list(
      B = result,
      G_list = G_list
    ))
  } else {
    return(result)
  }
}


#' Efficient Matrix Multiplication for \eqn{\textbf{A}^{T}\textbf{G}\textbf{A}}
#'
#' @param G List of G matrices (\eqn{\textbf{G}})
#' @param A Constraint matrix (\eqn{\textbf{A}})
#' @param K Number of partitions minus 1 (\eqn{K})
#' @param nc Number of columns per partition
#' @param nca Number of constraint columns
#' @param parallel Use parallel processing
#' @param cl Cluster object
#' @param chunk_size Chunk size for parallel
#' @param num_chunks Number of chunks
#' @param rem_chunks Remaining chunks
#'
#' @details
#' Computes \eqn{\textbf{A}^{T}\textbf{G}\textbf{A}} efficiently in parallel chunks using \code{AGAmult_chunk()}.
#'
#' @return Matrix product \eqn{\textbf{A}^{T}\textbf{G}\textbf{A}}
#'
#' @keywords internal
#' @export
AGAmult_wrapper <- function(G,
                            A,
                            K,
                            nc,
                            nca,
                            parallel,
                            cl,
                            chunk_size,
                            num_chunks,
                            rem_chunks) {
  if(parallel & !is.null(cl)) {
    # Handle remainder chunks
    if(rem_chunks > 0) {
      rem_start <- num_chunks * chunk_size
      G_rem <- G[(rem_start + 1):(rem_start + rem_chunks)]
      rem_result <- AGAmult_chunk(G_rem,
                                  A,
                                  rem_start,
                                  rem_start + rem_chunks - 1,
                                  nc)
    } else {
      rem_result <- matrix(0, nca, nca)
    }

    # Process main chunks in parallel
    chunk_results <- parallel::parLapply(cl, 1:num_chunks, function(chunk) {
      chunk_start <- (chunk - 1) * chunk_size
      G_chunk <- G[(chunk_start + 1):(chunk_start + chunk_size)]
      AGAmult_chunk(G_chunk, A, chunk_start, chunk_start + chunk_size - 1, nc)
    })

    # Sum all results
    return(Reduce('+', c(chunk_results, list(rem_result))))

  } else {
    return(AGAmult(G, A, K, nc, nca))
  }
}

#' Efficiently Construct U Matrix
#'
#' @param G List of G matrices (\eqn{\textbf{G}})
#' @param A Constraint matrix (\eqn{\textbf{A}})
#' @param K Number of partitions minus 1 (\eqn{K})
#' @param nc Number of columns per partition
#' @param nca Number of constraint columns
#'
#' @return \eqn{\textbf{U}} matrix for constraints
#'
#' @details
#' Computes \eqn{\textbf{U} = \textbf{I} - \textbf{G}\textbf{A}(\textbf{A}^{T}\textbf{G}\textbf{A})^{-1}\textbf{A}^{T}} efficiently, avoiding unnecessary
#' multiplication of blocks of \eqn{\textbf{G}} with all-0 elements.
#'
#' @keywords internal
#' @export
get_U <- function(G, A, K, nc, nca){
  AGAInv <- invert(AGAmult(G, A, K, nc, nca))
  ## [Change 2026-02-12] replaced %**% with crossprod/tcrossprod
  I_minus_U <- t(matmult_U(crossprod(t(A),
                                     crossprod(t(AGAInv), (-t(A)))), G, nc, K))
  return(I_minus_U + diag(nc*(K+1)))
}

#' Generate Grid Indices Without expand.grid()
#'
#' @param vec_list List of vectors to combine
#' @param indices Indices of combinations to return
#'
#' @details
#' Returns selected combinations from the cartesian product of \code{vec_list}
#' without constructing full \code{expand.grid()} for memory efficiency.
#'
#' @return Data frame of selected combinations
#'
#' @keywords internal
#' @export
expgrid <- function(vec_list, indices) {
  # Calculate the total number of combinations
  total_combinations <- prod(sapply(vec_list, length))

  # Check if any indices are out of bounds
  if (any(indices > total_combinations) || any(indices < 1)) {
    stop("Invalid indices: some are out of bounds")
  }

  # Initialize the result list
  result <- vector("list", length(vec_list))
  names(result) <- names(vec_list)

  # Calculate the "stride" for each vector
  strides <- cumprod(c(1, sapply(vec_list[-length(vec_list)], length)))

  # For each requested index
  for (i in indices) {
    # Convert to 0-based index for easier calculation
    idx <- i - 1

    # For each vector in the predictors list
    for (j in seq_along(vec_list)) {
      # Calculate which element of this vector
      # corresponds to the current combination
      element_idx <- (idx %/% strides[j]) %% length(vec_list[[j]]) + 1
      result[[j]] <- c(result[[j]], vec_list[[j]][element_idx])
    }
  }

  # Convert the result to a data frame
  as.data.frame(result)
}

#' Construct Smoothing Spline Penalty Matrix
#'
#' @description
#' Builds penalty matrix combining smoothing spline and ridge penalties with optional
#' predictor/partition-specific components. Handles custom penalties and scaling.
#'
#' @param custom_penalty_mat Matrix; optional custom ridge penalty structure
#' @param L1 Matrix; integrated squared second derivative penalty (\eqn{\textbf{L}_1})
#' @param wiggle_penalty,flat_ridge_penalty Numeric; smoothing and ridge penalty parameters
#' @param K Integer; number of interior knots (\eqn{K})
#' @param nc Integer; number of basis columns per partition
#' @param unique_penalty_per_predictor,unique_penalty_per_partition Logical; enable predictor/partition-specific penalties
#' @param penalty_vec Named numeric; custom penalty values for predictors/partitions
#' @param colnm_expansions Character; column names for linking penalties to predictors
#' @param just_Lambda Logical; return only combined penalty matrix (\eqn{\boldsymbol{\Lambda}})
#'
#' @return List containing:
#' \itemize{
#'   \item Lambda - Combined \eqn{nc \times nc} penalty matrix (\eqn{\boldsymbol{\Lambda}})
#'   \item L1 - Smoothing spline penalty matrix (\eqn{\textbf{L}_1})
#'   \item L2 - Ridge penalty matrix (\eqn{\textbf{L}_2})
#'   \item L_predictor_list - List of predictor-specific penalty matrices (\eqn{\textbf{L}_\text{predictor\_list}})
#'   \item L_partition_list - List of partition-specific penalty matrices (\eqn{\textbf{L}_\text{partition\_list}})
#' }
#'
#' If \code{just_Lambda=TRUE} and no partition penalties, returns only Lambda matrix \eqn{\boldsymbol{\Lambda}}.
#'
#' @keywords internal
#' @export
compute_Lambda <- function(custom_penalty_mat,
                           L1,
                           wiggle_penalty,
                           flat_ridge_penalty,
                           K,
                           nc,
                           unique_penalty_per_predictor,
                           unique_penalty_per_partition,
                           penalty_vec,
                           colnm_expansions,
                           just_Lambda = TRUE){

  ## Custom or flat ridge penalty
  if(any(!is.null(custom_penalty_mat))){
    L2 <- custom_penalty_mat*flat_ridge_penalty
  } else {
    ## By default, = 1 if diagonal and smoothing spline penalty is 0 for
    # corresponding index (adds a ridge to all terms not penalized by the
    # smoothing spline penalty)
    L2 <- diag(nc)*((diag(L1) == 0)*flat_ridge_penalty)
  }

  ## Smoothing penalty = wiggle * L1
  Lambda <- (L1 + L2)*wiggle_penalty
  L_predictor_list <- list()
  L_partition_list <- list()

  ## A unique penalty for each predictor
  if(unique_penalty_per_predictor){
    predictors <- names(penalty_vec)[grep('predictor', names(penalty_vec))]

    ## Scales the elements of the smoothing penalty uniquely by predictor
    for(j in 1:length(predictors)){
      inds <- grep(substr(predictors[j], 10, nchar(predictors[j])), colnm_expansions)
      L <- L1 * penalty_vec[predictors[j]] * wiggle_penalty

      ## All elements of smoothing penalty NOT touching this predictor are 0
      L[-inds,-inds] <- 0
      Lambda <- Lambda + L
      L_predictor_list <- c(L_predictor_list, list(L))
    }
  }

  ## A unique penalty for each partition
  if(unique_penalty_per_partition){
    ## Constructing a unique penalty matrix based on the smoothing penalty
    # scaled uniquely for each partition
    partition_names <- names(penalty_vec)[grep('partition', names(penalty_vec))]
    L_partition_list <- lapply(partition_names, function(part){
      L1*penalty_vec[part]*wiggle_penalty
    })
    names(L_partition_list) <- partition_names
  }
  if(just_Lambda & !unique_penalty_per_partition){
    return(Lambda)
  } else{
    return(list(Lambda = Lambda,
                L1 = L1,
                L2 = L2,
                L_predictor_list = L_predictor_list,
                L_partition_list = L_partition_list))
  }
}

## [Change 2026-02-15] Improved documentation and functionality
#  of compute_G_eigen

#' Compute Eigenvalues and Related Matrices for G
#'
#' @description
#' Computes partition-wise covariance matrices and their matrix square roots
#' from the penalized information matrix. For each partition \eqn{k}, forms
#' \eqn{G_k^{-1} = X_k^T W_k X_k + \Lambda + L_k + S_k}, then extracts
#' \eqn{G_k}, \eqn{G_k^{1/2}}, and (for non-Gaussian or non-identity link)
#' \eqn{G_k^{-1/2}} via eigendecomposition. Eigenvalues at or below zero
#' are clamped, yielding a pseudoinverse treatment of rank-deficient blocks.
#'
#' @param X_gram List of Gram matrices \eqn{X_k^T W_k X_k} by partition.
#' @param Lambda Penalty matrix \eqn{\Lambda} (shared across partitions).
#' @param K Integer; number of interior knots (partitions = \eqn{K+1}).
#' @param parallel Logical; use parallel processing across partitions.
#' @param cl Cluster object from \code{parallel::makeCluster}.
#' @param chunk_size,num_chunks,rem_chunks Partition distribution parameters
#'   for parallel workers.
#' @param family GLM family object. The link function determines whether
#'   \eqn{G^{-1/2}} is computed (needed for non-identity links where
#'   iterative projection requires transforming unconstrained estimates).
#' @param unique_penalty_per_partition Logical; if \code{TRUE}, add
#'   partition-specific penalties from \code{L_partition_list}.
#' @param L_partition_list List of partition-specific penalty matrices.
#' @param keep_G Logical; if \code{TRUE}, return the full \eqn{G} matrix
#'   in addition to the square roots. Set to \code{FALSE} during
#'   intermediate IRWLS iterations where only \eqn{G^{1/2}} and
#'   \eqn{G^{-1/2}} are needed.
#' @param schur_corrections List of Schur complement correction matrices
#'   by partition, accounting for nuisance parameter uncertainty.
#'
#' @details
#' The eigendecomposition of \eqn{G_k^{-1}} yields eigenvalues
#' \eqn{d_i} and eigenvectors \eqn{V}. The matrix powers are then:
#' \deqn{G_k = V \, \mathrm{diag}(1/d_i) \, V^T}
#' \deqn{G_k^{1/2} = V \, \mathrm{diag}(1/\sqrt{d_i}) \, V^T}
#' \deqn{G_k^{-1/2} = V \, \mathrm{diag}(\sqrt{d_i}) \, V^T}
#'
#' Eigenvalues \eqn{d_i \le 0} are clamped: the corresponding entries in
#' \eqn{1/d_i} and \eqn{1/\sqrt{d_i}} are set to zero (pseudoinverse),
#' while the corresponding entries in \eqn{\sqrt{d_i}} are also set to
#' zero.
#'
#' @return A list with components:
#' \describe{
#'   \item{G}{List of \eqn{G_k} matrices by partition (if
#'     \code{keep_G = TRUE}; \code{NULL} entries otherwise).}
#'   \item{Ghalf}{List of \eqn{G_k^{1/2}} matrices by partition.}
#'   \item{GhalfInv}{List of \eqn{G_k^{-1/2}} matrices by partition
#'     (only for non-Gaussian families or non-identity links).}
#' }
#'
#' @keywords internal
#' @export
compute_G_eigen <- function(X_gram,
                            Lambda,
                            K,
                            parallel,
                            cl,
                            chunk_size,
                            num_chunks,
                            rem_chunks,
                            family,
                            unique_penalty_per_partition,
                            L_partition_list,
                            keep_G = TRUE,
                            schur_corrections) {

  ## Determine once whether GhalfInv is needed. For Gaussian identity,
  #  the Lagrangian projection uses Ghalf directly and GhalfInv is not
  #  required. For all other families/links, the projection transforms
  #  unconstrained estimates via G^{-1/2} and GhalfInv must be returned.
  need_GhalfInv <- (paste0(family)[2] != 'identity') |
    (paste0(family)[1] != 'gaussian')

  ## Core computation for a single partition k. Shared by the sequential,
  #  parallel, and remainder code paths to avoid logic duplication.
  compute_one <- function(k) {
    ## Form G_k^{-1} = X^T W X + Lambda [+ L_k] + S_k
    if(unique_penalty_per_partition){
      eig <- tryCatch({
        eigen(X_gram[[k]] +
                Lambda +
                L_partition_list[[k]] +
                schur_corrections[[k]],
              symmetric = TRUE)
      }, error = function(e) NULL)
    } else {
      eig <- tryCatch({
        eigen(X_gram[[k]] +
                Lambda +
                schur_corrections[[k]],
              symmetric = TRUE)
      }, error = function(e) NULL)
    }
    if(is.null(eig)) {
      return(list(G = NULL, Ghalf = NULL, GhalfInv = NULL))
    }

    ## Clamp non-positive eigenvalues for pseudoinverse treatment:
    #  set 1/d_i = 0 for d_i <= 0, so those directions are projected out
    eigen_values <- eig$values
    eigen_values[eig$values <= 0] <- 1
    inv_eigen_values <- 1/eigen_values
    inv_eigen_values[eig$values <= 0] <- 0
    sqrt_inv_eigen_values <- sqrt(inv_eigen_values)

    ## G^{1/2} = V diag(1/sqrt(d_i)) V^T
    Ghalf <- eig$vectors %**% (t(eig$vectors) * sqrt_inv_eigen_values)

    ## G = V diag(1/d_i) V^T (only if requested)
    if(keep_G){
      G <- eig$vectors %**% (t(eig$vectors) * inv_eigen_values)
    } else {
      G <- NULL
    }

    ## G^{-1/2} = V diag(sqrt(d_i)) V^T (only for non-Gaussian/non-identity)
    #  [Change 2026-02-15] Replaced division by sqrt_inv_eigen_values with
    #  multiplication by sqrt_eigen_values, with explicit zero-clamping.
    #  The previous formulation divided by zero when eigenvalues were <= 0,
    #  producing Inf entries and triggering warnings.
    if(need_GhalfInv){
      sqrt_eigen_values <- sqrt(pmax(eig$values, 0))
      GhalfInv <- eig$vectors %**% (t(eig$vectors) * sqrt_eigen_values)
      return(list(G = G,
                  Ghalf = Ghalf,
                  GhalfInv = GhalfInv))
    } else {
      return(list(G = G,
                  Ghalf = Ghalf))
    }
  }

  ## [Change 2026-02-15] Consolidated three near-identical code paths
  # (sequential, parallel main chunks, parallel remainder chunks) into
  #  a single compute_one() helper. The previous implementation had an
  #  inconsistent family check in the remainder path: it tested only
  #  the link (paste0(family)[2] != 'identity') whereas the parallel
  #  and sequential paths tested both family and link. For non-Gaussian
  #  families with identity link, remainder chunks would omit GhalfInv
  #  while main chunks included it, causing downstream extraction to fail
  #  when (K+1) mod chunk_size != 0.

  if(parallel & !is.null(cl)) {
    ## Handle remainder partitions that don't fill a full chunk
    if(rem_chunks > 0) {
      rem_indices <- num_chunks * chunk_size + 1:rem_chunks
      rem <- lapply(rem_indices, compute_one)
    } else {
      rem <- list()
    }

    ## Process main chunks in parallel
    result <- c(
      Reduce("c", parallel::parLapply(cl, 1:num_chunks, function(chunk) {
        inds <- (chunk - 1)*chunk_size + 1:chunk_size
        lapply(inds, compute_one)
      })),
      rem
    )
  } else {
    ## Sequential computation
    result <- lapply(1:(K+1), compute_one)
  }

  ## Reorganize from list-of-partitions to list-of-matrix-types
  if(need_GhalfInv){
    result_processed <- list(
      G = lapply(result, `[[`, "G"),
      Ghalf = lapply(result, `[[`, "Ghalf"),
      GhalfInv = lapply(result, `[[`, "GhalfInv")
    )
  } else {
    result_processed <- list(
      G = lapply(result, `[[`, "G"),
      Ghalf = lapply(result, `[[`, "Ghalf")
    )
  }
  return(result_processed)
}

#' Compute Matrix Square Root Derivative
#'
#' @description
#' Calculates \eqn{d\textbf{G}^{1/2}/d\lambda} matrices for each partition using eigendecomposition.
#' Follows similar approach to \code{compute_G_eigen()} but for matrix derivatives.
#'
#' @param dG_dlambda List of \eqn{nc \times nc} \eqn{d\textbf{G}/d\lambda} matrices by partition
#' @param nc Integer; number of columns per partition
#' @param K Integer; number of interior knots (\eqn{K})
#' @param parallel,cl,chunk_size,num_chunks,rem_chunks Parallel computation parameters
#'
#' @return List of \eqn{nc \times nc} matrices containing \eqn{d\textbf{G}_k^{1/2}/d\lambda} for each partition k
#'
#' @keywords internal
#' @export
compute_dGhalf <- function(dG_dlambda,
                           nc,
                           K,
                           parallel,
                           cl,
                           chunk_size,
                           num_chunks,
                           rem_chunks) {
  if(parallel & !is.null(cl)) {
    ## Handle remainder chunks first
    if(rem_chunks > 0) {
      rem_indices <- num_chunks * chunk_size + 1:rem_chunks
      rem <- lapply(rem_indices, function(k) {
        mat_k <- dG_dlambda[[k]]
        mat_k[is.na(mat_k)] <- 0
        if(!any(!(mat_k == 0))){
          mat_k <- diag(nrow(mat_k))
        }
        eig <- eigen(mat_k, symmetric = TRUE)
        eigen_values <- eig$values
        eigen_values[eig$values <= 0] <- 1
        sqrt_eigen_values <- sqrt(eigen_values)
        sqrt_eigen_values[eig$values <= 0] <- 0
        eig$vectors %**% (t(eig$vectors) * sqrt_eigen_values)
      })
    } else {
      rem <- list()
    }

    ## Process main chunks in parallel
    result <- c(
      Reduce("c",parallel::parLapply(cl, 1:num_chunks, function(chunk) {
        inds <- (chunk - 1)*chunk_size + 1:chunk_size
        lapply(inds, function(k) {
          mat_k <- dG_dlambda[[k]]
          mat_k[is.na(mat_k)] <- 0
          if(!any(!(mat_k == 0))){
            mat_k <- diag(nrow(mat_k))
          }
          eig <- eigen(mat_k, symmetric = TRUE)
          eigen_values <- eig$values
          eigen_values[eig$values <= 0] <- 1
          sqrt_eigen_values <- sqrt(eigen_values)
          sqrt_eigen_values[eig$values <= 0] <- 0
          eig$vectors %**% (t(eig$vectors) * sqrt_eigen_values)
        })
      })),
      rem
    )

  } else {
    ## Sequential computation
    result <- lapply(dG_dlambda, function(mat_k) {
      mat_k[is.na(mat_k)] <- 0
      if(!any(!(mat_k == 0))){
        mat_k <- diag(nrow(mat_k))
      }
      eig <- eigen(mat_k, symmetric = TRUE)
      eigen_values <- eig$values
      eigen_values[eig$values <= 0] <- 1
      sqrt_eigen_values <- sqrt(eigen_values)
      sqrt_eigen_values[eig$values <= 0] <- 0
      eig$vectors %**% (t(eig$vectors) * sqrt_eigen_values)
    })
  }

  return(result)
}

#' Get Centers for Partitioning
#'
#' @param data Matrix of predictor data
#' @param K Number of partitions minus 1 (\eqn{K})
#' @param cluster_args List with custom centers and kmeans args
#' @param cluster_on_indicators Include binary predictors in clustering
#'
#' @details
#' Returns partition centers via:
#'
#' 1. Custom supplied centers if provided as a valid \eqn{K \times q} matrix
#'
#' 2. kmeans clustering on all non-spline variables if \code{cluster_on_indicators=TRUE}
#'
#' 3. kmeans clustering excluding binary variables if \code{cluster_on_indicators=FALSE}
#'
#' @return Matrix of cluster centers
#'
#' @keywords internal
#' @export
get_centers <- function(data, K, cluster_args, cluster_on_indicators) {

  ## If custom centers isn't null, return them
  if(any(!is.na(cluster_args[[1]]))){
    return(cluster_args[[1]])

    ## Partition clusters including 0/1 predictors
  } else if(cluster_on_indicators){
    km <- stats::kmeans(data,
                 K+1,
                 cluster_args[-1])

  } else {
    bin <- which(apply(data, 2, is_binary))
    data[,bin] <- 0
    km <- stats::kmeans(data,
                 K+1,
                 cluster_args[-1])


  }
  return(km$centers)
}

#' Find Neighboring Cluster Partitions Using Midpoint Distance Criterion
#'
#' @description
#' Identifies neighboring partitions by evaluating whether the midpoint between
#' cluster centers is closer to those centers than to any other center.
#'
#' @param centers Matrix; rows are cluster center coordinates
#' @param parallel Logical; use parallel processing
#' @param cl Cluster object for parallel execution
#' @param neighbor_tolerance Numeric; scaling factor for distance comparisons
#'
#' @return List where element i contains indices of centers neighboring center i
#'
#' @keywords internal
#' @export
find_neighbors <- function(centers, parallel, cl, neighbor_tolerance) {
  num_centers <- nrow(centers)
  neighbors <- vector("list", num_centers)

  if (!is.null(cl) & parallel) {
    ## Create index pairs for all comparisons
    pairs <- expand.grid(i = 1:(num_centers-1), j = 2:num_centers)
    pairs <- pairs[pairs$j > pairs$i, ]

    ## Split into chunks
    chunk_size <- max(1, ceiling(nrow(pairs) / length(cl)))
    chunks <- split(1:nrow(pairs), ceiling(seq_along(1:nrow(pairs))/chunk_size))

    ## Export data to cluster
    parallel::clusterExport(cl, "centers", envir = environment())

    ## Process chunks in parallel
    results <- parallel::parLapply(cl, chunks, function(chunk_indices) {
      chunk_results <- vector("list", length(chunk_indices))

      for (idx in seq_along(chunk_indices)) {
        pair_idx <- chunk_indices[idx]
        i <- pairs$i[pair_idx]
        j <- pairs$j[pair_idx]

        ## Check if midpoint is neighbor
        mid <- (centers[i,] + centers[j,])/2
        dist_ij <- sum((mid - centers[i,])^2) / neighbor_tolerance

        distances_to_others <- apply(centers[-c(i,j),,drop=FALSE], 1,
                                     function(k)sum((mid - k)^2))

        if(all(dist_ij < distances_to_others)) {
          chunk_results[[idx]] <- c(i, j)
        }
      }

      return(do.call(rbind, chunk_results[!sapply(chunk_results, is.null)]))
    })

    ## Combine results
    neighbor_pairs <- do.call(rbind, results)

    ## Convert to neighbor list
    if (!is.null(neighbor_pairs)) {
      for (k in 1:nrow(neighbor_pairs)) {
        i <- neighbor_pairs[k, 1]
        j <- neighbor_pairs[k, 2]
        neighbors[[i]] <- c(neighbors[[i]], j)
        neighbors[[j]] <- c(neighbors[[j]], i)
      }
    }

  } else {
    ## Sequential neighbor computation
    for (i in 1:(num_centers-1)) {
      for (j in (i+1):num_centers) {
        ## Get midpoint between centers i and j
        mid <- (centers[i,] + centers[j,])/2

        ## Compute scaled distance from midpoint to center i
        dist_ij <- sum((mid - centers[i,])^2) / neighbor_tolerance

        ## Get distances from midpoint to all other centers
        distances_to_others <- apply(centers[-c(i,j),,drop=FALSE], 1,
                                     function(k)sum((mid - k)^2))

        ## If midpoint closer to i,j than others, they are neighbors
        if(all(dist_ij < distances_to_others)) {
          neighbors[[i]] <- c(neighbors[[i]], j)
          neighbors[[j]] <- c(neighbors[[j]], i)
        }
      }
    }
  }

  return(neighbors)
}


#' Create Data Partitions Using Clustering
#'
#' @description
#' Partitions data support into clusters using Voronoi-like diagrams.
#'
#' @param data Numeric matrix of predictor variables
#' @param cluster_args Parameters for clustering
#' @param cluster_on_indicators Logical to include binary predictors
#' @param K Number of partitions minus 1 (\eqn{K})
#' @param parallel Logical to enable parallel processing
#' @param cl Cluster object for parallel computation
#' @param do_not_cluster_on_these Columns to exclude from clustering
#' @param neighbor_tolerance Scaling factor for neighbor detection
#'
#' @return
#' A list containing:
#' - centers: Cluster center coordinates
#' - knots: Knot points between centers
#' - assign_partition: Function to assign new data to partitions
#' - neighbors: List of neighboring partition indices
#'
#' @keywords internal
#' @export
make_partitions <- function(data,
                            cluster_args,
                            cluster_on_indicators,
                            K,
                            parallel,
                            cl,
                            do_not_cluster_on_these,
                            neighbor_tolerance) {

  q <- ncol(data)

  ## If K is 0, return as is
  if(K == 0){
    return(list(
      centers = matrix()[-1,,drop=FALSE],
      knots = matrix()[-1,,drop=FALSE],
      assign_partition = function(x)1,
      neighbors = list()
    ))
  }

  ## Do not cluster on indicator variables, if desired
  if(!cluster_on_indicators){
    bin <- which(apply(data, 2, is_binary))
    data[,bin] <- 0
  }

  ## Variables we do not want to cluster on
  if(length(do_not_cluster_on_these) > 0){
    data[,do_not_cluster_on_these] <- 0
  }

  ## Compute CVT centers
  initial_points <- get_centers(data,
                                K,
                                cluster_args,
                                cluster_on_indicators = cluster_on_indicators)

  ## Find neighboring partitions
  neighbors <- find_neighbors(initial_points,
                              parallel,
                              cl,
                              neighbor_tolerance)

  ## Function to find knots for neighboring partitions
  find_knots <- function(initial_points, neighbors) {
    additional_points <- list()

    for (i in 1:length(neighbors)) {
      for (j in neighbors[[i]]) {

        ## Do not double count
        if(j > i){
          next
        }

        ## Get centers
        c1 <- initial_points[i,]
        c2 <- initial_points[j,]

        ## Find midpoint
        mid <- (c1 + c2)/2

        ## Record as a constraint
        additional_points[[length(additional_points) + 1]] <-
          mid
        names(additional_points)[length(additional_points)] <-
          paste0(i,'_',j)
      }
    }

    if (length(additional_points) > 0) {
      return(do.call(rbind, additional_points))
    } else {
      return(matrix(nrow = 0, ncol = ncol(initial_points)))
    }
  }

  ## Knot points for smoothing constraints
  knots <- find_knots(initial_points, neighbors)

  ## Create partition assignment function
  assign_partition <- function(new_data) {

    if(!cluster_on_indicators){
      bin <- which(apply(data, 2, is_binary))
      new_data[,bin] <- 0
    }
    if(length(do_not_cluster_on_these) > 0){
      new_data[,do_not_cluster_on_these] <- 0
    }

    ## Calculate batch size (bytes_per_row = 8 bytes * num_centers)
    bytes_per_row <- 8 * nrow(initial_points)
    batch_size <- floor((0.5 * 1024^3) / bytes_per_row)

    ## Initialize results vector
    total_rows <- nrow(new_data)
    assignments <- numeric(total_rows)

    ## Process in batches for saving memory
    for(i in seq(1, total_rows, batch_size)) {
      nd_idx <- min(i + batch_size - 1, total_rows)
      batch <- new_data[i:nd_idx, ]

      ## Assign cluster = nearest-neighbor of the centers
      nn <- FNN::get.knnx(initial_points, batch, k=1)
      assignments[i:nd_idx] <- nn$nn.index
    }

    return(assignments - 0.5)
  }

  ## Return the components: centers, knots, and partition assign fxn
  rownames(initial_points) <- paste0("center_", 1:nrow(initial_points))
  return(list(
    centers = initial_points,
    knots = knots,
    assign_partition = assign_partition,
    neighbors = neighbors
  ))
}

#' Compute Integrated Squared Second Derivative Penalty Matrix for Smoothing Splines
#'
#' @description
#' Generates a penalty matrix representing the integrated squared second derivatives
#' for smoothing spline basis functions, which controls the smoothness of the fitted curve.
#'
#' @param colnm_expansions Character vector of column names for basis expansions
#' @param C Numeric matrix of basis expansions
#' @param power1_cols Indices of linear term columns
#' @param power2_cols Indices of quadratic term columns
#' @param power3_cols Indices of cubic term columns
#' @param power4_cols Indices of quartic term columns
#' @param interaction_single_cols Indices of single interaction columns
#' @param interaction_quad_cols Indices of quadratic interaction columns
#' @param triplet_cols Indices of triplet interaction columns
#' @param nc Number of cubic expansions
#' @param select_cols Optional vector of column indices to penalize (default: all linear terms)
#'
#' @return
#' A symmetric \eqn{p \times p} penalty matrix \eqn{\textbf{P}} representing integrated squared second derivatives
#' for basis expansions in a single partition of the smoothing spline.
#'
#' @details
#' This function computes the analytic form of the traditional integrated, squared, second-derivative evaluated over the bounds of the input data.
#' If \eqn{f(x) = \textbf{X}\boldsymbol{\beta}}, then the penalty is based on \eqn{\int \{ f''(x) \}^2 dx = \boldsymbol{\beta}^{T}(\int \textbf{X}''^{T}\textbf{X}'' dx)\boldsymbol{\beta}}.
#' This function computes the matrix \eqn{\textbf{P} = \int \textbf{X}''^{T}\textbf{X}'' dx}.
#' When scaled by a non-negative scalar (wiggle penalty, predictor penalties and/or partition penalties), this becomes the smoothing spline penalty.
#'
#' @keywords internal
#' @export
get_2ndDerivPenalty <- function(colnm_expansions,
                                C,
                                power1_cols,
                                power2_cols,
                                power3_cols,
                                power4_cols,
                                interaction_single_cols,
                                interaction_quad_cols,
                                triplet_cols,
                                nc,
                                select_cols = NULL){

  mat <- matrix(0, nrow = nc, ncol = nc)
  if(!any(!is.null(select_cols))){
    select_cols <- 1:length(power1_cols)
  }

  ## 2nd derivative penalty matrix,
  # note f(x) = Xb,
  # f''(x) = X''b
  # { f''(x) }^2 = (X''b)^{T}X''b = b^{T}(X''^{T}X'')b = b^{T}[ P ]b
  # below, we are computing the [ P ] matrix entries
  if(length(power1_cols) == 0){
    return(mat)
  } else {
    output <- mat
    for(v in select_cols){
      maxv <- max(C[,power1_cols[v]])
      minv <- min(C[,power1_cols[v]])
      diff1 <- maxv - minv
      diff2 <- maxv^2 - minv^2
      diff3 <- maxv^3 - minv^3
      diff4 <- maxv^4 - minv^4
      diff5 <- maxv^5 - minv^5

      ## univariate penalties
      if(length(power2_cols) > 0){
        mat[power2_cols[v], power2_cols[v]] <- 4*diff1
        # (2)^2 => 4v
      }
      if(length(power3_cols) > 0){
        mat[power3_cols[v], power3_cols[v]] <- 12*diff3
        # (6v)^2 => 12v^3
      }
      if(length(power4_cols) > 0){
        mat[power4_cols[v], power4_cols[v]] <- (144/5)*diff5
        # (12v^2)^2 => (144/5)v^4
      }
      if(length(power2_cols) > 0 & length(power3_cols) > 0){
        mat[power2_cols[v], power3_cols[v]] <- 6*diff2 # (2 * 6v) => 6v^2
        mat[power3_cols[v], power2_cols[v]] <- mat[power2_cols[v],
                                                   power3_cols[v]]
      }
      if(length(power2_cols) > 0 & length(power4_cols) > 0){
        mat[power2_cols[v], power4_cols[v]] <- 8*diff3 # (2 * 12v^2) => 8v^3
        mat[power4_cols[v], power2_cols[v]] <- mat[power2_cols[v],
                                                   power4_cols[v]]
      }
      if(length(power3_cols) > 0 & length(power4_cols) > 0){
        mat[power3_cols[v], power4_cols[v]] <- 18*diff4 # (6v * 12v^2) => 18v^4
        mat[power4_cols[v], power3_cols[v]] <- mat[power3_cols[v],
                                                   power4_cols[v]]
      }

      ## linear by linear interaction penalties if available
      if(length(power1_cols) > 1){
        if(length(interaction_single_cols) > 0){
          interaction_singles <- interaction_single_cols[
            grep(paste0("_",v,"_"),  colnm_expansions[interaction_single_cols])]
          if(length(interaction_singles) > 0){

            mat[interaction_singles, interaction_singles] <-
              diff1 # (1)^2 => v

            if(length(power2_cols) > 0){
              mat[interaction_singles, power2_cols[v]] <-
                2*diff1 # (2 * 1) => 2v
              mat[power2_cols[v], interaction_singles] <-
                mat[interaction_singles, power2_cols[v]]
            }
            if(length(power3_cols) > 0){
              mat[interaction_singles, power3_cols[v]] <-
                3*diff2 # (6v * 1) => 3v^2
              mat[power3_cols[v], interaction_singles] <-
                mat[interaction_singles, power3_cols[v]]
            }
            if(length(power4_cols) > 0){
              mat[interaction_singles, power4_cols[v]] <-
                4*diff3 # (12v^2 * 1) => 4v^3
              mat[power4_cols[v], interaction_singles] <-
                mat[interaction_singles, power4_cols[v]]
            }
          }
        }

        ## linear by quadratic interaction penalties if available
        if(length(interaction_quad_cols) > 0){

          ## interaction quadratic terms, for this variable
          interaction_quads <-
            interaction_quad_cols[grep(paste0("_",v,"_"),
                                       colnm_expansions[interaction_quad_cols])]
          if(length(interaction_quads) > 0){
            for(w in 1:length(power1_cols[-v])){

              ## the other variable, with interactions affecting quadratic terms
              wvar <- c(power1_cols[-v])[w]
              maxw <- max(C[,wvar])
              minw <- min(C[,wvar])
              diffw <- maxw - minw
              diffw2 <- maxw^2 - minw^2

              ## quadratic interaction indices
              interq <- interaction_quads[grep(colnm_expansions[wvar],
                                               colnm_expansions[interaction_quads])]
              if(length(interq) > 0){
                if(length(power2_cols) > 0){
                  ## this is the _w_x_v_^2 term
                  nchv <- nchar(colnm_expansions[power2_cols[v]])
                  interqv2 <- interq[substr(colnm_expansions[interq],
                                            nchar(colnm_expansions[interq]) - nchv + 1,
                                            nchar(colnm_expansions[interq])) ==
                                       colnm_expansions[power2_cols[v]]]
                  ## this is the _v_x_w_^2 term
                  interqv1 <- interq[-which(interq == interqv2)]
                } else {
                  interqv2 <- c()
                  interqv1 <- interq
                }

                ## Some variables might not have the expansions as others
                if(length(interqv2) == 0 | length(interqv1) == 0){
                  next
                }

                ## Compute integrated squared second derivative
                base_val1 <- 4*diffw2*diff1
                base_val2 <- 2*diffw*diff2
                mat[interqv1, interqv1] <- base_val1 # (2w)^2 => 4w^2 * v
                mat[interqv1, interqv2] <- base_val1 +
                  base_val2 # (2w * [2w + 2v]) => 4w^2*v + 2wv^2
                mat[interqv2, interqv1] <- base_val1 +
                  base_val2
                mat[interqv2, interqv2] <- base_val1 +
                  base_val2 +
                  (4/3)*diff3# ([2w + 2v])^2 => 4w^2*v + 2wv^2 + (4/3)v^3


                ## Compute linear-quadratic x quadratic term
                # integrated squared second derivative
                if(length(power2_cols) > 0){
                  base_val1 <- 4*diffw*diff1
                  base_val2 <- 2*diff2
                  mat[power2_cols[v], interqv1] <- base_val1
                  mat[interqv1, power2_cols[v]] <- base_val1
                  mat[power2_cols[v], interqv2] <- base_val1 + base_val2
                  mat[interqv2, power2_cols[v]] <- base_val1 + base_val2
                }
                if(length(power3_cols) > 0){
                  base_val1 <- 6*diffw*diff2
                  base_val2 <- 4*diff3
                  mat[power3_cols[v], interqv1] <- base_val1
                  mat[interqv1, power3_cols[v]] <- base_val1
                  mat[power3_cols[v], interqv2] <- base_val1 + base_val2
                  mat[interqv2, power3_cols[v]] <- base_val1 + base_val2
                }
                if(length(power4_cols) > 0){
                  base_val1 <- 8*diffw*diff3
                  base_val2 <- 6*diff4
                  mat[power4_cols[v], interqv1] <- base_val1
                  mat[interqv1, power4_cols[v]] <- base_val1
                  mat[power4_cols[v], interqv2] <- base_val1 + base_val2
                  mat[interqv2, power4_cols[v]] <- base_val1 + base_val2
                }

                ## Compute linear-linear interaction x linear
                # integrated squared second derivative
                if(length(interaction_single_cols) > 0){
                  interaction_singles <-
                    interaction_single_cols[grep(paste0("_",v,"_"),
                                                 colnm_expansions[interaction_single_cols])]
                  if(length(interaction_singles) > 0){
                    for(j in 1:length(interaction_singles)){
                      base_val1 <- 2*diffw*diff1
                      base_val2 <- (2/3)*diff3
                      mat[interaction_singles[j], interqv1] <- base_val1
                      mat[interqv1, interaction_singles[j]] <- base_val1
                      mat[interaction_singles[j], interqv2] <- base_val1 +
                        base_val2
                      mat[interqv2, interaction_singles[j]] <- base_val1 +
                        base_val2

                    }
                  }
                }

                ## Compute 3-way interaction x interaction
                # integrated squared second derivatives
                if(length(triplet_cols) > 0){
                  triplets <- triplet_cols[grep(paste0("_",v,"_"),
                                                colnm_expansions[triplet_cols])]
                  if(length(triplets) > 0){
                    for(tr in 1:length(triplets)){
                      other2_vars <- unlist(strsplit(colnm_expansions[triplets[tr]],
                                                     'x'))
                      other2_vars <- other2_vars[other2_vars !=
                                                   colnm_expansions[power1_cols[v]]]
                      v1 <- C[,other2_vars[1]]
                      v2 <- C[,other2_vars[2]]
                      diffv1 <- max(v1) - min(v1)
                      diffv2 <- max(v2) - min(v2)
                      if(other2_vars[1] == colnm_expansions[wvar]){
                        base_val1 <- 2*(diffw2 + 2*diffw*diffv2)*diff1
                        base_val2 <- 2*(diffw + diffv2)*diff2
                        mat[triplets[tr], interqv1] <- base_val1
                        mat[interqv1, triplets[tr]] <- base_val1
                        mat[triplets[tr], interqv2] <- base_val1 + base_val2
                        mat[interqv2, triplets[tr]] <- base_val1 + base_val2
                      } else if(other2_vars[2] == colnm_expansions[wvar]){
                        base_val1 <- 2*(diffw2 + 2*diffw*diffv1)*diff1
                        base_val2 <- 2*(diffw + diffv1)*diff2
                        mat[triplets[tr], interqv1] <- base_val1
                        mat[interqv1, triplets[tr]] <- base_val1
                        mat[triplets[tr], interqv2] <- base_val1 + base_val2
                        mat[interqv2, triplets[tr]] <- base_val1 + base_val2
                      } else {
                        base_val1 <- 2*diffw*(diffv1 + diffv2)*diff1
                        base_val2 <- 2*(diffv1 + diffv2)*diff2
                        mat[triplets[tr], interqv1] <- base_val1
                        mat[interqv1, triplets[tr]] <- base_val1
                        mat[triplets[tr], interqv2] <- base_val1 + base_val2
                        mat[interqv2, triplets[tr]] <- base_val1 + base_val2
                      }
                    }
                  }
                }
              }
            }
          }
        }

        ### Three-way Interaction Terms (vwu) ###
        # Second derivative is w + u
        # Cross terms with all other basis functions
        # Cases depend on which variable is w
        if(length(triplet_cols) > 0) {
          triplets <-
            triplet_cols[grep(paste0("_", v, "_"), colnm_expansions[triplet_cols])]
          if (length(triplets) > 0) {
            other2_vars <- lapply(triplets, function(tr) {
              vars <- unlist(strsplit(colnm_expansions[tr], 'x'))
              vars[vars != colnm_expansions[power1_cols[v]]]
            })
            for (tr in 1:length(other2_vars)) {

              ## The first of two "other" variables, of 3-way interaction
              maxw <- max(C[,other2_vars[[tr]][1]])
              minw <- min(C[,other2_vars[[tr]][1]])
              diffw <- maxw - minw
              diffw2 <- maxw^2 - minw^2

              ## The second of two "other" variables, of 3-way interaction
              maxu <- max(C[,other2_vars[[tr]][2]])
              minu <- min(C[,other2_vars[[tr]][2]])
              diffu <- maxu - minu
              diffu2 <- maxu^2 - minu^2

              ## Adapt this code to handle 3-way terms, i.e.
              ## the second derivative of vwu with respect to v is (w + u)
              ## integral for diagonal term =
              # int^{v = maxv}_{v = minv} w + u dv du dw => (w + u)*v
              trip_inter <- intersect(intersect(
                grep(other2_vars[[tr]][1], colnm_expansions),
                grep(other2_vars[[tr]][2], colnm_expansions)),
                triplets
              )
              mat[trip_inter, trip_inter] <-
                (diffw2 + diffu2 + 2*diffw*diffu)*diff1

              if(length(power2_cols) > 0){
                mat[power2_cols[v], trip_inter] <- 2*diff1*(diffu + diffw)
                mat[trip_inter, power2_cols[v]] <-
                  mat[power2_cols[v], trip_inter]
              }
              if(length(power3_cols) > 0){
                mat[power3_cols[v], trip_inter] <- 3*diff2*(diffu + diffw)
                mat[trip_inter, power3_cols[v]] <-
                  mat[power3_cols[v], trip_inter]
              }
              if(length(power4_cols) > 0){
                mat[power4_cols[v], trip_inter] <- 4*diff3*(diffu + diffw)
                mat[trip_inter, power4_cols[v]] <-
                  mat[power4_cols[v], trip_inter]
              }
              if(length(interaction_single_cols) > 0){
                interaction_singles <-
                  interaction_single_cols[grep(paste0("_",v,"_"),
                                               colnm_expansions[interaction_single_cols])]
                if(length(interaction_singles) > 0){
                  for(j in 1:length(interaction_singles)){
                    mat[interaction_singles[j], trip_inter] <-
                      (diffu + diffw)*diff1  # (u+w)*1 => (u+w)v
                    mat[trip_inter, interaction_singles[j]] <-
                      mat[interaction_singles[j], trip_inter]
                  }
                }
              }
            }
          }
        }
      }
      ## Elementwise-sum the output for other predictors to the matrix here
      output <- output + mat
    }
    return(output)
  }
}

#' Wrapper for Smoothing Spline Penalty Computation
#'
#' @description
#' Computes smoothing spline penalty matrix with optional parallel processing.
#' Calls  \code{\link{get_2ndDerivPenalty}} after
#' processing spline vs. nonspline terms and preparing for parallel if desired.
#'
#' @param K Number of partitions (\eqn{K+1})
#' @param colnm_expansions Column names of basis expansions
#' @param C Basis expansion matrix of two rows, first of all maximums, second of all minimums, for all variables of interest = \code{rbind(apply(C, 2, max)), rbind(apply(C, 2, min)))} for cubic expansions "C"
#' @param power1_cols Linear term columns
#' @param power2_cols Quadratic term columns
#' @param power3_cols Cubic term columns
#' @param power4_cols Quartic term columns
#' @param interaction_single_cols Single interaction columns
#' @param interaction_quad_cols Quadratic interaction columns
#' @param triplet_cols Triplet interaction columns
#' @param nonspline_cols Predictors not treated as spline effects
#' @param nc Number of cubic expansions
#' @param parallel Logical to enable parallel processing
#' @param cl Cluster object for parallel computation
#'
#' @return
#' A \eqn{p \times p} penalty matrix for smoothing spline regularization containing the
#' elementwise sum of the integrated squared second derivative of the fitted
#' function with respect to predictors of interest.
#'
#' Function is exported for reference purposes - use at your own risk!
#'
#' @keywords internal
#' @export
get_2ndDerivPenalty_wrapper <- function(K,
                                        colnm_expansions,
                                        C,
                                        power1_cols,
                                        power2_cols,
                                        power3_cols,
                                        power4_cols,
                                        interaction_single_cols,
                                        interaction_quad_cols,
                                        triplet_cols,
                                        nonspline_cols,
                                        nc,
                                        parallel,
                                        cl) {
  ## Modification such that we can get the same operations performed for
  # nonspline terms too without affecting the rest of lgspline
  colnm_expansions_og <- colnm_expansions
  if(length(nonspline_cols) > 0){
    for(jj in 1:length(nonspline_cols)){
      power1_cols <- c(power1_cols, nonspline_cols[jj])
      ## For each power, detect if already present for spline effects
      # If so, append a 0-column for the categorical variable
      # Else, skip
      if(length(power2_cols) > 0){
        colnm_expansions <- c(colnm_expansions, paste0(colnm_expansions[nonspline_cols[jj]],'^2'))
        C <- cbind(C, 0)
        power2_cols <- c(power2_cols, ncol(C))
      }
      if(length(power3_cols) > 0){
        colnm_expansions <- c(colnm_expansions, paste0(colnm_expansions[nonspline_cols[jj]],'^3'))
        C <- cbind(C, 0)
        power3_cols <- c(power3_cols, ncol(C))
      }
      if(length(power4_cols) > 0){
        colnm_expansions <- c(colnm_expansions, paste0(colnm_expansions[nonspline_cols[jj]],'^4'))
        C <- cbind(C, 0)
        power4_cols <- c(power4_cols, ncol(C))
      }
    }
    ## Update colnames and number of columns of expansions in C with
    # new nonspline power terms
    colnames(C) <- colnm_expansions
    nc <- ncol(C)
  }

  ## If parallel processing
  if(parallel & (K > 1)){
    ## Determine chunk size based on cluster length
    chunk_size <- max(1, floor(2 * length(cl)))

    ## Total number of columns to process
    total_cols <- length(power1_cols)

    ## Initialize result matrix
    result <- matrix(0, nrow = nc, ncol = nc)

    ## Process in chunks
    for(start in seq(1, total_cols, by = chunk_size)) {
      ## Determine end of current chunk
      end <- min(start + chunk_size - 1, total_cols)

      ## Process current chunk in parallel
      chunk_result <- Reduce("+",
                             parallel::parLapply(cl,
                                                 start:end,
                                                 function(select_col) {
                                                   get_2ndDerivPenalty(colnm_expansions,
                                                                       C,
                                                                       power1_cols,
                                                                       power2_cols,
                                                                       power3_cols,
                                                                       power4_cols,
                                                                       interaction_single_cols,
                                                                       interaction_quad_cols,
                                                                       triplet_cols,
                                                                       nc,
                                                                       select_col)
                                                 }))
      ## Add chunk result to overall result
      result <- result + chunk_result
    }
  } else {
    ## Otherwise, compute serial
    result <- get_2ndDerivPenalty(colnm_expansions,
                                  C,
                                  power1_cols,
                                  power2_cols,
                                  power3_cols,
                                  power4_cols,
                                  interaction_single_cols,
                                  interaction_quad_cols,
                                  triplet_cols,
                                  nc)
  }
  colnames(result) <- colnm_expansions
  rownames(result) <- colnm_expansions
  ## Isolate the entries excluding appended
  result <- result[colnm_expansions_og, colnm_expansions_og]
  return(result)
}

#' Compute Log-Likelihood for Weibull Accelerated Failure Time Model
#'
#' @description
#' Calculates the log-likelihood for a Weibull accelerated failure time (AFT)
#' survival model, supporting right-censored survival data.
#'
#' @param log_y Numeric vector of logarithmic response/survival times
#' @param log_mu Numeric vector of logarithmic predicted survival times
#' @param status Numeric vector of censoring indicators
#'   (1 = event, 0 = censored)
#'   Indicates whether an event of interest occurred (1) or the observation was
#'   right-censored (0). In survival analysis, right-censoring occurs when the
#'   full survival time is unknown, typically because the study ended or the
#'   subject was lost to follow-up before the event of interest occurred.
#' @param scale Numeric scalar representing the Weibull scale parameter
#'   (sigma), equivalent to \code{survreg$scale}. This is the square root of
#'   the dispersion stored in \code{lgspline$sigmasq_tilde}.
#' @param weights Optional numeric vector of observation weights (default = 1)
#'
#' @return
#' A numeric scalar representing the total log-likelihood of the model
#'
#' @details
#' The function computes log-likelihood contributions for a Weibull AFT model,
#' explicitly accounting for right-censored observations. It supports optional
#' observation weighting to accommodate complex sampling designs.
#'
#' This both provides a tool for actually fitting Weibull AFT models, and
#' boilerplate code for users who wish to incorporate Lagrangian multiplier
#' smoothing splines into their own custom models.
#'
#' @examples
#'
#' ## Minimal example of fitting a Weibull Accelerated Failure Time model
#' # Simulating survival data with right-censoring
#' set.seed(1234)
#' x1 <- rnorm(1000)
#' x2 <- rbinom(1000, 1, 0.5)
#' yraw <- rexp(exp(0.01*x1 + 0.01*x2))
#' # status: 1 = event occurred, 0 = right-censored
#' status <- rbinom(1000, 1, 0.25)
#' yobs <- ifelse(status, runif(1, 0, yraw), yraw)
#' df <- data.frame(
#'   y = yobs,
#'   x1 = x1,
#'   x2 = x2
#' )
#'
#' ## Fit model using lgspline with Weibull AFT specifics
#' model_fit <- lgspline(y ~ spl(x1) + x2,
#'                       df,
#'                       unconstrained_fit_fxn = unconstrained_fit_weibull,
#'                       family = weibull_family(),
#'                       need_dispersion_for_estimation = TRUE,
#'                       dispersion_function = weibull_dispersion_function,
#'                       glm_weight_function = weibull_glm_weight_function,
#'                       schur_correction_function = weibull_schur_correction,
#'                       status = status,
#'                       opt = FALSE,
#'                       K = 1)
#'
#' loglik_weibull(log(model_fit$y), log(model_fit$ytilde), status,
#'   sqrt(model_fit$sigmasq_tilde))
#'
#' @export
loglik_weibull <- function(log_y, log_mu, status, scale, weights = 1) {

  ## Log-likelihood contributions
  ## scale = sigma (Weibull scale), NOT sigma^2 (dispersion)
  z <- (log_y - log_mu) / scale
  logL <- status * (-log(scale) +
                      z -
                      log_y) -
    exp(z)

  ## Return sum of log likelihood
  return(sum(logL * weights))
}


#' Compute Gradient of Log-Likelihood of Weibull Accelerated Failure Model
#'
#' @description
#' Calculates the gradient of log-likelihood for a Weibull accelerated failure
#' time (AFT) survival model, supporting right-censored survival data.
#'
#' @param X Design matrix
#' @param y Response vector
#' @param mu Predicted mean vector
#' @param order_list List of observation indices per partition
#' @param dispersion Dispersion parameter (sigma^2 = scale^2). The lgspline
#'   framework stores and passes dispersion (sigma^2); the Weibull scale
#'   (sigma) matching \code{survreg$scale} is \code{sqrt(dispersion)}.
#' @param VhalfInv Inverse square root of correlation matrix (if applicable)
#' @param observation_weights Observation weights
#' @param status Censoring indicator (1 = event, 0 = censored)
#'
#' @return
#' A numeric vector representing the gradient with respect to coefficients.
#'
#' @details
#' Needed if using "blockfit", correlation structures, or quadratic programming
#' with Weibull AFT models.
#'
#' The gradient is computed on a scale that omits the 1/sigma prefactor.
#' Specifically, the true score is (1/sigma) * X^T diag(w) (exp(z) - status),
#' but both this function and the corresponding information matrix used
#' internally omit 1/sigma and 1/sigma^2 respectively, so the Newton-Raphson
#' step G*u remains correct. This matches the convention in
#' \code{\link{unconstrained_fit_weibull}}.
#'
#' @examples
#'
#' set.seed(1234)
#' t1 <- rnorm(1000)
#' t2 <- rbinom(1000, 1, 0.5)
#' yraw <- rexp(exp(0.01*t1 + 0.01*t2))
#' status <- rbinom(1000, 1, 0.25)
#' yobs <- ifelse(status, runif(1, 0, yraw), yraw)
#' df <- data.frame(
#'   y = yobs,
#'   t1 = t1,
#'   t2 = t2
#' )
#'
#' ## Example using blockfit for t2 as a linear term - output does not look
#' # different, but internal methods used for fitting change
#' model_fit <- lgspline(y ~ spl(t1) + t2,
#'                       df,
#'                       unconstrained_fit_fxn = unconstrained_fit_weibull,
#'                       family = weibull_family(),
#'                       need_dispersion_for_estimation = TRUE,
#'                       qp_score_function = weibull_qp_score_function,
#'                       dispersion_function = weibull_dispersion_function,
#'                       glm_weight_function = weibull_glm_weight_function,
#'                       schur_correction_function = weibull_schur_correction,
#'                       K = 1,
#'                       blockfit = TRUE,
#'                       opt = FALSE,
#'                       status = status,
#'                       verbose = TRUE)
#'
#' print(summary(model_fit))
#'
#' @export
weibull_qp_score_function = function(X,
                                     y,
                                     mu,
                                     order_list,
                                     dispersion,
                                     VhalfInv,
                                     observation_weights,
                                     status){
  ## dispersion = sigma^2; scale = sigma = sqrt(dispersion)
  scale <- sqrt(dispersion)
  order_indices <- unlist(order_list)
  ## Score of beta (scaled by sigma^2 for numerical convenience).
  # True score: (1/sigma) * X^T diag(w) (exp(z) - status)
  #  Scaled score (x sigma^2): sigma * X^T diag(w) (exp(z) - status)
  #  This matches the convention in unconstrained_fit_weibull, where
  #  the information is also unscaled (missing 1/sigma^2), so the
  #  Newton-Raphson step G*u remains correct.
  crossprod(X , cbind(
    observation_weights *
    (exp((log(y) - log(mu)) / scale) - status[order_indices]) *
    scale)
  )
}

#' Correction for the Variance-Covariance Matrix for Uncertainty in Scale
#'
#' @description
#' Computes the Schur complement \eqn{\textbf{S}} such that
#' \eqn{\textbf{G}^* = (\textbf{G}^{-1} + \textbf{S})^{-1}} properly
#' accounts for uncertainty in estimating the Weibull scale parameter when
#' estimating the variance-covariance matrix. Otherwise, the
#' variance-covariance matrix is optimistic and assumes the scale is known,
#' when it was in fact estimated. Note that the parameterization adds the
#' output of this function elementwise (not subtract) so for most cases, the
#' output of this function will be negative or a negative
#' definite/semi-definite matrix.
#'
#' @param X Block-diagonal matrices of spline expansions
#' @param y Block-vector of response
#' @param B Block-vector of coefficient estimates
#' @param dispersion Scalar, estimate of dispersion (sigma^2 = scale^2). The
#'   lgspline framework stores and passes dispersion (sigma^2); the Weibull
#'   scale (sigma) matching \code{survreg$scale} is \code{sqrt(dispersion)}.
#' @param order_list List of partition orders
#' @param K Number of partitions minus 1 (\eqn{K})
#' @param family Distribution family
#' @param observation_weights Optional observation weights (default = 1)
#' @param status Censoring indicator (1 = event, 0 = censored)
#'   Indicates whether an event of interest occurred (1) or the observation was
#'   right-censored (0). In survival analysis, right-censoring occurs when the
#'   full survival time is unknown, typically because the study ended or the
#'   subject was lost to follow-up before the event of interest occurred.
#'
#' @return
#' List of \eqn{p \times p} matrices representing the Schur-complement
#' corrections \eqn{\textbf{S}_k} to be elementwise added to each block of the
#' information matrix, before inversion.
#'
#' @details
#' Adjusts the variance-covariance matrix unscaled for coefficients to account
#' for uncertainty in estimating the Weibull scale parameter, that otherwise
#' would be lost if simply using
#' \eqn{\textbf{G}=(\textbf{X}^{T}\textbf{W}\textbf{X} + \textbf{L})^{-1}}.
#' This is accomplished using a correction based on the Schur complement so we
#' avoid having to construct the entire variance-covariance matrix, or
#' modifying the procedure for \code{\link{lgspline}} substantially.
#' For any model with nuisance parameters that must have uncertainty accounted
#' for, this tool will be helpful.
#'
#' This both provides a tool for actually fitting Weibull accelerated failure
#' time (AFT) models, and boilerplate code for users who wish to incorporate
#' Lagrangian multiplier smoothing splines into their own custom models.
#'
#' @examples
#'
#' ## Minimal example of fitting a Weibull Accelerated Failure Time model
#' # Simulating survival data with right-censoring
#' set.seed(1234)
#' t1 <- rnorm(1000)
#' t2 <- rbinom(1000, 1, 0.5)
#' yraw <- rexp(exp(0.01*t1 + 0.01*t2))
#' # status: 1 = event occurred, 0 = right-censored
#' status <- rbinom(1000, 1, 0.25)
#' yobs <- ifelse(status, runif(1, 0, yraw), yraw)
#' df <- data.frame(
#'   y = yobs,
#'   t1 = t1,
#'   t2 = t2
#' )
#'
#' ## Fit model using lgspline with Weibull Schur correction
#' model_fit <- lgspline(y ~ spl(t1) + t2,
#'                       df,
#'                       unconstrained_fit_fxn = unconstrained_fit_weibull,
#'                       family = weibull_family(),
#'                       need_dispersion_for_estimation = TRUE,
#'                       dispersion_function = weibull_dispersion_function,
#'                       glm_weight_function = weibull_glm_weight_function,
#'                       schur_correction_function = weibull_schur_correction,
#'                       status = status,
#'                       opt = FALSE,
#'                       K = 1)
#'
#' print(summary(model_fit))
#'
#' ## Fit model using lgspline without Weibull Schur correction
#' naive_fit <- lgspline(y ~ spl(t1) + t2,
#'                       df,
#'                       unconstrained_fit_fxn = unconstrained_fit_weibull,
#'                       family = weibull_family(),
#'                       need_dispersion_for_estimation = TRUE,
#'                       dispersion_function = weibull_dispersion_function,
#'                       glm_weight_function = weibull_glm_weight_function,
#'                       status = status,
#'                       opt = FALSE,
#'                       K = 1)
#'
#' print(summary(naive_fit))
#'
#' @export
weibull_schur_correction <- function(X,
                                    y,
                                    B,
                                    dispersion,
                                    order_list,
                                    K,
                                    family,
                                    observation_weights,
                                    status){
  ## dispersion = sigma^2; scale = sigma = sqrt(dispersion)
  scale <- sqrt(dispersion)
  lapply(1:(K+1), function(k){
    if(nrow(X[[k]]) < 1){
      return(0)
    } else {
      mu <- family$linkinv(c(X[[k]] %**% B[[k]]))
      s <- status[order_list[[k]]]
      obs <- y[[k]]
      z <- (log(obs) - log(mu)) / scale
      exp_z <- exp(z)
      zexp_z <- z * exp_z
      weights <- c(observation_weights[[k]])

      ## Correction via Schur complement
      # Extract true conditional variance-covariance of beta coefficients
      # conditional upon estimate of scale (sigma).
      # I = ( I_bb I_bs^{T} )
      #     ( I_bs I_ss     )
      # for b = beta, s = scale (sigma)
      # Note that: I_bb = invert(G[[k]]) is incorrect, I_bb is part of
      # Fisher info
      I_bs <- t(X[[k]]) %**% cbind(weights * zexp_z * scale)
      I_ss <- -sum(
        weights * (
          (s + 2*s*z + zexp_z + exp_z * z^2)
        )
      )
      # compl gets elementwise added to G[[k]] for all k = 1...K+1
      compl <- I_bs %**% matrix(-1/I_ss) %**% t(I_bs)
      # Schur complement correction to pass on to compute_G_eigen()
      return(compl)
    }
  })
}

#' @rdname weibull_schur_correction
#' @export
weibull_shur_correction <- weibull_schur_correction

#' Estimate Scale for Weibull Accelerated Failure Time Model
#'
#' @description
#' Computes maximum log-likelihood scale estimate (sigma) for a Weibull
#' accelerated failure time (AFT) survival model.
#'
#' This both provides a tool for actually fitting Weibull AFT Models, and
#' boilerplate code for users who wish to incorporate Lagrangian multiplier
#' smoothing splines into their own custom models.
#'
#' @param log_y Logarithm of response/survival times
#' @param log_mu Logarithm of predicted survival times
#' @param status Censoring indicator (1 = event, 0 = censored)
#'   Indicates whether an event of interest occurred (1) or the observation was
#'   right-censored (0). In survival analysis, right-censoring occurs when the
#'   full survival time is unknown, typically because the study ended or the
#'   subject was lost to follow-up before the event of interest occurred.
#' @param weights Optional observation weights (default = 1)
#'
#' @return
#' Scalar representing the estimated Weibull scale (sigma), equivalent to
#' \code{survreg$scale}. The dispersion (as stored in
#' \code{lgspline$sigmasq_tilde}) is sigma^2.
#'
#' @details
#' Calculates maximum log-likelihood estimate of scale (sigma) for Weibull AFT
#' model accounting for right-censored observations using Brent's method for
#' optimization.
#'
#' @examples
#'
#' ## Simulate exponential data with censoring
#' set.seed(1234)
#' mu <- 2  # mean of exponential distribution
#' n <- 500
#' y <- rexp(n, rate = 1/mu)
#'
#' ## Introduce censoring (25% of observations)
#' status <- rbinom(n, 1, 0.75)
#' y_obs <- ifelse(status, y, NA)
#'
#' ## Compute scale estimate
#' scale_est <- weibull_scale(
#'   log_y = log(y_obs[!is.na(y_obs)]),
#'   log_mu = log(mu),
#'   status = status[!is.na(y_obs)]
#' )
#'
#' print(scale_est)
#'
#'
#' @export
weibull_scale <- function(log_y, log_mu, status, weights = 1){
  optim(
    1,
    fn = function(par){
      ## par is scale (sigma), passed directly to loglik_weibull
      -loglik_weibull(log_y, log_mu, status, par, weights)
    },
    method = 'Brent',
    lower = 1e-64,
    upper = 100
  )$par
}

## [Change 2026-02-14]: weibull_family -- add $aic, $loglik, $dev.resids methods
# for logLik.lgspline compatibility

#' Weibull Family for Survival Model Specification
#'
#' @description
#' Creates a compatible family object for Weibull accelerated failure time (AFT)
#' models with customizable tuning options.
#'
#' This both provides a tool for actually fitting Weibull AFT Models, and
#' boilerplate code for users who wish to incorporate Lagrangian multiplier
#' smoothing splines into their own custom models.
#'
#' @return
#' A list containing family-specific components for survival model estimation
#'
#' @details
#' Provides a comprehensive family specification for Weibull AFT models,
#' including family name, link function, inverse link function, custom loss
#' function for model tuning, and methods for AIC and log-likelihood
#' computation compatible with \code{logLik.lgspline}.
#'
#' Supports right-censored survival data with flexible parameter estimation.
#'
#' Note on scale vs. dispersion: throughout this package, the lgspline object
#' stores \code{sigmasq_tilde} which equals sigma^2 (dispersion), where
#' sigma is the Weibull scale parameter matching \code{survreg$scale}.
#' Functions that accept a \code{dispersion} argument receive sigma^2;
#' functions that accept a \code{scale} argument receive sigma.
#'
#' @examples
#'
#' ## Simulate survival data with covariates
#' set.seed(1234)
#' n <- 1000
#' t1 <- rnorm(n)
#' t2 <- rbinom(n, 1, 0.5)
#'
#' ## Generate survival times with Weibull-like structure
#' lambda <- exp(0.5 * t1 + 0.3 * t2)
#' yraw <- rexp(n, rate = 1/lambda)
#'
#' ## Introduce right-censoring
#' status <- rbinom(n, 1, 0.75)
#' y <- ifelse(status, yraw, runif(1, 0, yraw))
#'
#' ## Prepare data
#' df <- data.frame(y = y, t1 = t1, t2 = t2, status = status)
#'
#' ## Fit model using custom Weibull family
#' model_fit <- lgspline(y ~ spl(t1) + t2,
#'                       df,
#'                       unconstrained_fit_fxn = unconstrained_fit_weibull,
#'                       family = weibull_family(),
#'                       need_dispersion_for_estimation = TRUE,
#'                       dispersion_function = weibull_dispersion_function,
#'                       glm_weight_function = weibull_glm_weight_function,
#'                       schur_correction_function = weibull_schur_correction,
#'                       status = status,
#'                       opt = FALSE,
#'                       K = 1)
#'
#' summary(model_fit)
#'
#' ## Log-likelihood now works via logLik.lgspline:
#' # logLik(model_fit)
#'
#' @export
weibull_family <- function() list(
  family = "weibull",
  link = "log",
  linkfun = log,
  linkinv = exp,

  ## [Change 2026-02-12] Added dev.resids for base R
  #  compatibility when logLik.lgspline checks for family$dev.resids
  dev.resids = function(y, mu, wt) {
    ## Deviance residuals for Weibull AFT: -2 * log-likelihood
    #  This is a simplified version; full version uses status via
    #  custom_dev.resids below
    wt * (y/mu - log(y/mu) - 1) * 2
  },

  ## Custom loss used in place of MSE for computing GCV
  custom_dev.resids =
    function(y, mu, order_indices, family, observation_weights, status){
      log_mu <- log(mu)
      log_y <- log(y)
      status <- status[order_indices]

      ## Initialize scale (sigma) using null model
      init_scale <-
        weibull_scale(log_y, mean(log_y), status,
                      observation_weights)
      ## Find scale (sigma) given current mu
      scale <- optim(
        init_scale,
        fn = function(par){
          ## par is scale (sigma)
          -loglik_weibull(log_y, log_mu, status, par, observation_weights)
        },
        lower = init_scale/10, # [Change 2026-02-17] Increased from 5-10
        upper = init_scale*10,
        method = 'Brent'
      )$par

      ## -2 * log-likelihood
      dev <- -2*(
        status * (-log(scale) +
                    (1/scale - 1)*log_y -
                    log_mu/scale) -
          (exp((log_y - log_mu)/scale))
      )
      return(dev * observation_weights)
    },

  ## [Change 2026-02-12] Added $loglik method for
  #  logLik.lgspline compatibility
  #  Accepts a fitted lgspline model object and returns the scalar
  #  log-likelihood value
  loglik = function(model_fit) {
    ## Extract necessary components from the fitted model
    log_y <- log(model_fit$y)
    log_mu <- log(model_fit$ytilde)
    ## sigmasq_tilde is dispersion (sigma^2); take sqrt to get scale (sigma)
    scale <- sqrt(model_fit$sigmasq_tilde)

    ## Status must be stored in model_fit$extra or model_fit$call_args
    #  Try common locations
    status <- NULL
    if(!is.null(model_fit$status)) {
      status <- model_fit$status
    } else if(!is.null(model_fit$extra_args$status)) {
      status <- model_fit$extra_args$status
    } else {
      ## If status not found, assume all events observed
      warning("Censoring status not found in model object; ",
              "assuming all events observed (status = 1).")
      status <- rep(1, length(log_y))
    }

    weights <- if(!is.null(model_fit$observation_weights)){
      model_fit$observation_weights
    } else {
      rep(1, length(log_y))
    }

    loglik_weibull(log_y, log_mu, status, scale, weights)
  },

  ##  [Change 2026-02-12] Added $aic method for
  #  logLik.lgspline compatibility
  #  Returns AIC = -2*loglik + 2*edf
  #  edf = effective degrees of freedom (trace of hat matrix + 1 for scale)
  aic = function(model_fit) {
    ll <- weibull_family()$loglik(model_fit)
    ## edf: trace of smoother matrix + 1 for scale parameter
    edf <- model_fit$edf
    if(is.null(edf)){
      ## Fallback: use trace_XUGX if available
      edf <- model_fit$trace_XUGX
    }
    if(is.null(edf)){
      warning("Cannot compute AIC: effective degrees of freedom not found.")
      return(NA)
    }
    ## +1 for the scale parameter
    -2 * ll + 2 * (edf + 1)
  }
)


#' Estimate Weibull Dispersion for Accelerated Failure Time Model
#'
#' @description
#' Computes the dispersion parameter (sigma^2 = scale^2) for a Weibull
#' accelerated failure time (AFT) model, supporting right-censored survival
#' data. The returned value is sigma^2, where sigma is the Weibull scale
#' parameter matching \code{survreg$scale}.
#'
#' This both provides a tool for actually fitting Weibull AFT Models, and
#' boilerplate code for users who wish to incorporate Lagrangian multiplier
#' smoothing splines into their own custom models.
#'
#' @param mu Predicted survival times
#' @param y Observed response/survival times
#' @param order_indices Indices to align status with response
#' @param family Weibull AFT model family specification
#' @param observation_weights Optional observation weights
#' @param status Censoring indicator (1 = event, 0 = censored)
#'   Indicates whether an event of interest occurred (1) or the observation was
#'   right-censored (0). In survival analysis, right-censoring occurs when the
#'   full survival time is unknown, typically because the study ended or the
#'   subject was lost to follow-up before the event of interest occurred.
#'
#' @return
#' Dispersion estimate (sigma^2) for the Weibull AFT model, i.e., the squared
#' scale parameter. The Weibull scale (sigma) matching \code{survreg$scale} is
#' \code{sqrt()} of this value.
#'
#' @seealso \code{\link{weibull_scale}} for the underlying scale estimation
#'   function
#'
#' @examples
#'
#' ## Simulate survival data with covariates
#' set.seed(1234)
#' n <- 1000
#' t1 <- rnorm(n)
#' t2 <- rbinom(n, 1, 0.5)
#'
#' ## Generate survival times with Weibull-like structure
#' lambda <- exp(0.5 * t1 + 0.3 * t2)
#' yraw <- rexp(n, rate = 1/lambda)
#'
#' ## Introduce right-censoring
#' status <- rbinom(n, 1, 0.75)
#' y <- ifelse(status, yraw, runif(1, 0, yraw))
#'
#' ## Example of using dispersion function
#' mu <- mean(y)
#' order_indices <- seq_along(y)
#' weights <- rep(1, n)
#'
#' ## Estimate dispersion (= scale^2 = sigma^2)
#' dispersion_est <- weibull_dispersion_function(
#'   mu = mu,
#'   y = y,
#'   order_indices = order_indices,
#'   family = weibull_family(),
#'   observation_weights = weights,
#'   VhalfInv = NULL,
#'   status = status
#' )
#'
#' print(dispersion_est)          # sigma^2
#' print(sqrt(dispersion_est))    # sigma (comparable to survreg$scale)
#'
#' @export
weibull_dispersion_function <- function(mu,
                                        y,
                                        order_indices,
                                        family,
                                        observation_weights,
                                        VhalfInv,
                                        status){

  ## Maximizes log-likelihood of right-censored data
  log_mu <- log(mu)
  log_y <- log(y)
  observation_weights <- c(observation_weights)
  status <- status[order_indices]

  ## Initialize scale (sigma) using null model
  init_scale <-
    weibull_scale(log_y,
                  mean(log_y),
                  status,
                  observation_weights)
  ## Find scale (sigma) given current mu
  scale <- optim(
    init_scale,
    fn = function(par){
      ## par is scale (sigma)
      -loglik_weibull(log_y,
                      log_mu,
                      status,
                      par,
                      observation_weights)
    },
    lower = init_scale/10, # [Change 2026-02-17] Increased from 5 to 10
    upper = init_scale*10,
    method = 'Brent'
  )$par

  ## Return dispersion = sigma^2 = scale^2
  return(scale^2)
}

#' Weibull GLM Weight Function for Constructing Information Matrix
#'
#' @description
#' Computes diagonal weight matrix \eqn{\textbf{W}} for the information matrix
#' \eqn{\textbf{G} = (\textbf{X}^{T}\textbf{W}\textbf{X} + \textbf{L})^{-1}}
#' in Weibull accelerated failure time (AFT) models.
#'
#' @param mu Predicted survival times
#' @param y Observed response/survival times
#' @param order_indices Order of observations when partitioned to match
#'   "status" to "response"
#' @param family Weibull AFT family
#' @param dispersion Estimated dispersion parameter (sigma^2 = scale^2). The
#'   lgspline framework stores and passes dispersion (sigma^2); the Weibull
#'   scale (sigma) matching \code{survreg$scale} is \code{sqrt(dispersion)}.
#' @param observation_weights Weights of observations submitted to function
#' @param status Censoring indicator (1 = event, 0 = censored)
#'   Indicates whether an event of interest occurred (1) or the observation was
#'   right-censored (0). In survival analysis, right-censoring occurs when the
#'   full survival time is unknown, typically because the study ended or the
#'   subject was lost to follow-up before the event of interest occurred.
#'
#' @return
#' Vector of weights for constructing the diagonal weight matrix
#' \eqn{\textbf{W}} in the information matrix
#' \eqn{\textbf{G} = (\textbf{X}^{T}\textbf{W}\textbf{X} + \textbf{L})^{-1}}.
#'
#' @details
#' This function generates weights used in constructing the information matrix
#' after unconstrained estimates have been found. Specifically, it is used in
#' the construction of the \eqn{\textbf{U}} and \eqn{\textbf{G}} matrices
#' following initial unconstrained parameter estimation.
#'
#' These weights are analogous to the variance terms in generalized linear
#' models (GLMs). Like logistic regression uses \eqn{\mu(1-\mu)}, Poisson
#' regression uses \eqn{e^{\mu}}, and Linear regression uses constant weights,
#' Weibull AFT models use \eqn{\exp((\log y - \log \mu)/\sigma)} where
#' \eqn{\sigma = \sqrt{\text{dispersion}}} is the scale parameter.
#'
#' @examples
#'
#' ## Demonstration of glm weight function in constrained model estimation
#' set.seed(1234)
#' n <- 1000
#' t1 <- rnorm(n)
#' t2 <- rbinom(n, 1, 0.5)
#'
#' ## Generate survival times
#' lambda <- exp(0.5 * t1 + 0.3 * t2)
#' yraw <- rexp(n, rate = 1/lambda)
#'
#' ## Introduce right-censoring
#' status <- rbinom(n, 1, 0.75)
#' y <- ifelse(status, yraw, runif(1, 0, yraw))
#'
#' ## Fit model demonstrating use of custom glm weight function
#' model_fit <- lgspline(y ~ spl(t1) + t2,
#'                       data.frame(y = y, t1 = t1, t2 = t2),
#'                       unconstrained_fit_fxn = unconstrained_fit_weibull,
#'                       family = weibull_family(),
#'                       need_dispersion_for_estimation = TRUE,
#'                       dispersion_function = weibull_dispersion_function,
#'                       glm_weight_function = weibull_glm_weight_function,
#'                       schur_correction_function = weibull_schur_correction,
#'                       status = status,
#'                       opt = FALSE,
#'                       K = 1)
#'
#' print(summary(model_fit))
#'
#'
#' @export
weibull_glm_weight_function <- function(mu,
                                        y,
                                        order_indices,
                                        family,
                                        dispersion,
                                        observation_weights,
                                        status){
  ## dispersion = sigma^2; scale = sigma = sqrt(dispersion)
  scale <- sqrt(dispersion)
  val <- exp((log(y) - log(mu)) / scale)
  if(any(!is.finite(val))){
    return(rep(1, length(val)))
  }
  newval <- val * c(observation_weights)
  return(newval)
}

#' Unconstrained Weibull Accelerated Failure Time Model Estimation
#'
#' @description
#' Estimates parameters for an unconstrained Weibull accelerated failure time
#' (AFT) model supporting right-censored survival data.
#'
#' This both provides a tool for actually fitting Weibull AFT Models, and
#' boilerplate code for users who wish to incorporate Lagrangian multiplier
#' smoothing splines into their own custom models.
#'
#' @param X Design matrix of predictors
#' @param y Survival/response times
#' @param LambdaHalf Square root of penalty matrix
#'   (\eqn{\boldsymbol{\Lambda}^{1/2}})
#' @param Lambda Penalty matrix (\eqn{\boldsymbol{\Lambda}})
#' @param keep_weighted_Lambda Flag to retain weighted penalties
#' @param family Distribution family specification
#' @param tol Convergence tolerance (default 1e-8)
#' @param K Number of partitions minus one (\eqn{K})
#' @param parallel Flag for parallel processing
#' @param cl Cluster object for parallel computation
#' @param chunk_size Processing chunk size
#' @param num_chunks Number of computational chunks
#' @param rem_chunks Remaining chunks
#' @param order_indices Observation ordering indices
#' @param weights Optional observation weights
#' @param status Censoring status indicator (1 = event, 0 = censored)
#'   Indicates whether an event of interest occurred (1) or the observation was
#'   right-censored (0). In survival analysis, right-censoring occurs when the
#'   full survival time is unknown, typically because the study ended or the
#'   subject was lost to follow-up before the event of interest occurred.
#'
#' @return
#' Optimized beta parameter estimates (\eqn{\boldsymbol{\beta}}) for Weibull
#' AFT model
#'
#' @details
#' Estimation Approach:
#' The function employs a two-stage optimization strategy for fitting
#' accelerated failure time models via maximum likelihood:
#'
#' 1. Outer Loop: Estimate Scale Parameter (sigma) using Brent's method
#'
#' 2. Inner Loop: Estimate Regression Coefficients (given sigma) using
#'    damped Newton-Raphson.
#'
#' Note: the score and information inside the Newton-Raphson are both scaled
#' by sigma^2 (i.e., both omit the 1/sigma and 1/sigma^2 prefactors
#' respectively). Since the Newton-Raphson step is
#' G*u = (X^T W X)^{-1} X^T v, the sigma^2 factors cancel and the step
#' remains correct.
#'
#' @examples
#'
#' ## Simulate survival data with covariates
#' set.seed(1234)
#' n <- 1000
#' t1 <- rnorm(n)
#' t2 <- rbinom(n, 1, 0.5)
#'
#' ## Generate survival times with Weibull-like structure
#' lambda <- exp(0.5 * t1 + 0.3 * t2)
#' yraw <- rexp(n, rate = 1/lambda)
#'
#' ## Introduce right-censoring
#' status <- rbinom(n, 1, 0.75)
#' y <- ifelse(status, yraw, runif(1, 0, yraw))
#' df <- data.frame(y = y, t1 = t1, t2 = t2)
#'
#' ## Fit model using lgspline with Weibull AFT unconstrained estimation
#' model_fit <- lgspline(y ~ spl(t1) + t2,
#'                       df,
#'                       unconstrained_fit_fxn = unconstrained_fit_weibull,
#'                       family = weibull_family(),
#'                       need_dispersion_for_estimation = TRUE,
#'                       dispersion_function = weibull_dispersion_function,
#'                       glm_weight_function = weibull_glm_weight_function,
#'                       schur_correction_function = weibull_schur_correction,
#'                       status = status,
#'                       opt = FALSE,
#'                       K = 1)
#'
#' ## Print model summary
#' summary(model_fit)
#'
#' @keywords internal
#' @export
unconstrained_fit_weibull <- function(X,
                                      y,
                                      LambdaHalf,
                                      Lambda,
                                      keep_weighted_Lambda,
                                      family,
                                      tol = 1e-8,
                                      K,
                                      parallel,
                                      cl,
                                      chunk_size,
                                      num_chunks,
                                      rem_chunks,
                                      order_indices,
                                      weights,
                                      status # status goes in the ellipsis arg
) {

  ## Weight if non-null
  if(any(!is.null(weights))){
    if(length(weights) > 0){
      weights <- c(weights)
    }
  } else {
    weights <- rep(1, length(y))
  }

  log_y <- log(y)

  ## Initialize scale (sigma)
  init_scale <- weibull_scale(log_y,
                              mean(log_y),
                              status[order_indices],
                              weights)

  ## First, use outer-loop to optimize scale (sigma)
  #  Then given scale, optimize beta.
  #  Note: the score and information inside damped_newton_r are both
  #  scaled by sigma^2 (omitting 1/sigma and 1/sigma^2 prefactors),
  #  so the Newton-Raphson step remains correct.
  scale <- optim(init_scale,
                 fn = function(par){
                   ## par is scale (sigma)
                   scale <- par
                   beta <- cbind(damped_newton_r(
                     c(mean(log_y), rep(0, ncol(X)-1)),
                     function(par){
                       beta <- cbind(par)
                       log_mu <- c(X %**% beta)
                       loglik_weibull(log_y,
                                      log_mu,
                                      status[order_indices],
                                      scale,
                                      weights) +
                         -0.5*c(crossprod(beta, Lambda) %**% beta)
                     },
                     function(par){
                       beta <- cbind(par)
                       eta <- c(X %**% beta)
                       z <- (log_y - eta) / scale
                       zeta <- exp(z)
                       ## Scaled score (x sigma^2): true score is
                       #  (1/sigma) * X^T diag(w)(exp(z) - status),
                       #  multiplied through by sigma^2 gives
                       #  sigma * X^T diag(w)(exp(z) - status)
                       grad_beta_sc <- crossprod(
                            X,
                            (weights*(zeta - status[order_indices])) * scale
                          ) +
                          Lambda %**% beta
                       cbind(grad_beta_sc)
                     },
                     function(par){
                       beta <- cbind(par)
                       eta <- X %**% beta
                       z <- (log_y - eta) / scale
                       zeta <- c(exp(z))
                       ## Scaled information (x sigma^2): true info is
                       #  (1/sigma^2) * X^T diag(w*exp(z)) X,
                       #  multiplied through by sigma^2 gives
                       #  X^T diag(w*exp(z)) X
                       info_sc <- crossprod(X, weights * zeta * X) + Lambda
                       info_sc
                     },
                     tol
                   ))
                   ## Return negated penalized log-likelihood
                   #  as the objective to minimize
                   log_mu <- X %**% beta
                   -loglik_weibull(log_y,
                                   log_mu,
                                   status[order_indices],
                                   par,
                                   weights) +
                     0.5*c(crossprod(beta, Lambda) %**% beta)
                 },
                 lower = init_scale/10, # [Change 2026-02-17] Increased to 10
                 upper = init_scale*10, # (was 5)
                 method = 'Brent')$par

  ## Now optimize beta, given optimal scale (sigma)
  beta <- cbind(damped_newton_r(
    c(mean(log_y), rep(0, ncol(X)-1)),
    function(par){
      beta <- cbind(par)
      log_mu <- c(X %**% beta)
      loglik_weibull(log_y,
                     log_mu,
                     status[order_indices],
                     scale,
                     weights) +
        -0.5*c(crossprod(beta, Lambda) %**% beta)
    },
    function(par){
      beta <- cbind(par)
      eta <- c(X %**% beta)
      z <- (log_y - eta) / scale
      zeta <- exp(z)
      ## Scaled score (x sigma^2); see note above
      grad_beta_sc <- crossprod(
          X,
          (weights*(zeta - status[order_indices])) * scale
         ) +
        Lambda %**% beta
      cbind(grad_beta_sc)
    },
    function(par){
      beta <- cbind(par)
      eta <- X %**% beta
      z <- (log_y - eta) / scale
      zeta <- c(exp(z))
      ## Scaled information (x sigma^2); see note above
      info_sc <- crossprod(X, weights * zeta * X) + Lambda
      info_sc
    },
    tol
  ))

  return(beta)
}

#' Unconstrained Generalized Linear Model Estimation
#'
#' @description
#' Fits generalized linear models without smoothing constraints
#' using penalized maximum likelihood estimation. This is applied to each
#' partition to obtain the unconstrained estimates, prior to imposing the
#' smoothing constraints.
#'
#' @param X Design matrix of predictors
#' @param y Response variable vector
#' @param LambdaHalf Square root of penalty matrix
#'   (\eqn{\boldsymbol{\Lambda}^{1/2}})
#' @param Lambda Penalty matrix (\eqn{\boldsymbol{\Lambda}})
#' @param keep_weighted_Lambda Logical flag to control penalty matrix handling:
#'   - \code{TRUE}: Return coefficients directly from weighted penalty fitting
#'   - \code{FALSE}: Apply damped Newton-Raphson optimization to refine
#'     estimates
#' @param family Distribution family specification
#' @param tol Convergence tolerance
#' @param K Number of partitions minus one (\eqn{K})
#' @param parallel Flag for parallel processing
#' @param cl Cluster object for parallel computation
#' @param chunk_size Processing chunk size
#' @param num_chunks Number of computational chunks
#' @param rem_chunks Remaining chunks
#' @param order_indices Observation ordering indices
#' @param weights Optional observation weights
#' @param ... Additional arguments passed to \code{glm.fit}
#'
#' @return
#' Optimized parameter estimates for canonical generalized linear models.
#'
#' For fitting non-canonical GLMs, use \code{keep_weighted_Lambda = TRUE}
#' since the score and hessian equations below are no longer valid.
#'
#' For Gamma(link='log') using \code{keep_weighted_Lambda = TRUE} is
#' misleading. The information is weighted by a constant (shape parameter)
#' rather than some mean-variance relationship. So
#' \code{keep_weighted_Lambda = TRUE} is highly recommended for log-link Gamma
#' models. This constant flushes into the penalty terms, and so the
#' formulation of the information matrix is valid.
#'
#' For other scenarios, like probit regression, there will be diagonal weights
#' incorporated into the penalty matrix for providing initial MLE estimates,
#' which technically imposes a prior distribution on beta coefficients that
#' isn't by intent.
#'
#' Heuristically, it shouldn't affect much, as these will be updated to their
#' proper form when providing estimates under constraint; lgspline otherwise
#' does use the correct form of score and information afterwards,
#' regardless of canonical/non-canonical status,
#' as long as 'glm_weight_function' and 'qp_score_function' are properly
#' specified.
#'
#'
#' @keywords internal
#' @export
unconstrained_fit_default <- function(X,
                                      y,
                                      LambdaHalf,
                                      Lambda,
                                      keep_weighted_Lambda,
                                      family,
                                      tol,
                                      K,
                                      parallel,
                                      cl,
                                      chunk_size,
                                      num_chunks,
                                      rem_chunks,
                                      order_indices,
                                      weights,
                                      ...){


  if(nrow(X) == 0){
    return(cbind(rep(0, ncol(X))))
  }

  ## Weight if non-null
  # Yields first NR updated as
  # = (X^{T}V^{-1}X + Lambda)^{-1}X^{T}V^{-1}y for V^{-1} = diag(weights)
  if(any(!is.null(weights))){
    if(length(weights) == length(y)){
      weights <- c(weights)
    } else {
      weights <- rep(1, length(y))
    }
  } else {
    weights <- rep(1, length(y))
  }

  ## Ordinary fit using Tikhinov parameterization
  mod <- try({glm.fit(x = rbind(X, LambdaHalf),
                      y = cbind(c(y, rep(family$linkinv(0),
                                         nrow(LambdaHalf)))),
                      family = family,
                      weights = c(weights, rep(mean(weights),
                                               nrow(LambdaHalf))),
                      ...)}, silent = TRUE)
  if(keep_weighted_Lambda & any(!(inherits(mod, 'try-error')))){
    return(cbind(mod$coefficients))
  }

  if(any(inherits(mod, 'try-error'))){
    init <- c(family$linkfun(mean(y)), rep(0, ncol(X)-1))
  } else {
    init <- c(coef(mod))
  }
  if(any(is.na(init))){
    init <- rep(0, length(init))
  }
  if(any(!is.finite(init))){
    init <- rep(0, length(init))
  }

  ## Remove weights from Tikhinov penalties using damped nr
  res <- cbind(damped_newton_r(
    ## initial guess
    init,
    ## proportional to log-likelihood
    function(par){
      -sum(weights*family$dev.resids(
        y,
        family$linkinv(c(X %**% cbind(par))),
        wt = 1))*0.5 -
        0.5*c(crossprod(par, Lambda) %**% cbind(par))
    },
    ## score
    function(par){
      c(crossprod(X, weights*cbind(y - family$linkinv(X %**% cbind(par)))) -
           Lambda %**% cbind(par))
    },
    ## information
    function(par){
      crossprod(X, weights*c(family$variance(X %**% cbind(par))) * X) +
        Lambda
    },
    tol))
  return(res)
}

#' Collapse Matrix List into a Single Block-Diagonal Matrix
#'
#' @description
#' Transforms a list of matrices into a single block-diagonal matrix. This is
#' useful for quadratic programming problems especially, where the
#' block-diagonal operations may not be plausible.
#'
#' @param matlist List of input matrices
#'
#' @return
#' Block-diagonal matrix combining input matrices
#'
#' @keywords internal
#' @export
collapse_block_diagonal <- function(matlist){
  nrows <- sapply(matlist, nrow)
  ncols <- sapply(matlist, ncol)
  Reduce('rbind', lapply(1:length(matlist), function(k){
    mat <- matrix(0,
                  nrow = nrows[k],
                  ncol = sum(ncols))
    mat[,sum(ncols[-c(k:length(nrows))]) +
          1:ncols[k]] <-
      matlist[[k]]
    mat
  }))
}

#' Compute Newton-Raphson Parameter Update with Numerical Stabilization
#'
#' @description
#' Performs parameter update in iterative optimization.
#'
#' Called by \code{\link{damped_newton_r}} in the update step
#'
#' @param gradient_val Numeric vector of gradient values (\eqn{\textbf{u}})
#' @param neghessian_val Negative Hessian matrix (\eqn{\textbf{G}^{-1}}
#'   approximately)
#'
#' @return
#' Numeric vector of parameter updates (\eqn{\textbf{G}\textbf{u}})
#'
#' @details
#' This helper function is a core component of Newton-Raphson optimization.
#' It provides a computationally-stable approach to computing
#' \eqn{\textbf{G}\textbf{u}}, for information matrix \eqn{\textbf{G}} and
#' score vector \eqn{\textbf{u}}, where the Newton-Raphson update can be
#' expressed as
#' \eqn{\boldsymbol{\beta}^{(m+1)} = \boldsymbol{\beta}^{(m)} + \textbf{G}\textbf{u}}.
#'
#' @seealso \code{\link{damped_newton_r}} for the full optimization routine
#'
#' @keywords internal
#' @export
nr_iterate <- function(gradient_val, neghessian_val){
  sc <- sqrt(mean(abs(neghessian_val))) # for computational stability
  ## [Change 2026-02-14] used crossprod instead of %**% here
  crossprod(t(invert(neghessian_val / sc)), cbind(gradient_val / sc))
}


#' Damped Newton-Raphson Parameter Optimization
#'
#' @description
#' Performs iterative parameter estimation with adaptive step-size dampening
#'
#' Internal function for fitting unconstrained GLM models using damped
#' Newton-Raphson optimization technique.
#'
#' @param parameters Initial parameter vector to be optimized
#' @param loglikelihood Function computing log-likelihood for current parameters
#' @param gradient Function computing parameter gradients
#' @param neghessian Function computing negative Hessian matrix
#' @param tol Numeric convergence tolerance (default 1e-7)
#' @param max_cnt Maximum number of optimization iterations (default 64)
#' @param max_dmp_steps Maximum damping step attempts (default 16)
#'
#' @return
#' Optimized parameter estimates after convergence or reaching iteration limit
#'
#' @details
#' Implements a robust damped Newton-Raphson optimization algorithm.
#'
#' @seealso
#' - \code{\link{nr_iterate}} for parameter update computation
#'
#' @keywords internal
#' @export
damped_newton_r <- function(parameters,
                            loglikelihood,
                            gradient,
                            neghessian,
                            tol = 1e-7,
                            max_cnt = 64,
                            max_dmp_steps = 16){

  ## Initialize convergence checker
  converge_eps <- 100
  eps <- 1000

  ## Initialize old parameters and new parameters to update
  new_param <- parameters
  old_param <- parameters

  ## Do not go past max_cnt
  master_count <- 0
  while(eps > tol & master_count < max_cnt){

    ## Reset components
    new_param <- c(new_param)
    old_param <- new_param
    prev_objective <- loglikelihood(old_param)
    new_objective <- prev_objective - 1
    count <- 0
    if(is.na(prev_objective) |
       is.nan(prev_objective) |
       !is.finite(prev_objective)){
      cat('\n \t Error Encountered, Number of N.R. steps so far: ',
          master_count, '\n')
      stop('\n \t NA/NaN/non-finite value detected when running',
           ' unconstrained damped Newton-Raphson.',
           ' \n \t Try re-fitting a simpler model, using',
           ' greater/smaller penalties, ',
           ' experimenting with different knot locations, or reducing the',
           ' number of knots.')
    }

    ## Damp iterations, only updates if performance improves
    while((new_objective <= prev_objective) & count < max_dmp_steps){
      new_param <- old_param +
        (2^(-count)) * nr_iterate(gradient(old_param),
                                  neghessian(old_param))
      new_objective <- loglikelihood(new_param)
      if(is.na(new_objective) |
         is.nan(new_objective) |
         !is.finite(new_objective)){
        new_objective <- -Inf
      }
      count <- count + 1
    }

    ## Check for change in old vs. new parameters
    eps <- max(abs(old_param - new_param))
    # Break if NAs/NaNs/Infs occur
    if(any(is.na(eps) | is.nan(eps) | !is.finite(eps))){
      new_param <- old_param
      eps <- 0
    }
    master_count <- master_count + 1
  }
  return(new_param)
}


#' Generate Interaction Variable Patterns
#'
#' @description
#' Generates all possible interaction patterns for 2 or 3 variables. This is
#' used in part for identifying which interactions and expansions to exclude
#' (provided to "exclude_these_expansions" argument of lgspline) based on
#' formulas provided.
#'
#' @param vars Character vector of variable names
#'
#' @return
#' Character vector of interaction pattern strings
#'
#' @keywords internal
#' @export
get_interaction_patterns <- function(vars) {
  if(length(vars) == 2) {
    ## Linear: both orderings
    lin <- c(paste0(vars[1], "x", vars[2]),
             paste0(vars[2], "x", vars[1]))

    ## Quadratic: both orderings for each squared term
    quad <- c(paste0(vars[1], "x", vars[2], "^2"),
              paste0(vars[2], "x", vars[1], "^2"))

    return(c(lin, quad))
  } else if(length(vars) == 3) {
    ## All possible 3-way orderings
    return(c(paste0(vars[1], "x", vars[2], "x", vars[3]),
             paste0(vars[1], "x", vars[3], "x", vars[2]),
             paste0(vars[2], "x", vars[1], "x", vars[3]),
             paste0(vars[2], "x", vars[3], "x", vars[1]),
             paste0(vars[3], "x", vars[1], "x", vars[2]),
             paste0(vars[3], "x", vars[2], "x", vars[1])))
  }
}


#' BFGS Implementation for REML Parameter Estimation
#'
#' @description
#' BFGS optimizer designed for REML optimization of correlation parameters. Combines
#' function evaluation and gradient computation into single call to avoid redundant
#' model refitting.
#'
#' @param par Numeric vector of initial parameter values.
#' @param fn Function returning list(objective, gradient). Must return both objective
#' value and gradient vector matching length(par).
#' @param control List of control parameters:
#' \describe{
#'   \item{maxit}{Maximum iterations, default 100}
#'   \item{abstol}{Absolute convergence tolerance, default sqrt(.Machine$double.eps)}
#'   \item{reltol}{Relative convergence tolerance, default sqrt(.Machine$double.eps)}
#'   \item{initial_damp}{Initial damping factor, default 1}
#'   \item{min_damp}{Minimum damping before termination, default 2^-10}
#'   \item{trace}{Print iteration progress, default FALSE}
#' }
#'
#' @return List containing:
#' \describe{
#'   \item{par}{Parameter vector minimizing objective}
#'   \item{value}{Minimum objective value}
#'   \item{counts}{Number of iterations}
#'   \item{convergence}{TRUE if converged within maxit}
#'   \item{message}{Description of termination status}
#'   \item{vcov}{Final approximation of inverse-Hessian, useful for inference}
#' }
#'
#' @details
#' Implements BFGS, used internally by \code{lgspline()} for optimizing correlation parameters via REML
#' when argument for computing gradient \code{VhalfInv_grad} is not NULL.
#'
#' This is more efficient than native BFGS, since gradient and loss can be computed simultaneously,
#' avoiding re-computing components in "fn" and "gr" separately.
#'
#' @examples
#' \donttest{
#'
#' ## Minimize Rosenbrock function
#' fn <- function(x) {
#'   # Objective
#'   f <- 100*(x[2] - x[1]^2)^2 + (1-x[1])^2
#'   # Gradient
#'   g <- c(-400*x[1]*(x[2] - x[1]^2) - 2*(1-x[1]),
#'          200*(x[2] - x[1]^2))
#'   list(f, g)
#' }
#' (res <- efficient_bfgs(c(0.5, 2.5), fn))
#'
#' ## Compare to
#' (res0 <- stats::optim(c(0.5, 2.5), function(x)fn(x)[[1]], hessian = TRUE))
#' solve(res0$hessian)
#' }
#'
#' @keywords internal
#' @export
efficient_bfgs <- function(par, fn, control = list()) {
  ctrl <- list(
    maxit = 50,
    abstol = sqrt(.Machine$double.eps),
    reltol = sqrt(.Machine$double.eps),
    initial_damp = 1,
    min_damp = 2^-16,
    trace = FALSE
  )
  ctrl[names(control)] <- control

  n_params <- length(par)
  x <- par
  Inv <- diag(n_params)
  damp <- ctrl$initial_damp
  best_x <- x
  best_f <- Inf
  best_Inv <- Inv

  result <- fn(c(x))
  if(length(result) != 2) stop("fn must return list of (objective, gradient)")
  f <- result[[1]]
  grad <- result[[2]]
  if(is.null(grad) | any(is.na(grad))) grad <- approx_grad(x, fn)
  if(length(grad) != n_params) stop("gradient must match parameter length")

  ## [Change 2026-02-12] Replaced %**% with crossprod/tcrossprod
  for(iter in 1:ctrl$maxit) {
    prev_x <- x
    prev_f <- f
    prev_grad <- grad

    p <- -crossprod(t(Inv), cbind(grad))
    x_new <- x + damp * p

    result <- fn(c(x_new))
    f_new <- result[[1]]
    grad_new <- result[[2]]
    if(is.null(grad_new) | any(is.na(grad_new))) grad_new <-
      approx_grad(x_new, fn)

    if(is.na(f_new) || is.nan(f_new) || !is.finite(f_new)) {
      damp <- damp/2
      if(damp < ctrl$min_damp) break
      next
    }

    if(f_new < f || iter <= 2) {
      if(f_new < best_f) {
        best_f <- f_new
        best_x <- x_new
        best_Inv <- Inv
      }

      s0 <- cbind(x_new - x)
      y0 <- cbind(grad_new - grad)
      denom <- sum(y0 * s0)

      if(abs(denom) > 1e-16) {
        rho0 <- 1/denom
        term1 <- diag(n_params) - rho0 * tcrossprod(s0, y0)
        term2 <- diag(n_params) - rho0 * tcrossprod(y0, s0)
        Inv <- crossprod(t(term1), crossprod(t(Inv), term2)) +
          rho0 * tcrossprod(s0)
      }

      x <- x_new
      f <- f_new
      grad <- grad_new
      damp <- 1

      if(iter > 2) {
        if(abs(f - prev_f) < ctrl$abstol * (abs(prev_f) + ctrl$reltol)) break
        if(max(abs(x - prev_x)) < ctrl$abstol) break
      }

    } else {
      x <- best_x*(1-damp) + x*damp
      f <- best_f*(1-damp) + f*damp
      Inv <- best_Inv*(1-damp) + Inv*damp
      damp <- damp/2
      if(damp < ctrl$min_damp) break
    }

    if(ctrl$trace) {
      cat(sprintf("Iter %d: f = %f, |grad| = %f, damp = %f\n",
                  iter, f, sqrt(sum(grad^2)), damp))
    }
  }

  list(
    par = c(best_x),
    value = best_f,
    counts = iter,
    convergence = (iter < ctrl$maxit),
    message = if(iter == ctrl$maxit)"Maximum iterations reached" else "Converged",
    vcov = best_Inv
  )
}

#' Finite-difference Gradient Computer
#'
#' @description
#' Computes finite-difference approximation of gradient given input of arguments
#' x and function fn
#'
#' @param x Numeric vector of function arguments
#' @param fn Function returning list(objective, gradient)
#' @param eps Numeric scalar, finite difference tolerance
#'
#' @return Numeric vector of finite-difference approximated gradient
#'
#' @details
#' Used within \code{efficient_bfgs} if needed externally, but internally, this function
#' is actually ignored since when \code{VhalfInv_grad} is not supplied, \code{stats::optim()}
#' is used instead.
#'
#' @keywords internal
#' @export
approx_grad <- function(x, fn, eps = sqrt(.Machine$double.eps)) {
  grad <- numeric(length(x))
  for(i in 1:length(x)) {
    ## Scale base epsilon by parameter magnitude
    h1 <- eps * max(1, abs(x[i]))
    h2 <- h1/2
    x_eps <- c(x)

    ## First gradient at h1
    x_eps[i] <- x[i] + 0.5*h1
    f1_plus <- fn(x_eps)[[1]]
    x_eps[i] <- x[i] - 0.5*h1
    f1_minus <- fn(x_eps)[[1]]
    g1 <- (f1_plus - f1_minus)/h1

    ## Second gradient at h2
    x_eps[i] <- x[i] + 0.5*h2
    f2_plus <- fn(x_eps)[[1]]
    x_eps[i] <- x[i] - 0.5*h2
    f2_minus <- fn(x_eps)[[1]]
    g2 <- (f2_plus - f2_minus)/h2

    ## Richardson extrapolation
    grad[i] <- (4*g2 - g1)/3
  }
  -grad
}

#' Calculate Matrix Square Root for Symmetric Matrices
#'
#' @param mat A symmetric, positive-definite matrix \eqn{\textbf{M}}
#'
#' @return A matrix \eqn{\textbf{B}} such that \eqn{\textbf{B}\textbf{B} = \textbf{M}}
#'
#' @details
#' For matrix \eqn{\textbf{M}}, computes \eqn{\textbf{B}} where \eqn{\textbf{B}\textbf{B} = \textbf{M}} using eigenvalue decomposition:
#'
#' 1. Compute eigendecomposition \eqn{\textbf{M} = \textbf{V}\textbf{D}\textbf{V}^T}
#'
#' 2. Set eigenvalues below \code{sqrt(.Machine$double.eps)} to 0 for stability
#'
#' 3. Take elementwise square root of eigenvalues: \eqn{\textbf{D}^{1/2}}
#'
#' 4. Reconstruct as \eqn{\textbf{B} = \textbf{V} \textbf{D}^{1/2} \textbf{V}^T}
#'
#' This provides the unique symmetric positive-definite square root.
#'
#' You can use this to help construct a custom \code{Vhalf_fxn} for fitting
#' correlation structures, see \code{\link{lgspline}}.
#'
#' @examples
#' ## Identity matrix
#' m1 <- diag(2)
#' matsqrt(m1)  # Returns identity matrix
#'
#' ## Compound symmetry correlation matrix
#' rho <- 0.5
#' m2 <- matrix(rho, 3, 3) + diag(1-rho, 3)
#' B <- matsqrt(m2)
#' # Verify: B %**% B approximately equals m2
#' all.equal(B %**% B, m2)
#'
#' ## Example for correlation structure
#' n_blocks <- 2  # Number of subjects
#' block_size <- 3  # Measurements per subject
#' rho <- 0.7  # Within-subject correlation
#' # Correlation matrix for one subject
#' R <- matrix(rho, block_size, block_size) +
#'      diag(1-rho, block_size)
#' # Full correlation matrix for all subjects
#' V <- kronecker(diag(n_blocks), R)
#' Vhalf <- matsqrt(V)
#'
#' @export
matsqrt <- function(mat) {
  ## [Change 2026-02-18] use symmetric = TRUE
  eig <- eigen(mat, symmetric = TRUE)
  sqrtv <- suppressWarnings(
    suppressMessages(
      sqrt(
        pmax(
          eig$values, .Machine$double.eps
          )
        )
      )
    )
  ## [Change 2026-02-12]: replaced %**% with tcrossprod
  tcrossprod(eig$vectors * rep(sqrtv, each = nrow(eig$vectors)), eig$vectors)
}

#' Calculate Matrix Inverse Square Root for Symmetric Matrices
#'
#' @param mat A symmetric, positive-definite matrix \eqn{\textbf{M}}
#'
#' @return A matrix \eqn{\textbf{B}} such that \eqn{\textbf{B}\textbf{B} = \textbf{M}^{-1}}
#'
#' @details
#' For matrix \eqn{\textbf{M}}, computes \eqn{\textbf{B}} where \eqn{\textbf{B}\textbf{B} = \textbf{M}^{-1}} using eigenvalue decomposition:
#'
#' 1. Compute eigendecomposition \eqn{\textbf{M} = \textbf{V}\textbf{D}\textbf{V}^T}
#'
#' 2. Set eigenvalues below \code{sqrt(.Machine$double.eps)} to 0
#'
#' 3. Take elementwise reciprocal square root: \eqn{\textbf{D}^{-1/2}}
#'
#' 4. Reconstruct as \eqn{\textbf{B} = \textbf{V} \textbf{D}^{-1/2} \textbf{V}^T}
#'
#' For nearly singular matrices, eigenvalues below the numerical threshold
#' are set to 0, and their reciprocals in \eqn{\textbf{D}^{-1/2}} are also set to 0.
#'
#' This implementation is particularly useful for whitening procedures in GLMs
#' with correlation structures and for computing variance-covariance matrices
#' under constraints.
#'
#' You can use this to help construct a custom \code{VhalfInv_fxn} for fitting
#' correlation structures, see \code{\link{lgspline}}.
#'
#' @examples
#' ## Identity matrix
#' m1 <- diag(2)
#' matinvsqrt(m1)  # Returns identity matrix
#'
#' ## Compound symmetry correlation matrix
#' rho <- 0.5
#' m2 <- matrix(rho, 3, 3) + diag(1-rho, 3)
#' B <- matinvsqrt(m2)
#' # Verify: B %**% B approximately equals solve(m2)
#' all.equal(B %**% B, solve(m2))
#'
#' ## Example for GLM correlation structure
#' n_blocks <- 2  # Number of subjects
#' block_size <- 3  # Measurements per subject
#' rho <- 0.7  # Within-subject correlation
#' # Correlation matrix for one subject
#' R <- matrix(rho, block_size, block_size) +
#'      diag(1-rho, block_size)
#' ## Full correlation matrix for all subjects
#' V <- kronecker(diag(n_blocks), R)
#' ## Create whitening matrix
#' VhalfInv <- matinvsqrt(V)
#'
#' # Example construction of VhalfInv_fxn for lgspline
#' VhalfInv_fxn <- function(par) {
#'   rho <- tanh(par)  # Transform parameter to (-1, 1)
#'   R <- matrix(rho, block_size, block_size) +
#'        diag(1-rho, block_size)
#'   kronecker(diag(n_blocks), matinvsqrt(R))
#' }
#'
#' @export
matinvsqrt <- function(mat) {
  ## [Change 2026-02-18] use symmetric = TRUE
  eig <- eigen(mat, symmetric = TRUE)
  sqrtv <- suppressWarnings(
    suppressMessages(
      sqrt(
        pmax(
          eig$values, .Machine$double.eps
          )
        )
      )
    )
  ## [Change 2026-02-12]: replaced %**% with tcrossprod
  tcrossprod(eig$vectors / rep(sqrtv, each = nrow(eig$vectors)), eig$vectors)
}

## ====================================================================== ##
## [Change 2026-02-14] Include internal helpers for correlation str fitting
## ====================================================================== ##

#' Compute Euclidean distance matrix for a cluster block
#' @keywords internal
.compute_dist_block <- function(spacetime, inds) {
  spacetime <- cbind(spacetime)
  diffs <- outer(spacetime[inds, 1], spacetime[inds, 1],
                 function(x, y) (x - y)^2)
  if(ncol(spacetime) > 1){
    diffs <- diffs + Reduce("+", lapply(
      2:ncol(spacetime),
      function(v) outer(spacetime[inds, v], spacetime[inds, v],
                        function(x, y) (x - y)^2)
    ))
  }
  sqrt(diffs / ncol(spacetime))
}

#' Compute rank-based distance matrix for AR(1) structures
#' @keywords internal
.rank_dists <- function(spacetime, inds) {
  diffs <- .compute_dist_block(spacetime, inds)
  block_size <- length(inds)
  unique_diffs <- unique(as.vector(diffs))
  matrix(match(diffs, sort(unique_diffs)) - 1, block_size, block_size)
}

#' Evaluate the REML gradient with respect to a single correlation parameter
#'
#' Computes the derivative of the negative REML objective with respect to a
#' scalar correlation parameter, given the matrix derivative
#' \eqn{\partial \mathbf{V}/\partial \rho}.
#'
#' @details
#' \subsection{Notation}{
#' \eqn{\mathbf{D} = \mathrm{diag}(d_1,\ldots,d_N)} is the diagonal matrix of
#' observation weights (\code{observation_weights}), and
#' \eqn{\mathbf{W} = \mathrm{diag}(w_1,\ldots,w_N)} is the diagonal matrix of
#' GLM working weights evaluated at the current fitted values via
#' \code{glm_weight_function}.  For a GLM with canonical link,
#' \eqn{w_i = \mathrm{Var}(Y_i \mid \mu_i) \cdot d_i} encodes the curvature
#' of the log-likelihood with respect to the linear predictor; for example,
#' for logistic regression \eqn{w_i = \mu_i(1-\mu_i)}.  In Gaussian identity
#' models both reduce to scalar multiples of the identity.
#'
#' \eqn{\mathbf{V}} is the \eqn{N \times N} correlation matrix implied by
#' \code{VhalfInv}, with
#' \eqn{\mathbf{V}^{-1} = (\mathbf{V}^{-1/2})^\top \mathbf{V}^{-1/2}}.
#'
#' The penalized observed information at the current iterate is
#' \deqn{\mathbf{M} = (\mathbf{X}^*)^\top
#'   \mathbf{V}^{-1}\mathbf{W}\mathbf{D}\,\mathbf{X}^*
#'   + \mathbf{U}^\top\boldsymbol{\Lambda}\mathbf{U},}
#' where \eqn{\mathbf{X}^* = \mathbf{X}\mathbf{U}} is the constrained design
#' (\eqn{N \times P}) and the first term is the quadratic approximation to the
#' penalized log-likelihood Hessian at the current \eqn{\boldsymbol{\mu}}.
#' For non-Gaussian families this is a local approximation (the IRWLS/Fisher
#' scoring Hessian), not an exact GLS information matrix.
#'
#' The constraint projection is
#' \eqn{\mathbf{U} = \mathbf{I} -
#'   \mathbf{G}\mathbf{A}(\mathbf{A}^\top\mathbf{G}\mathbf{A})^{-1}
#'   \mathbf{A}^\top}
#' with \eqn{\mathbf{G} = \mathbf{M}^{-1}}.
#' \eqn{\mathbf{U}} is idempotent (\eqn{\mathbf{U}^2 = \mathbf{U}}) but
#' \emph{not} symmetric, so
#' \eqn{\mathbf{U}^\top\boldsymbol{\Lambda}\mathbf{U} \neq
#'   \mathbf{U}\boldsymbol{\Lambda}\mathbf{U}^\top}.
#'
#' The \code{U} stored in \code{model_fit$U} is on the expansion- and
#' response-standardised scale.  The rescaled version used here is
#' \deqn{\tilde{\mathbf{U}} =
#'   \mathbf{U} \cdot \mathrm{diag}(1, s_1, \ldots, s_{p-1},
#'     1, s_1, \ldots) / \hat{\sigma}_y,}
#' where \eqn{s_j} are the expansion scales and \eqn{\hat{\sigma}_y}
#' standardises the response.  All quantities below use
#' \eqn{\tilde{\mathbf{U}}} in place of \eqn{\mathbf{U}}.
#' }
#'
#' \subsection{REML objective and its gradient}{
#' The REML objective is constructed by integrating out the fixed effects from
#' the penalized log-likelihood, using a Laplace approximation to the
#' marginal likelihood for non-Gaussian families.  This approximation is exact
#' for Gaussian identity models and is the standard extension used in
#' restricted maximum likelihood estimation for GLMMs.
#'
#' The REML correction term is
#' \eqn{-\tfrac{1}{2}\log|\mathbf{M}|}, where \eqn{\mathbf{M}} is the
#' penalized observed information defined above.  Differentiating with respect
#' to \eqn{\rho} (noting only \eqn{\mathbf{V}} depends on \eqn{\rho})
#' gives the REML correction gradient
#' \deqn{-\frac{1}{2}\mathrm{tr}\!\Bigl(
#'   \mathbf{M}^{-1}
#'   (\mathbf{X}^*)^\top\mathbf{V}^{-1}
#'   \frac{\partial\mathbf{V}}{\partial\rho}
#'   \mathbf{V}^{-1}\mathbf{W}\mathbf{D}\,\mathbf{X}^*
#'   \Bigr).}
#' }
#'
#' \subsection{Full gradient}{
#' \deqn{\frac{\partial(-\ell_R)}{\partial\rho}
#'   = \frac{1}{N}\Biggl[
#'     \underbrace{
#'       \frac{1}{2}\mathrm{tr}\!\Bigl(
#'         \mathbf{V}^{-1}\frac{\partial\mathbf{V}}{\partial\rho}
#'       \Bigr)
#'     }_{\text{(i) log-det of }\mathbf{V}}
#'     \underbrace{
#'       -\frac{1}{2\tilde{\sigma}^2}
#'       \mathbf{r}^\top\frac{\partial\mathbf{V}}{\partial\rho}\mathbf{r}
#'     }_{\text{(ii) residual quadratic form}}
#'     \underbrace{
#'       -\frac{1}{2}\mathrm{tr}\!\Bigl(
#'         \mathbf{M}^{-1}(\mathbf{X}^*)^\top\mathbf{V}^{-1}
#'         \frac{\partial\mathbf{V}}{\partial\rho}
#'         \mathbf{V}^{-1}\mathbf{W}\mathbf{D}\,\mathbf{X}^*
#'       \Bigr)
#'     }_{\text{(iii) REML correction}}
#'   \Biggr],}
#' where the whitened residual is
#' \deqn{\mathbf{r} =
#'   \mathbf{V}^{-1/2}\mathbf{D}^{1/2}\mathbf{W}^{-1/2}
#'   (\mathbf{y} - \boldsymbol{\mu}),}
#'   and \eqn{r_i =
#'   [\mathbf{V}^{-1/2}(\mathbf{y}-\boldsymbol{\mu})]_i
#'   \sqrt{d_i}/\sqrt{w_i}}.
#'
#' Term (i) does not involve \eqn{\mathbf{D}} or \eqn{\mathbf{W}}; the
#' log-determinant of \eqn{\mathbf{V}} depends only on the correlation
#' structure.  For non-Gaussian families terms (ii) and (iii) are evaluated
#' at the current IRWLS iterate and constitute a local approximation.
#' }
#'
#' @param dV \eqn{N \times N} numeric matrix giving
#'   \eqn{\partial\mathbf{V}/\partial\rho} for one scalar parameter
#'   \eqn{\rho}, with the chain rule through any reparameterization already
#'   applied.
#' @param model_fit A fitted \code{lgspline} object; see
#'   \code{\link{lgspline}}.
#' @param glm_weight_function The GLM weight function used during model
#'   fitting; see the \code{glm_weight_function} argument of
#'   \code{\link{lgspline}}.
#' @param ... Additional arguments forwarded to \code{glm_weight_function}.
#' @return A scalar: the gradient of the negative REML objective with respect
#'   to \eqn{\rho}, divided by \eqn{N}.
#' @seealso \code{\link{lgspline}}, \code{\link{lgspline.fit}}
#' @keywords internal
reml_grad_from_dV <- function(dV, model_fit, glm_weight_function, ...) {

  ## (ii) Residual quadratic form  -0.5 * r' (dV) r / sigma^2
  #
  #  r_i = [V^{-1/2} (y - mu)]_i * sqrt(d_i) / sqrt(w_i)
  #
  #  rep(1, N) is passed for observation_weights inside glm_weight_function
  #  to avoid double-counting D; D enters via sqrt(model_fit$weights) instead.
  #  For logistic: w_i = mu_i(1-mu_i), so r downweights residuals in regions
  #  of high curvature and upweights by observation weight d_i.
  glm_wt_inv <- sqrt(c(glm_weight_function(model_fit$ytilde,
                                           model_fit$y,
                                           1:model_fit$N,
                                           model_fit$family,
                                           model_fit$sigmasq_tilde,
                                           rep(1, model_fit$N),
                                           ...))) /
    sqrt(unlist(model_fit$weights)[model_fit$og_order])

  resid     <- model_fit$y - model_fit$ytilde
  VinvResid <- model_fit$VhalfInv %**% cbind(resid) / glm_wt_inv

  quad_term <- -0.5 * ((t(VinvResid) %**% dV) %**% VinvResid) /
    model_fit$sigmasq_tilde

  ## (i) Log-determinant trace  0.5 * tr(V^{-1} dV/dtheta)
  #
  #  tr(V^{-1} dV) = sum(VhalfInv * (VhalfInv %*% dV))  element-wise
  #  no D or W involvement -- log|V| depends only on the correlation structure
  trace_term <- 0.5 * sum((model_fit$VhalfInv %**% dV) * model_fit$VhalfInv)

  ## (iii) REML correction  -0.5 * tr( M^{-1} (X*)' V^{-1} dV V^{-1} W D X* )
  #
  #  Rescale U to match the unstandardised X in model_fit$X:
  #    U_tilde = U * diag(1, s_1,...,s_{p-1}, 1, s_1,...) / sd_y
  #
  #  VhalfInvX = V^{-1/2} X* = V^{-1/2} X U_tilde   (N x P)
  U <- t(t(model_fit$U) *
           rep(c(1, model_fit$expansion_scales), model_fit$K + 1)) /
    model_fit$sd_y

  VhalfInvX <- model_fit$VhalfInv %**%
    collapse_block_diagonal(model_fit$X)[unlist(model_fit$og_order), ] %**%
    U

  ## Raw block-diagonal penalty on the standardised scale
  if (length(model_fit$penalties$L_partition_list) != (model_fit$K + 1)) {
    model_fit$penalties$L_partition_list <-
      lapply(1:(model_fit$K + 1), function(k) 0)
  }

  Lambda_raw <- collapse_block_diagonal(
    lapply(1:(model_fit$K + 1), function(k)
      c(1, model_fit$expansion_scales) *
        (model_fit$penalties$L_partition_list[[k]] +
           model_fit$penalties$Lambda) %**%
        diag(c(1, model_fit$expansion_scales)) /
        model_fit$sd_y^2
    )
  )

  ## Penalty projected onto constrained subspace: U' Lambda_raw U
  #  U is NOT symmetric so  U' L U  !=  U L U'
  ULambdaU <- t(U) %**% Lambda_raw %**% U

  ## Gram block of M:  (X*)' V^{-1} W D X*
  #    = VhalfInvX' diag(w_i * d_i) VhalfInvX
  #
  #  glm_wt_i = sqrt(w_i * d_i); for logistic sqrt(mu_i(1-mu_i) * d_i).
  #  Scaling rows of VhalfInvX by glm_wt then taking crossprod gives
  #  VhalfInvX' diag(w_i d_i) VhalfInvX = (X*)' V^{-1} W D X* exactly.
  #  rep(1, N) passed again to avoid double-counting D.
  glm_wt <- sqrt(c(glm_weight_function(model_fit$ytilde,
                                       model_fit$y,
                                       1:model_fit$N,
                                       model_fit$family,
                                       model_fit$sigmasq_tilde,
                                       rep(1, model_fit$N),
                                       ...))) *
    sqrt(unlist(model_fit$weights)[model_fit$og_order])

  XVinvX_inv <- invert(crossprod(t(t(VhalfInvX) * c(glm_wt))) + ULambdaU)

  ## Trace computation, numerically stabilised by spectral norm of VInvX
  #
  #  VInvX    = V^{-1} X* = V^{-1/2} VhalfInvX   (N x P)
  #  VInvX_sc = VInvX / sc,  sc = ||VInvX||_2
  #
  #  tr( M^{-1} VInvX' dV VInvX )
  #    = sc^2 * sum( diag( (M^{-1} VInvX_sc') (dV VInvX_sc) ) )
  VInvX <- model_fit$VhalfInv %**% VhalfInvX
  sc    <- sqrt(norm(VInvX, "2"))
  VInvX <- VInvX / sc

  dXVinvX     <- (XVinvX_inv %**% t(VInvX)) %**% (dV %**% VInvX)
  XVinvX_term <- -0.5 * colSums(cbind(c(diag(dXVinvX) * sc))) * sc

  as.numeric(quad_term + trace_term + XVinvX_term) / model_fit$N
}
