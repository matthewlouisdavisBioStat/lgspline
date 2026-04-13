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
#' @param drop_first Logical; if \code{TRUE} and more than one dummy column is
#'   created, drop the first column.
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
#'
#' @export
create_onehot <- function(x, drop_first = FALSE) {
  ## For convenience, coerce to factor to ensure stable level ordering
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
#'
#' x <- c(1, 2, 3, 4, 5)
#' std(x)
#' print(mean(x))
#' print(sd(x))
#'
#'
#' @export
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
#'
#' x <- runif(5)
#' softplus(x)
#'
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
#'
#' M1 <- matrix(1:4, 2, 2)
#' M2 <- matrix(5:8, 2, 2)
#' M1 %**% M2
#'
#'
#' @keywords internal
#' @export
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
#' @return Inverted matrix, or the identity matrix if all inversion attempts fail
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
#' 4. Returns identity matrix if all methods fail, with optional warning
#'
#' For eigendecomposition, uses a small ridge penalty (\code{1e-16}) for stability and
#' zeroes eigenvalues below machine precision.
#'
#' @examples
#'
#' ## Well-conditioned matrix
#' A <- matrix(c(4,2,2,4), 2, 2)
#' invert(A) %*% A
#'
#' ## Singular matrix falls back to M.P. generalized inverse
#' B <- matrix(c(1,1,1,1), 2, 2)
#' invert(B) %*% B
#'
#'
#' @keywords internal
#' @export
invert <- function(mat, include_warnings = FALSE){

  ## Try Cholesky first for SPD matrices
  #  chol2inv(chol(M)) is faster than solve(M) for SPD M
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

  ## generalized inverse with small ridge penalty
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
#'
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
#' as this function, in the same order, excluding the argument
#' "custom_basis_fxn" itself.
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


    ## Pre-allocate once, then fill left-to-right in the same order used
    #  to build the column-name bookkeeping.
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

  ## Walk the expansion columns and recover the derivative from the column
  #  naming convention rather than rebuilding the basis from scratch.
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
#' @param N Integer; number of rows in input matrix
#' @param mat Numeric matrix; data to be partitioned
#' @param K Integer; number of interior knots (resulting in \eqn{K+1} partitions)
#'
#' @return List of length \eqn{K+1}, each element containing the submatrix for that partition
#'
#' @keywords internal
#' @export
knot_expand_list <- function(partition_codes,
                             partition_bounds,
                             N,
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
#' Matrices should be square and have compatible dimensions.
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
#' @param p_expansions Number of columns in the basis expansion per partition
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
#' @return List containing first-derivative matrices and, unless
#'   \code{just_first_derivatives = TRUE}, second-derivative matrices
#'
#' @keywords internal
#' @export
make_derivative_matrix  <-  function(
    p_expansions,
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

  ## Build derivative matrices one variable at a time, keyed by the raw
  #  expansion names used throughout the fitting code.
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
#' @param p_expansions Integer; number of columns in basis expansion per partition
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
make_constraint_matrix <- function(p_expansions,
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
  checkered_fitted_expand = checkered[, rep(1:(K+1), each = p_expansions), drop = FALSE]

  ## Expand out the constraints, located at knots
  constrain_fitted <- CKnots[,rep(1:p_expansions, K+1), drop = FALSE] *
    checkered_fitted_expand

  ## When non-spline and spline present, repeat fitted constraint for
  # spline-only
  if(length(nonspline_cols) > 0 & length(power1_cols) > 0){
    CKnots0 <- CKnots
    CKnots0[,nonspline_cols] <- 0
    constrain_fitted0 <- CKnots0[,rep(1:p_expansions, K+1), drop = FALSE] *
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
    constrain_first_deriv <- first_deriv[, rep(1:p_expansions, K+1),
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
    new_constr <- first_derivs[, rep(1:p_expansions, K+1), drop = FALSE] *
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
    constrain_second_deriv <- second_deriv[, rep(1:p_expansions, K+1),
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
      second_deriv_interaction[, rep(1:p_expansions, K+1), drop = FALSE] *
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
    return(cbind(rep(0, (K+1)*p_expansions)))
  } else {
    ## Otherwise, return A
    return(A)
  }
}

#' Test if Vector is Binary
#'
#' @param x Vector to test
#' @return Logical indicating if x has at most 2 unique values
#'
#' @keywords internal
#' @export
is_binary <- function(x){
  if(length(unique(x)) > 2){
    return(FALSE)
  }
  TRUE
}


#' Compute Derivative of Inverse-Information Matrix G with Respect to Lambda
#'
#' @description
#' Calculates the derivative of the inverse-information matrix
#' \eqn{\textbf{G}} with respect to the smoothing parameter \eqn{\lambda},
#' supporting both shared and partition-specific penalties.
#'
#' @param G A list of inverse-information matrices \eqn{\textbf{G}} for each
#'   partition
#' @param dPenalty_dlambda Derivative of the shared penalty matrix with
#'   respect to \eqn{\lambda}.
#' @param K Number of partitions minus 1 (\eqn{K})
#' @param unique_penalty_per_partition Logical indicating partition-specific penalties
#' @param dPenalty_partition_list Optional list of derivatives of the
#'   partition-specific penalty matrices with respect to \eqn{\lambda}.
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
                               dPenalty_dlambda,
                               K,
                               unique_penalty_per_partition,
                               dPenalty_partition_list,
                               parallel,
                               cl,
                               chunk_size,
                               num_chunks,
                               rem_chunks) {
  unique_penalty_per_partition <- isTRUE(unique_penalty_per_partition)

  if(!unique_penalty_per_partition){
    penalty_derivative <- dPenalty_dlambda
  } else {
    penalty_derivative <- lapply(dPenalty_partition_list, function(dP_k){
      dPenalty_dlambda + dP_k
    })
  }

  ## [2026-02-12] replaced %**% with crossprod/tcrossprod
  if(parallel & !is.null(cl)) {
    ## Handle remainder chunks first
    if(rem_chunks > 0) {
      rem_indices <- num_chunks*chunk_size + 1:rem_chunks
      rem <- lapply(rem_indices, function(k) {
        if(unique_penalty_per_partition){
          -crossprod(t(G[[k]]), crossprod(t(penalty_derivative[[k]]), G[[k]]))
        } else {
          -crossprod(t(G[[k]]), crossprod(t(penalty_derivative), G[[k]]))
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
            -crossprod(t(G[[k]]), crossprod(t(penalty_derivative[[k]]), G[[k]]))
          } else {
            -crossprod(t(G[[k]]), crossprod(t(penalty_derivative), G[[k]]))
          }
        })
      })),
      rem
    )

  } else {
    ## Sequential computation
    result <- lapply(1:(K+1), function(k){
      if(unique_penalty_per_partition){
        -crossprod(t(G[[k]]), crossprod(t(penalty_derivative[[k]]), G[[k]]))
      } else {
        -crossprod(t(G[[k]]), crossprod(t(penalty_derivative), G[[k]]))
      }
    })
  }

  return(result)
}

## The tuning hat-diagonal helpers below are documented separately.
#' Compute the Diagonal of the Tuning Hat Matrix
#'
#' @description
#' Computes the per-observation diagonal of the constrained hat matrix used
#' during tuning, without forming the full \eqn{N \times N} matrix.
#'
#' @param X List of partition-specific design matrices.
#' @param G List of partition-specific \eqn{\mathbf{G}_k} matrices.
#' @param A Constraint matrix.
#' @param AGAInv Inverse of \eqn{\mathbf{A}^{\top}\mathbf{G}\mathbf{A}}.
#' @param K Number of interior knots.
#' @param p_expansions Number of columns per partition.
#' @param eps Numeric floor used to keep \eqn{1 - h_{ii}} away from zero.
#'
#' @return List with \code{h_diag} and \code{C_list}, both partition-wise.
#'
#' @keywords internal
.compute_hat_diag <- function(X,
                              G,
                              A,
                              AGAInv,
                              K,
                              p_expansions,
                              eps = 1e-7) {

  n_partitions <- K + 1
  h_diag <- vector("list", n_partitions)
  C_list <- vector("list", n_partitions)

  for (k in seq_len(n_partitions)) {
    rows <- (k - 1) * p_expansions + seq_len(p_expansions)
    A_k <- A[rows, , drop = FALSE]
    V_k <- G[[k]] %**% A_k
    C_k <- G[[k]] - V_k %**% AGAInv %**% t(V_k)
    h_k <- rowSums((X[[k]] %**% C_k) * X[[k]])
    h_diag[[k]] <- pmin(pmax(h_k, 0), 1 - eps)
    C_list[[k]] <- 0.5 * (C_k + t(C_k))
  }

  list(h_diag = h_diag, C_list = C_list)
}


#' Compute the Derivative of the Tuning Hat Matrix Diagonal
#'
#' @description
#' Computes \eqn{d h_{ii} / d\lambda_1} partition-wise for the exact
#' leave-one-out tuning criterion.
#'
#' @param X List of partition-specific design matrices.
#' @param G List of partition-specific \eqn{\mathbf{G}_k} matrices.
#' @param dG_dlambda List of derivatives \eqn{d\mathbf{G}_k / d\lambda_1}.
#' @param A Constraint matrix.
#' @param AGAInv Inverse of \eqn{\mathbf{A}^{\top}\mathbf{G}\mathbf{A}}.
#' @param K Number of interior knots.
#' @param p_expansions Number of columns per partition.
#'
#' @return List with \code{dh_diag} and \code{dC_list}, both partition-wise.
#'
#' @keywords internal
.compute_dhat_diag <- function(X,
                               G,
                               dG_dlambda,
                               A,
                               AGAInv,
                               K,
                               p_expansions) {

  n_partitions <- K + 1
  V_list <- vector("list", n_partitions)
  dV_list <- vector("list", n_partitions)
  D_theta <- matrix(0, ncol(A), ncol(A))

  for (k in seq_len(n_partitions)) {
    rows <- (k - 1) * p_expansions + seq_len(p_expansions)
    A_k <- A[rows, , drop = FALSE]
    V_k <- G[[k]] %**% A_k
    dV_k <- dG_dlambda[[k]] %**% A_k
    V_list[[k]] <- V_k
    dV_list[[k]] <- dV_k
    D_theta <- D_theta + crossprod(A_k, dV_k)
  }

  W <- AGAInv %**% D_theta %**% AGAInv
  dh_diag <- vector("list", n_partitions)
  dC_list <- vector("list", n_partitions)

  for (k in seq_len(n_partitions)) {
    X_k <- X[[k]]
    V_k <- V_list[[k]]
    dV_k <- dV_list[[k]]
    dC_k <- dG_dlambda[[k]] -
      dV_k %**% AGAInv %**% t(V_k) -
      V_k %**% AGAInv %**% t(dV_k) +
      V_k %**% W %**% t(V_k)
    dC_k <- 0.5 * (dC_k + t(dC_k))

    ## Compute the leverage derivative observation-wise rather than taking
    # diag(X dC X') from the full matrix product. This is algebraically
    # equivalent, but is more numerically stable for the local terms that are
    # heavily reweighted inside the exact LOO gradient.
    q_k <- X_k %**% V_k
    dq_k <- X_k %**% dV_k
    dh_diag[[k]] <- rowSums((X_k %**% dG_dlambda[[k]]) * X_k) -
      2 * rowSums((dq_k %**% AGAInv) * q_k) +
      rowSums((q_k %**% W) * q_k)

    dC_list[[k]] <- dC_k
  }

  list(dh_diag = dh_diag, dC_list = dC_list)
}


#' Evaluate Exact Leave-One-Out Criterion at a Given Penalty Configuration
#'
#' @description
#' Computes the exact leave-one-out tuning criterion for the same transformed
#' linearized problem used elsewhere in the tuning machinery.
#'
#' @inheritParams .compute_gcvu
#'
#' @return List containing the fitted objects plus \code{criterion_value},
#'   \code{LOO_u}, \code{h_diag}, and \code{C_list}.
#'
#' @keywords internal
.compute_loocv <- function(par,
                           log_penalty_vec,
                           env,
                           ...) {
  verbose <- env$verbose
  eps <- 1e-7

  if (verbose) cat("        loocv_fxn start\n")

  wiggle_penalty <- exp(par[1])
  flat_ridge_penalty <- exp(par[2])
  if (env$unique_penalty_per_predictor || env$unique_penalty_per_partition) {
    penalty_vec <- exp(c(par[-c(1:2)]))
  } else {
    penalty_vec <- c()
  }

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
  L1 <- Lambda_list[[2]]
  L2 <- Lambda_list[[3]]

  if (verbose) cat("        compute_G_eigen\n")
  schur_corrections <- lapply(seq_len(env$K + 1), function(k) 0)
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

  if (verbose) cat("        loocv_fxn fit coefficients\n")
  return_G_getB <- TRUE
  B_list <- .fit_coefficients(G_list, Lambda,
                              Lambda_list$L_partition_list,
                              env, return_G_getB, ...)
  G_list <- B_list$G_list
  B <- B_list$B

  if (verbose) cat("        loocv_fxn AGAmult_wrapper\n")
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

  if (verbose) cat("        loocv_fxn get predictions\n")
  preds <- .compute_tuning_predictions(env$X, B, env$K,
                                       env$parallel & env$parallel_matmult,
                                       env$cl,
                                       env$chunk_size,
                                       env$num_chunks,
                                       env$rem_chunks)

  if (verbose) cat("        loocv_fxn residuals\n")
  residuals <- .compute_tuning_residuals(env$y, preds, env$delta,
                                         env$family,
                                         env$observation_weights,
                                         env$K, env$order_list,
                                         ...)

  if (verbose) cat("        loocv_fxn hat diagonal\n")
  hat_result <- .compute_hat_diag(env$X,
                                  G_list$G,
                                  env$A,
                                  AGAInv,
                                  env$K,
                                  env$p_expansions,
                                  eps = eps)

  denom_list <- lapply(hat_result$h_diag, function(h_k) pmax(1 - h_k, eps))
  loo_resid_sq <- mapply(function(r_k, denom_k) {
    (r_k / denom_k)^2
  }, residuals, denom_list, SIMPLIFY = FALSE)
  LOO_u <- mean(unlist(loo_resid_sq))

  if (verbose) cat("        loocv_fxn penalization operations\n")
  mp <- .compute_meta_penalty(wiggle_penalty, penalty_vec,
                              env$meta_penalty,
                              env$unique_penalty_per_predictor,
                              env$unique_penalty_per_partition)

  criterion_value <- LOO_u + mp

  if (verbose) cat("        done LOO,", LOO_u, "\n")

  return(list(criterion_value = criterion_value,
              LOO_u          = LOO_u,
              B              = B,
              G_list         = G_list,
              Lambda         = Lambda,
              L1             = L1,
              L2             = L2,
              L_predictor_list = Lambda_list$L_predictor_list,
              L_partition_list = Lambda_list$L_partition_list,
              residuals      = residuals,
              AGAInv         = AGAInv,
              h_diag         = hat_result$h_diag,
              C_list         = hat_result$C_list))
}


#' Compute Closed-Form Gradient of Exact Leave-One-Out Criterion
#'
#' @description
#' Computes the analytical gradient of the exact leave-one-out tuning
#' criterion with respect to the log penalty parameters.
#'
#' @param par Numeric vector; log-scale penalty parameters.
#' @param log_penalty_vec Numeric vector; log-scale predictor/partition
#'   penalties.
#' @param outlist List or NULL; cached results from \code{.compute_loocv()}.
#' @param env List; tuning environment.
#' @param ... Additional arguments passed to fitting functions.
#'
#' @return List with \code{criterion_value}, \code{gradient}, and
#'   \code{outlist}.
#'
#' @details
#' The implementation below keeps the exact LOO gradient fully analytic,
#' including the leverage derivative. In empirical diagnostics, the
#' observation-wise leverage component can be numerically delicate even when
#' the aggregate gradient remains informative; see the developer diagnostics in
#' \code{inst/diagnostics/gradient_diagnostics.R} for inspection tools.
#'
#' @keywords internal
.compute_loocv_gradient <- function(par,
                                    log_penalty_vec,
                                    outlist = NULL,
                                    env,
                                    ...) {
  verbose <- env$verbose
  eps <- 1e-7

  if (verbose) cat("        loo gradient start\n")

  wiggle_penalty <- exp(par[1])
  flat_ridge_penalty <- exp(par[2])
  if (env$unique_penalty_per_predictor || env$unique_penalty_per_partition) {
    penalty_vec <- exp(c(par[-c(1:2)]))
  } else {
    penalty_vec <- c()
  }

  lambda_1 <- wiggle_penalty
  lambda_2 <- flat_ridge_penalty

  if (any(is.null(outlist))) {
    outlist <- .compute_loocv(par, log_penalty_vec, env, ...)
  }

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
  AGAInvAGXy <- GhalfXy_temp_list[[2]]

  if (verbose) cat("        compute_dG_dlambda \n")
  dPenalty_dlambda1 <- outlist$Lambda / lambda_1
  dPenalty_partition_dlambda1 <- if (env$unique_penalty_per_partition) {
    lapply(outlist$L_partition_list, function(L_k) L_k / lambda_1)
  } else {
    list()
  }
  dG_dlambda <- compute_dG_dlambda(outlist$G_list$G,
                                   dPenalty_dlambda1,
                                   env$K,
                                   env$unique_penalty_per_partition,
                                   dPenalty_partition_dlambda1,
                                   env$parallel & env$parallel_matmult,
                                   env$cl,
                                   env$chunk_size,
                                   env$num_chunks,
                                   env$rem_chunks)

  if (verbose) cat("        compute_dGhalf \n")
  dGhalf <- compute_dGhalf(dG_dlambda,
                           env$p_expansions,
                           env$K,
                           env$parallel & env$parallel_eigen,
                           env$cl,
                           env$chunk_size,
                           env$num_chunks,
                           env$rem_chunks)

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

  dPenalty_dlambda2 <- (lambda_1 / lambda_2) * outlist$L2
  dPenalty_partition_dlambda2 <- if (env$unique_penalty_per_partition) {
    lapply(outlist$L_partition_list, function(L_k) 0 * L_k)
  } else {
    list()
  }
  dG_dlambda2 <- compute_dG_dlambda(outlist$G_list$G,
                                    dPenalty_dlambda2,
                                    env$K,
                                    env$unique_penalty_per_partition,
                                    dPenalty_partition_dlambda2,
                                    env$parallel & env$parallel_matmult,
                                    env$cl,
                                    env$chunk_size,
                                    env$num_chunks,
                                    env$rem_chunks)
  dGhalf2 <- compute_dGhalf(dG_dlambda2,
                            env$p_expansions,
                            env$K,
                            env$parallel & env$parallel_eigen,
                            env$cl,
                            env$chunk_size,
                            env$num_chunks,
                            env$rem_chunks)
  dG_u_dlambda2_Xyr <- compute_dG_u_dlambda_xy(
    AGAInvAGXy,
    outlist$AGAInv,
    outlist$G_list$G,
    env$A,
    dG_dlambda2,
    env$p_expansions,
    env$R_constraints,
    env$K,
    env$Xy,
    outlist$G_list$Ghalf,
    dGhalf2,
    GhalfXy_temp,
    env$parallel & env$parallel_matmult,
    env$cl,
    env$chunk_size,
    env$num_chunks,
    env$rem_chunks)

  if (verbose) cat("        compute_dhat_diag \n")
  dhat_result <- .compute_dhat_diag(env$X,
                                    outlist$G_list$G,
                                    dG_dlambda,
                                    env$A,
                                    outlist$AGAInv,
                                    env$K,
                                    env$p_expansions)
  dhat_result2 <- .compute_dhat_diag(env$X,
                                     outlist$G_list$G,
                                     dG_dlambda2,
                                     env$A,
                                     outlist$AGAInv,
                                     env$K,
                                     env$p_expansions)

  ## .compute_loocv() clips h_diag into [0, 1 - eps]. The derivative used in
  ## the LOO gradient must respect that same clipping; once an observation is
  ## on a boundary, the derivative of the clipped value is zero.
  dh_diag_clipped1 <- mapply(
    function(h_k, dh_k) dh_k * ((h_k > 0) & (h_k < 1 - eps)),
    outlist$h_diag,
    dhat_result$dh_diag,
    SIMPLIFY = FALSE
  )
  dh_diag_clipped2 <- mapply(
    function(h_k, dh_k) dh_k * ((h_k > 0) & (h_k < 1 - eps)),
    outlist$h_diag,
    dhat_result2$dh_diag,
    SIMPLIFY = FALSE
  )

  dr_list <- vector("list", env$K + 1)
  for (k in seq_len(env$K + 1)) {
    rows <- (k - 1) * env$p_expansions + seq_len(env$p_expansions)
    q_k <- dG_u_dlambda1_Xyr[rows, , drop = FALSE]
    dr_list[[k]] <- -c(env$X[[k]] %**% q_k)
  }

  dLOO_dlambda1 <- (2 / env$N_obs) * sum(unlist(mapply(
    function(r_k, h_k, dr_k, dh_k) {
      denom_k <- pmax(1 - h_k, eps)
      w_k <- r_k / denom_k
      (w_k / denom_k) * (dr_k + w_k * dh_k)
    },
    outlist$residuals,
    outlist$h_diag,
    dr_list,
    dh_diag_clipped1,
    SIMPLIFY = FALSE
  )))

  dr_list2 <- vector("list", env$K + 1)
  for (k in seq_len(env$K + 1)) {
    rows <- (k - 1) * env$p_expansions + seq_len(env$p_expansions)
    q_k <- dG_u_dlambda2_Xyr[rows, , drop = FALSE]
    dr_list2[[k]] <- -c(env$X[[k]] %**% q_k)
  }

  dLOO_dlambda2 <- (2 / env$N_obs) * sum(unlist(mapply(
    function(r_k, h_k, dr_k, dh_k) {
      denom_k <- pmax(1 - h_k, eps)
      w_k <- r_k / denom_k
      (w_k / denom_k) * (dr_k + w_k * dh_k)
    },
    outlist$residuals,
    outlist$h_diag,
    dr_list2,
    dh_diag_clipped2,
    SIMPLIFY = FALSE
  )))

  ## The exact LOO gradient is kept fully analytic here. Empirically, the
  # observation-wise leverage derivative can be numerically delicate even when
  # the aggregated tuning gradient remains useful; see diagnostics for caveats.
  log_wiggle_gradient <- dLOO_dlambda1 * lambda_1
  log_flat_gradient <- dLOO_dlambda2 * lambda_2

  gradient <- cbind(c(log_wiggle_gradient,
                      log_flat_gradient))

  if (env$unique_penalty_per_predictor) {
    predictor_penalties <- penalty_vec[grep("predictor", names(penalty_vec))]
    predictor_penalty_gradient <- sapply(
      seq_along(predictor_penalties), function(j) {
        mean(diag(outlist$L_predictor_list[[j]])) /
          mean(diag(outlist$Lambda)) *
          log_wiggle_gradient
      })
    gradient <- cbind(c(c(gradient), predictor_penalty_gradient))
  }

  if (env$unique_penalty_per_partition) {
    partition_penalties <- penalty_vec[grep("partition", names(penalty_vec))]
    partition_penalty_gradient <- sapply(
      seq_along(partition_penalties), function(j) {
        mean(diag(outlist$L_partition_list[[j]])) /
          mean(diag(outlist$Lambda +
                      outlist$L_partition_list[[j]])) *
          log_wiggle_gradient
      })
    gradient <- cbind(c(gradient, partition_penalty_gradient))
  }

  regularizer <- .compute_meta_penalty_gradient(
    wiggle_penalty, penalty_vec, env$meta_penalty,
    env$unique_penalty_per_predictor,
    env$unique_penalty_per_partition)

  if (verbose) cat("        loo gradient end \n")

  return(list(criterion_value = outlist$criterion_value,
              LOO_u = outlist$LOO_u,
              gradient = c(gradient) + regularizer,
              outlist = outlist))
}


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

    ## resids1: residuals from regressing G^{1/2}X'y onto
    #           G^{1/2}A -- the constrained least-squares component
    resids1 <- do.call('.lm.fit', list(x = GhalfA / comp_stab_sc,
                                       y = GhalfXy_temp / comp_stab_sc)
    )$residuals * comp_stab_sc

    ## resids2: residuals from regressing (dG^{1/2}/dl)X'y onto
    #           (dG^{1/2}/dl)A -- derivative-side correction term
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


#' Compute Eigenvalues and Related Matrices for G
#'
#' @description
#' Computes partition-wise inverse-information matrices and their matrix square
#' roots from the penalized information matrix via eigendecomposition.
#'
#' @param X_gram List of Gram matrices \eqn{X_k^T W_k X_k} by partition.
#' @param Lambda Penalty matrix \eqn{\Lambda} (shared across partitions).
#' @param K Integer; number of interior knots (partitions = \eqn{K+1}).
#' @param parallel Logical; use parallel processing across partitions.
#' @param cl Cluster object from \code{parallel::makeCluster}.
#' @param chunk_size,num_chunks,rem_chunks Partition distribution parameters.
#' @param family GLM family object.
#' @param unique_penalty_per_partition Logical; if \code{TRUE}, add
#'   partition-specific penalties from \code{L_partition_list}.
#' @param L_partition_list List of partition-specific penalty matrices.
#' @param keep_G Logical; if \code{TRUE}, return the full \eqn{G} matrix.
#' @param schur_corrections List of Schur complement correction matrices.
#'
#' @return A list with components \code{G} (or \code{NULL} blocks when
#'   \code{keep_G = FALSE}), \code{Ghalf}, and optionally \code{GhalfInv}.
#'
#' @keywords internal
#' @export
compute_G_eigen <- function(X_gram, Lambda, K, parallel, cl,
                            chunk_size, num_chunks, rem_chunks,
                            family, unique_penalty_per_partition,
                            L_partition_list, keep_G = TRUE,
                            schur_corrections) {

  ## Helper to compute eigen-based matrix powers efficiently
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

  ## Gaussian identity only needs G and G^(1/2); the other families also
  #  need G^(-1/2) for transformed-score and variance calculations.
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

  ## Repackage the per-partition eigendecomposition output into the list
  #  structure expected downstream by get_B(), blockfit, and inference.
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
#' Computes the least-squares projection component
#' \eqn{\mathbf{G}^{1/2}\mathbf{A}(\mathbf{A}^{T}\mathbf{G}\mathbf{A})^{-1}
#' \mathbf{A}^{T}\mathbf{G}\mathbf{X}^{T}\mathbf{y}} together with the
#' intermediate product
#' \eqn{(\mathbf{A}^{T}\mathbf{G}\mathbf{A})^{-1}\mathbf{A}^{T}\mathbf{G}\mathbf{X}^{T}\mathbf{y}}
#' for reuse downstream.
#'
#' @return Unnamed two-element list containing the projected result vector and
#'   the intermediate
#'   \eqn{(\mathbf{A}^{T}\mathbf{G}\mathbf{A})^{-1}\mathbf{A}^{T}\mathbf{G}\mathbf{X}^{T}\mathbf{y}}
#'   product.
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


#' Efficient Matrix Multiplication of G and A Matrices
#'
#' @param G List of G matrices
#' @param A Constraint matrix
#' @param K Number of partitions minus 1
#' @param p_expansions Number of columns per partition
#' @param R_constraints Number of constraint columns
#' @param parallel Use parallel processing
#' @param cl Cluster object
#' @param chunk_size Size of parallel chunks
#' @param num_chunks Number of chunks
#' @param rem_chunks Remaining chunks
#'
#' @details
#' Computes G %**% A when G has block diagonal structure and A is a matrix.
#' Processes in parallel chunks if enabled.
#'
#' @return List of matrix products, one per partition
#'
#' @keywords internal
#' @export
GAmult_wrapper <- function(G,
                           A,
                           K,
                           p_expansions,
                           R_constraints,
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
      rem <- GAmult(G_rem, A[(num_chunks*chunk_size)*p_expansions + 1:((rem_chunks)*p_expansions),],
                    rem_chunks-1, p_expansions, R_constraints)
    } else {
      rem <- list()
    }

    ## Process main chunks in parallel
    result <- c(
      Reduce("c",parallel::parLapply(cl, 1:num_chunks, function(chunk) {
        start_idx <- (chunk - 1) * chunk_size
        G_chunk <- G[(start_idx + 1):min(start_idx + chunk_size, length(G))]
        A_chunk <- A[(start_idx*p_expansions + 1):min((start_idx + chunk_size)*p_expansions,
                                                      nrow(A)), ]
        GAmult(G_chunk,
               A_chunk,
               chunk_size-1,
               p_expansions,
               R_constraints)
      })),
      rem
    )

  } else {
    ## Sequential computation using original C++ function
    result <- GAmult(G, A, K, p_expansions, R_constraints)
  }

  return(result)
}



#' Efficient Matrix Multiplication for \eqn{\textbf{A}^{T}\textbf{G}\textbf{A}}
#'
#' @param G List of G matrices (\eqn{\textbf{G}})
#' @param A Constraint matrix (\eqn{\textbf{A}})
#' @param K Number of partitions minus 1 (\eqn{K})
#' @param p_expansions Number of columns per partition
#' @param R_constraints Number of constraint columns
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
                            p_expansions,
                            R_constraints,
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
                                  p_expansions)
    } else {
      rem_result <- matrix(0, R_constraints, R_constraints)
    }

    # Process main chunks in parallel
    chunk_results <- parallel::parLapply(cl, 1:num_chunks, function(chunk) {
      chunk_start <- (chunk - 1) * chunk_size
      G_chunk <- G[(chunk_start + 1):(chunk_start + chunk_size)]
      AGAmult_chunk(G_chunk, A, chunk_start, chunk_start + chunk_size - 1, p_expansions)
    })

    # Sum all results
    return(Reduce('+', c(chunk_results, list(rem_result))))

  } else {
    return(AGAmult(G, A, K, p_expansions, R_constraints))
  }
}


#' Efficiently Construct U Matrix
#'
#' @param G List of G matrices (\eqn{\textbf{G}})
#' @param A Constraint matrix (\eqn{\textbf{A}})
#' @param K Number of partitions minus 1 (\eqn{K})
#' @param p_expansions Number of columns per partition
#' @param R_constraints Number of constraint columns
#'
#' @return \eqn{\textbf{U}} matrix for constraints
#'
#' @details
#' Computes \eqn{\textbf{U} = \textbf{I} - \textbf{G}\textbf{A}(\textbf{A}^{T}\textbf{G}\textbf{A})^{-1}\textbf{A}^{T}} efficiently, avoiding unnecessary
#' multiplication of blocks of \eqn{\textbf{G}} with all-0 elements.
#'
#' @keywords internal
#' @export
get_U <- function(G, A, K, p_expansions, R_constraints){
  AGAInv <- invert(AGAmult(G, A, K, p_expansions, R_constraints))
  I_minus_U <- t(matmult_U(crossprod(t(A),
                          crossprod(t(AGAInv), (-t(A)))), G, p_expansions, K))
  return(I_minus_U + diag(p_expansions*(K+1)))
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
#' predictor/partition-specific components.
#'
#' @param custom_penalty_mat Matrix; optional custom ridge penalty structure
#' @param L1 Matrix; integrated squared second derivative penalty (\eqn{\textbf{L}_1})
#' @param wiggle_penalty,flat_ridge_penalty Numeric; smoothing and ridge penalty parameters
#' @param K Integer; number of interior knots (\eqn{K})
#' @param p_expansions Integer; number of basis columns per partition
#' @param unique_penalty_per_predictor,unique_penalty_per_partition Logical; enable predictor/partition-specific penalties
#' @param penalty_vec Named numeric; custom penalty values for predictors/partitions
#' @param colnm_expansions Character; column names for linking penalties to predictors
#' @param just_Lambda Logical; return only combined penalty matrix (\eqn{\boldsymbol{\Lambda}})
#'
#' @return List containing Lambda, L1, L2, L_predictor_list, L_partition_list;
#'   or just Lambda if \code{just_Lambda=TRUE} and no partition penalties.
#'
#' @keywords internal
#' @export
compute_Lambda <- function(custom_penalty_mat,
                           L1,
                           wiggle_penalty,
                           flat_ridge_penalty,
                           K,
                           p_expansions,
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
    L2 <- diag(p_expansions)*((diag(L1) == 0)*flat_ridge_penalty)
  }

  ## Combined shared penalty = wiggle * (L1 + L2)
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


#' Compute Eigenvalues and Related Matrices for G
#'
#' @description
#' Computes partition-wise inverse-information matrices and their matrix square roots
#' from the penalized information matrix via eigendecomposition.
#'
#' @param X_gram List of Gram matrices \eqn{X_k^T W_k X_k} by partition.
#' @param Lambda Penalty matrix \eqn{\Lambda} (shared across partitions).
#' @param K Integer; number of interior knots (partitions = \eqn{K+1}).
#' @param parallel Logical; use parallel processing across partitions.
#' @param cl Cluster object from \code{parallel::makeCluster}.
#' @param chunk_size,num_chunks,rem_chunks Partition distribution parameters.
#' @param family GLM family object.
#' @param unique_penalty_per_partition Logical; if \code{TRUE}, add
#'   partition-specific penalties from \code{L_partition_list}.
#' @param L_partition_list List of partition-specific penalty matrices.
#' @param keep_G Logical; if \code{TRUE}, return the full \eqn{G} matrix.
#' @param schur_corrections List of Schur complement correction matrices.
#'
#' @return A list with components \code{G} (or \code{NULL} blocks when
#'   \code{keep_G = FALSE}), \code{Ghalf}, and optionally \code{GhalfInv}.
#'
#' @keywords internal
#' @export
compute_G_eigen <- function(X_gram, Lambda, K, parallel, cl,
                            chunk_size, num_chunks, rem_chunks,
                            family, unique_penalty_per_partition,
                            L_partition_list, keep_G = TRUE,
                            schur_corrections) {

  ## Determine once whether GhalfInv is needed
  need_GhalfInv <- (paste0(family)[2] != 'identity') |
    (paste0(family)[1] != 'gaussian')

  ## Core computation for a single partition k
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

    ## Clamp non-positive eigenvalues for pseudoinverse treatment
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

  if(parallel & !is.null(cl)) {
    if(rem_chunks > 0) {
      rem_indices <- num_chunks * chunk_size + 1:rem_chunks
      rem <- lapply(rem_indices, compute_one)
    } else {
      rem <- list()
    }

    result <- c(
      Reduce("c", parallel::parLapply(cl, 1:num_chunks, function(chunk) {
        inds <- (chunk - 1)*chunk_size + 1:chunk_size
        lapply(inds, compute_one)
      })),
      rem
    )
  } else {
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


#' Compute Eigen-Based Square-Root Factors for \eqn{d\textbf{G}/d\lambda}
#'
#' @description
#' Applies an eigendecomposition-based matrix square root to each partition-wise
#' \eqn{d\textbf{G}/d\lambda} matrix after replacing \code{NA} entries with 0
#' and treating all-0 inputs with an identity fallback.
#'
#' @param dG_dlambda List of \eqn{p \times p} \eqn{d\textbf{G}/d\lambda} matrices by partition
#' @param p_expansions Integer; number of columns per partition
#' @param K Integer; number of interior knots (\eqn{K})
#' @param parallel,cl,chunk_size,num_chunks,rem_chunks Parallel computation parameters
#'
#' @return List of \eqn{p \times p} eigen-based square-root factors derived from
#'   the partition-wise \eqn{d\textbf{G}/d\lambda} matrices
#'
#' @keywords internal
#' @export
compute_dGhalf <- function(dG_dlambda,
                           p_expansions,
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


#' Compute Integrated Squared Second Derivative Penalty Matrix
#'
#' @description
#' Computes the \eqn{p \times p} integrated squared second-derivative penalty matrix
#' \eqn{\boldsymbol{\Lambda}_s} for one partition of a monomial spline,
#' such that
#' \deqn{\boldsymbol{\beta}_k^\top \boldsymbol{\Lambda}_s
#'        \boldsymbol{\beta}_k
#'   = \int_{\mathbf{a}}^{\mathbf{b}}
#'     \|\tilde{f}_k''(\mathbf{t})\|^2 \, d\mathbf{t},}
#' where \eqn{\mathbf{a}} and \eqn{\mathbf{b}} are the observed
#' predictor minimums and maximums, and
#' \eqn{\tilde{f}_k(\mathbf{t}) = \mathbf{x}_k^\top
#' \boldsymbol{\beta}_k} is the fitted function for partition
#' \eqn{\mathcal{P}_k}.
#'
#' The implementation uses a general monomial derivative rule that
#' handles all term types (marginal powers, two-way interactions,
#' quadratic interactions, three-way interactions) in a single unified
#' loop.
#'
#' @param colnm_expansions Character vector of length \eqn{p} giving
#'   column names for the basis expansion matrix \eqn{\mathbf{X}_k}.
#'   Each name encodes the monomial structure: predictor names appear
#'   literally for linear terms, with \code{^d} suffixes for degree
#'   \eqn{d}, and an \code{x} separator between factors in interaction
#'   terms (e.g.\ \code{"_1_alphax_2_beta^2"} for
#'   \eqn{t_1 t_2^2}).  Predictor names must not be substrings of one
#'   another.
#' @param C Numeric \eqn{N \times p} matrix of basis expansions
#'   \eqn{\mathbf{X}_k}.  Used only to read the observed range
#'   \eqn{[a_j, b_j]} of each predictor via \code{min}/\code{max} of
#'   the linear columns.
#' @param power1_cols Integer vector.  Column indices of linear terms
#'   \eqn{t_j}.
#' @param power2_cols Integer vector.  Column indices of quadratic terms
#'   \eqn{t_j^2}.
#' @param power3_cols Integer vector.  Column indices of cubic terms
#'   \eqn{t_j^3}.
#' @param power4_cols Integer vector.  Column indices of quartic terms
#'   \eqn{t_j^4}.
#' @param interaction_single_cols Integer vector.  Column indices of
#'   linear-by-linear interaction terms \eqn{t_j t_l}.
#' @param interaction_quad_cols Integer vector.  Column indices of
#'   linear-by-quadratic interaction terms \eqn{t_j t_l^2} and
#'   \eqn{t_j^2 t_l}.
#' @param triplet_cols Integer vector.  Column indices of three-way
#'   interaction terms \eqn{t_j t_l t_m}.
#' @param select_cols Optional integer vector of predictor indices
#'   (positions within \code{power1_cols}) whose curvature operators
#'   \eqn{D_v} are summed.  Defaults to all \eqn{q} predictors.
#'
#' @return
#' A symmetric positive semi-definite \eqn{p \times p} matrix
#' \eqn{\boldsymbol{\Lambda}_s} with
#' \deqn{[\boldsymbol{\Lambda}_s]_{ij}
#'   = \sum_{v=1}^{q} \int_{\mathbf{a}}^{\mathbf{b}}
#'     D_v(\phi_i) \, D_v(\phi_j) \, d\mathbf{t}.}
#'
#' @details
#' \subsection{Mathematical framework}{
#'
#' Let the basis expansion for partition \eqn{\mathcal{P}_k} be
#' \eqn{\mathbf{x}_k = (\phi_1(\mathbf{t}), \ldots,
#' \phi_p(\mathbf{t}))^\top} where each \eqn{\phi_i} is a
#' multivariate monomial
#' \deqn{\phi_i(\mathbf{t})
#'   = \prod_{j=1}^{q} t_j^{\alpha_{ij}}.}
#'
#' The second derivative of \eqn{\tilde{f}_k} decomposes into
#' \eqn{q} total curvature operators, one per predictor.  For
#' predictor \eqn{v}:
#' \deqn{D_v = \frac{\partial^2}{\partial t_v^2}
#'   + \sum_{s \neq v} \frac{\partial^2}{\partial t_v \, \partial t_s}.}
#'
#' The monomial derivative rule gives each second partial derivative
#' in closed form.  For the pure second derivative (\eqn{r = s = v}):
#' \deqn{\frac{\partial^2}{\partial t_v^2}
#'       \prod_{j} t_j^{\alpha_j}
#'   = \alpha_v(\alpha_v - 1) \;
#'     t_v^{\alpha_v - 2} \prod_{j \neq v} t_j^{\alpha_j}.}
#'
#' For a mixed second derivative (\eqn{s \neq v}):
#' \deqn{\frac{\partial^2}{\partial t_v \, \partial t_s}
#'       \prod_{j} t_j^{\alpha_j}
#'   = \alpha_v \alpha_s \;
#'     t_v^{\alpha_v - 1} t_s^{\alpha_s - 1}
#'     \prod_{j \neq v,s} t_j^{\alpha_j}.}
#'
#' In both cases a term is zero when the required exponent would be
#' negative (e.g.\ \eqn{\alpha_v < 2} for the pure case).  Applying
#' \eqn{D_v} to \eqn{\phi_i} produces a sum of monomials with known
#' coefficients and exponent vectors.
#' }
#'
#' \subsection{Integration}{
#'
#' Because every \eqn{D_v(\phi_i)} is polynomial, the product
#' \eqn{D_v(\phi_i) \, D_v(\phi_j)} is also polynomial and the
#' multivariate integral factorises over predictors:
#' \deqn{\int_{\mathbf{a}}^{\mathbf{b}}
#'       \prod_{j=1}^{q} t_j^{e_j} \, d\mathbf{t}
#'   = \prod_{j=1}^{q}
#'     \frac{b_j^{e_j+1} - a_j^{e_j+1}}{e_j + 1}.}
#'
#' Notably, this integral runs over \emph{all} \eqn{q} predictor
#' ranges, including predictors that do not appear in the integrand
#' (for which \eqn{e_j = 0} and the factor reduces to
#' \eqn{b_j - a_j}).
#' }
#'
#' \subsection{Single-predictor verification}{
#'
#' For \eqn{q = 1} with expansion
#' \eqn{\mathbf{x} = (1, t, t^2, t^3)^\top} on \eqn{[a, b]}, the
#' penalty matrix reduces to
#' \deqn{\boldsymbol{\Lambda}_s
#'   = \int_a^b \mathbf{x}'' \mathbf{x}''^\top \, dt
#'   = \begin{pmatrix}
#'       0 & 0 & 0 & 0 \\
#'       0 & 0 & 0 & 0 \\
#'       0 & 0 & 4(b - a) & 6(b^2 - a^2) \\
#'       0 & 0 & 6(b^2 - a^2) & 12(b^3 - a^3)
#'     \end{pmatrix},}
#' matching the formula in Section 2.3.
#' }
#'
#' @section Naming convention:
#' Column names in \code{colnm_expansions} must encode each monomial's
#' predictor content.  The parser checks each predictor name (taken
#' from the linear columns) against the column name:
#' \itemize{
#'   \item \code{predname^d} (literal caret) signals exponent \eqn{d}.
#'   \item \code{predname} without \code{^} signals exponent 1.
#'   \item Absence of \code{predname} signals exponent 0.
#' }
#' Multiple factors in an interaction are separated by \code{x}.
#' Predictor names must be unique and must not be substrings of one
#' another.
#'
#' @references
#' Reinsch, C. H. (1967). Smoothing by spline functions.
#' \emph{Numerische Mathematik}, 10(3), 177--183.
#'
#' @seealso \code{\link{lgspline}} for the full model fitting
#'   interface.
#'
#' @examples
#'
#' ## Verification example: 3-predictor model with all term types:
#' #  This example constructs the penalty matrix analytically and then
#' #  verifies selected entries against closed-form hand calculations.
#' #  Users can extend or adapt this example to audit new basis
#' #  expansions.
#'
#' set.seed(1234)
#' n <- 2000
#' t1 <- runif(n, 1, 4)     # predictor 1, support [1, 4]
#' t2 <- runif(n, -2, 3)    # predictor 2, support [-2, 3]
#' t3 <- runif(n, 0.5, 2)   # predictor 3, support [0.5, 2]
#'
#' ## Predictor names (must not be substrings of one another)
#' pn <- c("_1_aa", "_2_bb", "_3_cc")
#'
#' ## Build column names encoding the monomial structure
#' col_names <- c(
#'   pn,                                                # linear
#'   paste0(pn, "^2"),                                  # quadratic
#'   paste0(pn, "^3"),                                  # cubic
#'   paste0(pn[1], "x", pn[2]),                         # t1*t2
#'   paste0(pn[1], "x", pn[3]),                         # t1*t3
#'   paste0(pn[2], "x", pn[3]),                         # t2*t3
#'   paste0(pn[2], "x", pn[1], "^2"),                   # t1^2*t2
#'   paste0(pn[1], "x", pn[2], "^2"),                   # t1*t2^2
#'   paste0(pn[1], "x", pn[2], "x", pn[3])              # t1*t2*t3
#' )
#' p_expansions <- length(col_names)   # 14 basis functions
#'
#' ## Build the expansion matrix C
#' C <- cbind(
#'   t1, t2, t3,
#'   t1^2, t2^2, t3^2,
#'   t1^3, t2^3, t3^3,
#'   t1*t2, t1*t3, t2*t3,
#'   t2*t1^2, t1*t2^2,
#'   t1*t2*t3
#' )
#' colnames(C) <- col_names
#'
#' ## Compute the penalty matrix
#' Ls <- get_2ndDerivPenalty(
#'   colnm_expansions       = col_names,
#'   C                      = C,
#'   power1_cols            = 1:3,
#'   power2_cols            = 4:6,
#'   power3_cols            = 7:9,
#'   power4_cols            = integer(0),
#'   interaction_single_cols = 10:12,
#'   interaction_quad_cols  = 13:14,
#'   triplet_cols           = 15,
#'   p_expansions           = p_expansions,
#'   select_cols            = 1:3
#' )
#'
#' ## Hand-computed reference values (exact, using true bounds) ---
#' #  Notation: dt1 = 4-1 = 3, dt2 = 3-(-2) = 5, dt3 = 2-0.5 = 1.5
#' #
#' #  Entry [4,4]: t1^2 diagonal.
#' #    D_1(t1^2) = 2.  No other D_v contributes.
#' #    integral (2)^2 dt1 dt2 dt3 = 4 * 3 * 5 * 1.5 = 90
#' #
#' #  Entry [10,10]: t1*t2 diagonal.
#' #    D_1(t1*t2) = 1 (mixed d^2/dt1 dt2).
#' #    D_2(t1*t2) = 1 (mixed d^2/dt2 dt1).
#' #    integral 1 dt1 dt2 dt3 + integral 1 dt1 dt2 dt3 = 2*3*5*1.5 = 45
#' #
#' #  Entry [15,15]: t1*t2*t3 diagonal.
#' #    D_1(t1*t2*t3) = t3 + t2  (mixed partials d^2/dt1 dt2 and d^2/dt1 dt3)
#' #    D_2(t1*t2*t3) = t1 + t3  (similarly)
#' #    D_3(t1*t2*t3) = t1 + t2  (similarly)
#' #    Full integral = sum of 3 terms = 120 + 337.5 + 266.25 = 723.75
#'
#' ## Compare (allowing ~0.5% tolerance for data-derived bounds)
#' cat("Entry [4,4]:   analytical =", round(Ls[4,4], 2),
#'     "  exact = 90.00\n")
#' cat("Entry [10,10]: analytical =", round(Ls[10,10], 2),
#'     "  exact = 45.00\n")
#' cat("Entry [15,15]: analytical =", round(Ls[15,15], 2),
#'     "  exact = 723.75\n")
#' cat("Symmetric:", isSymmetric(Ls), "\n")
#'
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
                                p_expansions,
                                select_cols = NULL) {

  ## Guard: return zero matrix when there are no linear columns
  mat <- matrix(0, nrow = p_expansions, ncol = p_expansions)
  if (length(power1_cols) == 0L) return(mat)

  if (!any(!is.null(select_cols))) {
    select_cols <- seq_along(power1_cols)
  }

  ## 1. Cache observed predictor bounds a_j, b_j for j = 1, ..., q
  q_predictors <- length(power1_cols)
  a_lo <- vapply(power1_cols, function(j) min(C[, j]), numeric(1))
  b_hi <- vapply(power1_cols, function(j) max(C[, j]), numeric(1))

  ## 2. Monomial integration helpers
  mono_int <- function(k, a, b) {
    (b^(k + 1L) - a^(k + 1L)) / (k + 1L)
  }

  int_product <- function(e_vec) {
    val <- 1
    for (j in seq_len(q_predictors)) {
      val <- val * mono_int(e_vec[j], a_lo[j], b_hi[j])
    }
    val
  }

  ## 3. Parse exponent vectors alpha_i from column names
  pred_names <- colnm_expansions[power1_cols]

  parse_exponents <- function(col_idx) {
    nm <- colnm_expansions[col_idx]
    e  <- integer(q_predictors)
    for (j in seq_len(q_predictors)) {
      pnm <- pred_names[j]
      pat <- paste0(pnm, "\\^([0-9]+)")
      m   <- regmatches(nm, regexpr(pat, nm, perl = TRUE))
      if (length(m) == 1L && nchar(m) > 0L) {
        e[j] <- as.integer(sub(paste0(pnm, "\\^"), "", m))
      } else if (grepl(pnm, nm, fixed = TRUE)) {
        e[j] <- 1L
      }
    }
    e
  }

  ## Collect all basis columns that can carry a non-zero penalty
  penalised <- sort(unique(c(power1_cols, power2_cols, power3_cols,
                             power4_cols, interaction_single_cols,
                             interaction_quad_cols, triplet_cols)))
  penalised <- penalised[penalised >= 1L & penalised <= p_expansions]

  ## Exponent matrix Alpha: Alpha[i, j] = alpha_{ij}
  Alpha <- matrix(0L, nrow = p_expansions, ncol = q_predictors)
  for (idx in penalised) {
    Alpha[idx, ] <- parse_exponents(idx)
  }

  ## 4. Total curvature operator D_v applied to monomial phi_i
  Dv <- function(alpha, v) {
    out <- list()

    ## Pure second derivative d^2/dt_v^2
    if (alpha[v] >= 2L) {
      e    <- alpha
      e[v] <- e[v] - 2L
      out[[length(out) + 1L]] <- list(c = alpha[v] * (alpha[v] - 1L),
                                      e = e)
    }

    ## Mixed second derivatives d^2/(dt_v dt_s) for each s != v
    for (s in seq_len(q_predictors)) {
      if (s == v) next
      if (alpha[v] >= 1L && alpha[s] >= 1L) {
        e    <- alpha
        e[v] <- e[v] - 1L
        e[s] <- e[s] - 1L
        out[[length(out) + 1L]] <- list(c = alpha[v] * alpha[s],
                                        e = e)
      }
    }
    out
  }

  ## 5. Penalty matrix entry
  entry <- function(i, j) {
    total <- 0
    for (v in select_cols) {
      ti <- Dv(Alpha[i, ], v)
      tj <- Dv(Alpha[j, ], v)
      if (length(ti) == 0L || length(tj) == 0L) next
      for (a in ti) {
        for (b in tj) {
          total <- total + a$c * b$c * int_product(a$e + b$e)
        }
      }
    }
    total
  }

  ## 6. Populate Lambda_s (symmetric; only upper triangle computed)
  Ls <- matrix(0, nrow = p_expansions, ncol = p_expansions)

  for (ii in seq_along(penalised)) {
    i <- penalised[ii]
    for (jj in ii:length(penalised)) {
      j <- penalised[jj]
      v <- entry(i, j)
      if (v != 0) {
        Ls[i, j] <- v
        Ls[j, i] <- v
      }
    }
  }

  Ls
}

#' Wrapper for Integrated Second-Derivative Penalty Computation
#'
#' @description
#' Computes the integrated squared second-derivative penalty matrix with
#' optional parallel processing.
#'
#' @param K Number of partitions (\eqn{K+1})
#' @param colnm_expansions Column names of basis expansions
#' @param C Basis expansion matrix used to determine observed predictor ranges
#' @param power1_cols Linear term columns
#' @param power2_cols Quadratic term columns
#' @param power3_cols Cubic term columns
#' @param power4_cols Quartic term columns
#' @param interaction_single_cols Single interaction columns
#' @param interaction_quad_cols Quadratic interaction columns
#' @param triplet_cols Triplet interaction columns
#' @param nonspline_cols Predictors not treated as spline effects
#' @param p_expansions Number of columns in the basis expansion per partition
#' @param parallel Logical to enable parallel processing
#' @param cl Cluster object for parallel computation
#'
#' @return A \eqn{p \times p} integrated squared second-derivative penalty
#'   matrix.
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
                                        p_expansions,
                                        parallel,
                                        cl) {
  ## Temporarily append corresponding polynomial terms so nonspline columns
  #  are carried through the same analytic penalty bookkeeping.
  colnm_expansions_og <- colnm_expansions
  if(length(nonspline_cols) > 0){
    for(jj in 1:length(nonspline_cols)){
      power1_cols <- c(power1_cols, nonspline_cols[jj])
      if(length(power2_cols) > 0){
        colnm_expansions <- c(colnm_expansions,
                              paste0(colnm_expansions[nonspline_cols[jj]],'^2'))
        C <- cbind(C, 0)
        power2_cols <- c(power2_cols, ncol(C))
      }
      if(length(power3_cols) > 0){
        colnm_expansions <- c(colnm_expansions,
                              paste0(colnm_expansions[nonspline_cols[jj]],'^3'))
        C <- cbind(C, 0)
        power3_cols <- c(power3_cols, ncol(C))
      }
      if(length(power4_cols) > 0){
        colnm_expansions <- c(colnm_expansions,
                              paste0(colnm_expansions[nonspline_cols[jj]],'^4'))
        C <- cbind(C, 0)
        power4_cols <- c(power4_cols, ncol(C))
      }
    }
    ## Update colnames and number of columns with new nonspline power terms
    colnames(C) <- colnm_expansions
    p_expansions <- ncol(C)
  }

  ## If parallel processing is worthwhile, split the selected derivative
  if(parallel & (K > 1)){
    chunk_size <- max(1, floor(2 * length(cl)))
    total_cols <- length(power1_cols)
    result <- matrix(0, nrow = p_expansions, ncol = p_expansions)

    for(start in seq(1, total_cols, by = chunk_size)) {
      end <- min(start + chunk_size - 1, total_cols)
      chunk_result <- Reduce("+",
                             parallel::parLapply(cl,
                                                 start:end,
                                                 function(select_col) {
                                                   get_2ndDerivPenalty(
                                                     colnm_expansions,
                                                     C,
                                                     power1_cols,
                                                     power2_cols,
                                                     power3_cols,
                                                     power4_cols,
                                                     interaction_single_cols,
                                                     interaction_quad_cols,
                                                     triplet_cols,
                                                     p_expansions,
                                                     select_col)
                                                 }))
      result <- result + chunk_result
    }
  } else {
    result <- get_2ndDerivPenalty(colnm_expansions,
                                  C,
                                  power1_cols,
                                  power2_cols,
                                  power3_cols,
                                  power4_cols,
                                  interaction_single_cols,
                                  interaction_quad_cols,
                                  triplet_cols,
                                  p_expansions)
  }
  colnames(result) <- colnm_expansions
  rownames(result) <- colnm_expansions
  ## Isolate the entries excluding appended
  result <- result[colnm_expansions_og, colnm_expansions_og]
  return(result)
}


#' Unconstrained Generalized Linear Model Estimation
#'
#' @description
#' Fits generalized linear models without smoothing constraints
#' using penalized maximum likelihood estimation. This is applied to each
#' partition to obtain the unconstrained estimates, prior to imposing the
#' smoothing constraints.
#'
#' Hot-start estimates are initialized by treating the matrix square-root inverse
#' as pseudo-observations, appending to the rows of the design matrix, and
#' calling \code{\link{glm.fit}} replacing the response for pseudo-observations
#' with the inverse link function applied to the value of eta = XB = 0. For log
#' link, exp(0) = 1 is valid for most families. For cases like inverse link,
#' dividing by 0 is obviously not possible, so this is replaced with 1/tol
#' where tol is the convergence tolerance argument.
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
#' Numeric column vector of unconstrained coefficient estimates.
#'
#' For fitting non-canonical GLMs, use \code{keep_weighted_Lambda = TRUE}
#' since the score and Hessian equations below are no longer valid.
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

  ## Return 0 vector for empty partitions
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

  ## Determine if quasi family should be used for initialization
  family_init <- family

  ## Override with quasi when initialization approach,
  #  calling g(0) to get the vector of response corresponding to natural link
  #  eta = XB = 0,
  #  For example, 1/(1+exp(-0)) = 0.5 is not a valid value according to family=
  #  binomial(). Will only affect hot-start initialization.
  if(!is.null(family$family)){

    fam_name  <- family$family
    link_name <- family$link

    needs_quasi <- FALSE

    ## Binomial always needs quasi for augmented ridge trick
    if(fam_name == "binomial"){
      needs_quasi <- TRUE
    }

    ## Inverse-type links produce Inf at g(0)
    if(link_name %in% c("inverse", "1/mu^2")){
      needs_quasi <- TRUE
    }

    ## Other binomial-style links
    if(link_name %in% c("logit", "probit", "cloglog", "cauchit")){
      needs_quasi <- TRUE
    }

    if(needs_quasi){

      family_init <- switch(
        fam_name,
        "binomial" = quasibinomial(link = link_name),
        "Gamma" = quasi(link = link_name, variance = "mu^2"),
        "inverse.gaussian" = quasi(link = link_name, variance = "mu^3"),
        family
      )
    }
  }

  ## For cases like inverse link, we need finite values, so replace 0 with
  #  the "tol" value
  psuedoresponse <- family_init$linkinv(0)
  if(is.na(psuedoresponse)){
    psuedoresponse <- family_init$linkinv(tol)
  }
  if(!is.finite(psuedoresponse)){
    psuedoresponse <- family_init$linkinv(tol)
  }

  ## Ordinary fit using Tikhonov parameterization
  # This is not "correct" because the GLM weights are introduced into the
  # square-root penalty matrix
  # but makes good hot-start initialization
  mod <- try({
    glm.fit(
      x = rbind(X, LambdaHalf),
      y = cbind(c(y,
                  rep(psuedoresponse,
                      nrow(LambdaHalf)))),
      family = family_init,
      weights = c(weights,
                  rep(mean(weights),
                      nrow(LambdaHalf))),
      ...
    )
  }, silent = TRUE)

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

  ## Adjusts Tikhonov penalties using damped nr, with initial hot-start
  #  values from glm.fit using augmented data
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

#' Collapse Matrix List into a Single Dense Block-Layout Matrix
#'
#' @description
#' Transforms a list of matrices into a single dense block-layout matrix. This
#' is useful especially for quadratic programming problems, where operating on
#' lists of blocks is not convenient.
#'
#' @param matlist List of input matrices
#'
#' @return
#' Dense matrix with each input block placed in its own column range
#'
#' @keywords internal
#' @export
collapse_block_diagonal <- function(matlist){
  ## Lay each block into its own column range so dense routines such as
  #  solve.QP can work with the full system directly.
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
#' Numeric column vector of parameter updates (\eqn{\textbf{G}\textbf{u}})
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
  sc <- sqrt(max(abs(diag(neghessian_val))))
  if(sc < .Machine$double.eps) sc <- 1
  H <- neghessian_val / sc
  g <- cbind(gradient_val / sc)
  out <- try(solve(H, g), silent = TRUE)
  if(inherits(out, 'try-error')) {
    out <- invert(H) %**% g
  }
  out
}


#' Damped Newton-Raphson Parameter Optimization
#'
#' @description
#' Performs iterative parameter estimation with adaptive step-size damping
#'
#' Internal function for fitting unconstrained model components using damped
#' Newton-Raphson updates.
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
#' Final parameter vector returned at termination.
#'
#' @details
#' Implements a robust damped Newton-Raphson optimization algorithm.
#' The Newton direction is computed once per outer iteration and reused
#' across damping half-steps.
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

  ## Initialize
  new_param <- parameters
  old_param <- parameters
  eps <- Inf
  master_count <- 0

  ## Run the loop while change between coefficients > tol and
  #  max iterations have not been reached
  while(eps > tol & master_count < max_cnt){

    old_param <- c(new_param)
    prev_objective <- loglikelihood(old_param)

    if(!is.finite(prev_objective)){
      cat('\n \t Error Encountered, Number of N.R. steps so far: ',
          master_count, '\n')
      stop('\n \t NA/NaN/non-finite value detected when running',
           ' unconstrained damped Newton-Raphson.',
           ' \n \t Try re-fitting a simpler model, using',
           ' greater/smaller penalties, ',
           ' experimenting with different knot locations, or reducing the',
           ' number of knots.')
    }

    ## Compute Newton direction once, reuse across damping half-steps
    direction <- nr_iterate(gradient(old_param), neghessian(old_param))

    ## Damp: halve step size until objective improves
    accepted <- FALSE
    for(count in 0:(max_dmp_steps - 1)){
      candidate <- old_param + (2^(-count)) * c(direction)
      new_objective <- loglikelihood(candidate)
      if(!is.finite(new_objective)) new_objective <- -Inf
      if(new_objective > prev_objective){
        new_param <- candidate
        accepted <- TRUE
        break
      }
    }

    if(!accepted){
      ## No damping level improved; declare convergence at current point
      new_param <- old_param
      break
    }

    ## Otherwise re-compute difference between coefficient updates
    eps <- max(abs(old_param - new_param))
    if(!is.finite(eps)){
      new_param <- old_param
      break
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
#' Character vector of interaction pattern strings for 2- or 3-variable
#' inputs, or \code{NULL} otherwise.
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
#'   \item{maxit}{Maximum iterations, default 50}
#'   \item{abstol}{Absolute convergence tolerance, default sqrt(.Machine$double.eps)}
#'   \item{reltol}{Relative convergence tolerance, default sqrt(.Machine$double.eps)}
#'   \item{initial_damp}{Initial damping factor, default 1}
#'   \item{min_damp}{Minimum damping before termination, default 2^-16}
#'   \item{trace}{Print iteration progress, default FALSE}
#' }
#'
#' @return List containing:
#' \describe{
#'   \item{par}{Parameter vector minimizing objective}
#'   \item{value}{Minimum objective value}
#'   \item{counts}{Number of iterations}
#'   \item{convergence}{TRUE if termination occurred before \code{maxit};
#'     this reflects the current stopping rule rather than a separate
#'     post-hoc convergence check}
#'   \item{message}{Description of termination status}
#'   \item{vcov}{Final approximation of inverse-Hessian, useful for inference}
#' }
#'
#' @details
#' Implements BFGS, used internally by \code{lgspline()} for optimizing
#' correlation parameters via REML when an analytic REML gradient is supplied.
#'
#' This is more efficient than native BFGS, since gradient and loss can be computed simultaneously,
#' avoiding re-computing components in "fn" and "gr" separately.
#'
#' @examples
#'
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
#' Computes a finite-difference approximation derived from the objective in
#' \code{fn} at \code{x}, returned in the sign convention currently used by
#' \code{efficient_bfgs}.
#'
#' @param x Numeric vector of function arguments
#' @param fn Function returning list(objective, gradient)
#' @param eps Numeric scalar, finite difference tolerance
#'
#' @return Numeric vector of finite-difference approximated gradient values in
#'   the current optimizer sign convention
#'
#' @details
#' Used within \code{efficient_bfgs} when \code{fn} does not supply a usable
#' gradient. In the main \code{lgspline()} correlation-optimization path,
#' \code{stats::optim()} is used instead when no analytic REML gradient is
#' available. Internally this helper returns the negated central-difference
#' approximation rather than the raw derivative.
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
#' @param mat A symmetric matrix \eqn{\textbf{M}}
#'
#' @return A matrix \eqn{\textbf{B}} such that \eqn{\textbf{B}\textbf{B}}
#'   equals \eqn{\textbf{M}} on the positive-eigenvalue subspace, with
#'   non-positive components truncated to 0.
#'
#' @details
#'
#' Uses an eigenvalue-decomposition-based approach.
#'
#' Non-positive eigenvalues are set to 0 before taking fourth roots.
#'
#' This implementation is particularly useful for whitening procedures in GLMs
#' with correlation structures and for computing variance-covariance matrices
#' under constraints.
#'
#' You can use this to help construct a custom \code{Vhalf_fxn}, or more
#' directly to build the \eqn{\mathbf{V}^{1/2}} input supplied to
#' \code{\link{lgspline}} for correlation-aware fits.
#'
#' @examples
#'
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
  eig <- eigen(mat, symmetric = TRUE)
  wts <- ifelse(eig$values > 0,
                eig$values^0.25,
                0)
  crossprod(t(eig$vectors) * wts)
}

#' Calculate Matrix Inverse Square Root for Symmetric Matrices
#'
#' @param mat A symmetric matrix \eqn{\textbf{M}}
#'
#' @return A matrix \eqn{\textbf{B}} such that \eqn{\textbf{B}\textbf{B}}
#'   equals the Moore-Penrose-style inverse on the positive-eigenvalue
#'   subspace, with non-positive components truncated to 0.
#'
#' @details
#'
#' Uses an eigenvalue-decomposition-based approach.
#'
#' Non-positive eigenvalues are set to 0 before taking inverse fourth roots.
#'
#' This implementation is particularly useful for whitening procedures in GLMs
#' with correlation structures and for computing variance-covariance matrices
#' under constraints.
#'
#' You can use this to help construct a custom \code{VhalfInv_fxn} for
#' \code{\link{lgspline}}.  When only \code{VhalfInv} is supplied there,
#' the corresponding \code{Vhalf} is reconstructed internally by inversion
#' for the GEE code paths.
#'
#' @examples
#'
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
  eig <- eigen(mat, symmetric = TRUE)
  wts <- ifelse(eig$values > 0,
                eig$values^-0.25,
                0)
  crossprod(t(eig$vectors) * wts)
}

#' Compute Euclidean distance matrix for a cluster block
#'
#' Returns pairwise Euclidean distances within the cluster indexed by
#' \code{inds}. When \code{spacetime} has multiple columns, squared distances
#' are averaged across dimensions before taking the square root.
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
#'
#' Converts the unique within-block distances returned by
#' \code{.compute_dist_block()} to integer lags \code{0, 1, 2, \ldots} in
#' increasing order.
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
#' \code{glm_weight_function}, with the observation-weight contribution
#' carried separately by \eqn{\mathbf{D}}.  For canonical GLM families,
#' \eqn{w_i} is the usual IRWLS/Fisher-scoring weight on the mean--variance
#' scale; for example, logistic regression gives
#' \eqn{w_i = \mu_i(1-\mu_i)}.  The combined weighting entering the
#' information matrix is \eqn{\mathbf{W}\mathbf{D}}.  In Gaussian identity
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
#'   \mathrm{diag}\!\left(\sqrt{d_i}/\sqrt{w_i}\right)
#'   \mathbf{V}^{-1/2}(\mathbf{y} - \boldsymbol{\mu}),}
#' and \eqn{r_i =
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
#' @export
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


#' Safely Replace Variable Names in Printed Terms
#'
#' Replaces a variable name only when it appears as a standalone term,
#' avoiding accidental replacements inside other words.
#'
#' Useful for plot formula labels, for example replacing "t" with "that"
#' without changing "intercept" into "inthatercept".
#'
#' @param x Character vector to modify.
#' @param old Character scalar, original variable name.
#' @param new Character scalar, replacement variable name.
#'
#' @return Character vector with safe term-wise replacements.
#'
#' @keywords internal
safe_replace_var <- function(x, old, new) {
  old_esc <- gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", old)
  gsub(
    paste0("(?<![[:alnum:]_])", old_esc, "(?![[:alnum:]_])"),
    new,
    x,
    perl = TRUE
  )
}







