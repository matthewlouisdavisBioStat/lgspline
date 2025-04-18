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
create_onehot <- function(x) {
  mat <- model.matrix(~ x - 1)
  colnames(mat) <- sub("^x", "", colnames(mat))
  return(as.data.frame(mat))
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
#' @return Inverted matrix or identity matrix if all methods fail
#'
#' @details
#' Tries methods in order:
#'
#' 1. Direct inversion using \code{armaInv()}
#'
#' 2. Generalized inverse using eigendecomposition
#'
#' 3. Returns identity matrix with warning if both fail
#'
#' For eigendecomposition, uses a small ridge penalty (\code{1e-16}) for stability and
#' zeroes eigenvalues below machine precision.
#'
#' @examples
#' ## Well-conditioned matrix
#' A <- matrix(c(4,2,2,4), 2, 2)
#' invert(A) %**% A
#'
#' ## Singular matrix falls back to M.P. generalized inverse
#' B <- matrix(c(1,1,1,1), 2, 2)
#' invert(B) %**% B
#'
#' @keywords internal
#' @export
invert <- function(mat, include_warnings = FALSE){

  ## Try inversion
  t <- try({
    armaInv(mat)
  }, silent = TRUE)

  ## Try generalized inverse with small ridge penalty
  if(any(inherits(t, 'try-error'))){
    t <- try({
      eig <- eigen(mat + 1e-16*diag(nrow(mat)), symmetric = TRUE)
      d_inv <- ifelse(eig$values <= sqrt(.Machine$double.eps),
                      0,
                      1/eig$values)
      eig$vectors %**% (t(eig$vectors) * d_inv)
    },silent = TRUE)
  } else {
    return(t)
  }

  ## Return diagonal matrix with warning
  # Justification for identity is that for Newton-Raphson, this reduces
  # to approx. gradient-descent. But this is a terrible option for inference!
  if(any(inherits(t, 'try-error'))){
    if(include_warnings) warning('Matrix not inverted, returning identity: ',
                                 print(t))
    return(diag(nrow(mat)))
  } else {
    return(t)
  }
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
        gramMatrix(list_in[[k]])
      })
    } else {
      rem <- list()
    }

    # Process main chunks in parallel
    result <- c(
      Reduce("c", parallel::parLapply(cl, 1:num_chunks, function(chunk) {
        inds <- (chunk - 1)*chunk_size + 1:chunk_size
        lapply(inds, function(k) {
          gramMatrix(list_in[[k]])
        })
      })),
      rem
    )
  } else {
    result <- lapply(list_in, gramMatrix)
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

  if(parallel & !is.null(cl)) {
    ## Handle remainder chunks first
    if(rem_chunks > 0) {
      rem_indices <- num_chunks*chunk_size + 1:rem_chunks
      rem <- lapply(rem_indices, function(k) {
        if(unique_penalty_per_partition){
          G[[k]] %**% -negL_lambda[[k]] %**% G[[k]]
        } else {
          G[[k]] %**% negL_lambda %**% G[[k]]
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
            G[[k]] %**% -negL_lambda[[k]] %**% G[[k]]
          } else {
            G[[k]] %**% negL_lambda %**% G[[k]]
          }
        })
      })),
      rem
    )

  } else {
    ## Sequential computation
    result <- lapply(1:(K+1), function(k){
      if(unique_penalty_per_partition){
        G[[k]] %**% -negL_lambda[[k]]  %**% G[[k]]
      } else {
        G[[k]] %**% negL_lambda %**% G[[k]]
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
        sum(diag(dG_dlambda[[k]] %**% GXX[[k]]))))

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


#' Compute Derivative of \eqn{\textbf{U}\textbf{G}\textbf{X}^{T}\textbf{y}} with Respect to Lambda
#'
#' @param AGAInv_AGXy Product of \eqn{(\textbf{A}^{T}\textbf{G}\textbf{A})^{-1}} and \eqn{\textbf{A}^{T}\textbf{G}\textbf{X}^{T}\textbf{y}}
#' @param AGAInv Inverse of \eqn{\textbf{A}^{T}\textbf{G}\textbf{A}}
#' @param G List of \eqn{\textbf{G}} matrices
#' @param A Constraint matrix \eqn{\textbf{A}}
#' @param dG_dlambda List of \eqn{d\textbf{G}/d\lambda} matrices
#' @param nc Number of columns
#' @param nca Number of constraint columns
#' @param K Number of partitions minus 1 (\eqn{K})
#' @param Xy List of \eqn{\textbf{X}^{T}\textbf{y}} products
#' @param Ghalf List of \eqn{\textbf{G}^{1/2}} matrices
#' @param dGhalf List of \eqn{d\textbf{G}^{1/2}/d\lambda} matrices
#' @param GhalfXy_temp Temporary storage for \eqn{\textbf{G}^{1/2}\textbf{X}^{T}\textbf{y}}
#' @param parallel Use parallel processing
#' @param cl Cluster object
#' @param chunk_size Size of parallel chunks
#' @param num_chunks Number of chunks
#' @param rem_chunks Remaining chunks
#'
#' @details
#' Computes \eqn{d(\textbf{U}\textbf{G}\textbf{X}^{T}\textbf{y})/d\lambda}.
#' Uses efficient implementation avoiding large matrix construction.
#' For large problems (\eqn{K \ge 10}, \eqn{nc > 4}), uses chunked parallel processing.
#' For smaller problems, uses simpler least squares approach based on \eqn{\textbf{G}^{1/2}}.
#'
#' @return Vector of derivatives
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
    # Helper function for chunk processing
    process_chunk <- function(inds, G_chunk, dG_chunk, A_chunk, Xy_chunk) {
      # Process GXy and dGXy for chunk
      GXy_chunk <- lapply(seq_along(inds), function(i) {
        k <- inds[i]
        G_chunk[[i]] %**% Xy_chunk[[i]]
      })

      dGXy_chunk <- lapply(seq_along(inds), function(i) {
        k <- inds[i]
        dG_chunk[[i]] %**% Xy_chunk[[i]]
      })

      # Process AdGA for chunk
      AdGA_chunk <- Reduce("+", lapply(seq_along(inds), function(i) {
        k <- inds[i]
        start_idx <- (k-1)*nc + 1
        end_idx <- k*nc
        crossprod(A[start_idx:end_idx,], dG_chunk[[i]] %**%
                    A[start_idx:end_idx,])
      }))

      list(GXy = GXy_chunk, dGXy = dGXy_chunk, AdGA = AdGA_chunk)
    }

    if(parallel & !is.null(cl)) {
      # Handle remainder chunks
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

      # Process main chunks in parallel
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

      # Combine results
      GXy <- Reduce("c", c(lapply(chunk_results, `[[`, "GXy"),
                           rem_result$GXy))
      dGXy <- Reduce("c", c(lapply(chunk_results, `[[`, "dGXy"),
                            rem_result$dGXy))
      AdGA <- Reduce("+", c(lapply(chunk_results, `[[`, "AdGA"),
                            list(rem_result$AdGA)))

    } else {
      # Sequential computation
      GXy <- Reduce("c", lapply(1:(K+1), function(k) G[[k]] %**% Xy[[k]]))
      dGXy <- Reduce("c", lapply(1:(K+1), function(k) dG_dlambda[[k]] %**%
                                   Xy[[k]]))
      AdGA <- Reduce("+", lapply(1:(K+1), function(k) {
        start_idx <- (k-1)*nc + 1
        end_idx <- k*nc
        crossprod(A[start_idx:end_idx,], dG_dlambda[[k]] %**%
                    A[start_idx:end_idx,])
      }))

    }

    ## Common intermediate matrices that don't need recomputing in each chunk
    A_AGAInv_AGXy <- cbind(A %**% AGAInv_AGXy)
    A_term2b_mat <- cbind(A %**% (-AGAInv %**% AdGA %**% AGAInv_AGXy))
    A_term2c_mat <- A %**% cbind(AGAInv %**% crossprod(A, cbind(unlist(dGXy))))

    # Process terms in parallel
    if(parallel & !is.null(cl)) {
      # Handle remainder chunks
      if(rem_chunks > 0) {
        rem_indices <- num_chunks*chunk_size + 1:rem_chunks

        term2a_rem <- lapply(rem_indices, function(k) {
          dG_dlambda[[k]] %**% A_AGAInv_AGXy[(k-1)*nc + 1:nc,,drop=FALSE]
        })

        term2b_rem <- lapply(rem_indices, function(k) {
          start_idx <- (k-1)*nc + 1
          end_idx <- k*nc
          G[[k]] %**% A_term2b_mat[start_idx:end_idx,, drop=FALSE]
        })

        term2c_rem <- lapply(rem_indices, function(k) {
          start_idx <- (k-1)*nc + 1
          end_idx <- k*nc
          G[[k]] %**% A_term2c_mat[start_idx:end_idx,, drop=FALSE]
        })
      } else {
        term2a_rem <- term2b_rem <- term2c_rem <- list()
      }

      # Process main chunks in parallel
      chunk_results <- parallel::parLapply(cl, 1:num_chunks, function(chunk) {
        start_idx <- (chunk - 1) * chunk_size
        inds <- (start_idx + 1):min(start_idx + chunk_size, K+1)

        # Compute all three terms for this chunk
        terms <- lapply(inds, function(k){
          t2a <- dG_dlambda[[k]] %**%
            A_AGAInv_AGXy[(k-1)*nc + 1:nc,, drop=FALSE]
          idx <- (k-1)*nc + 1:nc
          t2b <- G[[k]] %**% A_term2b_mat[idx,, drop=FALSE]
          t2c <- G[[k]] %**% A_term2c_mat[idx,, drop=FALSE]
          list(t2a,t2b,t2c)
        })
        term2a_chunk <- lapply(terms, `[[`, 1)
        term2b_chunk <- lapply(terms, `[[`, 2)
        term2c_chunk <- lapply(terms, `[[`, 3)

        list(term2a = term2a_chunk,
             term2b = term2b_chunk,
             term2c = term2c_chunk)
      })

      # Combine results
      term2a <- Reduce("c", c(lapply(chunk_results, `[[`, "term2a"),
                              term2a_rem))
      term2b <- Reduce("c", c(lapply(chunk_results, `[[`, "term2b"),
                              term2b_rem))
      term2c <- Reduce("c", c(lapply(chunk_results, `[[`, "term2c"),
                              term2c_rem))

    } else {
      # Sequential computation
      terms <- lapply(1:(K+1), function(k){
        t2a <- dG_dlambda[[k]] %**% A_AGAInv_AGXy[(k-1)*nc + 1:nc,,drop=FALSE]
        idx <- (k-1)*nc + 1:nc
        t2b <- G[[k]] %**% A_term2b_mat[idx,, drop=FALSE]
        t2c <- G[[k]] %**% A_term2c_mat[idx,, drop=FALSE]
        list(t2a,t2b,t2c)
      })
      term2a <- lapply(terms, `[[`, 1)
      term2b <- lapply(terms, `[[`, 2)
      term2c <- lapply(terms, `[[`, 3)
    }

    # Combine final results
    term2a <- Reduce("c", term2a)
    term2b <- Reduce("c", term2b)
    term2c <- Reduce("c", term2c)
    return(cbind(unlist(dGXy) - (unlist(term2a) +
                                   unlist(term2b) +
                                   unlist(term2c))))
  } else {

    ## Transform to least squares problem using G^(1/2)
    GhalfA <- Reduce("rbind", GAmult_wrapper(Ghalf,
                                             A,
                                             K,
                                             nc,
                                             nca,
                                             parallel,
                                             cl,
                                             chunk_size,
                                             num_chunks,
                                             rem_chunks))

    ## Transform to least squares problem using d[G^(1/2)]/dlambda;
    dGhalfA <- Reduce("rbind", GAmult_wrapper(dGhalf,
                                              A,
                                              K,
                                              nc,
                                              nca,
                                              parallel,
                                              cl,
                                              chunk_size,
                                              num_chunks,
                                              rem_chunks))
    dGhalfXy <- cbind(unlist(
      matmult_block_diagonal(dGhalf,
                             Xy,
                             K,
                             parallel,
                             cl,
                             chunk_size,
                             num_chunks,
                             rem_chunks)))

    ## Get first residuals
    comp_stab_sc <- sqrt(mean(abs(GhalfA)))
    comp_stab_sc <- 1*(comp_stab_sc == 0) + comp_stab_sc
    resids1 <- do.call('.lm.fit', list(x = GhalfA / comp_stab_sc,
                                       y = GhalfXy_temp / comp_stab_sc)
    )$residuals * comp_stab_sc

    ## Get second residuals
    resids2 <- do.call('.lm.fit', list(x = dGhalfA / comp_stab_sc,
                                       y = dGhalfXy / comp_stab_sc)
    )$residuals * comp_stab_sc

    ## Compute b2 in parallel if requested
    if(parallel & !is.null(cl)) {
      # Handle remainder chunks for b2
      if(rem_chunks > 0) {
        rem_indices <- num_chunks * chunk_size + 1:rem_chunks
        b2_rem <- lapply(rem_indices, function(k) {
          dGhalf[[k]] %**% cbind(resids2[(k-1)*nc + 1:nc]) +
            Ghalf[[k]] %**% cbind(resids1[(k-1)*nc + 1:nc])
        })
      } else {
        b2_rem <- list()
      }

      ## Process main chunks for b2 in parallel
      b2_result <- c(
        Reduce("c",parallel::parLapply(cl, 1:num_chunks, function(chunk) {
          inds <- (chunk - 1) * chunk_size + 1:chunk_size
          lapply(inds, function(k) {
            dGhalf[[k]] %**% cbind(resids2[(k-1)*nc + 1:nc]) +
              Ghalf[[k]] %**% cbind(resids1[(k-1)*nc + 1:nc])
          })
        })),
        b2_rem
      )
      b2 <- unlist(b2_result)
    } else {
      b2 <- Reduce("c", lapply(1:(K+1), function(k) {
        dGhalf[[k]] %**% cbind(resids2[(k-1)*nc + 1:nc]) +
          Ghalf[[k]] %**% cbind(resids1[(k-1)*nc + 1:nc])
      }))
    }

    ## Add components to get final result
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
#' @param wiggle_penalty,flat_ridge_penalty Initial penalty values
#' @param invsoftplus_initial_wiggle,invsoftplus_initial_flat Initial grid search values (log scale)
#' @param unique_penalty_per_predictor,unique_penalty_per_partition Logical; allow predictor/partition-specific penalties
#' @param invsoftplus_penalty_vec Initial values for predictor/partition penalties (log scale)
#' @param meta_penalty The "meta" ridge penalty, a regularization for predictor/partition penalties to pull them on log-scale towards 0 (1 on raw scale)
#' @param keep_weighted_Lambda,iterate Logical controlling GLM fitting
#' @param qp_score_function,quadprog,qp_Amat,qp_bvec,qp_meq Quadratic programming parameters (see arguments of \code{\link[lgspline]{lgspline}})
#' @param tol Numeric; convergence tolerance
#' @param sd_y,delta Response standardization parameters
#' @param constraint_value_vectors List; constraint values
#' @param parallel Logical; enable parallel computation
#' @param cl,chunk_size,num_chunks,rem_chunks Parallel computation parameters
#' @param custom_penalty_mat Optional custom penalty matrix
#' @param order_list List; observation ordering by partition
#' @param glm_weight_function,shur_correction_function Functions for GLM weights and corrections
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
    invsoftplus_initial_wiggle,
    invsoftplus_initial_flat,
    unique_penalty_per_predictor,
    unique_penalty_per_partition,
    invsoftplus_penalty_vec,
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
    shur_correction_function,
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

  ## For compatibility without multiple predictors nor multiple partitions
  if(length(invsoftplus_penalty_vec) == 0){
    penalty_vec <- c()
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
    # 2 penalties of interest
    wiggle_penalty <- unlist(softplus(par[1])) # smoothing spline f''(x)^2
    flat_ridge_penalty <- unlist(softplus(par[2])) # flat ridge regression penalty
    if(unique_penalty_per_predictor | unique_penalty_per_partition){
      penalty_vec <- softplus(c(par[-c(1:2)]))
    } else {
      penalty_vec < c()
    }

    # Reparameterize
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

    ## Shur complements are not needed here, only in get_B
    # This is for initialization only
    shur_corrections <- lapply(1:(K+1), function(k)0)
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
                              shur_corrections)

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
               shur_correction_function,
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
    sum_W <- compute_trace_UGXX_wrapper(G_list$G,
                                        A,
                                        GXX,
                                        AGAInv,
                                        nc,
                                        nca,
                                        K,
                                        parallel & parallel_trace,
                                        cl,
                                        chunk_size,
                                        num_chunks,
                                        rem_chunks)

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
         (paste0(family)[2] == 'identity' | paste0(family)[1] == 'gaussian')){
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
    ## 2 penalties of interest
    # smoothing spline 2nd derivative squared penalty
    wiggle_penalty <- unlist(softplus(par[1])) # smoothing spline f''(x)^2
    flat_ridge_penalty <- unlist(softplus(par[2])) # flat ridge regression penalty
    if(unique_penalty_per_predictor | unique_penalty_per_partition){
      penalty_vec <- softplus(c(par[-c(1:2)]))
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
      shur_corrections <- lapply(1:(K+1), function(k)0)
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
                                shur_corrections)

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
                     shur_correction_function,
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
      sum_W <- compute_trace_UGXX_wrapper(G_list$G,
                                          A,
                                          GXX,
                                          AGAInv,
                                          nc,
                                          nca,
                                          K,
                                          parallel & parallel_trace,
                                          cl,
                                          chunk_size,
                                          num_chunks,
                                          rem_chunks)

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

    ## Chain rule for computing gradients of other elements of the penalty
    dGCV_u_dlambda2 <- mean(diag(outlist$L2))/
      mean(diag(outlist$Lambda)) *
      dGCV_u_dlambda1

    # Gradient
    if(verbose){
      cat('        Gradient start \n')
    }
    gradient <- cbind(c(dGCV_u_dlambda1 * lambda_1 / (1 + lambda_1),
                        dGCV_u_dlambda2 * lambda_2/ ( 1 + lambda_2)))

    ## Gradients for unique penalties per predictor
    if(unique_penalty_per_predictor){
      predictor_penalties <-
        penalty_vec[grep('predictor',names(penalty_vec))]

      predictor_penalty_gradient <-
        sapply(1:length(predictor_penalties),function(j){
        mean(diag(outlist$L_predictor_list[[j]]))/
        mean(diag(outlist$Lambda)) * # ratio of trace measures contribution
        dGCV_u_dlambda1 *
        predictor_penalties[j] /
        (1 + predictor_penalties[j]) # chain rule, the derivative being passed on
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
          partition_penalties[j] / # chain rule         # partition penalties
          (1 + partition_penalties[j])
        })                                            # have not been added to
                                                      # Lambda yet
      ## Combined gradient
      gradient <- cbind(c(gradient, partition_penalty_gradient))
    }

    ## Regularization penalty - pull penalty terms to 1
    if(unique_penalty_per_partition | unique_penalty_per_predictor){
      regulizer <- c(1e-32*(wiggle_penalty - 1),
                     0,
                     meta_penalty*(penalty_vec - 1))
      meta_penalty <- 0.5*meta_penalty*sum((penalty_vec - 1)^2) +
                      0.5*1e-32*((wiggle_penalty - 1)^2)

    } else {
      regulizer <- c(1e-32*(wiggle_penalty - 1),
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
        result <- gr_fxn(c(lambda[1],lambda[2], invsoftplus_penalty_vec), outlist)
        gradient <- result$gradient
        outlist <- result$outlist
      }

      ## For first two iterations, use steepest descent
      if(iter <= 2){
        new_lambda <- lambda - damp * gradient

        ## Reset to best solution if numerical issues
        if(any(!is.finite(softplus(new_lambda))) |
           any(is.nan(softplus(new_lambda))) |
           any(is.na(softplus(new_lambda)))){
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

        ## Reset to best solution if numerical issues
        if(any(!is.finite(softplus(new_lambda))) |
           any(is.nan(softplus(new_lambda))) |
           any(is.na(softplus(new_lambda)))){
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
                              invsoftplus_penalty_vec))
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
    initial_grid <- expand.grid(wiggle = invsoftplus_initial_wiggle,
                                flat = invsoftplus_initial_flat)

    ## Function to safely evaluate gcv_u
    safe_gcvu <- function(par) {
      tryCatch({
        result <-
          gcvu_fxn(c(unlist(par), invsoftplus_penalty_vec))$GCV_u
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
        try(quasi_nr_fxn(c(best_start, invsoftplus_penalty_vec)), silent = TRUE),
        warning = function(w) if (include_warnings) warning(w) else invokeRestart("muffleWarning"),
        message = function(m) if (include_warnings) message(m) else invokeRestart("muffleMessage")
      )
      if(any(inherits(res, 'try-error'))){
        if(include_warnings) print(res)
        if(include_warnings) warning('Custom BFGS implementation failed. Try use_custom_bfgs =',
        ' FALSE, or manual tuning. Resorting to best as selected from',
        ' grid search.')
        par <- c(best_start, invsoftplus_penalty_vec)
      } else {
        par <- res$par
      }
    } else {
      ## vs. base R (finite-difference approx.)
      res <- withCallingHandlers(
        try({optim(c(best_start, invsoftplus_penalty_vec),
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
        par <- c(best_start, invsoftplus_penalty_vec)
      } else {
        par <- res$par
      }
    }
    if(verbose){
      cat('    Finished tuning penalties\n')
    }

    ## Inflate the penalties, since there's an inherent bias towards
    # no penalization in-sample
    infl <- ((nr+2)/((nr-2)))^2
    wiggle_penalty <- softplus(par[1])*infl
    flat_ridge_penalty <- softplus(par[2])*infl

    ## Update penalty vec for predictor-and-partition specific penalties
    if(length(invsoftplus_penalty_vec) > 0){
      penalty_vec <- softplus(invsoftplus_penalty_vec)
      penalty_vec[1:length(penalty_vec)] <-
        softplus(c(par[-c(1,2)]))*infl
    }
  } else if(length(invsoftplus_penalty_vec) > 0){
    penalty_vec <- softplus(invsoftplus_penalty_vec)
  } else {
    penalty_vec <- c()
  }


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

#' Get Constrained GLM Coefficient Estimates
#'
#' @description
#' Core estimation function for Lagrangian multiplier smoothing splines. Computes
#' coefficient estimates under smoothness constraints and penalties for GLMs, handling
#' four distinct computational paths depending on model structure.
#'
#' @param X List of design matrices by partition
#' @param X_gram List of Gram matrices by partition
#' @param Lambda Combined penalty matrix (smoothing spline + ridge)
#' @param keep_weighted_Lambda Logical; retain GLM weights in smoothing spline penalty
#' @param unique_penalty_per_partition Logical; whether to use partition-specific penalties
#' @param L_partition_list List of partition-specific penalty matrices
#' @param A Matrix of smoothness constraints (continuity, differentiability)
#' @param Xy List of X^T y products by partition
#' @param y List of responses by partition
#' @param K Integer; number of knots (partitions = K+1)
#' @param nc Integer; number of coefficients per partition
#' @param nca Integer; number of columns in constraint matrix
#' @param Ghalf List of G^(1/2) matrices by partition
#' @param GhalfInv List of G^(-1/2) matrices by partition
#' @param parallel_eigen,parallel_aga,parallel_matmult,parallel_unconstrained Logical flags for parallel computation components
#' @param cl Cluster object for parallel processing
#' @param chunk_size,num_chunks,rem_chunks Parameters for chunking in parallel processing
#' @param family GLM family object (includes link function and variance functions)
#' @param unconstrained_fit_fxn Function for obtaining unconstrained estimates
#' @param iterate Logical; whether to use iterative optimization for non-linear links
#' @param qp_score_function Function; see description in \code{\link[lgspline]{lgspline}}. Accepts arguments in order "X, y, mu, order_list, dispersion, VhalfInv, ...".
#' @param quadprog Logical; whether to use quadratic programming for inequality constraints
#' @param qp_Amat,qp_bvec,qp_meq Quadratic programming constraint parameters
#' @param prevB List of previous coefficient estimates for warm starts
#' @param prevUnconB List of previous unconstrained coefficient estimates
#' @param iter_count,prev_diff Iteration tracking for convergence
#' @param tol Numeric; convergence tolerance
#' @param constraint_value_vectors List of additional constraint values
#' @param order_list List of observation orderings by partition
#' @param glm_weight_function Function for computing GLM weights
#' @param shur_correction_function Function for uncertainty corrections in information matrix
#' @param need_dispersion_for_estimation Logical; whether dispersion needs estimation
#' @param dispersion_function Function for estimating dispersion parameter
#' @param observation_weights Optional observation weights by partition
#' @param homogenous_weights Logical; whether weights are constant
#' @param return_G_getB Logical; whether to return G matrices with coefficient estimates
#' @param blockfit Logical; whether to use block-fitting approach for special structure
#' @param just_linear_without_interactions Vector of columns for non-spline effects
#' @param Vhalf,VhalfInv Square root and inverse square root correlation matrices for GEE fitting
#' @param ... Additional arguments passed to fitting functions
#'
#' @details
#' This function implements the method of Lagrangian multipliers for fitting constrained
#' generalized linear models with smoothing splines. The function follows one of four main
#' computational paths depending on the model structure:
#'
#' \bold{1. Pure GEE (No Blockfitting):}
#' When Vhalf and VhalfInv are provided without blockfitting, the function uses generalized
#' estimating equations to handle correlated data. The design matrix is arranged in
#' block-diagonal form, and the estimation explicitly accounts for the correlation structure
#' provided. Uses sequential quadratic programming (SQP) for optimization.
#'
#' \bold{2. Blockfitting (With or Without Correlation):}
#' When blockfit=TRUE and linear-only terms are specified, the function separates spline effects
#' from linear-only terms. The design matrix is restructured with spline terms in block-diagonal
#' form and linear terms stacked together. This structure is particularly efficient when there
#' are many non-spline effects. The function can handle both correlated (GEE) and uncorrelated
#' data in this path.
#'
#' \bold{3. Canonical Gaussian, No Correlation:}
#' For Gaussian family with identity link and no correlation structure, calculations are greatly
#' simplified. No unconstrained fit function is needed; estimation uses direct matrix operations.
#' This path takes advantage of the closed-form solution available for linear models with
#' Gaussian errors.
#'
#' \bold{4. GLMs, No Correlation:}
#' For non-Gaussian GLMs without correlation structure, the function requires an unconstrained_fit_fxn
#' to obtain initial estimates. It may use iterative fitting for non-canonical links. This path
#' first computes unconstrained estimates for each partition, then applies constraints using the
#' method of Lagrangian multipliers.
#'
#' All paths use efficient matrix decompositions and avoid explicitly constructing full matrices
#' when possible. For non-canonical links, iterative fitting is employed to converge to the optimal
#' solution.
#'
#' @return
#' For `return_G_getB = FALSE`: A list of coefficient vectors by partition.
#'
#' For `return_G_getB = TRUE`: A list with elements:
#' \describe{
#'   \item{B}{List of coefficient vectors by partition}
#'   \item{G_list}{List containing G matrices (covariance matrices), Ghalf (square root),
#'                 and GhalfInv (inverse square root)}
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
                  shur_correction_function,
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

  ## For generalized estimating equations only
  if((!is.null(Vhalf) &
     (!is.null(VhalfInv))) &
     !(blockfit & length(just_linear_without_interactions) > 0)){
     ## Re-order to match partition ordering
     Vhalf <- Vhalf[unlist(order_list), unlist(order_list)]
     VhalfInv <- VhalfInv[unlist(order_list), unlist(order_list)]
    if(is.null(observation_weights) | length(unlist(observation_weights)) == 0){
      observation_weights <- 1
    }

    ## New design matrix in block-diagonal form
    X_block <- collapse_block_diagonal(X)

    ## Initial beta coefficients is the 0 vector
    beta_block <- cbind(rep(0, ncol(X_block)))

    ## Penalties
    if(unique_penalty_per_partition){
      Lambda_block <- Reduce("rbind", lapply(1:(K+1),function(k){
        Reduce("cbind", lapply(1:(K+1), function(j){
          if(j == k) Lambda +
            L_partition_list[[k]] else 0*
            Lambda
        }))
      }))
    } else {
      Lambda_block <- Reduce("rbind", lapply(1:(K+1), function(k){
        Reduce("cbind", lapply(1:(K+1), function(j){
          if(j == k) Lambda else 0*
            Lambda
        }))
      }))
    }

    ## Response (recall this is VhalfInv %**% y right now)
    # The whitening transform was already applied
    y_block <- cbind(unlist(y))
    Vhalfy <- Vhalf %**% y_block

    ## Solve smoothing constraints and qp constraints simultaneously
    damp_cnt <- 0
    master_cnt <- 0
    err <- Inf

    ## Initial link-transformed predictions
    XB <- Vhalf %**% (X_block %**% beta_block)

    ## Dummy qp objects and constraint matrix A, if not available
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

    ## Invoke procedure for fitting SQP problems
    while(err > tol & damp_cnt < 10 & master_cnt < 100){
      ## Initialize counters
      master_cnt <- master_cnt + 1
      damp <- 2^(-(damp_cnt))

      ## Information matrix in block-diagonal form
      if(need_dispersion_for_estimation){
        dispersion_temp <- dispersion_function(
          family$linkinv(XB),
          Vhalfy,
          unlist(order_list),
          family,
          unlist(observation_weights),
          ...
        )
      } else {
        dispersion_temp <- 1
      }
      W <- c(glm_weight_function(family$linkinv(XB),
                                 Vhalfy,
                                 unlist(order_list),
                                 family,
                                 dispersion_temp,
                                 unlist(observation_weights),
                                 ...))

      ## Re-package into partition form for shur correction
      result <- lapply(1:(K+1),function(k){
        cbind(beta_block[(k-1)*nc + 1:nc])
      })
      shur_correction <-
        shur_correction_function(
          list(Vhalf %**% X_block),
          list(cbind(Vhalfy)),
          list(cbind(unlist(result))),
          dispersion_temp,
          list(unlist(order_list)),#order_list,
          0,#K,
          family,
          unlist(observation_weights),
          ...
        )
      if(!(any(unlist(shur_correction) != 0))){
        shur_correction_collapsed <- 0
        shur_correction <- 0
      } else {
        shur_correction_collapsed <- collapse_block_diagonal(shur_correction)
      }

      ## Information matrix (G^{-1})
      info <- t(X_block) %**% (W * X_block) +
        Lambda_block +
        shur_correction_collapsed

      ## Numerical stability
      sc <- sqrt(mean(abs(info)))

      ## Initialize without inequality constraint, single Newton-Raphson step
      if(master_cnt == 1){
        infoinv_block <- sc*invert(sc*info)
        U <- (diag(1, nrow(info)) -
                infoinv_block %**%
                A %**%
                invert(t(A) %**% infoinv_block %**% A) %**%
                t(A))
        ## = UG[X^{T}(y-g(XB))]
        # Lambda_block %**% beta_block not needed since beta_block initializes
        # to 0
        beta_new <- c(
            U %**%
            infoinv_block %**%
            t(X_block) %**%
            (y_block - VhalfInv %**% cbind(family$linkinv(c(XB))))
        )
        infoinv_block <- NULL
      } else {
        ## Else QP update with inequality constraints
        qp_score <- qp_score_function(
          X_block,
          y_block,
          VhalfInv %**% cbind(family$linkinv(XB)),
          unlist(order_list),
          dispersion_temp,
          VhalfInv,
          unlist(observation_weights),
          ...
        )
        beta_new <- try({quadprog::solve.QP(
          Dmat = info/sc,
          dvec = (qp_score -
                    Lambda_block %**% beta_block +
                    info %**% beta_block)/sc,
          Amat = qp_Amat,
          bvec = qp_bvec,
          meq = qp_meq
        )$solution}, silent = TRUE)
      }

      ## quadprog is very sensitive to positive definiteness
      if(any(inherits(beta_new, 'try-error'))){
        beta_new <- 0*beta_block
      }

      ## If no iteration, return solution as non-damped solution
      if((!iterate & master_cnt > 2) |
         (paste0(family)[1] == 'gaussian' &
          paste0(family)[2] == 'identity')){
        beta_block <- beta_new
        damp_cnt <- 11
        master_cnt <- 101
        err <- tol - 1

        ## Else, use damped SQP
      } else {
        ## Damped update
        beta_new <- (1-damp)*beta_block + damp*beta_new
        XB <- Vhalf %**% (X_block %**% cbind(beta_new))
        if(need_dispersion_for_estimation){
          dispersion_temp <- dispersion_function(
            family$linkinv(XB),
            Vhalfy,
            unlist(order_list),
            family,
            unlist(observation_weights),
            ...
          )
        } else {
          dispersion_temp <- 1
        }
        W <- c(glm_weight_function(family$linkinv(XB),
                                   Vhalfy,
                                   unlist(order_list),
                                   family,
                                   dispersion_temp,
                                   unlist(observation_weights),
                                   ...))

        ## Priority goes
        # 1) family$custom_dev.resids if available
        # 2) family$dev.resids if available
        # 3) MSE otherwise
        if(!is.null(family$custom_dev.resids)){
          raw <- family$custom_dev.resids(Vhalfy,
                                          family$linkinv(c(XB)),
                                          unlist(order_list),
                                          family,
                                          unlist(observation_weights),
                                          ...)
          err_new <- mean((
                        VhalfInv %**% cbind(sign(raw)*sqrt(abs(
                          raw
                        ))) / sqrt(c(W))
                      )^2)
        } else if(is.null(family$dev.resids)){
          err_new <- mean((unlist(observation_weights)*
                         (y_block - VhalfInv %**% cbind(family$linkinv(XB))))^2)
        } else {
          err_new <-
            mean(unlist(observation_weights)*
                   family$dev.resids(y_block,
                                     VhalfInv %**% cbind(family$linkinv(XB)),
                                     wt = 1/W))
        }
        if(is.null(err_new) | is.na(err_new) | !is.finite(err_new)){
            damp_cnt <- damp_cnt + 1
          } else if(err_new <= err){

          ## Update coefficients, set damp to 0
          prev_err <- err
          err <- err_new
          abs_diff <- max(abs(beta_new - beta_block))
          beta_block <- beta_new
          damp_cnt <- 0

          ## Convergence criteria met when improvement in performance is
          # negligible and change in beta coefficients is small and
          # more than 10 iterations have passed
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

      ## Re-package into partition form
      result <- lapply(1:(K+1),function(k){
        cbind(beta_block[(k-1)*nc + 1:nc])
      })
    }
    ## Return output, with G re-computed if desired
    if(return_G_getB){
      if(paste0(family)[1] == 'gaussian' &
         paste0(family)[2] == 'identity'){
        return(list(
          B = result,
          G_list = list(G = lapply(Ghalf, function(mat) mat %**% mat),
                        Ghalf = Ghalf)
        ))
      }
      ## Update estimates and variance components
      XB <- Vhalf %**% (X_block %**% cbind(unlist(result)))
      if(need_dispersion_for_estimation){
        dispersion_temp <- dispersion_function(
          family$linkinv(XB),
          Vhalfy,
          unlist(order_list),
          family,
          unlist(observation_weights),
          ...
        )
      } else {
        dispersion_temp <- 1
      }
      W <- c(glm_weight_function(family$linkinv(XB),
                                 Vhalfy,
                                 unlist(order_list),
                                 family,
                                 dispersion_temp,
                                 unlist(observation_weights),
                                 ...))

      shur_correction <-
        shur_correction_function(
          list(Vhalf %**% X_block),
          list(cbind(Vhalfy)),
          list(cbind(unlist(result))),
          dispersion_temp,
          list(unlist(order_list)),#order_list,
          0,#K,
          family,
          unlist(observation_weights),
          ...
        )
      if(!(any(unlist(shur_correction) != 0))){
        shur_correction_collapsed <- 0
        shur_correction <- 0
      } else {
        shur_correction_collapsed <- collapse_block_diagonal(shur_correction)
      }

      info <- t(X_block) %**% (W * X_block) +
        Lambda_block +
        shur_correction_collapsed

      G <- lapply(1:(K+1), function(k)NA)
      Ghalf <- lapply(1:(K+1), function(k)NA)
      GhalfInv <- lapply(1:(K+1), function(k)NA)
      for(k in 1:(K+1)){
        G[[k]] <- info[(k-1)*nc + 1:nc, (k-1)*nc + 1:nc]
        eig <- eigen(G[[k]], symmetric = TRUE)
        vals <- eig$values
        vals[eig$values <= sqrt(.Machine$double.eps)] <- 0
        sqrt_vals <- sqrt(vals)
        Ghalf[[k]] <- eig$vectors %**% (t(eig$vectors) * sqrt_vals)
        inv_sqvals <- ifelse(sqrt_vals > 0, 1/sqrt_vals, 0)
        GhalfInv[[k]] <- eig$vectors %**% (t(eig$vectors) * inv_sqvals)
      }

      G_list <- list(
        Ghalf = Ghalf,
        GhalfInv = GhalfInv,
        G = G
      )

      return(list(
        B = result,
        G_list = G_list
      ))
    } else {
      return(result)
    }
  }

  ## Fit in block-diagonal form if desired with many independent spline effects
  # Not parallelized nor memory efficient, except for case of many
  # non-interactive non-spline effects and few spline effects with many knots
  # Any linear constraints upon non-spline effects without interactions will be
  # lost.
  # Compatible with (in fact, assumes) GEE scenario
  if(blockfit &
     length(just_linear_without_interactions) > 0 &
     K > 0){
    if(is.null(observation_weights) | length(unlist(observation_weights)) == 0){
      observation_weights <- 1
    }
    ## Assumes GEE scenario
    # replaces correlation structure Vwith identity matrix if missing
    if(is.null(Vhalf) | is.null(VhalfInv)){
      cor_fl <- FALSE
      #Vhalf <- diag(nrow(cbind(unlist(y))))
      #VhalfInv <- Vhalf
    } else {
      cor_fl <- TRUE
      ## Re-order to match the block-partition ordering
      Vhalf <- Vhalf[unlist(order_list), unlist(order_list)]
      VhalfInv <- VhalfInv[unlist(order_list), unlist(order_list)]
    }

    ## Extract just-linear-without-interaction expansions
    # one partition
    flat_cols <- unlist(
                    c(sapply(paste0('_',just_linear_without_interactions,'_'),
                        function(col)which(colnames(X[[1]]) == col))))
    # all partitions
    flat_cols_collapsed <- Reduce("c",
                                  sapply(1:(K+1),function(kk){
                                    nc*(kk-1) + flat_cols
                                  }))

    ## New design matrix in block-diagonal form,
    # with linear-only terms (X2) binded together agnostic to partition
    # to the usual spline effects (X1)
    X1 <- Reduce("rbind", lapply(1:(K+1),function(k){
      Reduce("cbind",lapply(1:(K+1),function(j){
        if(nrow(X[[k]]) == 0){
          return(X[[k]][,-flat_cols,drop=FALSE])
        } else if(j == k) X[[k]][,-flat_cols,drop=FALSE] else 0*
                          X[[k]][,-flat_cols,drop=FALSE]
      }))
    }))
    X2 <- Reduce("rbind", lapply(1:(K+1),function(k){
      return(X[[k]][,flat_cols,drop=FALSE])
    }))
    X_block <- cbind(X1, X2)
    X1 <- NULL
    X2 <- NULL

    ## Initial beta coefficients is the 0 vector
    beta_block <- cbind(rep(0, ncol(X_block)))

    ## Penalties
    # Fix for Lambda_block construction
    # First extract components correctly
    nc_flat <- length(flat_cols)
    nc_minus_flat <- nc - nc_flat

    # Create correct non-flat block diagonal
    Lambda_non_flat <- Lambda[-flat_cols, -flat_cols]
    Lambda_flat <- Lambda[flat_cols, flat_cols]

    # Create block-diagonal matrix for non-flat penalties
    Lambda_block1 <- Reduce("rbind", lapply(1:(K+1), function(k){
      Reduce("cbind", lapply(1:(K+1), function(j){
        if(j == k) {
          if(unique_penalty_per_partition) {
            Lambda_non_flat + L_partition_list[[k]][-flat_cols, -flat_cols]
          } else {
            Lambda_non_flat
          }
        } else {
          matrix(0, nrow=nrow(Lambda_non_flat), ncol=ncol(Lambda_non_flat))
        }
      }))
    }))

    # Keep flat penalties as is
    Lambda_block2 <- Lambda_flat

    # Build correctly sized Lambda_block
    Lambda_block <- matrix(0, nrow=ncol(X_block), ncol=ncol(X_block))
    Lambda_block[1:nrow(Lambda_block1), 1:ncol(Lambda_block1)] <- Lambda_block1
    Lambda_block[(nrow(Lambda_block1)+1):(ncol(X_block)),
                 (ncol(Lambda_block1)+1):(ncol(X_block))] <- Lambda_block2

    y_block <- cbind(unlist(y))
    if(cor_fl){
      Vhalfy <- Vhalf %**% y_block
    } else {
      Vhalfy <- y_block
    }

    ## Solve smoothing constraints and qp constraints simultaneously
    damp_cnt <- 0
    master_cnt <- 0
    err <- Inf

    ## Initial link-transformed predictions
    if(cor_fl){
      XB <- Vhalf %**% (X_block %**% beta_block)
    } else {
      XB <- X_block %**% beta_block
    }

    ## Dummy qp objects and constraint matrix A, if not available
    # Below merges equality and inequality constraints together
    if(K == 0 | any(all.equal(unique(A), 0) == TRUE)){
      A <- cbind(rep(0, nc*(K+1)))
    }
    At <- t(A)
    A0 <- At[,-flat_cols_collapsed,drop=FALSE]
    A1 <- At[,flat_cols] / (K + 1)
    for(k in 2:(K+1)){
      A1 <- A1 + At[,flat_cols + (k-1)*nc] / (K + 1)
    }
    At <- cbind(A0, A1)
    A <- t(At)
    A1 <- NULL
    A0 <- NULL
    if(!is.null(qp_Amat)){
      qp_At <- t(qp_Amat)
      qp_A0 <- qp_At[,-flat_cols_collapsed,drop=FALSE]
      qp_A1 <- qp_At[,flat_cols] / (K + 1)
      for(k in 2:(K+1)){
        qp_A1 <- qp_A1 + qp_At[,flat_cols + (k-1)*nc] / (K + 1)
      }
      qp_At <- cbind(qp_A0, qp_A1)
      qp_Amat <- t(qp_At)
      qp_Amat <- cbind(A, qp_Amat)
    } else {
      qp_Amat <- A
    }
    qp_A1 <- NULL
    qp_A0 <- NULL

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

    ## Invoke procedure for fitting SQP problems
    nc_flat <- length(flat_cols)
    nc_minus_flat <- nc - nc_flat
    while(err > tol & damp_cnt < 10 & master_cnt < 100){
      ## Initialize counters
      master_cnt <- master_cnt + 1
      damp <- 2^(-(damp_cnt))

      ## Information matrix in block-diagonal form
      if(need_dispersion_for_estimation){
        dispersion_temp <- dispersion_function(
          family$linkinv(XB),
          Vhalfy,
          unlist(order_list),
          family,
          unlist(observation_weights),
          ...
        )
      } else {
        dispersion_temp <- 1
      }
      W <- c(glm_weight_function(family$linkinv(XB),
                                 Vhalfy,
                                 unlist(order_list),
                                 family,
                                 dispersion_temp,
                                 unlist(observation_weights),
                                 ...))

      ## Re-package into partition form for shur correction
      result <- lapply(1:(K+1),function(k){
        beta <- rep(0, nc)
        beta[-flat_cols] <- beta_block[(k-1)*nc_minus_flat +
                                         1:nc_minus_flat]
        beta[flat_cols] <- beta_block[(K+1)*nc_minus_flat +
                                        1:nc_flat]
        cbind(beta)
      })
      if(cor_fl){
        shur_correction <-
          shur_correction_function(
            list(Vhalf %**% collapse_block_diagonal(X)),
            list(cbind(Vhalfy)),
            list(cbind(unlist(result))),
            dispersion_temp,
            list(unlist(order_list)),#order_list,
            0,#K,
            family,
            unlist(observation_weights),
            ...
          )
      } else {
        shur_correction <-
          shur_correction_function(
            list(collapse_block_diagonal(X)),
            list(cbind(Vhalfy)),
            list(cbind(unlist(result))),
            dispersion_temp,
            list(unlist(order_list)),#order_list,
            0,#K,
            family,
            unlist(observation_weights),
            ...
          )
      }
      shur_correction <- lapply(1:(K + 1), function(k){
        shur_correction[[1]][(k-1)*nc + 1:nc, (k-1)*nc + 1:nc]
      })
      if(!(any(unlist(shur_correction) != 0))){
        shur_correction_collapsed <- 0
        shur_correction <- 0
      } else {
        ## Extract components from each partition's shur_correction
        shur_non_flat_list <- lapply(1:(K+1), function(k) {
          shur_correction[[k]][-flat_cols, -flat_cols]
        })

        ## Build block-diagonal for non-flat
        shur_correction1 <- Reduce("rbind", lapply(1:(K+1), function(k) {
          Reduce("cbind", lapply(1:(K+1), function(j) {
            if(j == k) shur_non_flat_list[[k]] else
              matrix(0, nrow=nrow(shur_non_flat_list[[k]]),
                     ncol=ncol(shur_non_flat_list[[k]]))
          }))
        }))

        ## Sum flat components
        shur_correction2 <- Reduce("+", lapply(1:(K+1), function(k) {
          shur_correction[[k]][flat_cols, flat_cols]
        }))
        if(K > 0) shur_correction2 <- shur_correction2 / (K+1)

        ## Build final matrix with same structure as Lambda_block
        shur_correction_collapsed <- matrix(0, nrow=ncol(X_block),
                                            ncol=ncol(X_block))
        shur_correction_collapsed[1:nrow(shur_correction1),
                                  1:ncol(shur_correction1)] <- shur_correction1
        shur_correction_collapsed[(nrow(shur_correction1)+1):ncol(X_block),
                                  (ncol(shur_correction1)+1):ncol(X_block)] <-
          shur_correction2
      }
      info <- t(X_block) %**% (c(W) * X_block) +
        Lambda_block +
        shur_correction_collapsed

      ## Numerical stability
      sc <- sqrt(mean(abs(info)))

      ## Initialize without inequality constraint, single Newton-Raphson step
      if(master_cnt == 1){
        infoinv_block <- sc*invert(sc*info)
        U <- (diag(1, nrow(info)) -
                infoinv_block %**%
                A %**%
                invert(t(A) %**% infoinv_block %**% A) %**%
                t(A))
        ## = UG[X^{T}(y-g(XB))]
        # Lambda_block %**% beta_block not needed since beta_block initializes
        # to 0
        if(cor_fl){
          beta_new <- c(
            U %**%
              infoinv_block %**%
              t(X_block) %**%
              (y_block - VhalfInv %**% cbind(family$linkinv(c(XB)))))
        } else {
          beta_new <- c(
            U %**%
              infoinv_block %**%
              t(X_block) %**%
              (y_block - family$linkinv(c(XB))))
        }
        infoinv_block <- NULL
      } else {
        ## Else QP update with inequality constraints
        if(cor_fl){
          mu <- VhalfInv %**% cbind(family$linkinv(XB))
        } else {
          mu <- cbind(family$linkinv(XB))
        }
        if(cor_fl){
          qp_score <- qp_score_function(
            X_block,
            y_block,
            mu,
            unlist(order_list),
            dispersion_temp,
            VhalfInv,
            unlist(observation_weights),
            ...
          )
        } else {
          qp_score <- qp_score_function(
            X_block,
            y_block,
            mu,
            unlist(order_list),
            dispersion_temp,
            NULL,
            unlist(observation_weights),
            ...
          )
        }
        beta_new <- try({quadprog::solve.QP(
          Dmat = info/sc,
          dvec = (qp_score -
                  Lambda_block %**% beta_block +
                          info %**% beta_block)/sc,
          Amat = qp_Amat,
          bvec = qp_bvec,
          meq = qp_meq
        )$solution}, silent = TRUE)
      }

      ## quadprog is very sensitive to positive definiteness
      if(any(inherits(beta_new, 'try-error'))){
        beta_new <- 0*beta_block
      }

      ## If no iteration, return solution as non-damped solution
      if((!iterate & master_cnt > 2) |
         (paste0(family)[1] == 'gaussian' &
          paste0(family)[2] == 'identity')){
        beta_block <- cbind(beta_new)
        damp_cnt <- 11
        master_cnt <- 101
        err <- tol - 1
        ## Else, use damped SQP
      } else {
        ## Damped update
        beta_new <- (1-damp)*beta_block + damp*beta_new
        if(cor_fl){
          XB <- Vhalf %**% (X_block %**% cbind(beta_new))
        } else {
          XB <- X_block %**% cbind(beta_new)
        }
        if(need_dispersion_for_estimation){
          dispersion_temp <- dispersion_function(
            family$linkinv(XB),
            Vhalfy,
            unlist(order_list),
            family,
            unlist(observation_weights),
            ...
          )
        } else {
          dispersion_temp <- 1
        }
        W <- c(glm_weight_function(family$linkinv(XB),
                                   Vhalfy,
                                   unlist(order_list),
                                   family,
                                   dispersion_temp,
                                   unlist(observation_weights),
                                   ...))

        ## Priority goes
        # 1) family$custom_dev.resids if available
        # 2) family$dev.resids if available
        # 3) MSE
        if(!is.null(family$custom_dev.resids)){
          raw <- family$custom_dev.resids(Vhalfy,
                                          family$linkinv(c(XB)),
                                          unlist(order_list),
                                          family,
                                          unlist(observation_weights),
                                          ...)
          if(cor_fl){
            err_new <- mean((
              VhalfInv %**% cbind(sign(raw)*sqrt(abs(
                raw
              ))) / sqrt(c(W))
            )^2)
          } else{
            err_new <- mean(raw)
          }
        } else if(is.null(family$dev.resids)){
          if(cor_fl){
            err_new <- mean(unlist(observation_weights)*
                          (y_block - VhalfInv %**% cbind(family$linkinv(XB))^2))
          } else {
            err_new <- mean(unlist(observation_weights)*
                              (y_block - cbind(family$linkinv(XB))^2))
          }
        } else {
          if(cor_fl){
            err_new <-
              mean(unlist(observation_weights)*
                     family$dev.resids(y_block,
                                       VhalfInv %**% family$linkinv(XB),
                                       wt = 1/W))
          } else {
            err_new <-
              mean(unlist(observation_weights)*
                     family$dev.resids(y_block,
                                       family$linkinv(XB),
                                       wt = 1/W))
          }
        }
        if(is.null(err_new) | is.na(err_new) | !is.finite(err_new)){
          damp_cnt <- damp_cnt + 1
        } else if(err_new <= err){

          ## Update coefficients, set damp to 0
          prev_err <- err
          err <- err_new
          abs_diff <- max(abs(beta_new - beta_block))
          beta_block <- beta_new
          damp_cnt <- 0

          ## Convergence criteria met when improvement in performance is
          # negligible and change in beta coefficients is small and
          # more than 10 iterations have passed
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

      ## Re-package into partition form
      nc_flat <- length(flat_cols)
      nc_minus_flat <- nc - nc_flat
      result <- lapply(1:(K+1),function(k){
        beta <- rep(0, nc)
        beta[-flat_cols] <- beta_block[(k-1)*nc_minus_flat +
                                           1:nc_minus_flat]
        beta[flat_cols] <- beta_block[(K+1)*nc_minus_flat +
                                      1:nc_flat]
        cbind(beta)
      })
    }
    ## Return output, with G re-computed if desired
    if(return_G_getB){
      if(paste0(family)[1] == 'gaussian' &
         paste0(family)[2] == 'identity'){
        return(list(
          B = result,
          G_list = list(G = lapply(Ghalf, function(mat) mat %**% mat),
                        Ghalf = Ghalf)
        ))
      }
      ## ## ## ##
      ## update
      if(cor_fl){
        XB <- Vhalf %**% (collapse_block_diagonal(X) %**% cbind(unlist(result)))
      } else {
        XB <- collapse_block_diagonal(X) %**% cbind(unlist(result))
      }
      if(need_dispersion_for_estimation){
        dispersion_temp <- dispersion_function(
          family$linkinv(XB),
          Vhalfy,
          unlist(order_list),
          family,
          unlist(observation_weights),
          ...
        )
      } else {
        dispersion_temp <- 1
      }
      W <- c(glm_weight_function(family$linkinv(XB),
                                 Vhalfy,
                                 unlist(order_list),
                                 family,
                                 dispersion_temp,
                                 unlist(observation_weights),
                                 ...))
      if(cor_fl){
        shur_correction <-
          shur_correction_function(
            list(Vhalf %**% collapse_block_diagonal(X)),
            list(cbind(Vhalfy)),
            list(cbind(unlist(result))),
            dispersion_temp,
            list(unlist(order_list)),#order_list,
            0,#K,
            family,
            unlist(observation_weights),
            ...
          )
      } else {
        shur_correction <-
          shur_correction_function(
            list(collapse_block_diagonal(X)),
            list(cbind(Vhalfy)),
            list(cbind(unlist(result))),
            dispersion_temp,
            list(unlist(order_list)),#order_list,
            0,#K,
            family,
            unlist(observation_weights),
            ...
          )
      }
      shur_correction <- lapply(1:(K + 1), function(k){
        shur_correction[[1]][(k-1)*nc + 1:nc, (k-1)*nc + 1:nc]
      })
      if(!(any(unlist(shur_correction) != 0))){
        shur_correction_collapsed <- 0
        shur_correction <- 0
      } else {
        ## Extract components from each partition's shur_correction
        shur_non_flat_list <- lapply(1:(K+1), function(k) {
          shur_correction[[k]][-flat_cols, -flat_cols]
        })

        ## Build block-diagonal for non-flat
        shur_correction1 <- Reduce("rbind", lapply(1:(K+1), function(k) {
          Reduce("cbind", lapply(1:(K+1), function(j) {
            if(j == k) shur_non_flat_list[[k]] else
              matrix(0, nrow=nrow(shur_non_flat_list[[k]]),
                     ncol=ncol(shur_non_flat_list[[k]]))
          }))
        }))

        ## Sum flat components
        shur_correction2 <- Reduce("+", lapply(1:(K+1), function(k) {
          shur_correction[[k]][flat_cols, flat_cols]
        }))
        if(K > 0) shur_correction2 <- shur_correction2 / (K+1)

        ## Build final matrix with same structure as Lambda_block
        shur_correction_collapsed <- matrix(0, nrow=ncol(X_block),
                                            ncol=ncol(X_block))
        shur_correction_collapsed[1:nrow(shur_correction1),
                                  1:ncol(shur_correction1)] <- shur_correction1
        shur_correction_collapsed[(nrow(shur_correction1)+1):ncol(X_block),
                                  (ncol(shur_correction1)+1):ncol(X_block)] <-
          shur_correction2
      }
      info <- t(X_block) %**% (c(W) * X_block) +
        Lambda_block +
        shur_correction_collapsed

      ## Replace the current G construction with this
      G <- lapply(1:(K+1), function(k) {NA})
      Ghalf <- lapply(1:(K+1), function(k) {NA})
      GhalfInv <- lapply(1:(K+1), function(k) {NA})

      ## Extract the sizes we need
      info_size <- dim(info)[1]
      non_flat_block_size <- (K+1) * nc_minus_flat
      flat_block_size <- nc_flat

      ## For each partition
      for(k in 1:(K+1)) {
        ## Initialize a matrix for this partition's G
        G[[k]] <- matrix(0, nrow=nc, ncol=nc)

        # Extract the non-flat portion for this partition
        non_flat_start <- (k-1) * nc_minus_flat + 1
        non_flat_end <- k * nc_minus_flat
        non_flat_indices <- non_flat_start:non_flat_end

        # Extract the flat portion (same for all partitions)
        flat_start <- non_flat_block_size + 1
        flat_end <- info_size
        flat_indices <- flat_start:flat_end

        # Copy the non-flat portion to G[[k]]
        G[[k]][-flat_cols, -flat_cols] <- info[non_flat_indices,
                                               non_flat_indices]

        # Copy the flat portion to G[[k]]
        G[[k]][flat_cols, flat_cols] <- info[flat_start:flat_end,
                                             flat_start:flat_end]

        # Copy the off-diagonal blocks (interactions between flat and non-flat)
        # = Upper-right: non-flat rows, flat columns
        G[[k]][-flat_cols, flat_cols] <- info[non_flat_indices,
                                              flat_indices]

        # = Lower-left: flat rows, non-flat columns (transpose of upper-right)
        G[[k]][flat_cols, -flat_cols] <- t(info[non_flat_indices,
                                                flat_indices])

        # Now compute Ghalf and GhalfInv as before
        eig <- eigen(G[[k]], symmetric = TRUE)
        vals <- eig$values
        vals[eig$values <= sqrt(.Machine$double.eps)] <- 0
        sqrt_vals <- sqrt(vals)
        Ghalf[[k]] <- eig$vectors %**% (t(eig$vectors) * sqrt_vals)
        inv_sqvals <- ifelse(sqrt_vals > 0, 1/sqrt_vals, 0)
        GhalfInv[[k]] <- eig$vectors %**% (t(eig$vectors) * inv_sqvals)
      }
      G_list <- list(
        Ghalf = Ghalf,
        GhalfInv = GhalfInv,
        G = G
      )
      return(list(
        B = result,
        G_list = G_list
      ))
    } else {
      return(result)
    }
  }

  ## For cases besides canonical Gaussian
  if(!(paste0(family)[2] == 'identity') |
     !(paste0(family)[1] == 'gaussian')){
    use_lm <- FALSE
    eig <- eigen(Lambda, symmetric = TRUE)
    LambdaHalf <- eig$vectors %**% (t(eig$vectors) * sqrt(ifelse(
      eig$values > 0,
      eig$values,
      0)))
    if(any(is.null(prevUnconB)) & parallel_unconstrained){
      ## Unconstrained estimate, extract from each partition
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
                    ## Adjust penalties if unique penalties are present for each partition
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

                    ## Fit ordinary model to each partition separately
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
    } else if(any(is.null(prevUnconB))){
      ## Starting unconstrained fits per partition
      unconstrained_estimate <- lapply(1:(K+1), function(k){

        ## Adjust penalties if unique penalties are present for each partition
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

        ## Fit ordinary model to each partition separately
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
    } else {
      unconstrained_estimate <- prevUnconB
    }

    ## No knots, we're done
    if(K == 0 & !quadprog & length(constraint_value_vectors) == 0){
      if(return_G_getB){

        ## ## ## ##
        ## update
        # dispersion
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
            mu,
            unlist(y),
            unlist(order_list),
            family,
            unlist(observation_weights),
            ...
          )
        } else {
          dispersion_temp <- 1
        }

        ## Weighted design matrix for convenience in computing G agnostic to
        # distribution
        Xw <- lapply(1:(K+1),
                     function(k){
                       var <- glm_weight_function(family$linkinv(X[[k]] %**%
                                        cbind(c(unconstrained_estimate[[k]]))),
                                                  y[[k]],
                                                  order_list[[k]],
                                                  family,
                                                  dispersion_temp,
                                                  observation_weights[[k]],
                                                  ...)
                       cbind(X[[k]] * c(sqrt(var)))
                     })
        ## X^{T}WX
        X_gram <- compute_gram_block_diagonal(Xw,
                                              parallel_matmult,
                                              cl,
                                              chunk_size,
                                              num_chunks,
                                              rem_chunks)

        ## Shur complements
        shur_corrections <- shur_correction_function(
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

        # Correlation matrix G, G^{1/2}, and G^{-1/2}
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
                                  shur_corrections)
        return(list(
          B = unconstrained_estimate,
          G_list = G_list
        ))
      } else {
        return(unconstrained_estimate)
      }
    }

    ## Else, iterate
    if(iter_count > 0){
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
          mu,
          unlist(y),
          unlist(order_list),
          family,
          unlist(observation_weights),
          ...
        )
      } else {
        dispersion_temp <- 1
      }
      shur_corrections <- shur_correction_function(
        X,
        y,
        prevB,
        dispersion_temp,
        order_list,
        K,
        family,
        observation_weights,
        ...
      )
      ## Update G^{-1/2} for left multiplying against unconstrained estimate
      GhalfInv <- compute_G_eigen(X_gram,
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
                                  shur_corrections)$GhalfInv
    }
    ## G^{1/2}X'y = G^{-1/2}GX'y
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

    ## Otherwise, with identity link, the computations greatly simplify
  } else {
    use_lm <- TRUE
    if(K == 0 & !quadprog & length(constraint_value_vectors) == 0){
      G <- list(Ghalf[[1]] %**% Ghalf[[1]])
      result <- list(G[[1]] %**% Xy[[1]])
      ## Return the coefficients directly
      if(return_G_getB){
        return(list(
          B = result,
          G_list = list(G = G,
                        Ghalf = Ghalf)
        ))
      } else {
        return(result)
      }
    } else {
      ## G^{1/2}X'y, the 'ystar' of the OLS problem
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

  ## \textbf{G}^{1/2}\textbf{A}, the 'xstar' of the ols problem
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
  )$residuals /  comp_stab_sc
  ## For imposing additional linear constraints homogenous across partitions
  if(length(constraint_value_vectors) > 0){
    if(any(unlist(constraint_value_vectors) != 0)){
      ## Compute the (I-U)b0 portion of Lagrngian solution,
      ## Ubhat + (I-U)b0
      comp_stab_sc <- 1/sqrt(K + 1)
      preds_star <- GhalfA %**%
        (invert(gramMatrix(GhalfA) * comp_stab_sc) %**%
           (t(A) %**%
              (Reduce("rbind", constraint_value_vectors) * comp_stab_sc)
           )
        )

      ## add these to the "residuals" of before
      resids_star <- resids_star + c(preds_star)
    }
  }

  ## Coefficients are the previous + update
  if(parallel_matmult & !is.null(cl)){
    # Handle remainder chunks first
    if(rem_chunks > 0) {
      rem_indices <- num_chunks * chunk_size + 1:rem_chunks
      rem <- lapply(rem_indices, function(k) {
        Ghalf[[k]] %**% cbind(resids_star[(k-1)*nc + 1:nc])
      })
    } else {
      rem <- list()
    }
    ## Process main chunks in parallel
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

  ## Sequential computation for U and B
  if(iterate & !use_lm & iter_count < 100){
    if(any(!is.null(prevB))){
      ## Change in previous estimate and current estimate of beta
      diff <- mean(
        unlist(
          sapply(1:(K+1),
                 function(k)mean(
                   abs(result[[k]] - prevB[[k]])))))
      if(diff < tol){
        iter_count <- Inf
      }
      ## if difference is getting bigger, stop - return previous
      if(c(prev_diff) <= c(diff)){
        result <- prevB
        iter_count <- Inf
      }
    } else {
      diff <- Inf
    }
    ## ## ## ##
    ## update
    # dispersion
    if(need_dispersion_for_estimation){
      dispersion_temp <- dispersion_function(
        family$linkinv(unlist(matmult_block_diagonal(X,
                                                     result,
                                                     K,
                                                     parallel_matmult,
                                                     cl,
                                                     chunk_size,
                                                     num_chunks,
                                                     rem_chunks))),
        unlist(y),
        unlist(order_list),
        family,
        unlist(observation_weights),
        ...
      )
    } else {
      dispersion_temp <- 1
    }
    # Weighted X for convenience
    Xw <- lapply(1:(K+1),
                 function(k){
                   if(nrow(X[[k]]) == 0){
                     return(X[[k]])
                   }
                   var <- glm_weight_function(family$linkinv(X[[k]] %**%
                                                cbind(c(result[[k]]))),
                                            y[[k]],
                                            order_list[[k]],
                                            family,
                                            dispersion_temp,
                                            observation_weights[[k]],
                                            ...)
                   if(length(var) == 1){
                     if(c(var) == 0){
                       return(t(t(X[[k]]) %**% rbind(0)))
                     } else {
                       return(t(t(X[[k]]) %**% rbind(sqrt(var))))
                     }
                   } else {
                     var <- diag(c(sqrt(var)))
                   }
                   t(t(X[[k]]) %**% var)
                 })
    # X^{T}WX
    X_gram <- compute_gram_block_diagonal(Xw,
                                          parallel_matmult,
                                          cl,
                                          chunk_size,
                                          num_chunks,
                                          rem_chunks)
    # Shur complements for accounting for uncertainty of estimating dispersion
    # or additional factors
    shur_corrections <- shur_correction_function(
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
                              keep_G = FALSE,
                              shur_corrections)

    ## Call recursively
    result <- get_B(X,
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
                    G_list$Ghalf,
                    G_list$GhalfInv,
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
                    prevB = result,
                    prevUnconB = unconstrained_estimate,
                    iter_count = iter_count + 1,
                    prev_diff = diff, # update arguments for each recursive call
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
                    blockfit,
                    just_linear_without_interactions,
                    Vhalf,
                    VhalfInv,
                    ...)
  } else {

    ## Impose quadratic programming constraints if desired
    # Not parallelized nor memory efficient
    if(quadprog){

      ## Big-matrix components (not memory efficient anymore)
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

      ## Solve smoothing constraints and qp constraints simultaneously
      damp_cnt <- 0
      master_cnt <- 0
      err <- Inf

      ## Initial link-transformed predictions
      XB <- X_block %**% beta_block

      ## Invoke procedure for fitting SQP problems
      while(err > tol & damp_cnt < 10 & master_cnt < 100){
        ## Initialize counters
        master_cnt <- master_cnt + 1
        damp <- 2^(-(damp_cnt))

        ## Information matrix in block-diagonal form
        if(need_dispersion_for_estimation){
          dispersion_temp <- dispersion_function(
            family$linkinv(XB),
            y_block,
            unlist(order_list),
            family,
            unlist(observation_weights),
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

        ## Re-package into partition form for shur correction
        result <- lapply(1:(K+1),function(k){
          cbind(beta_block[(k-1)*nc + 1:nc])
        })
        shur_correction <-
                    shur_correction_function(
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
        if(!(any(unlist(shur_correction) != 0))){
          shur_correction <- 0
        } else {
          shur_correction <- collapse_block_diagonal(shur_correction)
        }
        info <- (t(X_block) %**% (W * X_block)) + Lambda_block + shur_correction

        ## Numerical stability
        sc <- sqrt(mean(abs(info)))

        ## QP update
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

        ## quadprog is very sensitive to positive definiteness
        if(any(inherits(beta_new, 'try-error'))){
          beta_new <- 0*beta_block
        }

        ## If no iteration, return solution as non-damped solution
        if(!iterate & master_cnt > 1){
          beta_block <- beta_new
          damp_cnt <- 11
          master_cnt <- 101
          err <- tol - 1
          ## Else, use damped SQP
        } else {
          ## Damped update
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

            ## Update mean absolute value of score (err)
            # and coefficients, set damp to 0
            prev_err <- err
            err <- err_new
            abs_diff <- max(abs(beta_new - beta_block))
            beta_block <- beta_new
            damp_cnt <- 0

            ## Convergence criteria met when improvement in performance is
            # negligible and change in beta coefficients is small and
            # more than 10 iterations have passed
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

        ## Re-package into partition form
        result <- lapply(1:(K+1),function(k){
          cbind(beta_block[1:nc + (k-1)*nc])
        })
      }
    }
    if(return_G_getB){
      if(paste0(family)[1] == 'gaussian' &
         paste0(family)[2] == 'identity'){
          return(list(
            B = result,
            G_list = list(G = lapply(Ghalf, function(mat) mat %**% mat),
                          Ghalf = Ghalf)
          ))
      }
      ## ## ## ##
      ## update
      if(need_dispersion_for_estimation){
        dispersion_temp <- dispersion_function(
          family$linkinv(unlist(matmult_block_diagonal(X,
                                                       result,
                                                       K,
                                                       parallel_matmult,
                                                       cl,
                                                       chunk_size,
                                                       num_chunks,
                                                       rem_chunks))),
          unlist(y),
          unlist(order_list),
          family,
          unlist(observation_weights),
          ...
        )
      } else {
        dispersion_temp <- 1
      }
      ## Compute square root matrix of X^{T}WX => XW^{1/2}
      Xw <- lapply(1:(K+1),
                   function(k){
                     if(nrow(X[[k]]) == 0){
                       return(X[[k]])
                     }
                     var <- glm_weight_function(family$linkinv(X[[k]] %**%
                                                    cbind(c(result[[k]]))),
                                                y[[k]],
                                                order_list[[k]],
                                                family,
                                                dispersion_temp,
                                                observation_weights[[k]],
                                                ...)
                     if(length(var) == 1){
                       if(c(var) == 0){
                         return(t(t(X[[k]]) %**% rbind(0)))
                       } else {
                         return(t(t(X[[k]]) %**% rbind(sqrt(var))))
                       }
                     } else {
                       var <- diag(c(sqrt(var)))
                     }
                     t(t(X[[k]]) %**% var)
                   })
      X_gram <- compute_gram_block_diagonal(Xw,
                                            parallel_matmult,
                                            cl,
                                            chunk_size,
                                            num_chunks,
                                            rem_chunks)
      # Shur complements
      shur_corrections <- shur_correction_function(
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
                                shur_corrections)
      return(list(
        B = result,
        G_list = G_list
      ))
    } else {
      return(result)
    }
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
  I_minus_U <- t(matmult_U(A %**% (AGAInv %**% (-t(A))), G, nc, K))
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

#' Compute Eigenvalues and Related Matrices for G
#'
#' @param X_gram Gram matrix list (\eqn{\textbf{X}^{T}\textbf{X}})
#' @param Lambda Penalty matrix (\eqn{\boldsymbol{\Lambda}})
#' @param K Number of partitions minus 1 (\eqn{K})
#' @param parallel Use parallel processing
#' @param cl Cluster object
#' @param chunk_size Chunk size for parallel
#' @param num_chunks Number of chunks
#' @param rem_chunks Remaining chunks
#' @param family GLM family
#' @param unique_penalty_per_partition Use partition penalties
#' @param L_partition_list Partition penalty list (\eqn{\textbf{L}_\text{partition\_list}})
#' @param keep_G Return full G matrix (\eqn{\textbf{G}})
#' @param shur_corrections List of Shur complement corrections (\eqn{\textbf{S}})
#'
#' @details
#' Computes \eqn{\textbf{G}}, \eqn{\textbf{G}^{1/2}} and \eqn{\textbf{G}^{-1/2}} matrices via eigendecomposition of \eqn{\textbf{X}^{T}\textbf{X} + \boldsymbol{\Lambda}_\text{effective} + \textbf{S}}.
#' Handles partition-specific penalties and parallel processing.
#' For non-identity link functions, also returns \eqn{\textbf{G}^{-1/2}}.
#'
#' @return List containing combinations of:
#' \itemize{
#'   \item G - Full \eqn{\textbf{G}} matrix (if \code{keep_G=TRUE})
#'   \item Ghalf - \eqn{\textbf{G}^{1/2}} matrix
#'   \item GhalfInv - \eqn{\textbf{G}^{-1/2}} matrix (for non-identity links)
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
                            shur_corrections) {
  if(parallel & !is.null(cl)) {
    # Handle remainder chunks first
    if(rem_chunks > 0) {
      rem_indices <- num_chunks * chunk_size + 1: rem_chunks
      rem <- lapply(rem_indices, function(k) {
        ## Add partition-specific penalties if enabled
        if(unique_penalty_per_partition){
          eig <- tryCatch({
            eigen(X_gram[[k]] +
                    Lambda +
                    L_partition_list[[k]] +
                    shur_corrections[[k]],
                  symmetric = TRUE)
          }, error = function(e) NULL)
        } else {
          eig <- tryCatch({
            eigen(X_gram[[k]] +
                    Lambda +
                    shur_corrections[[k]],
                  symmetric = TRUE)
          }, error = function(e) NULL)
        }
        if(is.null(eig)) {
          return(list(G = NULL, Ghalf = NULL))
        }

        ## Handle numerical stability for eigenvalues
        eigen_values <- eig$values
        eigen_values[eig$values <= 0] <- 1
        inv_eigen_values <- 1/eigen_values
        inv_eigen_values[eig$values <= 0] <- 0
        sqrt_inv_eigen_values <- sqrt(inv_eigen_values)

        ## Compute matrix powers via eigendecomposition
        Ghalf <- eig$vectors %**% (t(eig$vectors) * sqrt_inv_eigen_values)
        if(keep_G){
          G <- eig$vectors %**% (t(eig$vectors) * inv_eigen_values)
        } else {
          G <- NULL
        }
        if((paste0(family)[2] != 'identity')){
          GhalfInv <- eig$vectors %**% (t(eig$vectors) / sqrt_inv_eigen_values)
          return(list(G = G,
                      Ghalf = Ghalf,
                      GhalfInv = GhalfInv))
        } else {
          return(list(G = G,
                      Ghalf = Ghalf))
        }
      })
    } else {
      rem <- list()
    }

    ## Process main chunks in parallel with same logic
    result <- c(
      Reduce("c",parallel::parLapply(cl, 1:num_chunks, function(chunk) {
        inds <- (chunk - 1)*chunk_size + 1:chunk_size
        lapply(inds, function(k) {
          if(unique_penalty_per_partition){
            eig <- tryCatch({
              eigen(X_gram[[k]] + Lambda +
                      L_partition_list[[k]] +
                      shur_corrections[[k]],
                    symmetric = TRUE)
            }, error = function(e) NULL)
          } else {
            eig <- tryCatch({
              eigen(X_gram[[k]] +
                      Lambda +
                      shur_corrections[[k]],
                    symmetric = TRUE)
            }, error = function(e) NULL)
          }
          if(is.null(eig)) {
            return(list(G = NULL, Ghalf = NULL))
          }
          eigen_values <- eig$values
          eigen_values[eig$values <= 0] <- 1
          inv_eigen_values <- 1/eigen_values
          inv_eigen_values[eig$values <= 0] <- 0
          sqrt_inv_eigen_values <- sqrt(inv_eigen_values)
          Ghalf <- eig$vectors %**% (t(eig$vectors) * sqrt_inv_eigen_values)
          if(keep_G){
            G <- eig$vectors %**% (t(eig$vectors) * inv_eigen_values)
          } else {
            G <- NULL
          }
          if((paste0(family)[2] != 'identity' |
              paste0(family)[1] != 'gaussian')){
            GhalfInv <- eig$vectors %**% (t(eig$vectors) /
                                            sqrt_inv_eigen_values)
            return(list(G = G,
                        Ghalf = Ghalf,
                        GhalfInv = GhalfInv))
          } else {
            return(list(G = G,
                        Ghalf = Ghalf))
          }
        })
      })),
      rem
    )
  } else {
    ## Sequential computation follows same logic as well
    result <- lapply(1:(K+1),function(k) {
      if(unique_penalty_per_partition){
        eig <- tryCatch({
          eigen(X_gram[[k]] +
                  Lambda +
                  L_partition_list[[k]] +
                  shur_corrections[[k]],
                symmetric = TRUE)
        }, error = function(e) NULL)
      } else {
        eig <- tryCatch({
          eigen(X_gram[[k]] +
                  Lambda +
                  shur_corrections[[k]],
                symmetric = TRUE)
        }, error = function(e) NULL)
      }
      if(is.null(eig)) {
        return(list(G = NULL, Ghalf = NULL))
      }
      eigen_values <- eig$values
      eigen_values[eig$values <= 0] <- 1
      inv_eigen_values <- 1/eigen_values
      inv_eigen_values[eig$values <= 0] <- 0
      sqrt_inv_eigen_values <- sqrt(inv_eigen_values)
      Ghalf <- eig$vectors %**% (t(eig$vectors) * sqrt_inv_eigen_values)
      if(keep_G){
        G <- eig$vectors %**% (t(eig$vectors) * inv_eigen_values)
      } else {
        G <- NULL
      }
      if((paste0(family)[2] != 'identity') |
         (paste0(family)[1] != 'gaussian')){
        GhalfInv <- eig$vectors %**% (t(eig$vectors) / sqrt_inv_eigen_values)
        return(list(G = G,
                    Ghalf = Ghalf,
                    GhalfInv = GhalfInv))
      } else {
        return(list(G = G,
                    Ghalf = Ghalf))
      }
    })
  }

  ## Reorganize results by matrix type (G, Ghalf, GhalfInv)
  if((paste0(family)[2] != 'identity') |
     (paste0(family)[1] != 'gaussian')){
    ## Case with GhalfInv
    result_processed <- list(
      G = lapply(result, `[[`, "G"),
      Ghalf = lapply(result, `[[`, "Ghalf"),
      GhalfInv = lapply(result, `[[`, "GhalfInv")
    )
  } else {
    ## Case without GhalfInv
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
#'                       shur_correction_function = weibull_shur_correction,
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
  z <- (log_y-log_mu)/scale
  logL <- status * (-log(scale) +
                      z -
                      log_y) -
    exp(z)

  ## Return sum of log likelihood
  return(sum(logL * weights))
}


#' Compute gradient of log-likelihood of Weibull accelerated failure model without penalization
#'
#' @description
#' Calculates the gradient of log-likelihood for a Weibull accelerated failure
#' time (AFT) survival model, supporting right-censored survival data.
#'
#' @param X Design matrix
#' @param y Response vector
#' @param mu Predicted mean vector
#' @param order_list List of observation indices per partition
#' @param dispersion Dispersion parameter (scale^2)
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
#'                       shur_correction_function = weibull_shur_correction,
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
  scale <- sqrt(dispersion)
  order_indices <- unlist(order_list)
  t(X) %**% cbind(exp((log(y) - log(mu))/scale) -
                    cbind(scale *
                            status[order_indices]) *
                    observation_weights)
}

#' Correction for the Variance-Covariance Matrix for Uncertainty in Scale
#'
#' @description
#' Computes the shur complement \eqn{\textbf{S}} such that \eqn{\textbf{G}^* = (\textbf{G}^{-1} + \textbf{S})^{-1}} properly
#' accounts for uncertainty in estimating dispersion when estimating
#' variance-covariance. Otherwise, the variance-covariance matrix is optimistic
#' and assumes the scale is known, when it was in fact estimated. Note that the
#' parameterization adds the output of this function elementwise (not subtract)
#' so for most cases, the output of this function will be negative or a
#' negative definite/semi-definite matrix.
#'
#' @param X Block-diagonal matrices of spline expansions
#' @param y Block-vector of response
#' @param B Block-vector of coefficient estimates
#' @param dispersion Scalar, estimate of dispersion, \eqn{ = \text{Weibull scale}^2}
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
#' List of \eqn{p \times p} matrices representing the shur-complement corrections \eqn{\textbf{S}_k} to be
#' elementwise added to each block of the information matrix, before inversion.
#'
#' @details
#' Adjusts the variance-covariance matrix unscaled for coefficients to account
#' for uncertainty in estimating the Weibull scale parameter, that otherwise
#' would be lost if simply using \eqn{\textbf{G}=(\textbf{X}^{T}\textbf{W}\textbf{X} + \textbf{L})^{-1}}. This is accomplished
#' using a correction based on the Shur complement so we avoid having to
#' construct the entire variance-covariance matrix, or modifying the procedure
#' for \code{\link{lgspline}} substantially.
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
#' ## Fit model using lgspline with Weibull shur correction
#' model_fit <- lgspline(y ~ spl(t1) + t2,
#'                       df,
#'                       unconstrained_fit_fxn = unconstrained_fit_weibull,
#'                       family = weibull_family(),
#'                       need_dispersion_for_estimation = TRUE,
#'                       dispersion_function = weibull_dispersion_function,
#'                       glm_weight_function = weibull_glm_weight_function,
#'                       shur_correction_function = weibull_shur_correction,
#'                       status = status,
#'                       opt = FALSE,
#'                       K = 1)
#'
#' print(summary(model_fit))
#'
#' ## Fit model using lgspline without Weibull shur correction
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
weibull_shur_correction <- function(X,
                                    y,
                                    B,
                                    dispersion,
                                    order_list,
                                    K,
                                    family,
                                    observation_weights,
                                    status){
  lapply(1:(K+1), function(k){
    if(nrow(X[[k]]) < 1){
      return(0)
    } else {
      mu <- family$linkinv(c(X[[k]] %**% B[[k]]))
      s <- status[order_list[[k]]]
      obs <- y[[k]]
      z <- (log(obs) - log(mu))/sqrt(dispersion)
      exp_z <- exp(z)
      zexp_z <- z*exp_z
      weights <- c(observation_weights[[k]])

      ## Correction via Shur complement
      # Extract true conditional variance-covariance of beta coefficients
      # conditional upon estimate of scale.
      # I = ( I_bb I_bs^{T} )
      #     ( I_bs I_ss     )
      # \eqn{\textbf{I} = \begin{pmatrix} \textbf{I}_{bb} & \textbf{I}_{bs}^{T} \\ \textbf{I}_{bs} & I_{ss} \end{pmatrix}}
      # for b = beta, s = dispersion (scale)
      # Note that: I_bb = invert(G[[k]]) is incorrect, I_bb is part of Fisher info
      I_bs <- t(X[[k]]) %**% cbind(weights * zexp_z * sqrt(dispersion))
      I_ss <- -sum(
        weights * (
          (s + 2*s*z + zexp_z + exp_z * z^2)
        )
      )
      # compl gets elementwise added to G[[k]] for all k = 1...K+1
      compl <- I_bs %**% matrix(-1/I_ss) %**% t(I_bs)
      # Shur complement correction to pass on to compute_G_eigen()
      return(compl)
    }
  })
}

#' Estimate Scale for Weibull Accelerated Failure Time Model
#'
#' @description
#' Computes maximum log-likelihood scale estimate of Weibull accelerated failure
#' time (AFT) survival model.
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
#' Scalar representing the estimated scale
#'
#' @details
#' Calculates maximum log-likelihood estimate of scale for Weibull AFT model
#' accounting for right-censored observations using Brent's method for
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
      -loglik_weibull(log_y, log_mu, status, par, weights)
    },
    method = 'Brent',
    lower = 1e-64,
    upper = 100
  )$par
}

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
#' Provides a comprehensive family specification for Weibull AFT models, including Family
#' name, link function, inverse link function, and custom loss function for model tuning
#'
#' Supports right-censored survival data with flexible parameter estimation.
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
#'                       shur_correction_function = weibull_shur_correction,
#'                       status = status,
#'                       opt = FALSE,
#'                       K = 1)
#'
#' summary(model_fit)
#'
#' @export
weibull_family <- function()list(family = "weibull",
                                 link = "log",
                                 linkfun = log,
                                 linkinv = exp,
                                 ## Custom loss used in place of MSE for computing GCV
                                 custom_dev.resids =
                                   function(y,
                                            mu,
                                            order_indices,
                                            family,
                                            observation_weights,
                                            status){
                                     log_mu <- log(mu)
                                     log_y <- log(y)
                                     status <- status[order_indices]

                                     ## Initialize scale
                                     init_scale <-
                                       weibull_scale(log_y,
                                                     mean(log_y),
                                                     status[order_indices],
                                                     observation_weights)
                                     ## Find scale
                                     scale <- optim(
                                       init_scale,
                                       fn = function(par){
                                         -loglik_weibull(log_y,
                                                         log_mu,
                                                         status,
                                                         par,
                                                         observation_weights)
                                       },
                                       lower = init_scale/5,
                                       upper = init_scale*5,
                                       method = 'Brent'
                                     )$par

                                     ## -2 * log-likelihood
                                     dev <- -2*(
                                       ## Log-likelihood contributions
                                       status * (-log(scale) +
                                                   (1/scale - 1)*log_y -
                                                   log_mu/scale) -
                                         (exp((log_y - log_mu)/scale))
                                     )
                                     return(dev * observation_weights)
                                   })


#' Estimate Weibull Dispersion for Accelerated Failure Time Model
#'
#' @description
#' Computes the scale parameter for a Weibull accelerated failure time (AFT)
#' model, supporting right-censored survival data.
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
#' Squared scale estimate for the Weibull AFT model (dispersion)
#'
#' @seealso \code{\link{weibull_scale}} for the underlying scale estimation function
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
#' ## Estimate dispersion
#' dispersion_est <- weibull_dispersion_function(
#'   mu = mu,
#'   y = y,
#'   order_indices = order_indices,
#'   family = weibull_family(),
#'   observation_weights = weights,
#'   status = status
#' )
#'
#' print(dispersion_est)
#'
#' @export
weibull_dispersion_function <- function(mu,
                                        y,
                                        order_indices,
                                        family,
                                        observation_weights,
                                        status){

  ## Maximizes log-likelihood of right-censored data
  log_mu <- log(mu)
  log_y <- log(y)
  observation_weights <- c(observation_weights)
  status <- status[order_indices]

  ## Initialize scale
  init_scale <-
    weibull_scale(log_y,
                  mean(log_y),
                  status[order_indices],
                  observation_weights)
  ## Find scale
  scale <- optim(
    init_scale,
    fn = function(par){
      -loglik_weibull(log_y,
                      log_mu,
                      status,
                      par,
                      observation_weights)
    },
    lower = init_scale/5,
    upper = init_scale*5,
    method = 'Brent'
  )$par

  return(scale^2)
}

#' Weibull GLM Weight Function for Constructing Information Matrix
#'
#' @description
#' Computes diagonal weight matrix \eqn{\textbf{W}} for the information matrix
#' \eqn{\textbf{G} = (\textbf{X}^{T}\textbf{W}\textbf{X} + \textbf{L})^{-1}} in Weibull accelerated failure time (AFT) models.
#'
#' @param mu Predicted survival times
#' @param y Observed response/survival times
#' @param order_indices Order of observations when partitioned to match "status" to "response"
#' @param family Weibull AFT family
#' @param dispersion Estimated dispersion parameter (\eqn{s^2})
#' @param observation_weights Weights of observations submitted to function
#' @param status Censoring indicator (1 = event, 0 = censored)
#'   Indicates whether an event of interest occurred (1) or the observation was
#'   right-censored (0). In survival analysis, right-censoring occurs when the
#'   full survival time is unknown, typically because the study ended or the
#'   subject was lost to follow-up before the event of interest occurred.
#'
#' @return
#' Vector of weights for constructing the diagonal weight matrix \eqn{\textbf{W}}
#' in the information matrix \eqn{\textbf{G} = (\textbf{X}^{T}\textbf{W}\textbf{X} + \textbf{L})^{-1}}.
#'
#' @details
#' This function generates weights used in constructing the information matrix
#' after unconstrained estimates have been found. Specifically, it is used in
#' the construction of the \eqn{\textbf{U}} and \eqn{\textbf{G}} matrices following initial unconstrained
#' parameter estimation.
#'
#' These weights are analogous to the variance terms in generalized linear
#' models (GLMs). Like logistic regression uses \eqn{\mu(1-\mu)}, Poisson regression uses
#' \eqn{e^{\mu}}, and Linear regression uses constant weights, Weibull AFT models use
#' \eqn{\exp((\log y - \log \mu)/s)} where \eqn{s} is the scale (= \eqn{\sqrt{\text{dispersion}}}) parameter.
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
#'                       shur_correction_function = weibull_shur_correction,
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
  val <- exp((log(y) - log(mu))/sqrt(dispersion))
  if(any(!is.finite(val))){
    return(rep(1, length(val)))
  }
  newval <- val * c(observation_weights)
  return(newval)
}

#' Compute Newton-Raphson Parameter Update with Numerical Stabilization
#'
#' @description
#' Performs parameter update in iterative optimization.
#'
#' Called by \code{\link{damped_newton_r}} in the update step
#'
#' @param gradient_val Numeric vector of gradient values (\eqn{\textbf{u}})
#' @param neghessian_val Negative Hessian matrix (\eqn{\textbf{G}^{-1}} approximately)
#'
#' @return
#' Numeric vector of parameter updates (\eqn{\textbf{G}\textbf{u}})
#'
#' @details
#' This helper function is a core component of Newton-Raphson optimization.
#' It provides a computationally-stable approach to computing \eqn{\textbf{G}\textbf{u}}, for
#' information matrix \eqn{\textbf{G}} and score vector \eqn{\textbf{u}}, where the Newton-Raphson update
#' can be expressed as \eqn{\boldsymbol{\beta}^{(m+1)} = \boldsymbol{\beta}^{(m)} + \textbf{G}\textbf{u}}.
#'
#' @seealso \code{\link{damped_newton_r}} for the full optimization routine
#'
#' @keywords internal
#' @export
nr_iterate <- function(gradient_val, neghessian_val){
  sc <- sqrt(mean(abs(neghessian_val))) # for computational stability
  invert(neghessian_val / sc) %**% cbind(gradient_val / sc)
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
      cat('\n \t Error Encountered, Number of N.R. steps so far: ', master_count, '\n')
      stop('\n \t NA/NaN/non-finite value detected when running unconstrained damped',
           ' Newton-Raphson.',
           ' \n \t Try re-fitting a simpler model, using greater/smaller penalties, ',
           ' experimenting with different knot locations, or reducing the',
           ' number of knots.')
    }

    ## Damp iterations, only updates if performance improves
    while((new_objective <= prev_objective) & count < max_dmp_steps){
      new_param <- old_param + (2^(-count))*nr_iterate(gradient(old_param),
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
    # Break if Nas/NaNs/Infs occur
    if(any(is.na(eps) | is.nan(eps) | !is.finite(eps))){
      new_param <- old_param
      eps <- 0
    }
    master_count <- master_count + 1
  }
  return(new_param)
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
#' @param LambdaHalf Square root of penalty matrix (\eqn{\boldsymbol{\Lambda}^{1/2}})
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
#' Optimized beta parameter estimates (\eqn{\boldsymbol{\beta}}) for Weibull AFT model
#'
#' @details
#' Estimation Approach:
#' The function employs a two-stage optimization strategy for fitting
#' accelerated failure time models via maximum likelihood:
#'
#' 1. Outer Loop: Estimate Scale Parameter using Brent's method
#'
#' 2. Inner Loop: Estimate Regression Coefficients (given scale) using
#'    damped Newton-Raphson.
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
#'                       shur_correction_function = weibull_shur_correction,
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
                                      status # status goes in the ellipse arg
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

  ## Initialize scale
  init_scale <- weibull_scale(log_y,
                              mean(log_y),
                              status[order_indices],
                              weights)

  ## First, use outer-loop to optimize scale
  # Then given scale, optimize beta
  scale <- optim(init_scale,
                 fn = function(par){
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
                         -0.5*c(t(beta) %**% Lambda %**% beta)
                     },
                     function(par){
                       beta <- cbind(par)
                       eta <- c(X %**% beta)
                       z <- (log_y - eta)/scale
                       zeta <- exp(z)
                       grad_beta <- t(X) %**%
                         (weights*(zeta - status[order_indices]))*scale +
                         Lambda %**% beta
                       cbind(grad_beta)
                     },
                     function(par){
                       beta <- cbind(par)
                       eta <- X %**% beta
                       z <- (log_y - eta)/scale
                       zeta <- c(exp(z))
                       info <- (t(X) %**% (weights * zeta * X) + Lambda)
                       info
                     },
                     tol
                   ))
                   log_mu <- X %**% beta
                   -loglik_weibull(log_y,
                                   log_mu,
                                   status[order_indices],
                                   par,
                                   weights) +
                     0.5*c(t(beta) %**% Lambda %**% beta)
                 },
                 lower = init_scale/5,
                 upper = init_scale*5,
                 method = 'Brent')$par

  ## Now optimize beta, given optimal scale
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
        -0.5*c(t(beta) %**% Lambda %**% beta)
    },
    function(par){
      beta <- cbind(par)
      eta <- c(X %**% beta)
      z <- (log_y - eta)/scale
      zeta <- exp(z)
      grad_beta <- t(X) %**%
        (weights*(zeta - status[order_indices]))*scale +
        Lambda %**% beta
      cbind(grad_beta)
    },
    function(par){
      beta <- cbind(par)
      eta <- X %**% beta
      z <- (log_y - eta)/scale
      zeta <- c(exp(z))
      info <- (t(X) %**% (weights * zeta * X) + Lambda)
      info
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
#' @param LambdaHalf Square root of penalty matrix (\eqn{\boldsymbol{\Lambda}^{1/2}})
#' @param Lambda Penalty matrix (\eqn{\boldsymbol{\Lambda}})
#' @param keep_weighted_Lambda Logical flag to control penalty matrix handling:
#'   - `TRUE`: Return coefficients directly from weighted penalty fitting
#'   - `FALSE`: Apply damped Newton-Raphson optimization to refine estimates
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
#' For fitting non-canonical GLMs, use \code{keep_weighted_Lambda = TRUE} since the
#' score and hessian equations below are no longer valid.
#'
#' For Gamma(link='log') using \code{keep_weighted_Lambda = TRUE} is misleading.
#' The information is weighted by a constant (shape parameter) rather than some
#' mean-variance relationship. So \code{keep_weighted_Lambda = TRUE} is highly
#' recommended for log-link Gamma models. This constant flushes into the
#' penalty terms, and so the formulation of the information matrix is valid.
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
#' as long as 'glm_weight_function' and 'qp_score_function' are properly specified.
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
                      y = cbind(c(y, rep(family$linkinv(0), nrow(LambdaHalf)))),
                      family = family,
                      weights = c(weights, rep(mean(weights), nrow(LambdaHalf))),
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
  tX <- t(X)
  res <- cbind(damped_newton_r(
    ## initial guess
    init,
    ## proportional to log-likelihood
    function(par){
      -sum(weights*family$dev.resids(
        y,
        family$linkinv(c(X %**% cbind(par))),
        wt = 1))*0.5 -
        0.5*c(t(par) %**% Lambda %**% cbind(par))
    },
    ## score
    function(par){
      c(tX %**% (weights*cbind(y - family$linkinv(X %**% cbind(par)))) -
          Lambda %**% cbind(par))
    },
    ## information
    function(par){
      tX %**% (weights*c(family$variance(X %**% cbind(par))) * X) +
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
  ## Basic control parameters shared from stats::optim
  ctrl <- list(
    maxit = 50,
    abstol = sqrt(.Machine$double.eps),
    reltol = sqrt(.Machine$double.eps),
    initial_damp = 1,
    min_damp = 2^-16,
    trace = FALSE
  )
  ctrl[names(control)] <- control

  ## Setup
  n_params <- length(par)
  x <- par
  Inv <- diag(n_params)
  damp <- ctrl$initial_damp
  best_x <- x
  best_f <- Inf
  best_Inv <- Inv

  ## Objective and gradient both evaluated by fn is key difference here
  result <- fn(c(x))
  if(length(result) != 2) stop("fn must return list of (objective, gradient)")
  f <- result[[1]]
  grad <- result[[2]]
  ## Use finite-difference automatically if grad is NULL or NA
  if(is.null(grad) |
     any(is.na(grad))) grad <- approx_grad(x, fn)
  if(length(grad) != n_params) stop("gradient must match parameter length")

  ## Iterate through, running damped BFGS
  for(iter in 1:ctrl$maxit) {
    prev_x <- x
    prev_f <- f
    prev_grad <- grad

    p <- -Inv %**% cbind(grad)
    x_new <- x + damp * p

    result <- fn(c(x_new))
    f_new <- result[[1]]
    grad_new <- result[[2]]
    ## Use finite-difference automatically if grad is NULL or NA
    if(is.null(grad_new) |
       any(is.na(grad_new))) grad_new <- approx_grad(x_new, fn)

    ## Error checking
    if(is.na(f_new) || is.nan(f_new) || !is.finite(f_new)) {
      damp <- damp/2
      if(damp < ctrl$min_damp) break
      next
    }

    ## If new objective is better or less than 3 iterations have occurred
    if(f_new < f || iter <= 2) {
      if(f_new < best_f) {
        best_f <- f_new
        best_x <- x_new
        best_Inv <- Inv
      }

      ## Classic BFGS - notation is NOT the same as general lgspline
      # i.e. y is NOT our response here
      s0 <- cbind(x_new - x)
      y0 <- cbind(grad_new - grad)
      denom <- sum(y0 * s0)

      if(abs(denom) > 1e-16) {
        rho0 <- 1/denom
        term1 <- diag(n_params) - rho0 * (s0 %**% t(y0))
        term2 <- diag(n_params) - rho0 * (y0 %**% t(s0))
        Inv <- term1 %**% Inv %**% term2 + rho0 * (s0 %**% t(s0))
      }

      x <- x_new
      f <- f_new
      grad <- grad_new
      damp <- 1

      if(iter > 2) {
        if(abs(f - prev_f) < ctrl$abstol * (abs(prev_f) + ctrl$reltol)) break
        if(max(abs(x - prev_x)) < ctrl$abstol) break
      }

      ## Otherwise, damp-step and try again
    } else {
      x <- best_x*(1-damp) + x*damp
      f <- best_f*(1-damp) + f*damp
      Inv <- best_Inv*(1-damp) + Inv*damp
      damp <- damp/2
      if(damp < ctrl$min_damp) break
    }

    ## Optional printout
    if(ctrl$trace) {
      cat(sprintf("Iter %d: f = %f, |grad| = %f, damp = %f\n",
                  iter, f, sqrt(sum(grad^2)), damp))
    }
  }

  ## Return output analogoous to stats::optim
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

#' Calculate Matrix Square Root
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
  eig <- eigen(mat)
  sqrtv <- suppressWarnings(suppressMessages(sqrt(eig$values)))
  eig$vectors %**% (t(eig$vectors) * sqrtv)
}

#' Calculate Matrix Inverse Square Root
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
  eig <- eigen(mat)
  sqrtv <- suppressWarnings(suppressMessages(sqrt(eig$values)))
  eig$vectors %**% (t(eig$vectors) / sqrtv)
}

