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
#' model.matrix() with dummy-intercept coding. Column names are cleaned by removing the
#' 'x' prefix added by model.matrix().
#'
#' @examples
#' x <- factor(c("A", "B", "A", "C"))
#' create_onehot(x)
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
#' Standardizes a numeric vector by centering and scaling to unit variance.
#'
#' @param x Numeric vector to standardize
#'
#' @return Standardized vector with mean 0 and standard deviation 1
#'
#' @examples
#' x <- c(1, 2, 3, 4, 5)
#' std(x)
#'
#' @noRd
std <- function(x){
  (x-mean(x))/sd(x)
}

#' Efficient Matrix Multiplication Operator
#'
#' @description
#' Operator wrapper around C++ efficient_matrix_mult() for matrix multiplication syntax.
#'
#' This is an internal helper function - use at your own risk.
#'
#' @param x Left matrix
#' @param y Right matrix
#'
#' @return Matrix product of x and y
#'
#' @examples
#' A <- matrix(1:4, 2, 2)
#' B <- matrix(5:8, 2, 2)
#' A %**% B
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
#'
#' @return Inverted matrix or identity matrix if all methods fail
#'
#' @details
#' Tries methods in order:
#' 1. Direct inversion using armaInv()
#' 2. Generalized inverse using eigendecomposition
#' 3. Returns identity matrix with warning if both fail
#'
#' For eigendecomposition, uses a small ridge penalty (1e-16) for stability and
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
#' @export
invert <- function(mat){

  ## Try inversion
  t <- try({
    armaInv(mat)
  },silent = TRUE)

  ## Try generalized inverse with small ridge penalty
  if(any(class(t) == 'try-error')){
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
  if(any(class(t) == 'try-error')){
    warning('Matrix not inverted, returning identity: ', print(t))
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
#' @param K Number of blocks minus 1
#' @param parallel Logical; whether to use parallel processing
#' @param cl Cluster object for parallel processing
#' @param chunk_size Number of blocks per chunk for parallel processing
#' @param num_chunks Number of chunks for parallel processing
#' @param rem_chunks Remaining blocks after chunking
#'
#' @return List containing products of corresponding blocks
#'
#' @details
#' When parallel=TRUE, splits computation into chunks processed in parallel.
#' Handles remainder chunks separately. Uses matmult_block_diagonal_cpp() for
#' actual multiplication.
#'
#' The function expects A and B to contain corresponding blocks that can be
#' matrix multiplied.
#'
#' @examples
#' A <- list(matrix(1:4,2,2), matrix(5:8,2,2))
#' B <- list(matrix(1:4,2,2), matrix(5:8,2,2))
#' matmult_block_diagonal(A, B, K=1, parallel=FALSE, cl=NULL,
#'                       chunk_size=1, num_chunks=1, rem_chunks=0)
#'
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
      Reduce("c",parLapply(cl, 1:num_chunks, function(k){
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
#' Creates a design matrix containing polynomial expansions and interaction terms
#' for predictor variables. Supports customizable term generation including
#' polynomial degrees, interaction types, and selective term exclusion.
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
#' @param ... Additional arguments passed to custom_basis_fxn
#'
#' @return Matrix with columns for intercept, polynomial terms, and specified interactions
#'
#' @noRd
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
                                      include_quadratic_interactions = TRUE,
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
       length(intersect(numerics, vars_that_interact)) > 1){
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
#' @details
#' * For first derivatives: linear terms → 1, quadratic → 2x, cubic → 3x², quartic → 4x³
#' * For second derivatives: linear terms → 0, quadratic → 2, cubic → 6x, quartic → 12x²
#' * Handles division by zero by adding small constant to denominator
#'
#' @noRd
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
#' Computes second derivatives for interaction terms including two-way linear,
#' quadratic, and three-way interactions. Handles special cases for each type.
#'
#' @param dat Numeric matrix; design matrix containing basis expansions
#' @param var Character; variable name to differentiate with respect to
#' @param interaction_single_cols Integer vector; column indices for linear-linear interactions
#' @param interaction_quad_cols Integer vector; column indices for linear-quadratic interactions
#' @param triplet_cols Integer vector; column indices for three-way interactions
#' @param colnm_C Character vector; column names of complete basis matrix
#' @param power1_cols Integer vector; column indices of linear terms
#' @param scale Numeric; scaling factor for normalization
#'
#' @details
#' Calculates second derivatives for:
#' * Linear-linear interactions (constant)
#' * Linear-quadratic interactions (linear in other variable)
#' * Three-way interactions (sum of other variables)
#'
#' @return Numeric matrix of second derivatives, same dimensions as input
#'
#' @noRd
take_interaction_2ndderivative <-
  function(dat,
           var,
           interaction_single_cols,
           interaction_quad_cols,
           triplet_cols,
           colnm_C,
           power1_cols,
           scale) {

    ## Initialize output matrix
    n_cols <- ncol(dat)
    dat_deriv <- matrix(0, nrow = nrow(dat), ncol = n_cols)
    colnames(dat_deriv) <- colnames(dat)
    variable_values <- dat[,var]

    ## Index of linear term
    v <- which(colnm_C[power1_cols] == var)

    ## Detect interactions, if relevant
    if (length(interaction_single_cols) > 0) {
      interaction_singles <-
        interaction_single_cols[grep(paste0("_", var, "_"),
                                     colnm_C[interaction_single_cols])]
      if (length(interaction_singles) > 0) {
        ## 2nd derivative of two-way linear-linear interactions is always 1
        dat_deriv[, interaction_singles] <- 1
      }
    }

    ## For two-way linear-quadratic interactions, more work is needed
    if (length(interaction_quad_cols) > 0) {
      interaction_quads <-
        interaction_quad_cols[grep(paste0("_", var, "_"),
                                   colnm_C[interaction_quad_cols])]
      if (length(interaction_quads) > 0) {
        for (w in 1:length(power1_cols[-v])) {
          ## The other variable, with interactions affecting quadratic terms
          wvar <- c(power1_cols[-v])[w]

          ## Quadratic interaction indices
          interq <-
            interaction_quads[grep(colnm_C[wvar], colnm_C[interaction_quads])]
          if (length(interq) > 0) {
            ## This is the _var^2_x_w term
            if(length(power2_cols) > 0){
              nchv <- nchar(colnm_C[power2_cols[v]])
              interqv2 <-
                interq[substr(colnm_C[interq], nchar(colnm_C[interq]) - nchv + 1,
                        nchar(colnm_C[interq])) == colnm_C[power2_cols[v]]]

              ## This is the _var_x_w_^2 term
              if(length(interqv2) > 0){
                interqv1 <- interq[-which(interq == interqv2)]
              }
            } else {
              ## This is the _var_x_w_^2 term
              interqv1 <- interq
            }

            ## 2nd derivative of each is 2*other, 2*other + 2*self respectively
            dat_deriv[, interqv1] <- 2 * dat[, colnm_C[wvar]]
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
        triplet_cols[grep(paste0("_", var, "_"), colnm_C[triplet_cols])]
      if (length(triplets) > 0) {
        other2_vars <- lapply(triplets, function(tr) {
          vars <- unlist(strsplit(colnm_C[tr], 'x'))
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
#' Divides a matrix into K+1 partitions based on knot locations, returning a list
#' of submatrices where each contains rows falling between adjacent knot boundaries.
#'
#' @param partition_codes Numeric vector; values determining partition assignment for each row
#' @param partition_bounds Numeric vector; ordered knot locations defining partition boundaries
#' @param nr Integer; number of rows in input matrix
#' @param mat Numeric matrix; data to be partitioned
#' @param K Integer; number of interior knots (resulting in K+1 partitions)
#'
#' @details
#' * For K=0, returns original matrix in a single-element list
#' * For K>0, creates K+1 partitions using (-Inf, knot_1, ..., knot_K, Inf) as boundaries
#' * Each row is assigned to partition i if partition_codes[row] falls between bounds[i] and bounds[i+1]
#' * Empty partitions return 0-row matrices
#'
#' @return List of length K+1, each element containing the submatrix for that partition
#'
#' @noRd
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
#' Creates block diagonal matrix by placing input matrices diagonally
#' with zeros elsewhere. Matrices must have compatible dimensions.
#'
#' @noRd
create_block_diagonal <- function(matrix_list) {

  # Calculate the total dimensions of the resulting matrix
  total_dim <- sum(sapply(matrix_list, nrow))

  # Create an empty matrix filled with zeros
  result <- matrix(0, nrow = total_dim, ncol = total_dim)

  # Fill the diagonal blocks
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
#' @return List of Gram matrices (X^{T}X) for each block
#'
#' @noRd
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
      Reduce("c", parLapply(cl, 1:num_chunks, function(chunk) {
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
#' @param power1_cols Indices of linear terms
#' @param interaction_single_cols Indices of first-order interactions
#' @param interaction_quad_cols Indices of quadratic interactions
#' @param triplet_cols Indices of three-way interactions
#' @param K Number of partitions minus 1
#' @param include_2way_interactions Include 2-way interactions
#' @param include_3way_interactions Include 3-way interactions
#' @param include_quadratic_interactions Include quadratic interactions
#' @param colnm_C Column names
#' @param C_scales Scale factors
#' @param just_first_derivatives Only compute first derivatives
#'
#' @return List containing first and second derivative matrices
#'
#' @noRd
make_derivative_matrix  <-  function(
    nc,
    Cpredictors,
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
    just_first_derivatives = FALSE){


  ## First derivative, for all numeric variables
  first_derivs <- lapply(colnm_C[power1_cols], function(v){
    take_derivative(dat = Cpredictors, var = v, scale = C_scales[v])
  })
  names(first_derivs) <- colnm_C[power1_cols]
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
      colnm_C[power1_cols], function(v) {
        take_interaction_2ndderivative(
          dat = Cpredictors,
          var = v,
          interaction_single_cols,
          interaction_quad_cols,
          triplet_cols,
          colnm_C,
          power1_cols,
          scale = C_scales[v]
        ) +
          take_derivative(dat = Cpredictors,
                          var = v,
                          second = TRUE,
                          scale = C_scales[v])
      }
    )
    names(second_derivs) <- colnm_C[power1_cols]
    ## No interactions present
  } else {
    second_derivs <- lapply(colnm_C[power1_cols], function(v){
      take_derivative(dat = Cpredictors,
                      var = v,
                      second = TRUE,
                      scale = C_scales[v])
    })
    names(second_derivs) <- colnm_C[power1_cols]
  }

  return(list('first_derivatives' = first_derivs,
              'second_derivatives' = second_derivs))
}

#' Create Smoothing Spline Constraint Matrix
#'
#' @description
#' Constructs constraint matrix A enforcing continuity and smoothness at knot boundaries
#' by constraining function values, derivatives, and interactions between partitions.
#'
#' @param nc Integer; number of columns in basis expansion
#' @param CKnots Matrix; basis expansions evaluated at knot points
#' @param power1_cols Integer vector; indices of linear terms
#' @param nonspline_cols Integer vector; indices of non-spline terms
#' @param interaction_single_cols Integer vector; indices of linear interaction terms
#' @param interaction_quad_cols Integer vector; indices of quadratic interaction terms
#' @param triplet_cols Integer vector; indices of three-way interaction terms
#' @param K Integer; number of interior knots (K+1 partitions)
#' @param include_constrain_fitted Logical; constrain function values at knots
#' @param include_constrain_first_deriv Logical; constrain first derivatives at knots
#' @param include_constrain_second_deriv Logical; constrain second derivatives at knots
#' @param include_constrain_interactions Logical; constrain interaction terms at knots
#' @param include_2way_interactions Logical; include two-way interactions
#' @param include_3way_interactions Logical; include three-way interactions
#' @param include_quadratic_interactions Logical; include quadratic interactions
#' @param colnm_C Character vector; column names for basis expansions
#' @param C_scales Numeric vector; scaling factors for standardization
#'
#' @details
#' Constraint matrix structure:
#' * Rows: Each constraint at knot locations
#' * Columns: Coefficients for each partition
#' * Values: +1/-1 pattern enforcing equality between adjacent partitions
#'
#' Constraints enforced (if enabled):
#' * Function value continuity at knots
#' * First derivative continuity at knots
#' * Second derivative continuity at knots
#' * Interaction term continuity at knots
#'
#' @return Matrix A of constraint coefficients. Rows correspond to constraints,
#' columns to coefficients across all partitions.
#'
#' @noRd
make_constraint_matrix <- function(nc,
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
                                   C_scales){

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
  # (note, knots are assigned to the "larger" partition b.)
  # If non-spline cols, set to 0. They shouldn't affect constr. on spline efx.
  CKnots_spline <- CKnots
  if(length(nonspline_cols) > 0 & length(power1_cols) > 0){
    CKnots_spline[,nonspline_cols] <-
      0*CKnots_spline[,nonspline_cols]
  } else if(length(nonspline_cols) > 0){
    ## If no spline effects present
    constrain_fitted <- CKnots[,rep(1:nc, K+1), drop = FALSE] *
                        checkered_fitted_expand
  }
  if(length(power1_cols) > 0){
    ## If spline effects present
    constrain_fitted = CKnots_spline[,rep(1:nc, K+1), drop = FALSE] *
                       checkered_fitted_expand
  }

  ## Zero-out the fitted constraint if not desired
  if(!include_constrain_fitted){
    constrain_fitted <- 0 * constrain_fitted
  }

  ## First derivative, for all numeric variables
  if(include_constrain_first_deriv){
    first_derivs <- lapply(colnames(CKnots)[power1_cols], function(v){
      take_derivative(dat = CKnots_spline, var = v, scale = C_scales[v])
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
                      scale = C_scales[v])
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
                            take_derivative(dat = CKnots_spline,
                                            var = v,
                                            second = TRUE,
                                            scale = C_scales[v])
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
                                      scale = C_scales[v],
                                      second = TRUE)
                          } else {
                            take_interaction_2ndderivative(
                              dat = CKnots_spline,
                              var = v,
                              interaction_single_cols,
                              interaction_quad_cols,
                              triplet_cols,
                              colnm_C,
                              power1_cols,
                              scale = C_scales[v]
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
#' @noRd
is_binary <- function(x){
  if(length(unique(x)) > 2){
    return(FALSE)
  }
  TRUE
}


#' Compute Derivative of Penalty Matrix G with Respect to Lambda
#'
#' @description
#' Calculates the derivative of the penalty matrix G with respect to the
#' smoothing parameter lambda, supporting both global and partition-specific penalties.
#' (i.e. the derivative of diagonal weight matrix 1/(1+x^{T}UGx) w.r.t. penalty)
#'
#' @param G A list of penalty matrices for each partition
#' @param L The base penalty matrix
#' @param K Number of partitions minus 1
#' @param lambda Smoothing parameter value
#' @param unique_penalty_per_partition Logical indicating partition-specific penalties
#' @param L_partition_list Optional list of partition-specific penalty matrices
#' @param parallel Logical to enable parallel processing
#' @param cl Cluster object for parallel computation
#' @param chunk_size Size of chunks for parallel processing
#' @param num_chunks Number of chunks
#' @param rem_chunks Remainder chunks
#'
#' @return
#' A list of derivative matrices for each partition
#'
#' @keywords internal
#' @noRd
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
      Reduce("c",parLapply(cl, 1:num_chunks, function(chunk) {
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

#' Compute Derivative of Penalty Matrix G with Respect to Lambda
#'
#' @description
#' Calculates the derivative of the penalty matrix G with respect to the
#' smoothing parameter lambda, supporting both global and partition-specific penalties.
#' (i.e. the derivative of diagonal weight matrix 1/(1+x'UGx) w.r.t. penalty)
#'
#' @param G A list of penalty matrices for each partition
#' @param L The base penalty matrix
#' @param K Number of partitions minus 1
#' @param lambda Smoothing parameter value
#' @param unique_penalty_per_partition Logical indicating partition-specific penalties
#' @param L_partition_list Optional list of partition-specific penalty matrices
#' @param parallel Logical to enable parallel processing
#' @param cl Cluster object for parallel computation
#' @param chunk_size Size of chunks for parallel processing
#' @param num_chunks Number of chunks
#' @param rem_chunks Remainder chunks
#'
#' @return
#' A list of derivative matrices for each partition
#'
#' @keywords internal
#' @noRd
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
    results <- parLapply(cl, chunk_indices, function(inds) {
      start_idx <- min(inds) - 1

      # Compute initial trace
      trace_part <- sum(sapply(inds, function(k)
        sum(diag(dG_dlambda[[k]] %**% GXX[[k]]))))

      # Get relevant chunk of A
      A_chunk <- A[(start_idx*nc + 1):min((max(inds))*nc, nrow(A)), ]
      len_inds <- length(inds) - 1

      # Compute both corrections at once
      correction1 <- compute_trace_correction(
        G[inds], A_chunk, GXX[inds], dGhalf_dlambda[inds],
        AGAInv, nc, len_inds)

      correction2 <- compute_trace_correction(
        G[inds], A_chunk, dGXX[inds], Ghalf[inds],
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


#' Calculate Trace of Matrix Product XUGX^{T}
#'
#' @param G List of G matrices
#' @param A Constraint matrix
#' @param GXX List of GX'X products
#' @param Ghalf List of G^(1/2) matrices
#' @param AGAInv Inverse of A^{T}GA
#' @param nc Number of columns
#' @param K Number of partitions minus 1
#' @param parallel Use parallel processing
#' @param cl Cluster object
#' @param chunk_size Size of parallel chunks
#' @param num_chunks Number of chunks
#' @param rem_chunks Remaining chunks
#'
#' @details
#' Computes trace(XUGX^{T}) = trace(GX^{T}X - PGX^{T}X) where P = GA(A^{T}GA)^(-1)A^{T}
#' Handles parallel computation by splitting into chunks.
#'
# = trace of UGX^{T}X = trace of (GX^{T}X - PGX^{T}X)
# for P = GA(A^{T}GA)^{-1}A^{T}
# trace of PGX^{T}X is the "correction" as described
#'
#' @return Trace value
#'
#' @noRd
compute_trace_UGXX_wrapper <- function(G,
                                       A,
                                       GXX,
                                       Ghalf,
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

    ## Compute traces in parallel with single parLapply call
    comp_stab_const <- sqrt(mean(abs(unlist(G))))
    results <- parLapply(cl, chunk_indices, function(inds) {
      start_idx <- min(inds) - 1

      ## Compute first trace
      trace_part <- sum(sapply(GXX[inds], function(x) sum(diag(x))))

      ## Compute correction
      A_chunk <- A[(start_idx*nc + 1):min((max(inds))*nc, nrow(A)), ]
      correction_part <- compute_trace_correction(G[inds],
                                                  A_chunk /
                                                    sqrt(comp_stab_const),
                                                  GXX[inds],
                                                  Ghalf[inds],
                                                  AGAInv,
                                                  nc,
                                                  length(inds)-1)

      c(trace_part, correction_part * comp_stab_const)
    })

    ## Efficiently sum results
    all_results <- do.call(rbind, results)
    trace <- mean(all_results[,1])
    correction <- mean(all_results[,2])

    return((trace - correction) * nrow(all_results))

  } else {

    ## Compute first trace
    trace_part <- sum(unlist(lapply(GXX, function(x) sum(diag(x)))))

    ## Compute correction
    comp_stab_const <- sqrt(mean(abs(unlist(G))))
    correction_part <- compute_trace_correction(G,
                                                A / sqrt(comp_stab_const),
                                                GXX,
                                                Ghalf,
                                                AGAInv,
                                                nc,
                                                K)

    return(c(trace_part - correction_part * comp_stab_const))
  }
}


#' Compute Derivative of UGX^{T}y with Respect to Lambda
#'
#' @param AGAInv_AGXy Product of (A^{T}GA)^{-1} and (A^{T}GX^{T}y)
#' @param AGAInv Inverse of A^{T}GA
#' @param G List of G matrices
#' @param A Constraint matrix
#' @param dG_dlambda List of dG/dλ matrices
#' @param nc Number of columns
#' @param nca Number of constraint columns
#' @param K Number of partitions minus 1
#' @param Xy List of X^{T}y products
#' @param Ghalf List of G^{1/2} matrices
#' @param dGhalf List of dG^{1/2}/dλ matrices
#' @param GhalfXy_temp Temporary storage for G^{1/2}X^{T}y
#' @param parallel Use parallel processing
#' @param cl Cluster object
#' @param chunk_size Size of parallel chunks
#' @param num_chunks Number of chunks
#' @param rem_chunks Remaining chunks
#'
#' @details
#' Computes d(UGX^{T}y)/dλ = d[(I-GA(A^{T}GA)^{-1}A^{T})GX^{T}y]/dλ
#' Uses efficient implementation avoiding large matrix construction.
#' For large problems (K >= 10, nc > 4), uses chunked parallel processing.
#' For smaller problems, uses simpler least squares approach.
#'
#' @return Vector of derivatives
#'
#' @noRd
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
      chunk_results <- parLapply(cl, 1:num_chunks, function(chunk) {
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
      chunk_results <- parLapply(cl, 1:num_chunks, function(chunk) {
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

    ## Transform to least squares problem using d[G^(1/2)]/dN;
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
    comp_stab_sc <- 1/sqrt(K + 1)
    resids1 <- do.call('.lm.fit', list(x = GhalfA * comp_stab_sc,
                                       y = GhalfXy_temp * comp_stab_sc)
    )$residuals / comp_stab_sc

    ## Get second residuals
    resids2 <- do.call('.lm.fit', list(x = dGhalfA * comp_stab_sc,
                                       y = dGhalfXy * comp_stab_sc)
      )$residuals / comp_stab_sc

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
        Reduce("c",parLapply(cl, 1:num_chunks, function(chunk) {
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


#' Compute Component G^{1/2}A(A^{T}dGA)^{-1}A^{T}GX^{T}y
#'
#' @param G List of G matrices
#' @param Ghalf List of G^{1/2} matrices
#' @param A Constraint matrix
#' @param AGAInv Inverse of A^{T}GA
#' @param Xy List of X^{T}y products
#' @param nc Number of columns
#' @param K Number of partitions minus 1
#' @param parallel Use parallel processing
#' @param cl Cluster object
#' @param chunk_size Size of parallel chunks
#' @param num_chunks Number of chunks
#' @param rem_chunks Remaining chunks
#'
#' @details
#' Computes G^{1/2}A(A^{T}dGA)^{-1}A^{T}GX^{T}y efficiently in parallel chunks.
#' Returns both the result and intermediate AGAInvAGXy product for reuse.
#'
#' @return List containing:
#' \itemize{
#'   \item Result vector
#'   \item AGAInvAGXy intermediate product
#' }
#'
#' @noRd
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
    AGXy_chunks <- parLapply(cl, 1:num_chunks, function(chunk) {
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
    result_chunks <- parLapply(cl, 1:num_chunks, function(chunk) {
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
#' @param use_custom_bfgs Logical; TRUE for analytic gradient BFGS, FALSE for finite differences
#' @param family GLM family with optional custom tuning loss
#' @param wiggle_penalty,flat_ridge_penalty Initial penalty values
#' @param log_initial_wiggle,log_initial_flat Initial grid search values (log scale)
#' @param unique_penalty_per_predictor,unique_penalty_per_partition Logical; allow predictor/partition-specific penalties
#' @param log_penalty_vec Initial values for predictor/partition penalties (log scale)
#' @param penalty_ridge Regularization for predictor/partition penalties
#' @param keep_weighted_Lambda,iterate Logical controlling GLM fitting
#' @param quadprog,qp_Amat,qp_bvec,qp_meq Quadratic programming parameters
#' @param tol Numeric; convergence tolerance
#' @param sd_y,delta Response standardization parameters
#' @param constraint_value_vectors List; constraint values
#' @param parallel,parallel_* Logical; enable parallel computation
#' @param cl,chunk_size,num_chunks,rem_chunks Parallel computation parameters
#' @param custom_penalty_mat Optional custom penalty matrix
#' @param order_list List; observation ordering by partition
#' @param glm_weight_function,shur_correction_function Functions for GLM weights and corrections
#' @param need_dispersion_for_estimation,dispersion_function Control dispersion estimation
#' @param observation_weights Optional observation weights
#' @param homogenous_weights Logical; TRUE if all weights equal
#' @param verbose Logical; print progress
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
#' @noRd
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
    iterate,
    quadprog,  qp_Amat, qp_bvec, qp_meq,
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
    verbose,
    ...
){

  if(verbose){
    cat('    Starting tuning\n')
  }

  ## For compatibility without multiple predictors nor multiple partitions
  if(length(log_penalty_vec) == 0){
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
  ## psuedocount
  if(is.null(delta) &
     any(is.null(family$custom_tuning_loss))){
    ## ignore psuedocount for y if not needed for link-fxn transform
    if((paste0(family)[2] == 'identity') |
       (paste0(family)[2] == 'log' & (min(unl_y) > 0 )) |
       ((any(paste0(family)[2] %in% c('inverse','1/mu^2'))) &
        (!(any(unl_y == 0)))) |
       (paste0(family)[2] == 'logit' & (!any(unl_y %in% c(0,1))))
    ){
      delta <- 0
      ## obtain psuedocount for y if needed for link-fxn transform
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
    wiggle_penalty <- unlist(exp(par[1])) # smoothing spline f''(x)^2
    flat_ridge_penalty <- unlist(exp(par[2])) # flat ridge regression penalty
    if(unique_penalty_per_predictor | unique_penalty_per_partition){
      penalty_vec <- exp(c(par[-c(1:2)]))
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
                                  colnm_C,
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
                                     rem_chunks))

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
      cat('        compute_trace_UGXX_wrapper\n')
    }
    sum_W <- compute_trace_UGXX_wrapper(G_list$G,
                                        A,
                                        GXX,
                                        G_list$Ghalf,
                                        AGAInv,
                                        nc,
                                        K,
                                        parallel & parallel_trace,
                                        cl,
                                        chunk_size,
                                        num_chunks,
                                        rem_chunks)

    if(verbose)cat('        Get predictions gcv_u fxn\n')
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
       is.null(family$custom_tuning_loss)){

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
        family$custom_tuning_loss(y[[k]],
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

    ## Regularization penalty
    if(verbose){
      cat('        gcvu_fxn penalization operations\n')
    }
    if(unique_penalty_per_partition | unique_penalty_per_predictor){
      regulizer <- 0.5*penalty_ridge*sum(penalty_vec^2) +
        -0.5*1e-32*((wiggle_penalty - 1))^2
    } else {
      regulizer <- -0.5*1e-32*((wiggle_penalty - 1))^2
    }

    if(verbose)cat('        done GCVu,', GCV_u, '\n')
    ## Output list, prevent from computing twice
    return(list(GCV_u = GCV_u + regulizer,
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
    wiggle_penalty <- unlist(exp(par[1])) # smoothing spline f''(x)^2
    flat_ridge_penalty <- unlist(exp(par[2])) # flat ridge regression penalty
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
                                    colnm_C,
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
                                       rem_chunks))


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
                                          G_list$Ghalf,
                                          AGAInv,
                                          nc,
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
         is.null(family$custom_tuning_loss)){

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
          family$custom_tuning_loss(y[[k]],
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

      ## Regularization penalty - pull wiggle/partition/predictor penalties to 1
      if(unique_penalty_per_partition | unique_penalty_per_predictor){
        regulizer <- 0.5*penalty_ridge*sum(penalty_vec^2) +
          -0.5*1e-32*((wiggle_penalty - 1))^2
      } else {
        regulizer <- -0.5*1e-32*((wiggle_penalty - 1))^2
      }

      if(verbose){
        cat('        gr fxn outlist\n')
      }
      outlist <- list(GCV_u = GCV_u + regulizer,
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
    gradient <- cbind(c(dGCV_u_dlambda1 * lambda_1,
                        dGCV_u_dlambda2 * lambda_2))

    ## Gradients for unique penalties per predictor
    if(unique_penalty_per_predictor){
      predictor_penalties <-
        penalty_vec[grep('predictor',names(penalty_vec))]

      predictor_penalty_gradient <-
        sapply(1:length(predictor_penalties),function(j){
        mean(diag(outlist$L_predictor_list[[j]]))/
        mean(diag(outlist$Lambda +
                  outlist$L_predictor_list[[j]])) * # ratio of trace measures contribution
        dGCV_u_dlambda1 *
        predictor_penalties[j] # chain rule, the derivative being passed on
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
          mean(diag(outlist$Lambda +
                    outlist$L_partition_list[[j]])) *
          dGCV_u_dlambda1 *
          partition_penalties[j] # chain rule
        })

      ## Combined gradient
      gradient <- cbind(c(gradient, partition_penalty_gradient))
    }

    ## Regularization penalty
    if(unique_penalty_per_partition | unique_penalty_per_predictor){
      regulizer <- c(-1e-32*(wiggle_penalty - 1),
                     0,
                     penalty_ridge*penalty_vec)

    } else {
      regulizer <- c(-1e-32*(wiggle_penalty - 1),
                     0)
    }
    if(verbose){
      cat('        Gradient end \n')
    }

    ## Return output
    return(list(GCV_u = outlist$GCV_u +
                  0.5*penalty_ridge*sum(penalty_vec^2) +
                  0.5*1e-32*((wiggle_penalty - 1)^2),
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

        ## Reset to best solution if numerical issues
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

        ## Reset to best solution if numerical issues
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
        result <- gcvu_fxn(c(unlist(par), log_penalty_vec))$GCV_u
        if(is.na(result) | is.nan(result)) {
          return(Inf)
        }
        return(result)
      }, error = function(e) {
        print(e)
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
      stop('All GCV criteria for the initial tuning grid were computed as NA, NaN, or non-finite: check your data for corrupt or missing values, try changing initial tuning grid, or try manual tuning instead. If you are setting no_intercept = TRUE, try experimenting with standardize_response = FALSE/TRUE.')
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
    ## My BFGS (closed-form gradient)
    if(use_custom_bfgs){
      res <- try(quasi_nr_fxn(c(best_start, log_penalty_vec)), silent = TRUE)
      if(any(class(res) == 'try-error')){
        print(res)
        warning('Custom BFGS implementation failed. Try use_custom_bfgs = FALSE, or manual tuning. Resorting to best as selected from grid search.')
        par <- c(best_start, log_penalty_vec)
      } else {
        par <- res$par
      }
    } else {
      ## vs. base R (finite-difference approx.)
      res <- try({optim(c(best_start, log_penalty_vec),
                   fn = function(par){
                     gcvu_fxn(par)$GCV_u
                    },
                   method = 'BFGS'
                  )}, silent = TRUE)
      if(any(class(res) == 'try-error')){
        print(res)
        warning('Base R BFGS failed. Try use_custom_bfgs = TRUE, or manual tuning. Resorting to best as selected from grid search.')
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
    infl <- ((nr+2)/((nr-2)))^2
    wiggle_penalty <- exp(par[1])*infl
    flat_ridge_penalty <- exp(par[2])*infl
    ## Update penalty vec for predictor-and-partition specific penalties
    if(length(log_penalty_vec) > 0){
      penalty_vec <- exp(log_penalty_vec)
      penalty_vec[1:length(penalty_vec)] <-
        exp(c(par[-c(1,2)]))*infl
    }
  } else if(length(log_penalty_vec) > 0){
    penalty_vec <- exp(log_penalty_vec)
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
                                colnm_C,
                                just_Lambda = FALSE)

  return(list("Lambda" = Lambda_list$Lambda,
              "flat_ridge_penalty" = flat_ridge_penalty,
              "wiggle_penalty" = wiggle_penalty,
              "other_penalties" = penalty_vec,
              "L_predictor_list" = Lambda_list$L_predictor_list,
              "L_partition_list" = Lambda_list$L_partition_list))
}

#' Matrix Multiplication for Block Diagonal Structure
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
#' Processes in parallel chunks if enabled.
#'
#' @return Matrix product
#'
#' @noRd
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
      Reduce("c",parLapply(cl, 1:num_chunks, function(chunk) {
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
#' Computes coefficient estimates under constraints and penalties for GLMs, handling
#' both regular and custom families, identifiability conditions, and quadratic
#' programming constraints. Uses efficient implementations avoiding explicit U matrix
#' construction.
#'
#' @param X List of design matrices by partition
#' @param X_gram List of Gram matrices by partition
#' @param Lambda Combined penalty matrix
#' @param keep_weighted_Lambda Logical; retain GLM weights in smoothing spline penalty as an artifact of Tikhinov parameterization
#' @param unique_penalty_per_partition Logical; allow partition-specific penalties
#' @param L_partition_list List of partition-specific penalty matrices
#' @param A Matrix of smoothness constraints
#' @param Xy List of X^{T}y products by partition
#' @param y List of responses by partition
#' @param K Integer; number of interior knots
#' @param nc Integer; columns per partition
#' @param nca Integer; constraint matrix columns
#' @param Ghalf,GhalfInv Lists of G^{1/2} and G^{-1/2} matrices
#' @param parallel_* Logical flags for parallel computation components
#' @param cl,chunk_size,num_chunks,rem_chunks Parallel processing parameters
#' @param family GLM family object with optional custom tuning loss
#' @param unconstrained_fit_fxn Function for unconstrained estimates
#' @param iterate Logical; use iterative optimization
#' @param quadprog,qp_Amat,qp_bvec,qp_meq Quadratic programming parameters
#' @param prevB,prevUnconB,iter_count,prev_diff Optimization history
#' @param tol Numeric; convergence tolerance
#' @param constraint_value_vectors List of additional constraint values
#' @param order_list List of observation orderings
#' @param glm_weight_function Function for GLM weights
#' @param shur_correction_function Function for uncertainty corrections (added elementwise to G^{-1})
#' @param need_dispersion_for_estimation,dispersion_function Control dispersion estimation
#' @param observation_weights Optional observation weights
#' @param homogenous_weights Logical; TRUE if weights are constant
#' @param return_G_getB Logical; return G matrices with coefficients
#' @param ... Additional arguments passed to fitting functions
#'
#' @details
#' Key features:
#' * Handles regular and custom GLM families
#' * Efficient implementation avoiding explicit U matrix
#' * Supports multiple types of constraints:
#'   - Smoothness at knots
#'   - Linear equality constraints
#'   - Quadratic programming constraints
#' * Iterative fitting for non-canonical links
#' * Parallel processing options for large problems
#' * Adjusts for uncertainty in dispersion estimation if needed
#' * Supports non-constant observation weights
#'
#' Implementation notes:
#' * Handles numerical stability via scaling
#' * Special case optimizations for Gaussian/identity
#' * Efficient parallel matrix operations (optional)
#' * Memory-efficient constraint handling
#'
#' @return List containing:
#' \itemize{
#'   \item B - List of coefficient vectors by partition
#'   \item G_list - Optional list of G matrices if return_G_getB=TRUE:
#'     \itemize{
#'       \item G - Correlation matrices
#'       \item Ghalf - Square root matrices
#'     }
#' }
#'
#' @noRd
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
                  ...){

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
        parLapply(cl,
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
                                            if(unique_penalty_per_partition){
                                              Lambda_temp <- Lambda +
                                                L_partition_list[[k]]
                                              eig <- eigen(Lambda_temp,
                                                           symmetric = TRUE)
                                              LambdaHalf_temp <-
                                                eig$vectors %**%
                                                (t(eig$vectors) * (sqrt(
                                                  ifelse(eig$values <= 0,0,
                                                         eig$values))))
                                            } else {
                                              LambdaHalf_temp <- LambdaHalf
                                              Lambda_temp <- Lambda
                                            }

                                            ## Fit ordinary model to each
                                            # partition seperately
                                            cbind(c(unconstrained_fit_fxn(
                                              X[[k]],
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

  ## G^{1/2}A, the 'xstar' of the ols problem
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
      Reduce("c",parLapply(cl, 1:num_chunks, function(chunk) {
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
                    ...)
  } else {

    ## Impose quadratic programming constraints if desired
    # Not parallelized nor memory efficient
    if(quadprog){

      ## Big-matrix components (not memory efficient anymore)
      X_block <- Reduce("rbind", lapply(1:(K+1),function(k){
        dummy <- 0*X[[k]]
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

      ## Solve smoothing constraints and monotonic constraints simultaneously
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
        info <- (t(X_block) %**% (X_block * W)) + Lambda_block + shur_correction

        ## Numerical stability
        sc <- sqrt(mean(abs(info)))

        ## QP update
        beta_new <- try({quadprog::solve.QP(
          Dmat = info/sc,
          dvec = (t(X_block) %**% (y_block - family$linkinv(XB)) -
                    Lambda_block %**% beta_block +
                    info %**% beta_block)/sc,
          Amat = cbind(A, qp_Amat),
          bvec = c(rep(0, ncol(A)), qp_bvec),
          meq = ncol(A) + qp_meq
        )$solution}, silent = TRUE)

        ## quadprog is very sensitive to positive definiteness
        if(any(class(beta_new) == 'try-error')){
          beta_new <- 0*beta_block
        }

        ## If no iteration, return solution as non-damped solution
        if(!iterate){
          damp_cnt <- 11
          master_cnt <- 101
          err <- tol - 1
          ## Else, use damped SQP
        } else {
          ## Damped update
          beta_new <- (1-damp)*beta_block + damp*beta_new
          XB <- X_block %**% beta_new
          err_new <- mean(abs(y_block - family$linkinv(XB)))
          if(err_new <= err){

            ## Update mean absolute value of score (err)
            # and coefficients, set damp to 0
            err <- err_new
            beta_block <- beta_new
            damp_cnt <- 0

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


#' Matrix Multiplication for A^{T}GA
#'
#' @param G List of G matrices
#' @param A Constraint matrix
#' @param K Number of partitions minus 1
#' @param nc Number of columns per partition
#' @param nca Number of constraint columns
#' @param parallel Use parallel processing
#' @param cl Cluster object
#' @param chunk_size Chunk size for parallel
#' @param num_chunks Number of chunks
#' @param rem_chunks Remaining chunks
#'
#' @details
#' Computes A^{T}GA efficiently in parallel chunks using AGAmult_chunk().
#'
#' @return Matrix product A^{T}GA
#'
#' @noRd
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
    chunk_results <- parLapply(cl, 1:num_chunks, function(chunk) {
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

#' Construct U Matrix
#'
#' @param G List of G matrices
#' @param A Constraint matrix
#' @param K Number of partitions minus 1
#' @param nc Number of columns per partition
#' @param nca Number of constraint columns
#'
#' @return U matrix for constraints
#'
#' @details
#' Computes U = I - GA(A^{T}GA)^{-1}A^{T} efficiently
#'
#' @noRd
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
#' Returns selected combinations from the cartesian product of vec_list
#' without constructing full expand.grid()
#'
#' @return Data frame of selected combinations
#'
#' @noRd
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
#' @param L1 Matrix; integrated squared second derivative penalty
#' @param wiggle_penalty,flat_ridge_penalty Numeric; smoothing and ridge penalty parameters
#' @param K Integer; number of interior knots
#' @param nc Integer; number of basis columns per partition
#' @param unique_penalty_per_predictor,unique_penalty_per_partition Logical; enable predictor/partition-specific penalties
#' @param penalty_vec Named numeric; custom penalty values for predictors/partitions
#' @param colnm_C Character; column names for linking penalties to predictors
#' @param just_Lambda Logical; return only combined penalty matrix
#'
#' @details
#' Penalty matrix structure:
#' * Base penalty = wiggle_penalty * (L1 + L2)
#' * L1 = integrated squared second derivative penalty
#' * L2 = ridge penalty (diagonal or custom)
#' * Predictor penalties scale L1 elements involving specific predictors
#' * Partition penalties create unique L1 scaling per partition
#'
#' Implementation features:
#' * Memory efficient for large problems
#' * Handles missing/zero penalties gracefully
#' * Maintains sparsity where possible
#' * Links penalties to predictors via column names
#'
#' @return List containing:
#' \itemize{
#'   \item Lambda - Combined nc x nc penalty matrix
#'   \item L1 - Smoothing spline penalty matrix
#'   \item L2 - Ridge penalty matrix
#'   \item L_predictor_list - List of predictor-specific penalty matrices
#'   \item L_partition_list - List of partition-specific penalty matrices
#' }
#'
#' If just_Lambda=TRUE and no partition penalties, returns only Lambda matrix.
#'
#' @noRd
compute_Lambda <- function(custom_penalty_mat,
                           L1,
                           wiggle_penalty,
                           flat_ridge_penalty,
                           K,
                           nc,
                           unique_penalty_per_predictor,
                           unique_penalty_per_partition,
                           penalty_vec,
                           colnm_C,
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
      inds <- grep(substr(predictors[j], 10, nchar(predictors[j])), colnm_C)
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
#' @param X_gram Gram matrix list
#' @param Lambda Penalty matrix
#' @param K Number of partitions minus 1
#' @param parallel Use parallel processing
#' @param cl Cluster object
#' @param chunk_size Chunk size for parallel
#' @param num_chunks Number of chunks
#' @param rem_chunks Remaining chunks
#' @param family GLM family
#' @param unique_penalty_per_partition Use partition penalties
#' @param L_partition_list Partition penalty list
#' @param keep_G Return full G matrix
#'
#' @details
#' Computes G, G^{1/2}, and G^{-1/2} matrices via eigendecomposition.
#' Handles partition-specific penalties and parallel processing.
#' For non-identity link functions, also returns G^{-1/2}.
#'
#' @return List containing combinations of:
#' \itemize{
#'   \item G - Full G matrix (if keep_G=TRUE)
#'   \item Ghalf - G^{1/2} matrix
#'   \item GhalfInv - G^{-1/2} matrix (for non-identity links)
#' }
#'
#' @noRd
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
          eig <- eigen(X_gram[[k]] +
                         Lambda +
                         L_partition_list[[k]] +
                         shur_corrections[[k]],
                       symmetric = TRUE)
        } else {
          eig <- eigen(X_gram[[k]] +
                       Lambda  +
                       shur_corrections[[k]],
                    symmetric = TRUE)
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
      Reduce("c",parLapply(cl, 1:num_chunks, function(chunk) {
        inds <- (chunk - 1)*chunk_size + 1:chunk_size
        lapply(inds, function(k) {
          if(unique_penalty_per_partition){
            eig <- eigen(X_gram[[k]] +
                           Lambda +
                           L_partition_list[[k]] +
                           shur_corrections[[k]],
                        symmetric = TRUE)
          } else {
            eig <- eigen(X_gram[[k]] +
                         Lambda   +
                         shur_corrections[[k]],
                       symmetric = TRUE)
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
            GhalfInv <- eig$vectors %**% (t(eig$vectors) / sqrt_inv_eigen_values)
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
        eig <- eigen(X_gram[[k]] +
                       Lambda +
                       L_partition_list[[k]] +
                       shur_corrections[[k]],
                     symmetric = TRUE)
      } else {
        eig <- eigen(X_gram[[k]] +
                       Lambda +
                       shur_corrections[[k]],
                     symmetric = TRUE)
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
#' Calculates dG^{1/2}/dλ matrices for each partition using eigendecomposition.
#' Follows similar approach to compute_G_eigen() but for matrix derivatives.
#'
#' @param dG_dlambda List of nc x nc dG/dλ matrices by partition
#' @param nc Integer; number of columns per partition
#' @param K Integer; number of interior knots
#' @param parallel,cl,chunk_size,num_chunks,rem_chunks Parallel computation parameters
#'
#' @details
#' Key steps:
#' * Handles NA/zero values in input matrices
#' * Uses eigendecomposition to compute matrix square root derivative
#' * Ensures numerical stability via eigenvalue thresholding
#' * Supports parallel processing with chunking
#'
#' Mathematical details:
#' * For each partition k:
#'   - Compute eigendecomposition of dG_k/dλ
#'   - Floor eigenvalues at 1, zero out negatives
#'   - Compute dG_k^{1/2}/dλ for all k = 1...(K+1) partitions
#'
#' @return List of nc x nc matrices containing dG_k^{1/2}/dλ for each partition k
#'
#' @noRd
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
      Reduce("c",parLapply(cl, 1:num_chunks, function(chunk) {
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
#' @param K Number of partitions minus 1
#' @param cluster_args List with custom centers and kmeans args
#' @param cluster_on_indicators Include binary predictors in clustering
#'
#' @details
#' Returns partition centers via:
#' 1. Custom supplied centers if provided
#' 2. kmeans clustering on all variables if cluster_on_indicators=TRUE
#' 3. kmeans clustering excluding binary variables if cluster_on_indicators=FALSE
#'
#' @return Matrix of cluster centers
#'
#' @noRd
get_centers <- function(data, K, cluster_args, cluster_on_indicators) {

  ## If custom centers isn't null, return them
  if(any(!is.na(cluster_args[[1]]))){
    return(cluster_args[[1]])

  ## Partition clusters including 0/1 predictors
  } else if(cluster_on_indicators){
    km <- kmeans(data,
                 K+1,
                 cluster_args[-1])

  } else {
    bin <- which(apply(data, 2, is_binary))
    data[,bin] <- 0
    km <- kmeans(data,
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
#' @details
#' Algorithm:
#' 1. For each pair of centers (i,j):
#'    - Compute midpoint m = (c_i + c_j)/2
#'    - Compare d(m,c_i) to d(m,c_k) for all other centers k
#'    - Centers are neighbors if midpoint closer to i,j than others
#' 2. Supports parallel processing for large numbers of centers:
#'    - Divides center pairs into chunks
#'    - Processes chunks in parallel
#'    - Combines results into neighbor list
#'
#' @return List where element i contains indices of centers neighboring center i
#'
#' @noRd
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
    clusterExport(cl, "centers", envir = environment())

    ## Process chunks in parallel
    results <- parLapply(cl, chunks, function(chunk_indices) {
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
#' Partitions data support into clusters using Voronoi-like diagrams
#'
#' @param data Numeric matrix of predictor variables
#' @param cluster_args Parameters for clustering
#' @param cluster_on_indicators Logical to include binary predictors
#' @param K Number of partitions minus 1
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
#' @noRd
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
#' @param colnm_C Character vector of column names for basis expansions
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
#' A symmetric p x p penalty matrix representing integrated squared second derivatives
#' for basis expansions in a single partition of the smoothing spline.
#'
#' @details
#' Computes penalty matrix entries through analytical antiderivative calculations
#' across predictor ranges. When scaled by lambda, this becomes the smoothing
#' spline penalty matrix that controls function complexity and smoothness.
#'
#' The penalty matrix captures second derivative characteristics across:
#' - Univariate polynomial terms (quadratic, cubic, quartic)
#' - Linear interactions
#' - Quadratic interactions
#' - Three-way interactions
#'
#' @noRd
get_2ndDerivPenalty <- function(colnm_C,
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
            grep(paste0("_",v,"_"),  colnm_C[interaction_single_cols])]
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
                                       colnm_C[interaction_quad_cols])]
          if(length(interaction_quads) > 0){
            for(w in 1:length(power1_cols[-v])){

              ## the other variable, with interactions affecting quadratic terms
              wvar <- c(power1_cols[-v])[w]
              maxw <- max(C[,wvar])
              minw <- min(C[,wvar])
              diffw <- maxw - minw
              diffw2 <- maxw^2 - minw^2

              ## quadratic interaction indices
              interq <- interaction_quads[grep(colnm_C[wvar],
                                               colnm_C[interaction_quads])]
              if(length(interq) > 0){
                if(length(power2_cols) > 0){
                  ## this is the _w_x_v_^2 term
                  nchv <- nchar(colnm_C[power2_cols[v]])
                  interqv2 <- interq[substr(colnm_C[interq],
                                            nchar(colnm_C[interq]) - nchv + 1,
                                            nchar(colnm_C[interq])) ==
                                       colnm_C[power2_cols[v]]]
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

                ## verified
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


                ## verified
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

                ## verified
                if(length(interaction_single_cols) > 0){
                  interaction_singles <-
                    interaction_single_cols[grep(paste0("_",v,"_"),
                                          colnm_C[interaction_single_cols])]
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

                ## verified
                if(length(triplet_cols) > 0){
                  triplets <- triplet_cols[grep(paste0("_",v,"_"),
                                                colnm_C[triplet_cols])]
                  if(length(triplets) > 0){
                    for(tr in 1:length(triplets)){
                      other2_vars <- unlist(strsplit(colnm_C[triplets[tr]],
                                                     'x'))
                      other2_vars <- other2_vars[other2_vars !=
                                                   colnm_C[power1_cols[v]]]
                      v1 <- C[,other2_vars[1]]
                      v2 <- C[,other2_vars[2]]
                      diffv1 <- max(v1) - min(v1)
                      diffv2 <- max(v2) - min(v2)
                      if(other2_vars[1] == colnm_C[wvar]){
                        base_val1 <- 2*(diffw2 + 2*diffw*diffv2)*diff1
                        base_val2 <- 2*(diffw + diffv2)*diff2
                        mat[triplets[tr], interqv1] <- base_val1
                        mat[interqv1, triplets[tr]] <- base_val1
                        mat[triplets[tr], interqv2] <- base_val1 + base_val2
                        mat[interqv2, triplets[tr]] <- base_val1 + base_val2
                      } else if(other2_vars[2] == colnm_C[wvar]){
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
            triplet_cols[grep(paste0("_", v, "_"), colnm_C[triplet_cols])]
          if (length(triplets) > 0) {
            other2_vars <- lapply(triplets, function(tr) {
              vars <- unlist(strsplit(colnm_C[tr], 'x'))
              vars[vars != colnm_C[power1_cols[v]]]
            })
            for (tr in 1:length(other2_vars)) {

              ## the first other variable, of 3-way interaction
              maxw <- max(C[,other2_vars[[tr]][1]])
              minw <- min(C[,other2_vars[[tr]][1]])
              diffw <- maxw - minw
              diffw2 <- maxw^2 - minw^2

              ## the second other variable, of 3-way interaction
              maxu <- max(C[,other2_vars[[tr]][2]])
              minu <- min(C[,other2_vars[[tr]][2]])
              diffu <- maxu - minu
              diffu2 <- maxu^2 - minu^2

              ## Adapt this code to handle 3-way terms, i.e.
              ## the second derivative of vwu with respect to v is (w + u)
              ## integral for diagonal term =
              # int^{v = maxv}_{v = minv} w + u dv du dw => (w + u)*v
              trip_inter <- intersect(intersect(
                grep(other2_vars[[tr]][1], colnm_C),
                grep(other2_vars[[tr]][2], colnm_C)),
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
                                           colnm_C[interaction_single_cols])]
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
#' Computes smoothing spline penalty matrix with optional parallel processing
#'
#' @param K Number of partitions
#' @param colnm_C Column names of basis expansions
#' @param C Basis expansion matrix
#' @param power1_cols Linear term columns
#' @param power2_cols Quadratic term columns
#' @param power3_cols Cubic term columns
#' @param power4_cols Quartic term columns
#' @param interaction_single_cols Single interaction columns
#' @param interaction_quad_cols Quadratic interaction columns
#' @param triplet_cols Triplet interaction columns
#' @param nc Number of cubic expansions
#' @param parallel Logical to enable parallel processing
#' @param cl Cluster object for parallel computation
#'
#' @return
#' A p x p penalty matrix for smoothing spline regularization
#'
#' @noRd
get_2ndDerivPenalty_wrapper <- function(K,
                                        colnm_C,
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

  ## Trick to get the same operations performed for nonspline terms too
  colnm_C_og <- colnm_C
  if(length(nonspline_cols) > 0){
    for(jj in 1:length(nonspline_cols)){
      power1_cols <- c(power1_cols, nonspline_cols[jj])
      ## For each power, detect if already present for spline effects
      # If so, append a 0-column for the categorical variable
      # Else, skip
      if(length(power2_cols) > 0){
        colnm_C <- c(colnm_C, paste0(colnm_C[nonspline_cols[jj]],'^2'))
        C <- cbind(C, 0)
        power2_cols <- c(power2_cols, ncol(C))
      }
      if(length(power3_cols) > 0){
        colnm_C <- c(colnm_C, paste0(colnm_C[nonspline_cols[jj]],'^3'))
        C <- cbind(C, 0)
        power3_cols <- c(power3_cols, ncol(C))
      }
      if(length(power4_cols) > 0){
        colnm_C <- c(colnm_C, paste0(colnm_C[nonspline_cols[jj]],'^4'))
        C <- cbind(C, 0)
        power4_cols <- c(power4_cols, ncol(C))
      }
    }
    ## Update colnames and number of columns of expansions in C with
    # new nonspline power terms
    colnames(C) <- colnm_C
    nc <- ncol(C)
  }

  if(parallel & (K > 1)){
    ## Compute penalties for each variable, sum elementwise
    result <- Reduce("+", parLapply(cl,
                                    1:length(power1_cols),
                                    function(select_col) {
      get_2ndDerivPenalty(colnm_C,
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

  } else {
    ## Otherwise, compute serial
    result <- get_2ndDerivPenalty(colnm_C,
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
  colnames(result) <- colnm_C
  rownames(result) <- colnm_C

  ## Isolate the entries excluding appended
  result <- result[colnm_C_og, colnm_C_og]
  return(result)
}

#' Compute Log-Likelihood for Weibull Accelerated Failure Time Model
#'
#' @description
#' Calculates the log-likelihood for a Weibull Accelerated Failure Time (AFT)
#' survival model, supporting right-censored survival data.
#'
#' @param log_y Numeric vector of logarithmic response/survival times
#' @param log_mu Numeric vector of logarithmic predicted survival times
#' @param status Numeric vector of censoring indicators (1 = event, 0 = censored)
#' @param scale Numeric scalar representing the Weibull distribution scale parameter
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
#' Specifically, the log-likelihood calculation incorporates:
#' - Censoring status
#' - Scale parameter
#' - Transformed survival times
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

#' Correction for the Variance-Covariance Matrix for Uncertainty in Scale
#'
#' @description
#' Computes the shur complement "S" such that G* = (G^{-1} - S)^{-1} properly
#' accounts for uncertainty in estimating dispersion when estimating
#' variance-covariance. Otherwise, the variance-covariance matrix is optimistic
#' and assumes the scale is known, when it was in fact estimated.
#'
#' @param X Block-diagonal matrices of spline expansions
#' @param y Block-vector of response
#' @param B Block-vector of coefficient estimates
#' @param dispersion Scalar, estimate of dispersion, = Weibull scale^2
#' @param order_list List of partition orders
#' @param K Number of partitions minus 1
#' @param family Distribution family
#' @param observation_weights Optional observation weights (default = 1)
#' @param status Censoring indicator (1 = event, 0 = censored)
#'
#' @return
#' List of p x p matrices representing the shur-complement corrections to be
#' elementwise added to each block of the information matrix, before inversion.
#'
#' @details
#' Adjusts the variance-covariance matrix unscaled (G) to account for uncertainty
#' in estimating the Weibull scale parameter, that otherwise would be lost
#' if naievely using just G=(X^{T}WX + L)^{-1}.
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
      # I = ( I_bb I_bs^{T} )
      #     ( I_bs I_ss     )
      # for b = beta, s = dispersion (scale)
      #I_bb <- t(X[[k]]) %**% cbind(weights * exp_z * X[[k]])
      I_bs <- t(X[[k]]) %**% cbind(weights * zexp_z * sqrt(dispersion))
      I_ss <- -sum(
        weights * (
          (s + 2*s*z + zexp_z + exp_z * z^2)
        )
      )
      compl <- I_bs %**% matrix(-1/I_ss) %**% t(I_bs)
      # Shur complement correction to pass on to compute_G_eigen()
      return(compl)
    }
  })
}

#' Estimate Scale for Weibull Accelerated Failure Time Model
#'
#' @description
#' Computes maximum log-likelihood scale estimate of Weibull AFT survival model
#'
#' @param log_y Logarithm of response/survival times
#' @param log_mu Logarithm of predicted survival times
#' @param status Censoring indicator (1 = event, 0 = censored)
#' @param weights Optional observation weights (default = 1)
#'
#' @return
#' Scalar representing the estimated scale
#'
#' @details
#' Calculates maximum log-likelihood estimate of scale for Weibull AFT model
#' accounting for right-censored observations using Brent
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
#' Creates a compatible family object for Weibull Accelerated Failure Time (AFT)
#' models with customizable tuning options.
#'
#' @return
#' A list containing family-specific components for survival model estimation
#'
#' @details
#' Provides a comprehensive family specification for Weibull AFT models, including:
#' - Family name
#' - Link function
#' - Inverse link function
#' - Custom loss function for model tuning
#'
#' Supports right-censored survival data with flexible parameter estimation
#'
#' @export
weibull_family <- function()list(family = "weibull",
     link = "log",
     linkfun = log,
     linkinv = exp,
     custom_tuning_loss =
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
#' Computes the scale parameter for a Weibull Accelerated Failure Time (AFT)
#' model, supporting right-censored survival data.
#'
#' @param mu Predicted survival times
#' @param y Observed response/survival times
#' @param order_indices Indices to align status with response
#' @param family Weibull AFT model family specification
#' @param observation_weights Optional observation weights
#' @param status Censoring indicator (1 = event, 0 = censored)
#'
#' @return
#' Squared scale estimate for the Weibull AFT model
#'
#' @details
#' Estimates model scale through:
#' - Initial scale estimation with intercept-only model
#' - Optimization of log-likelihood for right-censored data
#' - Handling of observation weights and censoring status
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


#' Weibull glm weight function for computing the diagonal W matrix of G = (X^{T}WX + L)^{-1}
#'
#' @description
#' Function for estimating dispersion, needed for estimating coefficients,
#' when using a Weibull AFT model
#'
#' @param y Response/survival times
#' @param Mu Predicted survival times
#' @param order_indices Order of observations when partitioned to match "status" to "response"
#' @param family Weibull AFT family
#' @param observation_weights Weights of observations submitted to function
#' @param status Censoring indicator (1 = event, 0 = censored)
#'
#' @return
#' Vector of weights such that W = diag(weights)
#'
#' @details
#' Not used for unconstrained fitting, but instead, for computing G under constraint afterwards
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
#' Performs parameter update in iterative optimization with scaled matrix inversion
#' to improve computational stability.
#'
#' @param gradient_val Numeric vector of gradient values
#' @param neghessian_val Negative Hessian matrix for parameter estimation
#'
#' @return
#' Numeric vector of parameter updates
#'
#' @details
#' Applies root mean absolute value scaling to:
#' - Mitigate numerical instability
#' - Improve matrix inversion conditioning
#' - Enhance optimization convergence
#'
#' @noRd
nr_iterate <- function(gradient_val, neghessian_val){
  sc <- sqrt(mean(abs(neghessian_val))) # for computational stability
  invert(neghessian_val / sc) %**% cbind(gradient_val / sc)
}

#' Damped Newton-Raphson Parameter Optimization
#'
#' @description
#' Performs iterative parameter estimation with adaptive step-size dampening
#'
#' @param parameters Initial parameter vector
#' @param loglikelihood Function computing log-likelihood
#' @param gradient Function computing parameter gradients
#' @param neghessian Function computing negative Hessian matrix
#' @param tol Convergence tolerance (default 1e-7)
#' @param max_cnt Maximum iteration limit (default 64)
#' @param max_dmp_steps Maximum damping step attempts (default 16)
#'
#' @return
#' Optimized parameter estimates
#'
#' @details
#' Implements damped Newton-Raphson optimization:
#' - Reduces step size if objective function does not improve
#' - Handles potential numerical instability
#' - Prevents over-parameterization issues
#'
#' @noRd
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
      cat('Number of N.R. steps so far: ', master_count, '\n')
      stop('NA/NaN/non-finite value detected when running unconstrained damped Newton-Raphson. While this is expected to sometimes occur in the damp-step inner-loop, this is not the case here. The error was detected in the outer loop, likely on the very first count (check above), and most often occurs with over-parameterized models. Try re-fitting a simpler model, using greater penalties, experimenting with different knot locations, or reducing the number of knots.')
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
#' Estimates parameters for an unconstrained Weibull Accelerated Failure Time model
#' supporting right-censored survival data
#'
#' @param X Design matrix of predictors
#' @param y Survival/response times
#' @param LambdaHalf Square root of penalty matrix
#' @param Lambda Penalty matrix
#' @param keep_weighted_Lambda Flag to retain weighted penalties
#' @param family Distribution family specification
#' @param tol Convergence tolerance (default 1e-8)
#' @param K Number of partitions minus one
#' @param parallel Flag for parallel processing
#' @param cl Cluster object for parallel computation
#' @param chunk_size Processing chunk size
#' @param num_chunks Number of computational chunks
#' @param rem_chunks Remaining chunks
#' @param order_indices Observation ordering indices
#' @param weights Optional observation weights
#' @param status Censoring status indicator
#'
#' @return
#' Optimized beta parameter estimates for Weibull AFT model
#'
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
#' using penalized maximum likelihood estimation
#'
#' @param X Design matrix of predictors
#' @param y Response variable vector
#' @param LambdaHalf Square root of penalty matrix
#' @param Lambda Penalty matrix
#' @param keep_weighted_Lambda Logical flag to control penalty matrix handling:
#'   - `TRUE`: Return coefficients directly from weighted penalty fitting
#'   - `FALSE`: Apply damped Newton-Raphson optimization to refine estimates
#' @param family Distribution family specification
#' @param tol Convergence tolerance
#' @param K Number of partitions minus one
#' @param parallel Flag for parallel processing
#' @param cl Cluster object for parallel computation
#' @param chunk_size Processing chunk size
#' @param num_chunks Number of computational chunks
#' @param rem_chunks Remaining chunks
#' @param order_indices Observation ordering indices
#' @param weights Optional observation weights
#'
#' @return
#' Optimized parameter estimates for generalized linear model
#'
#' @noRd
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
  # Yields first NR updated as (X^{T}VX + L)^{-1}X^{T}Vy for V = diag(weights)
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
  if(keep_weighted_Lambda & any(class(mod) != 'try-error')){
    return(cbind(mod$coefficients))
  }

  if(any(class(mod) == 'try-error')){
    init <- c(family$linkfun(mean(y)), rep(0, ncol(X)-1))
  } else {
    init <- c(coef(mod))
  }

  ## Remove weights from Tikhinov penalties using damped nr
  res <- cbind(damped_newton_r(
    ## initial guess
    init,
    ## log-likelihood
    function(par){
      -sum(weights*family$dev.resids(
        y,
        family$linkinv(c(X %**% cbind(par))),
        wt = 1)) -
      0.5*c(t(par) %**% Lambda %**% cbind(par))
    },
    ## score
    function(par){
      c(t(X) %**% (weights*cbind(y - family$linkinv(X %**% cbind(par)))) -
          Lambda %**% cbind(par))
    },
    ## information
    function(par){
      t(X) %**% (weights*c(family$variance(X %**% cbind(par))) * X) +
        Lambda
    },
    tol))
  return(res)
}

#' Collapse Matrix List into a Single Block-Diagonal Matrix
#'
#' @description
#' Transforms a list of matrices into a single block-diagonal matrix
#'
#' @param matlist List of input matrices
#'
#' @return
#' Block-diagonal matrix combining input matrices
#'
#' @noRd
collapse_block_diagonal <- function(matlist){
  nrows <- unlist(lapply(matlist, nrow))
  ncols <- unlist(lapply(matlist, ncol))
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
#' Generates all possible interaction patterns for 2 or 3 variables
#'
#' @param vars Character vector of variable names
#'
#' @return
#' Character vector of interaction pattern strings
#'
#' @details
#' Supports generating:
#' - Linear interactions for 2 variables
#' - Quadratic interactions for 2 variables
#' - Three-way interactions for 3 variables
#'
#' @noRd
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
