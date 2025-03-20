// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Efficient Matrix Multiplication
//'
//' Performs matrix multiplication using RcppArmadillo
//'
//' @param A First input matrix
//' @param B Second input matrix
//' @return Matrix product of A and B
//' @export
//' @keywords internal
// [[Rcpp::export]]
arma::mat efficient_matrix_mult(const arma::mat& A, const arma::mat& B) {
  return A * B;
}

//' Compute Gram Matrix
//'
//' Calculates X^T * X for the input matrix X
//'
//' @param X Input matrix
//' @return Gram matrix (X^T * X)
//' @export
//' @keywords internal
// [[Rcpp::export]]
arma::mat gramMatrix(const arma::mat& X) {
    return X.t() * X;
}

//' Matrix Inversion using Armadillo
//'
//' Computes the inverse of a matrix using Armadillo's inversion method
//'
//' @param x Input matrix to be inverted
//' @return Inverted matrix
//' @export
//' @keywords internal
// [[Rcpp::export]]
arma::mat armaInv(const arma::mat& x) {
    return arma::inv(x);
}

//' Block Matrix Multiplication
//'
//' Performs multiplication of a list of matrices with a block matrix
//'
//' @param G List of matrices G
//' @param A Input matrix A
//' @param K Number of blocks
//' @param nc Number of columns
//' @param nca Number of columns in A
//' @return List of multiplied matrices
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List GAmult(const Rcpp::List& G,
                  const arma::mat& A,
                  int K,
                  int nc,
                  int nca) {
  Rcpp::List result(K + 1);
  for (int k = 0; k < K + 1; ++k) {
    arma::mat G_k = Rcpp::as<arma::mat>(G[k]);
    int start_row = k * nc;
    int end_row = start_row + nc - 1;
    arma::mat A_slice = A.rows(start_row, end_row);
    result[k] = G_k * A_slice;
  }
  return result;
}

//' Compute Trace Correction
//'
//' Calculates a trace correction for matrix computations
//'
//' @param G List of matrices G
//' @param A Input matrix A
//' @param GXX List of GX^{T}X matrices
//' @param AGAInv Matrix (A^{T}GA)^{-1}
//' @param nc Number of columns of each partition of G
//' @param K Number of blocks
//' @return Trace correction value
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
double compute_trace_correction(const Rcpp::List& G,
                              const arma::mat& A,
                              const Rcpp::List& GXX,
                              const arma::mat& AGAInv,
                              int nc,
                              int K) {
    double correction = 0.0;
    
    // Kahan summation for enhanced numerical stability
    double c = 0.0;
    
    // Compute overall scaling to prevent over/underflow
    double max_norm = 0.0;
    for (int k = 0; k <= K; ++k) {
        arma::mat G_k = Rcpp::as<arma::mat>(G[k]);
        max_norm = std::max(max_norm, arma::norm(G_k, "fro"));
    }
    
    // Add small epsilon to prevent division by zero
    const double eps = std::numeric_limits<double>::epsilon();
    max_norm = std::max(max_norm, eps);
    
    for (int k = 0; k <= K; ++k) {
        arma::mat G_k = Rcpp::as<arma::mat>(G[k]);
        arma::mat A_k = A.rows(k*nc, (k+1)*nc-1);
        arma::mat GXX_k = Rcpp::as<arma::mat>(GXX[k]);
        
        // Normalize matrices relative to max_norm to prevent scaling issues
        G_k /= max_norm;
        A_k /= max_norm;
        GXX_k /= max_norm;
        
        // Check for zero/near-zero matrices
        if (arma::norm(G_k, "fro") > eps && 
            arma::norm(A_k, "fro") > eps && 
            arma::norm(GXX_k, "fro") > eps) {
            
            // Compute intermediate matrix products more carefully
            arma::mat temp1 = G_k * A_k;
            arma::mat temp2 = temp1 * AGAInv;
            arma::mat temp3 = A_k.t() * GXX_k;
            
            // Compute trace
            double trace_k = arma::trace(temp2 * temp3);
            
            // Kahan summation for enhanced numerical stability
            double y = trace_k - c;
            double t = correction + y;
            c = (t - correction) - y;
            correction = t;
        }
    }
    
    return correction;
}

//' AGAmult Computation Overall
//'
//' Performs A^{T}GA computation over all matrices in a list (G)
//'
//' @param G_chunk List of matrices for the chunk
//' @param A Input matrix A
//' @param K Number of partitions minus one
//' @param nc Number of columns
//' @return Resulting matrix
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
arma::mat AGAmult(const Rcpp::List& G,
  const arma::mat& A,
  int K,
  int nc,
  int nca) {
              arma::mat result(nca, nca, arma::fill::zeros);

              for (int k = 0; k < K + 1; ++k) {
                int start = k * nc;
                int finish = (k + 1) * nc - 1;

                arma::mat Asub = A.rows(start, finish);
                arma::uvec non0 = arma::find(arma::sum(arma::abs(Asub), 0) > 0);

                if (non0.n_elem > 0) {
                  Asub = Asub.cols(non0);
                  arma::mat G_k = Rcpp::as<arma::mat>(G[k]);
                  arma::mat GAsub = G_k * Asub;
                  arma::mat AGAsub = Asub.t() * GAsub;

                  result.submat(non0, non0) += AGAsub;
                }
              }

              return result;
}

//' Chunk-based AGAmult Computation
//'
//' Performs A^{T}GAmult computation on a specific chunk of matrices
//'
//' @param G_chunk List of matrices for the chunk
//' @param A Input matrix A
//' @param chunk_start Starting chunk index
//' @param chunk_end Ending chunk index
//' @param nc Number of columns
//' @return Resulting matrix
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
arma::mat AGAmult_chunk(const Rcpp::List& G_chunk,
                        const arma::mat& A,
                        int chunk_start,
                        int chunk_end,
                        int nc) {
  int nca = A.n_cols;
  arma::mat result(nca, nca, arma::fill::zeros);
  for (int i = 0; i < G_chunk.length(); ++i) {
    int k = chunk_start + i;
    arma::mat G_k = Rcpp::as<arma::mat>(G_chunk[i]);
    arma::mat A_k = A.rows(k * nc, (k + 1) * nc - 1);
    arma::mat GA_k = G_k * A_k;
    result += A_k.t() * GA_k;
  }
  return result;
}

//' Compute AGXy
//'
//' Computes A^T * G * Xy for specified range
//'
//' @param G List of G matrices
//' @param A Input matrix A
//' @param Xy List of Xy vectors
//' @param nc Number of columns
//' @param K Total number of blocks
//' @param start Start index
//' @param end End index
//' @return Resulting vector
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
arma::vec compute_AGXy(const Rcpp::List& G,
                       const arma::mat& A,
                       const Rcpp::List& Xy,
                       int nc,
                       int K,
                       int start,
                       int end) {
    int nca = A.n_cols;
    arma::vec AGXy(nca, arma::fill::zeros);
    for (int k = start; k <= end; ++k) {
        arma::mat G_k = Rcpp::as<arma::mat>(G[k-start]);
        arma::vec Xy_k = Rcpp::as<arma::vec>(Xy[k-start]);
        arma::mat A_k = A.rows(k*nc, (k+1)*nc-1);
        AGXy += A_k.t() * (G_k * Xy_k);
    }
    return AGXy;
}

//' Compute Result Blocks
//'
//' Computes blocks of results for specific computations
//'
//' @param G List of G matrices
//' @param Ghalf List of Ghalf matrices
//' @param A Input matrix A
//' @param AAGAInvAGXy Input vector
//' @param nc Number of columns
//' @param start Start index
//' @param end End index
//' @return Resulting vector of blocks
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
arma::vec compute_result_blocks(const Rcpp::List& G,
                                const Rcpp::List& Ghalf,
                                const arma::mat& A,
                                const arma::vec& AAGAInvAGXy,
                                int nc,
                                int start,
                                int end) {
  int block_count = end - start + 1;
  arma::vec result(nc * block_count);
  for (int k = 0; k < block_count; ++k) {
    arma::mat G_k = Rcpp::as<arma::mat>(G[k]);
    arma::mat Ghalf_k = Rcpp::as<arma::mat>(Ghalf[k]);
    arma::vec block_input = AAGAInvAGXy.subvec((k+start)*nc,
                                               (k+start+1)*nc-1);
    result.subvec(k*nc, (k+1)*nc-1) = Ghalf_k * (G_k * block_input);
  }
  return result;
}

//' Block Diagonal Matrix Multiplication
//'
//' Performs multiplication of block diagonal matrices
//'
//' @param A List of matrices A
//' @param B List of matrices B
//' @param K Number of blocks
//' @return List of multiplied matrices
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List matmult_block_diagonal_cpp(const Rcpp::List& A,
                                      const Rcpp::List& B,
                                      int K) {
    Rcpp::List result(K + 1);

    for (int k = 0; k <= K; ++k) {
        arma::mat A_k = Rcpp::as<arma::mat>(A[k]);
        arma::mat B_k = Rcpp::as<arma::mat>(B[k]);
        result[k] = A_k * B_k;
    }

    return result;
}

//' Vector-Matrix Multiplication for Block Diagonal Matrices
//'
//' Performs vector-matrix multiplication for block diagonal matrices
//'
//' @param A List of matrices A
//' @param b List of vectors b
//' @param K Number of blocks
//' @return List of resulting vectors
//' @export
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List vectorproduct_block_diagonal(const Rcpp::List& A,
                                        const Rcpp::List& b,
                                        int K) {
  Rcpp::List result(K + 1);

  for (int k = 0; k <= K; ++k) {
    arma::mat A_k = Rcpp::as<arma::mat>(A[k]);
    arma::vec b_k = Rcpp::as<arma::vec>(b[k]);
    result[k] = A_k.t() * b_k;
  }

  return result;
}

//' Matrix-Matrix Addition for Block Diagonal Matrices
//'
//' Performs matrix-matrix elementwise addition for block diagonal matrices
//'
//' @param A List of matrices A
//' @param B List of matrices B
//' @param K Number of blocks
//' @return List of resulting matrices
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List matadd_block_diagonal(const Rcpp::List& A,
                                  const Rcpp::List& B,
                                  int K) {
    Rcpp::List result(K + 1);

    for (int k = 0; k <= K; ++k) {
        arma::mat A_k = Rcpp::as<arma::mat>(A[k]);
        arma::mat B_k = Rcpp::as<arma::mat>(B[k]);
        result[k] = A_k + B_k;
    }

    return result;
}

//' Compute Derivative of W
//'
//' Computes derivative of diag(1/(1-XUGXt)) with respect to lambda, penalties of G
//'
//' @param A List of matrices A
//' @param GXX List of matrices of G times gram matrix of X transpose for each partition
//' @param Ghalf List of matrices of square-root of G for each partition
//' @param dG_dlambda List of derivatives of G matrix with respect to penalties lambda
//' @param dGhalf_dlambda List of derivatives of Ghalf matrix with respect to penalties lambda
//' @param AGAInv Matrix of A transpose times G times A within U = (I - GA(A^{T}GA)^{-1}A^{T})
//' @param nc Numeric numeric of basis expansions of predictors per partition
//' @param K Number of blocks
//' @return Matrix representing the derivative
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
double compute_dW_dlambda(const Rcpp::List& G,
                                    const arma::mat& A,
                                    const Rcpp::List& GXX,
                                    const Rcpp::List& Ghalf,
                                    const Rcpp::List& dG_dlambda,
                                    const Rcpp::List& dGhalf_dlambda,
                                    const arma::mat& AGAInv, int nc, int K) {
    // First compute tr(dG/dN; * X\'X)
    double trace = 0.0;
    for (int k = 0; k <= K; ++k) {
        arma::mat GXX_k = Rcpp::as<arma::mat>(GXX[k]);
        arma::mat dG_k = Rcpp::as<arma::mat>(dG_dlambda[k]);
        trace += arma::trace(dG_k * GXX_k);
    }

    // Compute correction terms efficiently
    double correction1 = 0.0;
    double correction2 = 0.0;

    // First correction term: tr(dG/dN; * A(A\'GA)^(-1)A\'GX\'X)
    for (int k = 0; k <= K; ++k) {
        arma::mat dGhalf_k = Rcpp::as<arma::mat>(dGhalf_dlambda[k]);
        arma::mat Ghalf_k = Rcpp::as<arma::mat>(Ghalf[k]);
        arma::mat A_k = A.rows(k*nc, (k+1)*nc-1);
        arma::mat GXX_k = Rcpp::as<arma::mat>(GXX[k]);

        arma::mat temp = dGhalf_k * A_k * AGAInv;
        arma::mat temp2 = A_k.t() * Ghalf_k.t() * GXX_k;
        correction1 += arma::trace(temp * temp2);
    }

    // Second correction term: tr(G * A(A\'GA)^(-1)A\'dG/dN; * X\'X)
    for (int k = 0; k <= K; ++k) {
        arma::mat Ghalf_k = Rcpp::as<arma::mat>(Ghalf[k]);
        arma::mat A_k = A.rows(k*nc, (k+1)*nc-1);
        arma::mat GXX_k = Rcpp::as<arma::mat>(GXX[k]);
        arma::mat dG_k = Rcpp::as<arma::mat>(dG_dlambda[k]);

        arma::mat temp = Ghalf_k * A_k * AGAInv;
        arma::mat temp2 = A_k.t() * dG_k * GXX_k;
        correction2 += arma::trace(temp * temp2);
    }

    // Return final result
    return trace - correction1 - correction2;
                                    }

//' Compute G^{1/2}(GA(A^{T}GA)^{-1}A^{T})GX^{T}y 
//'
//' Computes the equivalent vector of G^{1/2}X^{T}y in the proposed model-fitting procedure but for G^{1/2}(GA(A^{T}GA)^{-1}A^{T})GX^{T}y instead.
//'
//' @param G List of matrices of G times gram matrix of X transpose for each partition
//' @param Ghalf List of matrices of square-root of G for each partition
//' @param A List of matrices A
//' @param AGAInv Matrix of A transpose times G times A within U = (I - GA(A^{T}GA)^{-1}A^{T})
//' @param Xy List of vectors for partition-wise dot product of transpose(X) times vector y
//' @param nc Numeric numeric of basis expansions of predictors per partition
//' @param K Number of blocks
//' @return Matrix representing the operation
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
arma::vec compute_GhalfXy_temp(const Rcpp::List& G, const Rcpp::List& Ghalf,
                              const arma::mat& A, const arma::mat& AGAInv,
                              const Rcpp::List& Xy, int nc, int K) {
    int nca = A.n_cols;
    arma::vec result(nc * (K + 1));

    // Compute A^{T}GXy efficiently
    arma::vec AGXy(nca, arma::fill::zeros);
    for (int k = 0; k <= K; ++k) {
        arma::mat G_k = Rcpp::as<arma::mat>(G[k]);
        arma::vec Xy_k = Rcpp::as<arma::vec>(Xy[k]);
        arma::mat A_k = A.rows(k*nc, (k+1)*nc-1);

        AGXy += A_k.t() * (G_k * Xy_k);
    }

    // Compute A\'AGAInvAGXy once
    arma::vec AAGAInvAGXy = A * (AGAInv * AGXy);

    // Compute final result block by block
    for (int k = 0; k <= K; ++k) {
        arma::mat G_k = Rcpp::as<arma::mat>(G[k]);
        arma::mat Ghalf_k = Rcpp::as<arma::mat>(Ghalf[k]);
        arma::vec block_input = AAGAInvAGXy.subvec(k*nc, (k+1)*nc-1);
        result.subvec(k*nc, (k+1)*nc-1) = Ghalf_k * (G_k * block_input);
    }

    return result;
}

//' Left-Multiply a List of Block-Diagonal Matrices by U 
//'
//' Useful for computing UG where G is a list of block-diagonal matrices for each partition, and U is a square P by P dense matrix
//'
//' @param U Matrix of dimension P by P that projects onto the null space of A transpose
//' @param G List of G matrices for each partition
//' @param nc Numeric numeric of basis expansions of predictors per partition
//' @param K Number of blocks
//' @return Matrix of UG
//' @export
//' @keywords internal
// [[Rcpp::export]]
arma::mat matmult_U(const arma::mat& U, const Rcpp::List& G, int nc, int K) {
    int n_rows_U = U.n_rows;
    int n_cols_result = 0;

    // Calculate the total number of columns in the result
    for (int k = 0; k <= K; ++k) {
        arma::mat G_k = Rcpp::as<arma::mat>(G[k]);
        n_cols_result += G_k.n_cols;
    }

    arma::mat result(n_rows_U, n_cols_result);
    int col_offset = 0;

    for (int k = 0; k <= K; ++k) {
        arma::mat G_k = Rcpp::as<arma::mat>(G[k]);
        int start_col = k * nc;
        int end_col = start_col + nc - 1;

        arma::mat U_slice = U.cols(start_col, end_col);
        arma::mat product = U_slice * G_k;

        result.cols(col_offset, col_offset + G_k.n_cols - 1) = product;
        col_offset += G_k.n_cols;
    }

    return result;
}
