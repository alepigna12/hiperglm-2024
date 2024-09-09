//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Rcpp::IntegerVector;


// Compute matrix-vector product A[row_index, ] %*% v
// [[Rcpp::export]]
VectorXd row_subset_matvec(
    const Map<MatrixXd> A, const Map<VectorXd> v, IntegerVector row_index
) {
  int num_rows = row_index.size();
  VectorXd result = VectorXd::Zero(num_rows);
  for (int i = 0; i < num_rows; i++) {
    int row = row_index[i] - 1;
    result(i) = A.row(row).dot(v);
  }
  return result;
}


// Compute matrix-vector product A[row_index, ] %*% v, but using t(A) as an input.
// [[Rcpp::export]]
VectorXd row_subset_matvec_via_transpose(
  const Map<MatrixXd> tA, const Map<VectorXd> v, IntegerVector row_index
) {
  int num_rows = row_index.size();
  VectorXd result = VectorXd::Zero(num_rows);
  for (int i = 0; i < num_rows; i++) {
    int col = row_index[i] - 1;
    result(i) = tA.col(col).dot(v);
  }
  return result;
}


// Compute matrix-vector product tA[,col_index] %*% v, using tA as an input.
// [[Rcpp::export]]
VectorXd col_subset_matvec(
    const Map<MatrixXd> tA, const Map<VectorXd> v, IntegerVector col_index
) {
  int num_cols = col_index.size();
  int dim_result = tA.rows();
  VectorXd result = VectorXd::Zero(dim_result);
  for (int i = 0; i < num_cols; i++) {
    int col = col_index[i] - 1; 
    VectorXd col_vector = tA.col(col);
    
    result += col_vector * v(i);
  }
  return result;
}
