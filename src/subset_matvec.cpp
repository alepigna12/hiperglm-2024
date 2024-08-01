//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Rcpp::IntegerVector;


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


// Compute matrix-vector product t(A[row_index, ]) %*% v, but using A as an input.
// [[Rcpp::export]]
VectorXd transpose_row_subset_matvec(
    const Map<MatrixXd> A, const Map<VectorXd> v, IntegerVector row_index
) {
  int num_rows = row_index.size();
  int dim_result = A.cols();
  VectorXd result = VectorXd::Zero(dim_result);
  for (int i = 0; i < num_rows; i++) {
    int row = row_index[i] - 1; 
    VectorXd row_vector = A.row(row);
    
    result += row_vector * v(i);
  }
  return result;
}