#include "arap_single_iteration.h"

#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  Eigen::MatrixXd C = K.transpose() * U;
  Eigen::MatrixXd R(C.rows(), C.cols());

  for(int k = 0; k < data.n; k++) {
  	Eigen::Matrix3d C_k = C.block(k * 3, 0, 3, 3);
  	Eigen::Matrix3d R_k;
  	igl::polar_svd3x3(C_k, R_k);
  	R.block(k * 3, 0, 3, 3) = R_k;
  }

  igl::min_quad_with_fixed_solve(data, K*R, bc,  Eigen::VectorXd(), U);
}
