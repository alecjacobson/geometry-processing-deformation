#include "arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & V)
{
  // REPLACE WITH YOUR CODE
  Eigen::MatrixXd C = V.transpose() * K;
  Eigen::MatrixXd R(3*data.n, 3);
  Eigen::Matrix3d C_k, R_k;
  for (int k = 0; k < data.n; k++)
  {
    C_k = C.block(0,3*k,3,3);
    igl::polar_svd3x3(C_k, R_k);
    R.block(3*k,0,3,3) = R_k.transpose();
  }

  Eigen::MatrixXd B=K*R;
  Eigen::MatrixXd Beq;
  igl::min_quad_with_fixed_solve(data, B, bc, Beq, V);
}
