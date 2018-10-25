#include "arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  // REPLACE WITH YOUR CODE
  Eigen::MatrixXd C = (U.transpose() * K).transpose();
  Eigen::MatrixXd R(3 * U.rows(), 3);

  for (int i = 0; i < U.rows(); i++)
  {
    Eigen::Matrix3d R_block;
    Eigen::Matrix3d Ck = C.block(i * 3, 0, 3, 3);
    igl::polar_svd3x3(Ck, R_block);
    R.block(i * 3, 0, 3, 3) = R_block;
  }

  Eigen::MatrixXd Beq;
  Eigen::MatrixXd B = K * R;
  igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);
}
