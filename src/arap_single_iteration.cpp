#include "arap_single_iteration.h"
#include "igl/polar_svd3x3.h"
#include <igl/min_quad_with_fixed.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  // REPLACE WITH YOUR CODE

  // Construct B
  Eigen::MatrixXd R(3*U.rows(), 3);
  Eigen::MatrixXd C = K.transpose()*U;
  Eigen::Matrix3d Rk, Ck;
  for (int i = 0; i < U.rows(); ++i) {
    Ck = C.block<3,3>(3*i, 0);
    igl::polar_svd3x3(Ck, Rk);
    R.block<3,3>(3*i, 0) = Rk;
  }
  Eigen::MatrixXd B = K*R;

  // Solve for U
  bool success = min_quad_with_fixed_solve(
    data, B, bc, Eigen::MatrixXd(), U);

  if (!success)
    throw std::runtime_error("[arap_single_iteration] Solve returned false.");

}
