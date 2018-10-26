#include "arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{

  // local step: R
  Eigen::MatrixXd C, R;
  C = (U.transpose() * K).transpose();
  R.resize(3 * U.rows(), 3);

  // stack the chunks of R
  Eigen::Matrix3d Ck, Rk;
  for (int k = 0; k < U.rows(); k++) {
    Ck = C.block(k * 3, 0, 3, 3);
    igl::polar_svd3x3(Ck, Rk);
    R.block(k * 3, 0, 3, 3) = Rk;
  }

  // global step: U
  Eigen::MatrixXd B, Beq;
  B = K * R;
  igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);
}
