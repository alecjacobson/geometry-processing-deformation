#include <igl/min_quad_with_fixed.h>
#include "arap_single_iteration.h"
#include <igl/polar_svd3x3.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  Eigen::MatrixXd C = K.transpose() * U;
  Eigen::MatrixXd R(3 * U.rows(), 3);

  for (int k = 0; k < U.rows(); k++) {
    Eigen::Matrix3d Ck, Rk;
    for (int i = 0; i < 3; i++) {
      Ck.row(i) = C.row(3 * k + i);
    }
    igl::polar_svd3x3(Ck, Rk);
    for (int i = 0; i < 3; i++) {
      R.row(3 * k + i) = Rk.row(i);
    }
  }
  igl::min_quad_with_fixed_solve(data, K * R, bc, Eigen::MatrixXd(), U);
}
