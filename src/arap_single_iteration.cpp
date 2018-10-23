#include "arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>

using namespace Eigen;

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  // two level optimization: local and global
  // construct C first
  MatrixXd C = K.transpose() * U;
  // construct R: a stack of rotation matrices
  MatrixXd R(3 * data.n, 3); // a stack of rotation 3x3 matrix for each vertex
  Matrix3d Ck, Rk;
  for (int k = 0; k < data.n; k++) {
    Ck = C.block(3*k, 0, 3, 3);
    igl::polar_svd3x3(Ck, Rk);
    R.block(3 * k, 0, 3, 3) = Rk;
  }
  // solve for new U
  igl::min_quad_with_fixed_solve(data, K * R, bc, MatrixXd(), U);
}
