#include "arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  // Construct C: C' = V'K
  Eigen::MatrixXd C = U.transpose()* K;
  C = C.transpose();
  
  // Construct R
  int n = U.rows();
  Eigen::MatrixXd R(3*n, n);
  Eigen::Matrix3d Ck, Rk;
  for (int i = 0; i < n; i++) {
    Ck = C.block(0, 3 * i, 3, 3);
    igl::polar_svd3x3(Ck, Rk);
    R.block(0, 3 * i, 3, 3) = Rk;
  }
  
  // Construct B
  Eigen::MatrixXd B;
  B = K * R;
  
  /* Setup inputs to min_quad
   * data = data
   * W = empty vector
   * Y = bc
   * Beq = empty vector
   * Z = U
   */
  Eigen::VectorXd W, Beq, x, y ,z;
  igl::min_quad_with_fixed_solve(data, W, bc.col(0), Beq, x);
  igl::min_quad_with_fixed_solve(data, W, bc.col(1), Beq, y);
  igl::min_quad_with_fixed_solve(data, W, bc.col(2), Beq, z);
  U.col(0) = x;
  U.col(1) = y;
  U.col(2) = z;
}
