#include "arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>
#include <iostream>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  // Construct C: C' = V'K
  int n = U.rows();
  Eigen::MatrixXd C(3 * n, 3);
  C = (U.transpose()* K).transpose();
  
  // Construct R
  Eigen::MatrixXd R(3 * n, 3);
  Eigen::Matrix3d Ck, Rk;
  for (int i = 0; i < n; i++) {
    Ck = C.block(3 * i, 0, 3, 3);
    igl::polar_svd3x3(Ck, Rk);
    R.block(3 * i, 0, 3, 3) = Rk;
  }
  
  // Construct B
  Eigen::MatrixXd B(n, 3);
  B = (1.0/3) * K * R;
  
  /* Setup inputs to min_quad
   * data = data
   * B = B
   * Y = bc
   * Beq = empty vector
   */
  Eigen::VectorXd Beq(n);
  igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);
}
