#include "arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  // Get C
  Eigen::Matrix3d C = U.transpose() * K;
  // Solve for R
  Eigen::Matrix3d R;
  igl::polar_svd3x3(C, R);
  // Linear term
  Eigen::VectorXd B = K * R;
  // Empty constraints
  Eigen::VectorXd Beq;
  // Minimize energey trace( 0.5*U'*L*U + U'*B + constant )
  igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);
}
