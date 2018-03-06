#include "arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  // Get C - 3n x 3
  Eigen::MatrixXd C = K.transpose() * U;

  // Create R
  Eigen::MatrixXd R(C.rows(), C.cols());
  for(int k = 0; k < K.cols(); k+=3){
    // Solve for R_k
    Eigen::Matrix3d C_k = C.block(k,0, 3,3);
    C_k.normalize();
    Eigen::Matrix3d R_k;
    igl::polar_svd3x3(C_k, R_k);
    R.block(k,0, 3,3) = R_k;
  }
  // Linear term
  Eigen::MatrixXd B = K * R;
  // Empty constraints
  Eigen::VectorXd Beq;
  // Minimize energey trace( 0.5*U'*L*U + U'*B + constant )
  igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);
}
