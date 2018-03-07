#include "arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  int num_data_points = data.n;

  // Compute local weighted covariance matrices for each point
  // Transposing from the ReadMe because we will perform SVD on each 3x3 block
  Eigen::MatrixXd weighted_cov = K.transpose() * U;
  
  // We will solve for a closest rotation (3x3) for each of point (lives in point space)
  Eigen::MatrixXd R(3 * num_data_points, 3);

  for (int i = 0; i < num_data_points; i++) 
  {
    Eigen::Matrix3d cov_i = weighted_cov.block<3, 3>(3 * i, 0);

    Eigen::Matrix3d best_rotation;
    // Computes the closest rotation (best_rotation) to the input matrix cov_i 
    igl::polar_svd3x3(cov_i, best_rotation);
    
    R.block<3, 3>(3 * i, 0) = best_rotation;
  }

  Eigen::MatrixXd B = K * R;
  Eigen::VectorXd Beq;
  igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);
}
