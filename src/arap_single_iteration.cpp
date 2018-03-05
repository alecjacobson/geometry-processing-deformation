#include "arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  // REPLACE WITH YOUR CODE
	// Compute closest rotation R
	Eigen::MatrixXd C_T = U.transpose()*K;
	Eigen::MatrixXd R(3 * U.rows(), 3);
	for (int k = 0; k < U.rows(); k++) {
		Eigen::Matrix3d C_k = C_T.transpose().block<3, 3>(3 * k, 0);
		C_k = C_k / C_k.maxCoeff(); // Normalize C_k
		Eigen::Matrix3d R_k;
		igl::polar_svd3x3(C_k, R_k);
		R.block<3, 3>(3 * k, 0) = R_k;
	}
	// Compute B=KR
	Eigen::MatrixXd B = K * R;

	// Compute single arap minimization
	Eigen::VectorXd Beq;
	igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);
}
