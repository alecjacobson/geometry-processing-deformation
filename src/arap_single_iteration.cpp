#include "arap_single_iteration.h"
#include "igl/polar_svd3x3.h"
#include <igl/min_quad_with_fixed.h>
#include <iostream>

using namespace std;

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  // REPLACE WITH YOUR CODE
	// Local step
	Eigen::MatrixXd R(3*K.rows(), 3);
	const Eigen::MatrixXd C = K.transpose() * U;
	// Update each Rk
	Eigen::Matrix<double, 3, 3> Ck;
	Eigen::Matrix<double, 3, 3> Rk;
	for (int k = 0; k < K.rows(); k++) {
		Ck = C.block(3 * k, 0, 3, 3);
		igl::polar_svd3x3(Ck, Rk);
		R.block(3 * k, 0, 3, 3) = Rk;
	}

	// Global step
	Eigen::MatrixXd B = K * R;
	Eigen::MatrixXd Beq;
	igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);
}
