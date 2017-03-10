#include "arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/polar_svd3x3.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
	int n = U.rows();

	// Local optimization step
	Eigen::MatrixXd C = (U.transpose()*K).transpose();
	Eigen::Matrix3d C_i;
	Eigen::MatrixXd R(3 * n, 3);
	Eigen::Matrix3d R_i;
	for (int i = 0; i < 3*n; i += 3) {
		C_i = C.block(i, 0, 3, 3);
		igl::polar_svd3x3(C_i, R_i);
		R.block(i, 0, 3, 3) = R_i.transpose(); // might need to transpose this
	}

	// Global optimization step
	Eigen::MatrixXd B = K*R;
	Eigen::VectorXd Beq;
	igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);

}
