#include "arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>

void arap_single_iteration(
	const igl::min_quad_with_fixed_data<double> & data,
	const Eigen::SparseMatrix<double> & K,
	const Eigen::MatrixXd & bc,
	Eigen::MatrixXd & U)
{
	// local step
	Eigen::MatrixXd C = K.transpose() * U;
	Eigen::MatrixXd R;
	R.resizeLike(C);
	for (int i = 0; i < data.n; i++) {
		Eigen::Matrix3d Ck = C.block<3, 3>(i * 3, 0);
		Eigen::Matrix3d Rk;
		igl::polar_svd3x3(Ck, Rk);
		R.block<3, 3>(i * 3, 0) = Rk;
	}

	// global step
	Eigen::MatrixXd Beq;
	igl::min_quad_with_fixed_solve(data, K * R, bc, Beq, U);
}
