#include "arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>

void arap_single_iteration(
		const igl::min_quad_with_fixed_data<double> & data,
		const Eigen::SparseMatrix<double> & K,
		const Eigen::MatrixXd & bc,
		Eigen::MatrixXd & U) {
	int vertexCount = U.rows();

	// First we do local step. We have to solve with each individual rotation.
	Eigen::MatrixXd covariances = K.transpose() * U;
	Eigen::MatrixXd R(3 * vertexCount, 3);
	for (int index = 0; index < vertexCount; index++) {
		Eigen::Matrix3d thisCov = covariances.block<3, 3>(3 * index, 0);
		Eigen::Matrix3d thisR;
		igl::polar_svd3x3(thisCov, thisR);
		R.block<3, 3>(3 * index, 0) = thisR;
	}

	// Then we do global step
	Eigen::MatrixXd B = K * R;
	Eigen::VectorXd Beq = Eigen::VectorXd::Zero(3);
	igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);


}
