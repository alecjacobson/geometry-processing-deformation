#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>

#include <igl/min_quad_with_fixed.h>

void biharmonic_solve(const igl::min_quad_with_fixed_data<double> & data,
		const Eigen::MatrixXd & bc, Eigen::MatrixXd & D) {
	// REPLACE WITH YOUR CODE
	D = Eigen::MatrixXd::Zero(data.n, 3);

	Eigen::VectorXd coef = Eigen::VectorXd::Ones(data.n);
	Eigen::VectorXd coef2 = Eigen::VectorXd::Ones(data.n);
	min_quad_with_fixed(data, coef, bc, coef2, D);
}

