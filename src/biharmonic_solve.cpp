#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>

void biharmonic_solve(
	const igl::min_quad_with_fixed_data<double> & data,
	const Eigen::MatrixXd & bc,
	Eigen::MatrixXd & D)
{
  // REPLACE WITH YOUR CODE
	
	int n = data.n;
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(n, 1);
	Eigen::MatrixXd Beq = Eigen::MatrixXd::Zero(0,0);
	D.setZero(n, 3);
	for (int i = 0; i < 3; i++) {
		Eigen::VectorXd d;
		igl::min_quad_with_fixed_solve(data, B, bc.col(i), Beq, d);
		D.col(i) = d;
	}
}

