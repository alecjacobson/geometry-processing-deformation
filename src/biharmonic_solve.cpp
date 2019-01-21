#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
  // // REPLACE WITH YOUR CODE
  // D = Eigen::MatrixXd::Zero(data.n,3);

  	// Solve. The linear part of the cost function is just zero.
	Eigen::MatrixXd B (data.n, bc.cols());
	B.setZero(B.rows(), bc.cols());
	Eigen::MatrixXd Beq;

	igl::min_quad_with_fixed_solve(data, B, bc, Beq, D);

	return;
}

