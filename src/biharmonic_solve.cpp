#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>

using namespace Eigen;
using namespace igl;

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{	
	MatrixXd B = MatrixXd::Zero(data.n, bc.cols());
	VectorXd Beq(0);
	bool result = min_quad_with_fixed_solve(data, B, bc, Beq, D);
	assert(result);
}

