#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>

using namespace Eigen;

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
	igl::min_quad_with_fixed_solve(data, VectorXd::Zero(data.n), bc, MatrixXd(0, 0), D);
}

