#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>

using namespace Eigen;
using namespace igl;

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
  // REPLACE WITH YOUR CODE
  // D = Eigen::MatrixXd::Zero(data.n,3);
	
	MatrixXd Du = data.ldlt.solve

}

