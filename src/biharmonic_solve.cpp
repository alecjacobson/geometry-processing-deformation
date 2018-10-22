#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
  igl::min_quad_with_fixed_solve(data, Eigen::MatrixXd(data.n, 3), bc, Eigen::MatrixXd(), D);
}

