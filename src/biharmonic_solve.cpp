#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(data.n, 1);
  igl::min_quad_with_fixed_solve(data, B, bc, Eigen::MatrixXd(), D);
}

