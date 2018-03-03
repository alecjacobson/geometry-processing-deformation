#include "biharmonic_solve.h"
#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
  // Empty constraints
  Eigen::VectorXd Beq;
  // Linear term is 0 except at known values
  Eigen::VectorXd B  = Eigen::VectorXd::Zero(data.n);
  igl::min_quad_with_fixed_solve(data, B, bc, Beq, D);
}
