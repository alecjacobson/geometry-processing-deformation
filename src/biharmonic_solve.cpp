#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
  // REPLACE WITH YOUR CODE

  bool success = min_quad_with_fixed_solve(
    data, Eigen::VectorXd::Zero(data.n), bc, Eigen::MatrixXd(), D);

  if (!success)
    throw std::runtime_error("[biharmonic_solve] Solve returned false.");
}

