#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>

using namespace Eigen;

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
  // linear coefficients: dummy
  MatrixXd B = MatrixXd::Zero(data.n,3);
  // placeholder Beq
  MatrixXd Beq;
  min_quad_with_fixed_solve(data, B, bc, Beq, D);
}
