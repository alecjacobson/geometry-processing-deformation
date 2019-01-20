#include "biharmonic_solve.h"

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
  D.resize(data.n, 3);
  
  // from min_quad_with_fixed_solve code
  Eigen::VectorXd Beq;  
  Eigen::MatrixXd B = Eigen::VectorXd::Zero(data.n);

  igl::min_quad_with_fixed_solve(data, B, bc, Beq, D);
}