#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>
#include <fstream>

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
  printf("start bisol\n");
  const Eigen::VectorXd B_flat = Eigen::VectorXd::Zero(data.n);
  igl::min_quad_with_fixed_solve(data,B_flat,bc,Eigen::VectorXd(),D);
  // REPLACE WITH YOUR CODE
  //D = Eigen::MatrixXd::Zero(data.n,3);
  printf("end bisol\n");
}

