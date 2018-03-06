#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>
#include <iostream>

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
  // REPLACE WITH YOUR CODE
  Eigen::MatrixXd Dx;
  Eigen::MatrixXd Dy;
  Eigen::MatrixXd Dz;
  min_quad_with_fixed_solve(data,Eigen::MatrixXd::Zero(data.n,1),bc.col(0),Eigen::MatrixXd(), Dx);
  min_quad_with_fixed_solve(data,Eigen::MatrixXd::Zero(data.n,1),bc.col(1),Eigen::MatrixXd(), Dy);
  min_quad_with_fixed_solve(data,Eigen::MatrixXd::Zero(data.n,1),bc.col(2),Eigen::MatrixXd(), Dz);
  D = Eigen::MatrixXd::Zero(data.n, 3);
  D.col(0) = Dx;
  D.col(1) = Dy;
  D.col(2) = Dz;
}

