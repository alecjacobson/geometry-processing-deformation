#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>
#include <iostream>

using namespace std;

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
  // REPLACE WITH YOUR CODE
  // D = Eigen::MatrixXd::Zero(data.n,3);

  // My code
  Eigen::VectorXd B;
  B = Eigen::VectorXd::Zero(data.n); // need to init B size to pass the assertion in "slice" 
  Eigen::MatrixXd Beq;
  D.resize(data.n, 3);
  igl::min_quad_with_fixed_solve(data,B,bc,Beq,D);

}

