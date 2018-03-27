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
    int numV;
    numV = data.n;
    //cout << numV << endl;
    Eigen::MatrixXd Beq, B;
    Beq = Eigen::MatrixXd::Zero(1,3);
    B = Eigen::MatrixXd::Zero(numV,3);
    D.resize(numV, 3);
    min_quad_with_fixed_solve(data,B,bc,Beq,D);
}

