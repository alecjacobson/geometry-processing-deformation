#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{

    //get size
    int numV;
    numV = data.n;
    
    //Place holder matrices
    Eigen::MatrixXd Beq, B;
    Beq = Eigen::MatrixXd::Zero(1,3);
    B = Eigen::MatrixXd::Zero(numV,3);
    
    //Make D the correct size
    D.resize(numV, 3);
    
    //solve the optimization
    min_quad_with_fixed_solve(data,B,bc,Beq,D);
}

