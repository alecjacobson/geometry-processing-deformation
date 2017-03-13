#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>

#include <iostream>

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
    // REPLACE WITH YOUR CODE
    //std::cout<<"data n is: "<<data.n<<std::endl;
    D = Eigen::MatrixXd::Zero(data.n,3);
    //std::cout<<"bc is: "<<bc<<std::endl;
    
    Eigen::VectorXd B = Eigen::VectorXd::Zero(data.n, 1 );
    Eigen::VectorXd Beq, Dx, Dy, Dz; // empty equality constraints
    igl::min_quad_with_fixed_solve( data, B, bc, Beq, D );
}

