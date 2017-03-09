#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>
#include <iostream>

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
  /* Setup inputs to min_quad
   * data = data
   * B = empty vector
   * Y = bc
   * Beq = empty vector
   * Z = D
   */
  Eigen::VectorXd B, Beq, x, y ,z;
  igl::min_quad_with_fixed_solve(data, B, bc.col(0), Beq, x);
  std::cout << "solve here" << std::endl;
  igl::min_quad_with_fixed_solve(data, B, bc.col(1), Beq, y);
  std::cout << "solve here 2" << std::endl;
  igl::min_quad_with_fixed_solve(data, B, bc.col(2), Beq, z);
  std::cout << "solve here 3" << std::endl;
  
  D.resize(x.rows(), 3);
  D.col(0) = x;
  D.col(1) = y;
  D.col(2) = z;
  std::cout << "solve here 4" << std::endl;
}

