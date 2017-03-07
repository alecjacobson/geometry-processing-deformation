#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
  /* Setup inputs to min_quad
   * data = data
   * B = empty vector
   * Y = b by 1 list of constant fixed values
   * Beq = empty vector
   * Z = D
   */
  Eigen::VectorXd B, Beq, x, y ,z;
  igl::min_quad_with_fixed_solve(data, B, bc.col(0), Beq, x);
  igl::min_quad_with_fixed_solve(data, B, bc.col(1), Beq, y);
  igl::min_quad_with_fixed_solve(data, B, bc.col(2), Beq, z);
  
  D.resize(x.rows(), 3);
  D.col(0) = x;
  D.col(1) = y;
  D.col(2) = z;
}

