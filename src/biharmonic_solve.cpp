#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
  // REPLACE WITH YOUR CODE

  // Putting D^T Q D into the general quadratic minimization form min D^T A D + D^T B s.t. D(fixed_indices) = D_fixedvals, and A_eq x = Beq
  // A in data, so B = A_eq = B_eq = 0;
  D.resize(data.n, 3);
  Eigen::VectorXd B(data.n);
  B.setZero();
  Eigen::VectorXd Beq;
  Beq.setZero();

  // Solve for the displacements in each coordinate
  for (int i = 0; i < 3; i++) {
  	Eigen::VectorXd cur_col; 
	igl::min_quad_with_fixed_solve(data, B, bc.col(i), Beq, cur_col) ;
	D.col(i) = cur_col;
  }
}

