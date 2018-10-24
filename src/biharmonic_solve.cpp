#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
  // no linear term and linear constraints in biharmonic minimization problem, so set them to zero.
  D.resize(data.n, 3);
  Eigen::VectorXd B(data.n);
  B.setZero();
  Eigen::VectorXd Beq;
  Beq.setZero();

  // Solve for the displacements in each coordinate by minimizing D^T Q D with respect to boundary constraints
  for (int i = 0; i < 3; i++) {
  	Eigen::VectorXd cur_col; 
	igl::min_quad_with_fixed_solve(data, B, bc.col(i), Beq, cur_col) ;
	D.col(i) = cur_col;
  }
}

