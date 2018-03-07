#include "arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>
#include <iostream>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  Eigen::MatrixXd C = U.transpose() * K;
  Eigen::MatrixXd C_t = C.transpose();

  Eigen::MatrixXd R(3*data.n, 3);

  for (int i = 0; i < U.rows(); i+=3)
  {
	  Eigen::Matrix3d C_k;
	  Eigen::Matrix3d R_k;
	  C_k.row(0) = C_t.row(i);
	  C_k.row(1) = C_t.row(i+1);
	  C_k.row(2) = C_t.row(i+2);
	  igl::polar_svd3x3(C_k, R_k);
	  R.row(i) = R_k.row(0);
	  R.row(i+1) = R_k.row(1);
	  R.row(i+2) = R_k.row(2);
  }

  Eigen::VectorXd Beq;
  Eigen::MatrixXd B = K*R;

  igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);
}
