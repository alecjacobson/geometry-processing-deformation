#include "arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>
#include <iostream>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  // REPLACE WITH YOUR CODE
	int vetr_num = U.rows();
	Eigen::MatrixXd C;
	Eigen::MatrixXd R;
    Eigen::Matrix3d C2;
    Eigen::Matrix3d R2;
    Eigen::MatrixXd B;
    Eigen::MatrixXd Beq;

	C = K.transpose() * U;
	std::cout << C.rows() << std::endl;
    for (int i = 0; i < vetr_num; i++) {
		C2 = C.block(0, i*3, 3, 3);
		igl::polar_svd3x3(C2, R2);
		
		R.block(i*3, 0, 3, 3) = R2.transpose();
	}
	B = K * R;
    std::cout << B << std::endl;
	igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);
}
