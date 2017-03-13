#include "arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>
#include <iostream>

using namespace Eigen;

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
	int n = data.n;

	MatrixXd C = U.transpose()*K;

	C.transposeInPlace();

	//std::cout << C << std::endl;
	//std::cout << std::endl;

	MatrixXd R(3 * n, 3);

	Matrix3d r, c;

	for (int i = 0; i < n; ++i)
	{
		c = C.block(3 * i, 0, 3, 3);

		igl::polar_svd3x3(c, r);

		//r.setIdentity();

		R.block(3 * i, 0, 3, 3) = r;

		//std::cout << (r - Matrix3d::Identity()).norm() << ' ';
	}

	//std::cout << std::endl << std::endl;

	MatrixXd B = -K*R/3;

	//std::cout << U << std::endl;
	igl::min_quad_with_fixed_solve(data, B, bc, MatrixXd(0, 0), U);
	//std::cout << U << std::endl;
	//std::flush(std::cout);
}
