#include "arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>

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

	MatrixXd R(3 * n, 3);

	Matrix3d r, c;

	for (int i = 0; i < n; ++i)
	{
		c = C.block(3 * i, 0, 3, 3);

		igl::polar_svd3x3(c, r);

		R.block(3 * i, 0, 3, 3) = r;
	}

	MatrixXd B = -K*R/3;

	igl::min_quad_with_fixed_solve(data, B, bc, MatrixXd(0, 0), U);
}
