#include "arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/polar_svd3x3.h>

using namespace Eigen;
using namespace igl;

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
	//Local step
	int n = U.rows();
	MatrixXd C_t = K.transpose() * U;
	MatrixXd R(3 * n, 3);
		
	for (int i = 0; i < n; ++i) 
	{
		Matrix3d C_k = C_t.block(3 * i, 0, 3, 3);
		Matrix3d R_k;
		igl::polar_svd3x3(C_k, R_k);
		R.block(3 * i, 0, 3, 3) = R_k;		
	}
	
	// Global step
	MatrixXd B = K*R;	
	VectorXd Beq(0);
	igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);
}
