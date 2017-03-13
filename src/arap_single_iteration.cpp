#include "arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{	
	Eigen::MatrixXd  C = (U.transpose()*K).transpose();

	//Local Step
	Eigen::MatrixXd R(3 * U.rows(), 3);
	Eigen::Matrix3d cov, rot;
	for (int i = 0; i < U.rows(); i++) {
		cov = C.block(i*3, 0, 3, 3);		
		igl::polar_svd3x3<Eigen::Matrix3d>(cov, rot);
		R.block(i * 3, 0, 3, 3) = rot; //Eigen::Matrix3d::Identity();
	}

	//Global Step
	Eigen::MatrixXd B = K * R;
	igl::min_quad_with_fixed_solve(data, B, bc, Eigen::MatrixXd(), U);
	
}
