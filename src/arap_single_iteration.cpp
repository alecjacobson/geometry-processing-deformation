#include "arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include "closest_rotation.h"

void arap_single_iteration(
	const igl::min_quad_with_fixed_data<double> & data,
	const Eigen::SparseMatrix<double> & K,
	const Eigen::MatrixXd & bc,
	Eigen::MatrixXd & U)
{

	//For this, we first do a local step, finding the optimal rotations.
	int n = U.rows();
	Eigen::MatrixXd C = K.transpose()*U;
	Eigen::MatrixXd R;
	R.resizeLike(C);
	for (int i = 0; i < n; i++) {
		Eigen::Matrix3d rot;
		closest_rotation(C.block(3 * i, 0, 3, 3), rot);
		R.block(3 * i, 0, 3, 3) = rot;

		
	}

	auto B = K*R;

	//Finally, just solve for LV = B in each of the three dimensions.
	
	for (int i = 0; i < 3; i++) {
		Eigen::VectorXd u;
		igl::min_quad_with_fixed_solve(data, B.col(i), bc.col(i), Eigen::MatrixXd(0, 0), u);
		U.col(i) = u;
	}
}
