#include "arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include <Eigen/SVD>
#include <Eigen/Dense>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
	Eigen::MatrixXd C = U.transpose() * K;
	Eigen::MatrixXd R = Eigen::MatrixXd::Zero(3, 3 * U.cols());
	for (int i = 0; i < U.cols(); i++) {
		Eigen::MatrixXd c = C.block(3*i, 0, 3, 3);
		Eigen::MatrixXd r;
		closest_rotation(c, r);
		R.block(3*i, 0, 3, 3) = r;
	}
	Eigen::MatrixXd B = K * R;
	Eigen::MatrixXd Ux;
	Eigen::MatrixXd Uy;
	Eigen::MatrixXd Uz;
	min_quad_with_fixed_solve(data,B,bc.col(0),Eigen::MatrixXd(), Ux);
	min_quad_with_fixed_solve(data,B,bc.col(1),Eigen::MatrixXd(), Uy);
	min_quad_with_fixed_solve(data,B,bc.col(2),Eigen::MatrixXd(), Uz);
	U.col(0) = Ux;
	U.col(1) = Uy;
	U.col(2) = Uz;
}

void closest_rotation(
  const Eigen::MatrixXd & M,
  Eigen::MatrixXd & R)
{
	R = Eigen::Matrix3d::Identity();
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::Matrix3d omega;
	omega << 1, 0, 0,
	0, 1, 0,
	0, 0, (svd.matrixU() * svd.matrixV().transpose()).determinant();
	R = svd.matrixV() * omega * svd.matrixU().transpose();
}
