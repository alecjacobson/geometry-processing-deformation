#include "closest_rotation.h"
#include <Eigen/Dense>

void closest_rotation(
	const Eigen::Matrix3d & M,
	Eigen::Matrix3d & R)
{
	// Replace with your code
	//This is equivalent to using the SVD to write M=U.S.V^T
	//Then R = U.V^T

	Eigen::JacobiSVD<Eigen::Matrix3d> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix3d sigma = Eigen::Matrix3d::Identity();
	sigma(2, 2) = (svd.matrixU()*svd.matrixV().transpose()).determinant();

	R = svd.matrixU()*sigma*svd.matrixV().transpose();
}