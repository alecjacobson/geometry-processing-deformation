#include "arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>
#include <Eigen/Dense>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  Eigen::MatrixXd C = K; //bc.transpose() * K;

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(C, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::Matrix3d v = svd.matrixV();
  Eigen::Matrix3d u = svd.matrixU();
  Eigen::Matrix3d omega;
  Eigen::Matrix3d UVt = u * v.transpose();
  omega << 1, 0, 0,
		0, 1, 0,
		0, 0, UVt.determinant();
	
  Eigen::MatrixXd R = u * omega * v.transpose();

  Eigen::VectorXd Beq;
  Eigen::VectorXd B = K*R;

  igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);
}
