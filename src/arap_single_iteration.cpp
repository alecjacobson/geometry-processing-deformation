#include "arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>
void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  int n = U.rows();
  //Local step
  
  Eigen::MatrixXd C = (U.transpose() * K).transpose();
  Eigen::MatrixXd R(3 *n, 3);
  for (int i = 0; i < n; i++){
  	Eigen::Matrix3d temp = C.block(i *3, 0, 3, 3);
  	Eigen::Matrix3d temp2 = R.block(i * 3,0, 3, 3);
  	igl::polar_svd3x3(temp, temp2);
  	R.block(i * 3,0, 3, 3) = temp2;
  }
  
  //Global step
  Eigen::MatrixXd dummy2(0,0); //dummy for placeholding
  Eigen::MatrixXd B = K * R /3;
  igl::min_quad_with_fixed_solve(data, B, bc, dummy2, U);
}
