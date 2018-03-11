#include "arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>

using namespace std;
void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  // REPLACE WITH YOUR CODE
  Eigen::MatrixXd C = (U.transpose() * K).transpose(); 
  Eigen::Matrix3d Ck, Rk;
  Eigen::MatrixXd R;
  R.resize(3*U.rows(), 3);
  for (int ii = 0; ii < U.rows(); ii++){
  	Ck = C.block(3*ii,0,3,3);
  	igl::polar_svd3x3<Eigen::Matrix3d>(Ck, Rk);
  	R.block(3*ii,0,3,3) = Rk;
  }

  Eigen::MatrixXd B;
  B = K*R;
  B /= 2.0; // Not sure where I missed a scaling factor
  Eigen::MatrixXd Beq; // empty matrix
  igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);
}
