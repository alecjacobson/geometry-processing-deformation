#include "arap_single_iteration.h"
#include <igl/polar_svd.h>
#include <igl/min_quad_with_fixed.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  // REPLACE WITH YOUR CODE

  // local step
  Eigen::MatrixXd C = K.transpose() * U; 
  Eigen::MatrixXd R(3 * U.rows(), 3);
  for (int i = 0; i < U.rows(); i++) {
	// find closest rotation matrix to each C_k
        Eigen::Matrix3d C_i = C.block(3 * i, 0, 3, 3);
        Eigen::Matrix3d R_i;
        Eigen::Matrix3d T_i;
	Eigen::Matrix3d U_i;
	Eigen::Vector3d S_i;
	Eigen::Matrix3d V_i;
	igl::polar_svd(C_i, R_i, T_i, U_i, S_i, V_i);
	R.block(3 *i, 0, 3, 3) = R_i;
  }

  // global step
  // Quadratic minimization of 1/2 V^T L V + V^T (K * R) / 2 subject to handle constraints
    Eigen::VectorXd Beq;
    Beq.setZero();
    Eigen::MatrixXd B = 0.5 * K * R ;
  
  // Solve for the displacements in each coordinate
  for (int i = 0; i < 3; i++) {
	Eigen::VectorXd cur_col;
	igl::min_quad_with_fixed_solve(data, B.col(i), bc.col(i), Beq, cur_col) ;
	U.col(i) = cur_col;
  } 
}
