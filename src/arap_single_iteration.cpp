#include "arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U) {
    
    // Local:
    // compute the matrix of stacked weighted covariance matrices.
    Eigen::MatrixXd C = (U.transpose()*K).transpose();
    
    // compute the closest rotation matrix Rk for each Ck and stack them
    // to produce R.
    int n = U.rows();
    Eigen::MatrixXd R(3*n, 3);
    
    for (int i = 0; i < n; i++) {
        
        Eigen::Matrix3d Ck = C.block(i*3, 0, 3, 3);
        Eigen::Matrix3d Rk;
        igl::polar_svd3x3(Ck, Rk);
        
        R.block(i*3, 0, 3, 3) = Rk;
    }
    
    // Global:
    // // solve the minimization problem with bc constraints using the precomputed data.
    Eigen::MatrixXd B = K*R;
    Eigen::MatrixXd Beq;
    igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);
}
