#include "arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
    //local step
    Eigen::MatrixXd C, R, X;
    C = K.transpose() * U;
    R.resize(data.n * 3, 3);
    
    //compute R
    for (int i = 0; i < data.n; i++) {
        Eigen::Matrix3d Rk, Ck;
        Ck = C.block(i * 3, 0, 3, 3);
        igl::polar_svd3x3(Ck, Rk);
        R.block(i * 3, 0, 3, 3) = Rk;
    }
    
    igl::min_quad_with_fixed_solve(data, K * R, bc, X, U);
}
