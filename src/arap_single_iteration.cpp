#include "arap_single_iteration.h"

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
    //local step
    Eigen::MatrixXd C = K.transpose()*U;
    Eigen::MatrixXd R;
    R.resize(3*data.n,3);
    
    for(int k = 0; k < data.n; k++){
        Eigen::Matrix3d Ck = C.block<3,3>(3*k,0);
        Eigen::Matrix3d Rk;
        igl::polar_svd3x3(Ck,Rk);
        R.block<3,3>(3*k,0) = Rk;
    }
    
    //global step
    Eigen::VectorXd Beq;
    Eigen::MatrixXd B = K*R;
    igl::min_quad_with_fixed_solve(data,B,bc,Beq,U);
}
