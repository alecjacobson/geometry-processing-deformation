#include "arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>
#include <iostream>
using namespace std;

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  // REPLACE WITH YOUR CODE
    //Local Step
    
    int numV = data.n;
    Eigen::MatrixXd C,R, B, Beq;
    Eigen::Matrix<double, 3,3> Ck, Rk;
    C.resize(numV*3,3);
    R.resize(numV*3,3);
    B.resize(numV,3);
    C = K.transpose() * U;
    

    for (int vNo = 0; vNo < numV; vNo ++){
        Ck = C.middleRows(vNo*3,3);
        igl::polar_svd3x3(Ck, Rk);
        R.middleRows(vNo*3,3) = Rk;
    }
    
    //Global Step
    
    B = K * R;
    
    B.array() = B.array() / 6.0 ;
    Beq = Eigen::MatrixXd::Zero(1,3);
    
    min_quad_with_fixed_solve(data,B,bc,Beq,U);
    
}
