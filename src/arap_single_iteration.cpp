#include "arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  // REPLACE WITH YOUR CODE
    //Local Step
    
    //Extract number of vertices
    int numV = data.n;
    
    //C as defined in Readme
    //R stack of rotation matrices
    //B as defined in Readme
    //Beq placeholder matrices
    Eigen::MatrixXd C,R, B, Beq;
    
    //Ck is extracted 3x3 matrix from C
    //Rk is extracted 3x3 matrix from R
    
    Eigen::Matrix<double, 3,3> Ck, Rk;
    C.resize(numV*3,3);
    R.resize(numV*3,3);
    B.resize(numV,3);
    
    //compute C as in Readme
    C = K.transpose() * U;
    

    //Use polar_svd3x3 to compute closest rotation
    for (int vNo = 0; vNo < numV; vNo ++){
        Ck = C.middleRows(vNo*3,3);
        igl::polar_svd3x3(Ck, Rk);
        R.middleRows(vNo*3,3) = Rk;
    }
    
    //Global Step
    
    //Compute B as in ReadMe
    B = K * R;
    
    B.array() = B.array() / 6.0 ;
    Beq = Eigen::MatrixXd::Zero(1,3);
    
    //Solve the optimization
    min_quad_with_fixed_solve(data,B,bc,Beq,U);
    
}
