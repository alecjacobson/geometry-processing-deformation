#include "arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/polar_svd.h>
#include <igl/svd3x3.h>
#include <igl/min_quad_with_fixed.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  // REPLACE WITH YOUR CODE

  // First do the local step
  // For local step, we need to minimize tr(V'KR) w.r.t R.
  // Note the matrix property: tr(A'B) = <A, B>_Frobenius! So trace is related to Frobenious dot product!
  // Therefore, min_R tr(V'KR) is same as min_R <K'V, R>_Frobenius => same as HW02 (finding closest rotation).

  // Assume the initalization of the final vertices to be same as original vertices U
  // C = K'V

  Eigen::MatrixXd C = K.transpose() * U;
  
  // For the rotations
  Eigen::MatrixXd R;
  R.resize(3*U.rows(), 3);

  // Loop in steps of k
  for (int k=0; k<U.rows(); k++) {
    
    // polar_svd3x3 doesnt work with dynamic matrices of 3x3 size
    Eigen::Matrix3d R_k; 
    Eigen::Matrix3d C_k; 
    C_k = C.block(k*3, 0, 3, 3);

    // get the closest rotation
    igl::polar_svd3x3(C_k, R_k);

    // update R
    R.block(k*3, 0, 3, 3) = R_k;
  }

  // global step
  // minimize tr (V'LV + V'KR) w.r.t V
  // Compare this with min_quad_with_solve:
  // B = KR; Z = V; 
  Eigen::MatrixXd B = K * R;

  // dummy
  Eigen::MatrixXd Beq;

  // global step, replace U
  igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);
  
}
