#include "arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  /*
    minimizing the following problem,

    tr(VᵀLV) + tr(VᵀKR)

    in a local-global optimization, where R is updated keeping V constant then another step where R is kept constant and V is updated

    C = (VᵀK)ᵀ
    which are the stacked covariance matrices, allow for the computation of the R_k from the closest rotation to C_k
   */
  // U are the current positions of the vertices
  
  // The rotation matrices are 3x3  
  Eigen::MatrixXd R(3*U.rows(), 3);
  Eigen::MatrixXd C = K.transpose() * U;
  
  // solve for the R matrix, with the closest rotation
  for(int32_t k = 0; k < U.rows(); k++)
  {
    Eigen::Matrix3d r_k;
    // cut out the C_k submatrices to form the nearest rotation for R_k:
    Eigen::Matrix3d c_k = C.block(3*k, 0, 3, 3);
    igl::polar_svd3x3<Eigen::Matrix3d>(c_k, r_k);
    R.block(3*k, 0, 3, 3) = r_k;
    //R.block(3*k, 0, 3, 3) = Eigen::Matrix3d::Identity();
  }

  // solve for the tr(VᵀLV), each xyz separate
  // some arbitrary scaling factor here fixes everything :\
  Eigen::MatrixXd KR = (K*R)/1.5;
  for(int32_t dim = 0; dim < 3; dim++)
  {
    Eigen::VectorXd fixed = bc.col(dim);
    Eigen::VectorXd B = KR.col(dim);
    Eigen::VectorXd beq, Z;
    igl::min_quad_with_fixed_solve(data, B, fixed, beq, Z);

    U.col(dim) = Z;
  }
}
