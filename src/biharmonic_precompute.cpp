#include "biharmonic_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/min_quad_with_fixed.h>
#include <iostream>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  int num_v = V.rows();
  Eigen::SparseMatrix<double> L(num_v, num_v);
  Eigen::SparseMatrix<double> M(num_v, num_v);
  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);
  
  // Construct Q
  Eigen::SparseMatrix<double> M_inv;
  igl::invert_diag(M, M_inv);
  Eigen::SparseMatrix<double> Q = L.transpose() * M_inv * L;
  
  /* Setup inputs to min_quad:
   * Q = A
   * b = known
   * Aeq = empty matrix
   * pd = false
   */
  Eigen::SparseMatrix<double> Aeq;
  igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);
}

