#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <iostream>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  // Find L, M
  Eigen::SparseMatrix<double> L, M;
  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

  // Compute M_inverse
  Eigen::SparseMatrix<double> M_inverse(M.rows(), M.rows());
  Eigen::VectorXd M_inverse_dia;
  M_inverse_dia = Eigen::ArrayXd::Ones(M.rows()) / M.diagonal().array();
  for(int i = 0; i < M.rows(); i ++) {
    M_inverse.insert(i, i) = M_inverse_dia[i];
  }

  // Compute Q:
  Eigen::SparseMatrix<double> Q, Aeq;
  Q = L.transpose() * M_inverse * L;

  // Compute data:
  igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);
}

