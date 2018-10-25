#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  Eigen::SparseMatrix<double> L, M;
  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

  Eigen::SparseMatrix<double> M_i(M.rows(), M.cols());
  for (int i = 0; i < M.rows(); i++) {
    M_i.coeffRef(i, i) = 1.0 / M.coeff(i, i);
  }

  Eigen::SparseMatrix<double> Q = L.transpose() * M_i * L;
  igl::min_quad_with_fixed_precompute(Q, b, Eigen::SparseMatrix<double>(), false, data);
}

