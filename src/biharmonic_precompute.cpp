#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

using namespace Eigen;

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  SparseMatrix<double> L, M;
  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

  // compute M inverse: M is a diagonal matrix, easy to invert
  SparseMatrix<double> M_inv(V.rows(), V.rows());
  for (int i = 0; i < V.rows(); i++) {
    M_inv.insert(i, i) = 1.0 / M.coeff(i, i);
  }

  // compute Q
  SparseMatrix<double> Q, Aeq;
  Q = L.transpose() * M_inv * L;

  igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);
}
