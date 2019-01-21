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
  // REPLACE WITH YOUR CODE
  data.n = V.rows();

  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

  Eigen::SparseMatrix<double> Q = L.transpose() * M.cwiseInverse() * L;

  Eigen::SparseMatrix<double> Aeq;
  igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);
}

