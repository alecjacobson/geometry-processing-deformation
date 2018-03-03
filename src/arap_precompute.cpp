#include "arap_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix_entries.h>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  // Precompute data

  // Get Laplacian
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V,F,L);
  // Empty constraints
  Eigen::SparseMatrix<double> Aeq;
  // Minimize energey trace( 0.5*U'*L*U + U'*B + constant )
  igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);

  // TODO: Precompute K
}
