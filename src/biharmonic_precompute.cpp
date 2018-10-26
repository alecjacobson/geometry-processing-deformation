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

  // Compute L
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  // Compute M
  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MassMatrixType(), M);

  // Invert M (diagonal matrix)
  for (int i = 0; i < M.outerSize(); ++i)
    for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it) {
      it.valueRef() = 1.0/it.value();
    }

  // Compute Q
  Eigen::SparseMatrix<double> Q = L.transpose()*M*L;

  bool success = 
    igl::min_quad_with_fixed_precompute(Q, b, Eigen::SparseMatrix<double>(), true, data);

  if (!success)
    throw std::runtime_error("[biharmonic_precompute] Precompute returned false.");

}

