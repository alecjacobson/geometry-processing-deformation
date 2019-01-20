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
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(M);
  Eigen::SparseMatrix<double> M_inv_L = solver.solve(L).sparseView();
  Eigen::SparseMatrix<double> Q = L.transpose() * M_inv_L; // L^T M^-1 L

  igl::min_quad_with_fixed_precompute(Q, b, Eigen::SparseMatrix<double>(), false, data);
}

