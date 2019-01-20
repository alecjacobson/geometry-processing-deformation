#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  // Mass and cot matrix
  Eigen::SparseMatrix<double> M, L;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
  igl::cotmatrix(V, F, L);

  // Q = L.transpose() * M.inv() * L
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(M);

  Eigen::SparseMatrix<double> tmp = solver.solve(L).sparseView();
  // had to add sparseView here

  Eigen::SparseMatrix<double> Q = L.transpose() * tmp;

  //From source code for min_quad_with_fixed_precompute
  Eigen::SparseMatrix<double> Aeq;
  igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);
}