#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <Eigen/Sparse>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  Eigen::SparseMatrix<double> L;
  Eigen::SparseMatrix<double> M;
  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);
  
  // Compute ML'
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(M);
  Eigen::SparseMatrix<double> MinvL = solver.solve(L);
  
  Eigen::SparseMatrix<double> Q = L.transpose() * MinvL;
  
  igl::min_quad_with_fixed_precompute(Q, b, Eigen::SparseMatrix<double>(), false, data);
}

