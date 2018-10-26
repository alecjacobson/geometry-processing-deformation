#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include<Eigen/SparseCholesky>
#include <iostream>

using namespace Eigen;
using namespace std;

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  // REPLACE WITH YOUR CODE
  data.n = V.rows();
  Eigen::SparseMatrix<double> L(V.rows(), V.cols());
  igl::cotmatrix(V, F, L);
  Eigen::SparseMatrix<double> M(V.rows(), V.cols());
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

  // Calculate M^(-1) * L
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(M);
  Eigen::SparseMatrix<double> Q1 = solver.solve(L);

  // Calculate L^T * M^(-1) * L
  Eigen::SparseMatrix<double> Q = L.transpose() * Q1;

  // Precompute
  Eigen::SparseMatrix<double> B(Q.rows(), Q.cols());
  B.setZero();
  Eigen::SparseMatrix<double> Aeq;
  igl::min_quad_with_fixed_precompute(Q, b, Aeq, true, data);

}
