#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>


void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  // REPLACE WITH YOUR CODE
  //data.n = V.rows();

  Eigen::SparseMatrix<double> L, M, Q, M_in, Aeq, I;
  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_VORONOI, M);

  //compute M-1
  I.resize(M.rows(), M.cols());
  I.setIdentity();
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::UpLoType::Lower | Eigen::UpLoType::Upper> solver;
  solver.compute(M);
  M_in = solver.solve(I);

  //Q = Lt*M-1*L
  Q = L.transpose() * M_in * L;

  igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);
}

