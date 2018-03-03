#include "biharmonic_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/min_quad_with_fixed.h>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  // REPLACE WITH YOUR CODE
  //data.n = V.rows();

  // Get Laplacian
 Eigen::SparseMatrix<double> L;
 igl::cotmatrix(V,F,L);
 // Get Mass mastrix
 Eigen::SparseMatrix<double> M;
 igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
 // Get inverse of mass Matrix
 Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver;
 solver.compute(M);
 Eigen::SparseMatrix<double> I(M.rows(),M.rows());
 I.setIdentity();
 auto M_inv = solver.solve(I);

 // Get Q
 Eigen::MatrixXd Q;
 Q = L.transpose() * M_inv * L;
 Eigen::SparseMatrix<double> Q_sparse;
 Q_sparse = Q.sparseView();

 // Empty constraints
 Eigen::SparseMatrix<double> Aeq;

 // Minimize energey trace( 0.5*D'*Q*D + U'*B + constant )
 igl::min_quad_with_fixed_precompute(Q_sparse, b, Aeq, false, data);
}
