#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <iostream>

using namespace std;
void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  // REPLACE WITH YOUR CODE
  // data.n = V.rows();

  // My code
  Eigen::SparseMatrix<double> L; // laplacian
  igl::cotmatrix(V,F,L);

  Eigen::SparseMatrix<double> invM; // inverse mass matrix
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,invM);
  for (int ii = 0; ii < invM.rows(); ii++){
  	invM.coeffRef(ii,ii) = 1.0 / invM.coeffRef(ii,ii);
  }

  Eigen::SparseMatrix<double> Q; // bi-laplacian
  Q = Eigen::SparseMatrix<double>(L.transpose()) * invM * L;

  Eigen::SparseMatrix<double> Aeq; // empty matrix
  igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);
}

