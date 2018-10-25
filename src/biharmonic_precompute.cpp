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
  // data.n = V.rows();
  Eigen::SparseMatrix<double> L,M,A;
  igl::cotmatrix(V,F,L);
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);
  A = L.transpose() * M.cwiseInverse() * L;
  igl::min_quad_with_fixed_precompute(A,b,Eigen::SparseMatrix<double>(),false,data);
}

