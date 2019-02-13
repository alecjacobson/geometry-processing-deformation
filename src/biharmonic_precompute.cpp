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
  int n = V.rows();
  int k = b.rows();

  Eigen::SparseMatrix<double> laplacian(n, n);
  igl::cotmatrix(V, F, laplacian);
 
  Eigen::SparseMatrix<double> mass(n, n);
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, mass);

  Eigen::SparseMatrix<double> Q = laplacian.transpose() * mass.cwiseInverse() * laplacian;
  
  Eigen::SparseMatrix<double> Aeq;

  igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);
}

