#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data) {
    
    // construct the cotangent Laplacian.
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V,F,L);
    
    // construct the mass matrix.
    Eigen::SparseMatrix<double> M;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
    
    // construct Q. Because the mass matrix is diagonal, compute its inverse by
    // simply inverting each of the entries along the diagonal.
    Eigen::SparseMatrix<double> Q;
    Q = L.transpose()*(M.cwiseInverse())*L;
    
    // populate data.
    Eigen::SparseMatrix<double> Aeq;
    igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);

}

