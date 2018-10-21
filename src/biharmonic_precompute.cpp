#include "biharmonic_precompute.h"
#include <Eigen/SparseCholesky>
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

  // Construct mass and Laplacian Matrices
  Eigen::SparseMatrix<double> L;
  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_BARYCENTRIC, M);
  igl::cotmatrix(V, F, L);
  
  // Compute M^{-1}
  Eigen::SparseMatrix<double> I(V.rows(), V.rows());
  I.setIdentity();
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(M);
  Eigen::SparseMatrix<double> inv_M = ldlt.solve(M); 

  // Calculate Q = L^T * M^{-1} * L
  Eigen::SparseMatrix<double> Q = L.transpose() * inv_M * L; 

  // Precompute the matrix
  Eigen::SparseMatrix<double> A;
  A.setZero(); 
  igl::min_quad_with_fixed_precompute(Q, b, A, false, data);
}

