#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <iostream>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
	Eigen::SparseMatrix<double> L; 
	igl::cotmatrix(V, F, L);
	Eigen::SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
	Eigen::SparseMatrix<double> A = L.transpose() * M.cwiseInverse() * L; 
	min_quad_with_fixed_precompute(A, b, Eigen::SparseMatrix<double>(), false, data);
}

