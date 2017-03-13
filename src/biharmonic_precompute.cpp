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
	//Calculate Mass Matrix and the Cotangent Laplacian
	Eigen::SparseMatrix<double> M, L;
	igl::cotmatrix(V, F, L);
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
	
	//Mass Matrix is diagonal, so just invert the diagonal entries
	Eigen::SparseMatrix<double> M_inv(M.rows(), M.cols());
	for (int i = 0; i < M.rows();i++) {
		M_inv.coeffRef(i, i) = 1.0f / M.coeff(i, i);
	}

	//Calculate Q = L' x M^-1 x L
	Eigen::SparseMatrix<double> Q = L.transpose() * M_inv * L;
	igl::min_quad_with_fixed_precompute(Q, b, Eigen::SparseMatrix<double>(), false, data);
}

