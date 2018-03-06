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
  // // REPLACE WITH YOUR CODE
  // data.n = V.rows();

	// We want to factorize the the matrix that corresponds to Q in the readme. This requires
	// first computing the cotangent laplacian and mass matrices.

	// Get cotangent laplacian
	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);

	// Get mass matrix
	Eigen::SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

	// Compute the inverse of the mass matrix. Since we're using the diagonalized version,
	// all we have to do is compute reciprocals along the diagonal.
	Eigen::SparseMatrix<double> M_inv(M.cols(), M.cols());
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;

	for (int ii = 0; ii < M.cols(); ii++)
		triplets.push_back(T(ii, ii, 1.0/M.coeff(ii, ii)));

	M_inv.setFromTriplets(triplets.begin(), triplets.end());

	// Compute Q
	Eigen::SparseMatrix<double> Q (L.cols(), L.cols());
	Q = L.transpose()*M_inv*L;
	// std::cout << "Here" << std::endl;

	// Factorization

	// Create empty matrices for the equality constraints, of which we have none (because the known vertices are separate inputs)
	Eigen::SparseMatrix<double> Aeq;

	// Factorize
	igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);

	return;

}

