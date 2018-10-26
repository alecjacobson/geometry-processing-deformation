#include "arap_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/cotmatrix_entries.h>
#include <igl/massmatrix.h>
#include <iostream>

using namespace std;

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  // REPLACE WITH YOUR CODE
	// Compute L
	Eigen::SparseMatrix<double> L(V.rows(), V.rows());
	igl::cotmatrix(V, F, L);
	// Precompute for the global step
	Eigen::SparseMatrix<double> Aeq;
	igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);
	
	// Prepare matrix K
	K = Eigen::SparseMatrix<double>(V.rows(), 3 * V.rows());
	K.setZero();
	Eigen::MatrixXd C(F.rows(), 3);
	igl::cotmatrix_entries(V, F, C);
	std::vector<Eigen::Triplet<double>> triplets;
	for (int f = 0; f < F.rows(); f++) {
	// Looping over faces
		int i, j;
		// Looping over edges of the f-th face
		i = 1; j = 2;
		for (int k = 0; k < 3; k++) {
			// Looping over vertices
			for (int beta = 0; beta < 3; beta++) {
				triplets.push_back({ F(f,i), 3 * F(f,k) + beta,
					 1./6.*(V(F(f,i), beta) - V(F(f,j), beta))*L.coeff(F(f, i), F(f, j) ) });
				triplets.push_back({ F(f,j), 3 * F(f,k) + beta,
					-1./6.*(V(F(f,i), beta) - V(F(f,j), beta))*L.coeff(F(f, i), F(f, j) ) });
			}
		}
		i = 2; j = 0;
		for (int k = 0; k < 3; k++) {
			for (int beta = 0; beta < 3; beta++) {
				triplets.push_back({ F(f,i), 3 * F(f,k) + beta,
					 1./6.*(V(F(f,i), beta) - V(F(f,j), beta))*L.coeff(F(f, i), F(f, j)) });
				triplets.push_back({ F(f,j), 3 * F(f,k) + beta,
					-1./6.*(V(F(f,i), beta) - V(F(f,j), beta))*L.coeff(F(f, i), F(f, j)) });
			}
		}
		i = 0; j = 1;
		for (int k = 0; k < 3; k++) {
			for (int beta = 0; beta < 3; beta++) {
				triplets.push_back({ F(f,i), 3 * F(f,k) + beta,
					 1./6.*(V(F(f,i), beta) - V(F(f,j), beta))*L.coeff(F(f, i), F(f, j)) });
				triplets.push_back({ F(f,j), 3 * F(f,k) + beta,
					-1./6.*(V(F(f,i), beta) - V(F(f,j), beta))*L.coeff(F(f, i), F(f, j)) });
			}
		}
	}
	K.setFromTriplets(triplets.begin(), triplets.end());
}
