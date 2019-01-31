#include "arap_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix_entries.h>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  // REPLACE WITH YOUR CODE	

	// Compute cotangent laplacian
	typedef Eigen::SparseMatrix<double> SpMat;
	SpMat L;
	igl::cotmatrix(V, F, L);

	// Precompute ARAP minimizer
	SpMat Aeq;
	igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);

	// Compute cotangent of each angle in mesh accessible by C(face_idx, edge_idx)
	Eigen::MatrixXd C;
	igl::cotmatrix_entries(V, F, C);

	// Compute K
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;
	triplets.reserve(F.rows() * 3 * 3 * 6);
	for (int f = 0; f < F.rows(); f++) {
		for (int k = 0; k < 3; k++) { // Loop over three edges
			int i, j;
			switch (k) { // Identify corresponding vertices
			case 0:
				i = F(f, 1);
				j = F(f, 2);
				break;
			case 1:
				i = F(f, 2);
				j = F(f, 0);
				break;
			case 2:
				i = F(f, 0);
				j = F(f, 1);
				break;
			}
			Eigen::Vector3d e = C(f, k)*(V.row(i) - V.row(j));
			for (int beta = 0; beta < 3; beta++) { // Loop over three dimension of vector3d e
				// Push back for each vertex
				triplets.push_back(T(i, 3 * F(f, 0) + beta, e(beta)));
				triplets.push_back(T(i, 3 * F(f, 1) + beta, e(beta)));
				triplets.push_back(T(i, 3 * F(f, 2) + beta, e(beta)));
				triplets.push_back(T(j, 3 * F(f, 0) + beta, -e(beta)));
				triplets.push_back(T(j, 3 * F(f, 1) + beta, -e(beta)));
				triplets.push_back(T(j, 3 * F(f, 2) + beta, -e(beta)));
			}
		}
	}
	K.resize(V.rows(), V.rows() * 3);
	K.setFromTriplets(triplets.begin(), triplets.end());
	K = K / 3.0; // Compensate for adding Kij three times (each time per vertex)
}
