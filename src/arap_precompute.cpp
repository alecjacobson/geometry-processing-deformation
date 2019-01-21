#include "arap_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/cotmatrix_entries.h>

void arap_precompute(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::VectorXi & b,
	igl::min_quad_with_fixed_data<double> & data,
	Eigen::SparseMatrix<double> & K)
{
	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);
	Eigen::SparseMatrix<double> Aeq;
	igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);

	// compute K
	Eigen::MatrixXd C(F.rows(), 3);
	igl::cotmatrix_entries(V, F, C);

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(F.rows() * 54);
	for (int f = 0; f < F.rows(); f++) {
		for (int i = 0; i < 3; i++) {
			// half edge ij
			int j = (i + 1) % 3;
			Eigen::RowVector3d e = C(f, (i + 2) % 3) * (V.row(i) - V.row(j));
			for (int k = 0; k < 3; k++) {
				for (int alpha = 0; alpha < 3; alpha++) {
					for (int beta = 0; beta < 3; beta++) {
						int col = 3 * F(f, k) + beta;
						tripletList.push_back(T(F(f, i), col, e(beta)/6.0));
						tripletList.push_back(T(F(f, j), col, -e(beta)/6.0));
					}
				}
			}
		}
	}
	K.resize(V.rows(), 3 * V.rows());
	K.setFromTriplets(tripletList.begin(), tripletList.end());
}
