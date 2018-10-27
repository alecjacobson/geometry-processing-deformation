#include "arap_precompute.h"
typedef Eigen::Triplet<double> T;

#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>

void arap_precompute(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F,
		const Eigen::VectorXi & b, igl::min_quad_with_fixed_data<double> & data,
		Eigen::SparseMatrix<double> & K) {
	// REPLACE WITH YOUR CODE
	// reconstruct K
	int vnum = V.rows();
	int fnum = F.rows();

	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);

	Eigen::SparseMatrix<double> Aeq;
	igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);

	K.resize(vnum, vnum * 3);
	std::vector<T> coef;

	for (int i = 0; i < fnum; i++) {
		for (int j = 0; j < 3; j++) {
			int v1 = F(i, j);
			int v2 = F(i, (j + 1) % 3);
			double c = L(v1, v2);
			Eigen::VectorXd beta = c * (V.row(v1) - V.row(v2));

			for (int k = 0; k < 3; k++) {
				int v3 = F(i, k);

				for (int b = 0; b < 3; b++) {
					T tmp1(v1, 3 * v3 + b, beta(b));
					T tmp2(v2, 3 * v3 + b, -beta(b));
					coef.push_back(tmp1);
					coef.push_back(tmp2);
				}
			}
		}
	}
	K.setFromTriplets(coef.begin(), coef.end());

}
