#include "arap_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/cotmatrix_entries.h>
#include <igl/min_quad_with_fixed.h>
#include <vector>

void arap_precompute(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::VectorXi & b,
	igl::min_quad_with_fixed_data<double> & data,
	Eigen::SparseMatrix<double> & K)
{

	//Start by doing the precomputation to generate data.

	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);
	L *= (1 / 6.0);

	Eigen::SparseMatrix<double> Aeq(0, 0);
	igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);
	
	//Next, construct K.
	int num_faces = F.rows();

	//Fingers crossed that this is correct. Super unsure.
	std::vector<Eigen::Triplet<double>> entries;
	for (int f = 0; f < num_faces; f++) {
		for (int v = 0; v < 3; v++) {
			int i = F(f, v);
			int j = F(f, (v + 1) % 3);


			Eigen::RowVector3d edge = V.row(i) - V.row(j);
			float c = L.coeff(i, j);

			for (int k = 0; k < 3; k++) {
				int ki = F(v, k);
				entries.emplace_back(i, 3 * ki + 0, c*edge(0));
				entries.emplace_back(i, 3 * ki + 1, c*edge(1));
				entries.emplace_back(i, 3 * ki + 2, c*edge(2));

				entries.emplace_back(j, 3 * ki + 0, -c*edge(0));
				entries.emplace_back(j, 3 * ki + 1, -c*edge(1));
				entries.emplace_back(j, 3 * ki + 2, -c*edge(2));
			}
		}
	}
	K.resize(V.rows(), 3 * V.rows());
	K.setFromTriplets(entries.begin(), entries.end());

	K *= (1 / 6.0);

}
