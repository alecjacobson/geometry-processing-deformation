#include "arap_precompute.h"
#include "igl/cotmatrix.h"
#include <igl/min_quad_with_fixed.h>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
	int N = V.rows();

	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);
	Eigen::MatrixXd cot(F.rows(), 3);
	igl::cotmatrix_entries(V, F, cot);

	K.resize(N, 3 * N);
	std::vector<Eigen::Triplet<double>> triplets;	

	for (int f = 0; f < F.rows(); f++)
	{
		for (int e = 0; e < 3; e++)
		{
			int k = F(f, e);

			for (int d = 0; d < 3;d++) {
				triplets.emplace_back(F(f, 1), 3 * k + d, (cot(f, 0) / 3.0) * (V(F(f,1), d) - V(F(f,2), d)));
				triplets.emplace_back(F(f, 2), 3 * k + d, (cot(f, 0) / 3.0) * (V(F(f, 2), d) - V(F(f, 1), d)));
			}

			for (int d = 0; d < 3;d++) {
				triplets.emplace_back(F(f, 2), 3 * k + d, (cot(f, 1) / 3.0) * (V(F(f, 2), d) - V(F(f, 0), d)));
				triplets.emplace_back(F(f, 0), 3 * k + d, (cot(f, 1) / 3.0) * (V(F(f, 0), d) - V(F(f, 2), d)));
			}

			for (int d = 0; d < 3;d++) {
				triplets.emplace_back(F(f, 0), 3 * k + d, (cot(f, 2) / 3.0) * (V(F(f, 0), d) - V(F(f, 1), d)));
				triplets.emplace_back(F(f, 1), 3 * k + d, (cot(f, 2) / 3.0) * (V(F(f, 1), d) - V(F(f, 0), d)));
			}						
		}
	}

	K.setFromTriplets(triplets.begin(), triplets.end());
		
	igl::min_quad_with_fixed_precompute(L, b, Eigen::SparseMatrix<double>(), false, data);
}
