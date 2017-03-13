#include "arap_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/cotmatrix_entries.h>
#include <igl/min_quad_with_fixed.h>

using namespace Eigen;

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
	int n = V.rows(), f = F.rows();

	SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);
	L = -L;

	igl::min_quad_with_fixed_precompute(L, b, SparseMatrix<double>(0, 0), true, data);

	MatrixXd C(f, 3);
	igl::cotmatrix_entries(V, F, C);

	K.resize(n, 3*n);
	std::vector<Triplet<double>> K_val;

	K_val.reserve(54 * f);

	for (int face = 0; face < f; ++face)
	{
		for (int edge = 0; edge < 3; ++edge)
		{
			auto vi = F(face, (edge + 1) % 3);
			auto vj = F(face, (edge + 2) % 3);

			for (int k = 0; k < 3; ++k)
			{
				auto vk = F(face, k);

				for (int dim = 0; dim < 3; ++dim)
				{
					K_val.push_back({ vi, 3 * vk + dim, C(face, edge) * (V(vi, dim) - V(vj, dim)) });
					K_val.push_back({ vj, 3 * vk + dim, -C(face, edge) * (V(vi, dim) - V(vj, dim)) });
				}
			}
		}
	}

	K.setFromTriplets(K_val.begin(), K_val.end());
}
