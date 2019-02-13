#include "arap_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
	int n = V.rows();
	int m = F.rows();
	int k = b.rows();

	Eigen::SparseMatrix<double> L(n, n);
	igl::cotmatrix(V, F, L);

	Eigen::SparseMatrix<double> Aeq;

	igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(m*3*3*6);

	K.resize(n, 3*n);
	for (int f = 0; f < m; f++)
	{
		auto face = F.row(f);
		for (int e = 0; e < 3; e++)
		{
			int i = face(e % 3);
			int j = face((e+1) % 3);
			int k = face((e+2) % 3);

			Eigen::Vector3d v_i = V.row(i);
			Eigen::Vector3d v_j = V.row(j);
			Eigen::Vector3d e_ij = L.coeff(i,j)*(v_i-v_j);
			for (int beta = 0; beta < 3; beta++)
			{
				tripletList.push_back(T(i, 3 * i + beta, e_ij(beta)/2));
				tripletList.push_back(T(j, 3 * i + beta, -e_ij(beta)/2));
				tripletList.push_back(T(i, 3 * j + beta, e_ij(beta)/2));
				tripletList.push_back(T(j, 3 * j + beta, -e_ij(beta)/2));
				tripletList.push_back(T(i, 3 * k + beta, e_ij(beta)/2));
				tripletList.push_back(T(j, 3 * k + beta, -e_ij(beta)/2));
			}
			
		}
	}

	K.setFromTriplets(tripletList.begin(), tripletList.end());
	K = K / 3.0;
}
