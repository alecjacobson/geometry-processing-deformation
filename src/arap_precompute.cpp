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

	/*
	Eigen::SparseMatrix<double> M(n, n);
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);

	Eigen::SparseMatrix<double> Q = L.transpose() * M.cwiseInverse() * L;
	*/

	Eigen::SparseMatrix<double> Aeq;

	igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(n);

	K.resize(n, 3*n);
	for (int f = 0; f < m; f++)
	{
		auto face = F.row(f);
		for (int k = 0; k < 3; k++)
		{
			int i = face(k);
			int j = 0;
			if (k + 1 < 3)
				j = face(k + 1);
			else
				j = face(0);

			Eigen::Vector3d v_i = V.row(i);
			Eigen::Vector3d v_j = V.row(j);
			Eigen::Vector3d e_ij = L.coeff(i,j)*(v_i-v_j);
			for (int beta = 0; beta < 3; beta++)
			{
				tripletList.push_back(T(i, 3*k + beta, e_ij(beta)));
				tripletList.push_back(T(j, 3*k + beta, -e_ij(beta)));
			}
		}
	}

	K.setFromTriplets(tripletList.begin(), tripletList.end());
}
