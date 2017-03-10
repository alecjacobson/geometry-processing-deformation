#include "arap_precompute.h"
#include "igl/cotmatrix_entries.h"
#include "igl/cotmatrix.h"
#include <vector>
#include <igl/min_quad_with_fixed.h>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
 
	int n = V.rows();
	int n_f = F.rows();

	// Precompute min_quad_with_fixed data
	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);
	L = -L;
	Eigen::SparseMatrix<double> Aeq;
	igl::min_quad_with_fixed_precompute(L, b, Aeq, true, data);

	//Compute K
	K = Eigen::SparseMatrix<double>(n, 3 * n);
	Eigen::MatrixXd CL;
	igl::cotmatrix_entries(V, F, CL);
	Eigen::VectorXd v0;
	Eigen::VectorXd v1;
	Eigen::VectorXd v2;
	double c12 = 0;
	double c20 = 0;
	double c01 = 0;
	int i = 0;
	int j = 0;
	int l = 0;
	Eigen::VectorXd row_jl;
	Eigen::VectorXd row_li;
	Eigen::VectorXd row_ij;
	Eigen::VectorXd row;
	std::vector<Eigen::Triplet<double>> triplets;
	for (int k = 0; k < n_f; ++k) {
		i = F(k,0);
		j = F(k,1);
		l = F(k,2);

		v0 = V.row(i);
		v1 = V.row(j);
		v2 = V.row(l);

		c12 = CL(k,0);
		c20 = CL(k,1);
		c01 = CL(k,2);

		row_jl = c12*(v1 - v2);
		row_li = c20*(v2 - v0);
		row_ij = c01*(v0 - v1);

		row = row_jl + row_li + row_ij;

		triplets.push_back(Eigen::Triplet<double>(i, 3 * i, row(0)));
		triplets.push_back(Eigen::Triplet<double>(i, 3 * i + 1, row(1)));
		triplets.push_back(Eigen::Triplet<double>(i, 3 * i + 2, row(2)));

		triplets.push_back(Eigen::Triplet<double>(j, 3 * j, row(0)));
		triplets.push_back(Eigen::Triplet<double>(j, 3 * j + 1, row(1)));
		triplets.push_back(Eigen::Triplet<double>(j, 3 * j + 2, row(2)));

		triplets.push_back(Eigen::Triplet<double>(l, 3 * l, row(0)));
		triplets.push_back(Eigen::Triplet<double>(l, 3 * l + 1, row(1)));
		triplets.push_back(Eigen::Triplet<double>(l, 3 * l + 2, row(2)));

	}
	K.setFromTriplets(triplets.begin(), triplets.end());
}
