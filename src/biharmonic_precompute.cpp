#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>

using namespace Eigen;

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
	SparseMatrix<double> M, L;

	int n = V.rows(), f = F.rows();

	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
	igl::cotmatrix(V, F, L);

	SparseMatrix<double> M_inv(n, n);
	std::vector <Triplet<double>> M_inv_val;
	M_inv_val.reserve(n);

	for (int i = 0; i < M.rows(); ++i)
		M_inv_val.push_back({ i, i, 1 / M.coeff(i, i) });

	M_inv.setFromTriplets(M_inv_val.begin(), M_inv_val.end());

	Eigen::SparseMatrix<double> Q(L.transpose()*M_inv*L);

	Q.eval();

	igl::min_quad_with_fixed_precompute(Q, b, SparseMatrix<double>(0, 0), true, data);
}

