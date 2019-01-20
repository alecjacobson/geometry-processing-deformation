#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>

#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

typedef Eigen::Triplet<double> T;

void biharmonic_precompute(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F,
		const Eigen::VectorXi & b,
		igl::min_quad_with_fixed_data<double> & data) {
	// REPLACE WITH YOUR CODE
	data.n = V.rows();
	int vnum = data.n;

	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);

	Eigen::SparseMatrix<double> M;
	massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);

	Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
	solver.compute(M);
	Eigen::SparseMatrix<double> I(vnum, vnum);
	I.setIdentity();
	auto M_inv = solver.solve(I);

	// Q
	Eigen::SparseMatrix<double> Q = L.transpose() * M_inv * L;

	// build know function
	// Aeq * D = Beq
	Eigen::SparseMatrix<double> Aeq;
	std::vector<T> coef;
	for (int p = 0; p < b.rows(); p++) {
		T tmp(p, b(p), 1.0);
		coef.push_back(tmp);
	}
	Aeq.setFromTriplets(coef.begin(), coef.end());

	min_quad_with_fixed_precompute(Q, b, Aeq, true, data);
}

