#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

void biharmonic_precompute(
		const Eigen::MatrixXd & V,
		const Eigen::MatrixXi & F,
		const Eigen::VectorXi & b,
		igl::min_quad_with_fixed_data<double> & data) {

	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);

	Eigen::SparseMatrix<double> massMatrix;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, massMatrix);

	// Can't invert the mass matrix, but with some linear algebra, we can
	// express this as a linear system and solve. Works out because the mass
	// matrix is symmetric.
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	solver.compute(massMatrix);
	Eigen::SparseMatrix<double> x = solver.solve(L);

	// Then Q = transpose(L) * invert(mass) * L = transpose(x) * L
	Eigen::SparseMatrix<double> Q = x.transpose() * L;

	// More or less like in the smoothing assignment...
	Eigen::SparseMatrix<double> Aeq(1, V.rows());
	igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);
}

