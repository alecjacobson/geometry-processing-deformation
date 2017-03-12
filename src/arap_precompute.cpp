#include "arap_precompute.h"

#include <igl/cotmatrix.h>
#include <igl/cotmatrix_entries.h>
#include <igl/min_quad_with_fixed.h>

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

	Eigen::SparseMatrix<double> Aeq(0, 0);
	igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);
	
	//Next, construct K.
	Eigen::MatrixX3d cotangents;
	igl::cotmatrix_entries(V, F, cotangents);
	int num_faces = F.rows();
	for (int f = 0; f < num_faces; f++) {
		for (int v = 0; v < 3; v++) {
			Eigen::RowVector3f edge = V.row(F(f, v)) - V.row(F(f, (v + 1) % 3));


		}
	}

}
