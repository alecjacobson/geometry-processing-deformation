#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

using namespace Eigen;
using namespace igl;

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  // REPLACE WITH YOUR CODE
  //data.n = V.rows();
	
	SparseMatrix<double> L;
	cotmatrix(V, F, L);

	SparseMatrix<double> M;	
	massmatrix(V, F, MASSMATRIX_TYPE_DEFAULT, M);

	//Computer M inverse
	SimplicialLDLT <SparseMatrix<double> > solver;
	solver.compute(M);
	SparseMatrix<double> I(M.rows(), M.cols());
	I.setIdentity();
	SparseMatrix<double> M_inv = solver.solve(I);

	SparseMatrix<double> Q = L.transpose() * M_inv * L;

	Eigen::SparseMatrix<double> Aeq;
	min_quad_with_fixed_data<double> data;
	min_quad_with_fixed_precompute(Q, b, Aeq, false, data);
}

