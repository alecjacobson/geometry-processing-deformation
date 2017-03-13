#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  // REPLACE WITH YOUR CODE
	int n = V.rows();
	data.n = n;
	
	Eigen::SparseMatrix<double> Aeq;

	Eigen::SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solverM;
	solverM.compute(M);
	Eigen::SparseMatrix<double> I(M.rows(),M.cols());
	I.setIdentity();
	Eigen::SparseMatrix<double> M_inv = solverM.solve(I);

	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);

	Eigen::SparseMatrix<double> Q = L.transpose()*M_inv*L;

	igl::min_quad_with_fixed_precompute(Q, b, Aeq, true, data);
}

