#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

void biharmonic_precompute(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::VectorXi & b,
	igl::min_quad_with_fixed_data<double> & data)
{
	// get mass matrix and laplacian matrix
	Eigen::SparseMatrix<double> M, L;
	igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M); // default ?
	igl::cotmatrix(V, F, L);
	
	// Q = L^T * M^{-1} * L
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	solver.compute(M);
	Eigen::SparseMatrix<double> ML = solver.solve(L);
	Eigen::SparseMatrix<double> Q = L.transpose() * ML;

	Eigen::SparseMatrix<double> Aeq;
	// default is to have 0 linear equality constraints - note from min_quad_with_fixed.cpp
	igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);
}

