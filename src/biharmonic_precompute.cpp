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
  // REPLACE WITH YOUR CODE
	// Compute cotangent laplacian for the mesh
	typedef Eigen::SparseMatrix<double> SpMat;
	SpMat L;
	igl::cotmatrix(V, F, L);

	// Compute mass matrix for the mesh
	SpMat M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

	// Compute M^-1
	// Method 1: since M is a diagonal matrix, M^-1 = (1.0/M.diagonal()).asdiagonal
	Eigen::VectorXd diagonal = 1.0/(M.diagonal().array());
	SpMat M_inv;
	M_inv.resize(M.rows(), M.cols());
	for (int i = 0; i < M_inv.rows(); i++) {
		M_inv.insert(i, i) = diagonal[i];
	}
	// Method 2: solve M^1-M = I (didn't see much difference in performance with small meshes but it doesn't seem to
	// work with large meshes like knight.off, don't know why...)
	/*
	SpMat I(M.rows(), M.cols());
	I.setIdentity();
	Eigen::SimplicialCholesky<SpMat> chol(M);
	SpMat M_inv = chol.solve(I);
	*/
	
	// Compute Q = L^TM^-1L
	SpMat Q = L.transpose()*M_inv*L;

	// Precompute for biharmonic deformation
	SpMat Aeq;
	igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);
}

