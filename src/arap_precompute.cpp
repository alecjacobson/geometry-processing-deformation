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
  // REPLACE WITH YOUR CODE
  Eigen::SparseMatrix<double> A; 
  igl::cotmatrix(V, F, A);
  min_quad_with_fixed_precompute(A, b, Eigen::SparseMatrix<double>(), false, data);

  Eigen::MatrixXd C;
  igl::cotmatrix_entries(V, F, C);

  K.resize(V.rows(), 3 * V.rows());

  for (int f = 0; f < F.rows(); f++) { //each face
  	for (int edge_num = 0; edge_num < 3; edge_num++) { // each edge 
  		Eigen::RowVector3d edge = V.row(F(f, edge_num % 3)) - V.row(F(f, (edge_num + 1) % 3));
  		for (int k = 0; k < 3; k++) { // each vertex
  			for (int beta = 0; beta < 3; beta++) { // xyz
  				K.coeffRef(F(f, edge_num % 3), 3 * F(f, k) + beta) += C(f, k) * edge(beta);
  				K.coeffRef(F(f, (edge_num + 1) % 3), 3 * F(f, k) + beta) -= C(f, k) * edge(beta);
  			}
  		}
	}
  }
}
