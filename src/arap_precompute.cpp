#include "arap_precompute.h"
#include <igl/all_edges.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  // REPLACE WITH YOUR CODE

  // Construct Laplacian Matrix
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  // Compute the matrix K 
  typedef Eigen::Triplet<double> T;
  std::vector<T> triplets;
  triplets.reserve(27 * V.rows());
  for (int i = 0; i < F.rows(); i++) {
	Eigen::Vector3d weights_1 = L.coeff(F(i,0), F(i,1)) * (V.row(F(i,0)) - V.row(F(i,1))); 
	for (int j = 0; j <= 2; j++) {
		for (int k = 0; k <= 2; k++) {
			triplets.push_back(T(F(i,0), 3 * F(i,j) + k, weights_1(k))); 	
			triplets.push_back(T(F(i,1), 3 * F(i,j) + k, -weights_1(k))); 	
 		}
	}

	Eigen::Vector3d weights_2 = L.coeff(F(i,1), F(i,2)) * (V.row(F(i,1)) - V.row(F(i,2)));
	for (int j = 0; j <= 2; j++) {
		for (int k = 0; k <= 2; k++) {
			triplets.push_back(T(F(i,1), 3 * F(i,j) + k, weights_2(k))); 	
			triplets.push_back(T(F(i,2), 3 * F(i,j) + k, -weights_2(k))); 	
 		}
	}

	Eigen::Vector3d weights_3 = L.coeff(F(i,2), F(i,0)) * (V.row(F(i,2)) - V.row(F(i,0)));
	for (int j = 0; j <= 2; j++) {
		for (int k = 0; k <= 2; k++) {
			triplets.push_back(T(F(i,2), 3 * F(i,j) + k, weights_3(k))); 	
			triplets.push_back(T(F(i,0), 3 * F(i,j) + k, -weights_3(k))); 	
 		}
	}
  }
  K.resize(V.rows(), 3 * V.rows());
  K.setFromTriplets(triplets.begin(), triplets.end()); 

  Eigen::SparseMatrix<double> A;
  A.setZero();
  igl::min_quad_with_fixed_precompute(L, b, A, false, data);
}
