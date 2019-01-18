#include "arap_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>
#include <iostream>
typedef Eigen::Triplet<double> T;

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  // Compute L:
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  // Compute K:
  std::vector<T> triplet_list;
  for(int i = 0; i < F.rows(); i ++) {
    for(int j = 0; j < 3; j ++) {
      // Find vi, and vj:
      int vi = F(i, j);
      int vj = F(i, (j+1)%3);
      Eigen::Vector3d Vi = V.row(vi);
      Eigen::Vector3d Vj = V.row(vj);

      // Compute eij:
      Eigen::Vector3d eij = L.coeff(vi, vj) * (Vi - Vj);

      // K could be three possible vertex:
      for (int k = 0; k < 3; k ++) {
        int vk = F(i, k);
        // Loop over all 3 indexes:
        for(int l = 0; l < 3; l ++) {
          triplet_list.push_back(T(vi, 3*vk+l, eij[l]/6.0));
          triplet_list.push_back(T(vj, 3*vk+l, -eij[l]/6.0));
        }
      }
    }
  }

  K.resize(V.rows(), 3 * V.rows());
  K.setFromTriplets(triplet_list.begin(),triplet_list.end());
  igl::min_quad_with_fixed_precompute(L, b, Eigen::SparseMatrix<double>(), false, data);
}
