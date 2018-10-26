#include "arap_precompute.h"
#include <vector>
#include <igl/cotmatrix_entries.h>
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  // prepare K
  int i, j, k;
  double c_ij;
  Eigen::RowVectorXd e_ij;
  Eigen::MatrixXd half_cotangents;
  igl::cotmatrix_entries(V, F, half_cotangents);
  std::vector<Eigen::Triplet<double>> triplets;

  for (int face_i = 0; face_i < F.rows(); face_i++) { // faces
    for (int edge_i = 0; edge_i < 3; edge_i++) { // half-edges

      i = F(face_i, (edge_i + 1) % 3);
      j = F(face_i, (edge_i + 2) % 3);
      c_ij = half_cotangents(face_i, edge_i);
      e_ij = c_ij * (V.row(i) - V.row(j)) / 3.0;

      for(int vertex_i = 0; vertex_i < 3; vertex_i++) { // loop over each vertex on the face

        k = F(face_i, (edge_i + vertex_i) % 3);

        for (int beta = 0; beta < 3; beta++) { // beta is the component (x, y, or z)

          triplets.push_back(Eigen::Triplet<double>(i, 3 * k + beta, e_ij(beta)));
          triplets.push_back(Eigen::Triplet<double>(j, 3 * k + beta, -e_ij(beta)));

        }
      }
    }
  }

  K.resize(V.rows(), 3 * V.rows());
  K.setFromTriplets(triplets.begin(), triplets.end());

  // prepare L, precompute data
  Eigen::SparseMatrix<double> L, Aeq;
  igl::cotmatrix(V, F, L);
  igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);

}
