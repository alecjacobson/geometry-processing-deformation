#include "arap_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix_entries.h>
#include <iostream>

using namespace std;

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  int num_points = V.rows();
  int num_faces = F.rows();

  K.resize(num_points, num_points * 3);

  // Used in the first term of the Global step
  Eigen::SparseMatrix<double> Laplacian;
  igl::cotmatrix(V, F, Laplacian);

  // Update Data with pre computation of the Cotangent Laplacian
  Eigen::SparseMatrix<double> Aeq;
  igl::min_quad_with_fixed_precompute(Laplacian, b, Aeq, false, data);

  // Shape of C = (Num Faces, 3)
  // For each face, it contains 1/2 * cotangent weights (from the Cotangent Laplacian)
  // This will be used to construct K
  Eigen::MatrixXd C;
  igl::cotmatrix_entries(V, F, C);

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  tripletList.reserve(num_faces * 3 * 3 * 2);

  // For each face
  for (int f = 0; f < num_faces; f++) 
  {

    // For each half-edge (i, j) that opposite_vertex is opposite to
    for (int opposite_vertex = 0; opposite_vertex < 3; opposite_vertex++) 
    {
      int i = (opposite_vertex + 1) % 3;
      int j = (opposite_vertex + 2) % 3;
      int V_i = F(f, i);
      int V_j = F(f, j);

      // Weighted edge difference vector
      // c_ij = cotangent of the angle, opposite to the half-edge ij
      Eigen::RowVector3d e_ij = C(f, opposite_vertex) * (V.row(V_i) - V.row(V_j)) / 3.0; 

      // sum is over all rotations k, such that the half-edge ij 
      // belongs to the half-edges of the faces incident on the k-th vertex
      // that is: k| i j ∈ F(k)
      // Well, this means k can either be i or j, or the third vertex ("k") of the f-th face

      // Basically... for every vertex of the face, we have a valid k
      for (int k = 0; k < 3; k++)
      { 
        int V_k = F(f, k);

        // For each component Beta of the edge difference vector
        for (int beta = 0; beta < 3; beta++)
        { 
          tripletList.push_back(T(V_i, 3 * V_k + beta, e_ij(beta)));        
          tripletList.push_back(T(V_j, 3 * V_k + beta, -1.0 * e_ij(beta)));
        }
      }
    }
  }
  K.setFromTriplets(tripletList.begin(), tripletList.end());
}
