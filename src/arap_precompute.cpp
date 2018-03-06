#include "arap_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix_entries.h>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  // Precompute data
  // Get Laplacian
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V,F,L);
  // Empty constraints
  Eigen::SparseMatrix<double> Aeq;
  // Minimize energey trace( 0.5*U'*L*U + U'*B + constant )
  igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);

  //Precompute K
  K.resize(V.rows(), 3 * V.rows()); // n x 3n
  typedef Eigen::Triplet< double > Triplet;
  std::vector< Triplet > triplets;

  // gather cotangent weights
  Eigen::MatrixXd C;                      // 0     1     2
  igl::cotmatrix_entries(V,F,C); //cols are [1,2],[2,0],[0,1]

  // Loop over half edges
  for(int f_idx = 0; f_idx < F.rows(); f_idx++){
    int i,j,k; // indicies into V
    for(int edge = 0; edge < 3; edge++){
      i = F(f_idx, (1 + edge) % 3);
      j = F(f_idx, (2 + edge) % 3);
      k = F(f_idx, (0 + edge) % 3);

      Eigen::Vector3d v_i_minus_j = V.row(i) - V.row(j);
      // loop over beta
      for(int beta = 0; beta < 3; beta++){
        triplets.push_back(Triplet(i, 3*k + beta,  v_i_minus_j[beta] * C(f_idx,edge)));
        triplets.push_back(Triplet(j, 3*k + beta, -v_i_minus_j[beta] * C(f_idx,edge)));
      }
    }
  }
  K.setFromTriplets(triplets.begin(), triplets.end());
}
