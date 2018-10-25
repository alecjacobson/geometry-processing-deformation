#include "arap_precompute.h"
#include <vector>
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>


void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  // REPLACE WITH YOUR CODE

  // Construct K
  std::vector<Eigen::Triplet<double> > triplets;

  // get the cotmatrix
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  // iterate over faces
  for (int f=0; f<F.rows(); f++) {

    // iterate over 3 edges of a face using modulo
    for (int e=0; e<3; e++) {
      int v_i_idx = F(f, e);
      int v_j_idx = F(f, (e + 1) % 3);

      Eigen::Vector3d v_i = V.row(v_i_idx);
      Eigen::Vector3d v_j = V.row(v_j_idx);

      // v_i -> v_j represents edge e_ij
      Eigen::Vector3d e_ij = (v_i - v_j) * L.coeff(v_i_idx, v_j_idx);

      // loop over 3 dimensions x, y, z (beta)
      for (int beta=0; beta<3; beta++) {
        // index k can be index v_i or v_j or the third vertex
        
        // case 1: k is i 
        int k_idx = v_i_idx; 
        triplets.push_back(Eigen::Triplet<double>(v_i_idx, 3*k_idx + beta, (1.0/6.0)*e_ij(beta)));
        triplets.push_back(Eigen::Triplet<double>(v_j_idx, 3*k_idx + beta, (-1.0/6.0)*e_ij(beta)));

        // case 2: k is j
        k_idx = v_j_idx; 
        triplets.push_back(Eigen::Triplet<double>(v_i_idx, 3*k_idx + beta, (1.0/6.0)*e_ij(beta)));
        triplets.push_back(Eigen::Triplet<double>(v_j_idx, 3*k_idx + beta, (-1.0/6.0)*e_ij(beta)));

        // case 3: k i the third vertex
        k_idx = F(f, (e + 2) % 3);
        triplets.push_back(Eigen::Triplet<double>(v_i_idx, 3*k_idx + beta, (1.0/6.0)*e_ij(beta)));
        triplets.push_back(Eigen::Triplet<double>(v_j_idx, 3*k_idx + beta, (-1.0/6.0)*e_ij(beta)));

      }

    }
  }

  K.resize(V.rows(), 3 * V.rows());
  K.setFromTriplets(triplets.begin(), triplets.end());

  // data matrix precomputed similar to biharmonic_precompute
  // need to minimize_D tr(D'LD + D'KR)

  // not using this 
  Eigen::SparseMatrix<double> Aeq;

  Eigen::SparseMatrix<double> L_times_2 = 2.0 * L;

  // pd should be false else numerical issue as L is not positive definite
  igl::min_quad_with_fixed_precompute(L_times_2, b, Aeq, false, data);

}
