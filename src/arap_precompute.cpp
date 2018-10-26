#include "arap_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/cotmatrix_entries.h>
#include <igl/massmatrix.h>
#include <igl/min_quad_with_fixed.h>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  // REPLACE WITH YOUR CODE
  typedef Eigen::Triplet<double> T;

  // Compute cotangents
  // #F by 3 list of 1/2*cotangents corresponding angles
  // for triangles, columns correspond to edges [1,2],[2,0],[0,1]
  Eigen::MatrixXd C;
  igl::cotmatrix_entries(V, F, C);

  // Construct K
  std::vector<T> tripletList;

  // loop through each face (n)
  for (int n = 0; n < F.rows(); ++n) {
    // reference to current face
    const Eigen::RowVector3i& face = F.row(n);

    // loop through each half-edge ij
    for (int ij = 0; ij < 3; ++ij) {
      int i = face(ij);
      int j = face((ij + 1)%3);
      int l = (ij + 2)%3;

      Eigen::RowVector3d e_ij = C(n, l)*(V.row(i) - V.row(j))/3.0;

      // for all possible k
      for (int a = 0; a < 3; ++a) {
        int k = face(a);

        // loop through beta
        for (int beta = 0; beta < 3; ++beta) {
          tripletList.push_back(T(i, 3*k + beta, e_ij(beta)));
          tripletList.push_back(T(j, 3*k + beta, -e_ij(beta)));
        } // end loop beta
      } // end loop a (k)
    } // end loop v
  } // end loop n

  K.resize(V.rows(), V.rows()*3);
  K.setFromTriplets(tripletList.begin(), tripletList.end());

  // Compute L
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  bool success = 
    igl::min_quad_with_fixed_precompute(L, b, Eigen::SparseMatrix<double>(), false, data);

  if (!success)
    throw std::runtime_error("[arap_precompute] Precompute returned false.");

}
