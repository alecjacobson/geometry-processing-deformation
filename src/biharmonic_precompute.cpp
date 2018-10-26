#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <vector>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{

  Eigen::SparseMatrix<double> L, M, Q, Aeq; // just leave Aeq empty...?

  // find L and M
  igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);
  igl::cotmatrix(V, F, L);

  // invert M (diagonal matrix, EZ!)
  Eigen::SparseMatrix<double> M_inverse(M.rows(), M.rows());
  std::vector<Eigen::Triplet<double>> inverted_coeff_triplets;
  for (int i = 0; i < M.rows(); i++) {
    inverted_coeff_triplets.push_back(Eigen::Triplet<double>(i, i, 1.0 / M.coeff(i, i)));
  }
  M_inverse.setFromTriplets(inverted_coeff_triplets.begin(), inverted_coeff_triplets.end());

  // compute Q
  Q = L.transpose() * M_inverse * L;

  // now we can precompute
  igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);
}

