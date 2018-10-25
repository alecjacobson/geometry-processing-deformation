#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <vector>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  // REPLACE WITH YOUR CODE
  data.n = V.rows();

  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);

  // Find the inverse of the Mass matrix. 
  // Since it is a diagonal sparse,
  // manually compute element by element and create M_inv
  Eigen::SparseMatrix<double> M_inv;
  M_inv.resize(M.rows(), M.cols());
  std::vector<Eigen::Triplet<double> > triplets;

  for (int j=0; j<M.rows(); j++) {
    triplets.push_back(Eigen::Triplet<double>(j, j, 1.0/M.coeff(j,j)));
  }

  M_inv.setFromTriplets(triplets.begin(), triplets.end());

  // multiply by 2.0 since the igl library uses a 1/2 in front of the
  // quadratic term
  Eigen::SparseMatrix<double> Q = 2.0 * L.transpose() * M_inv * L;

  Eigen::SparseMatrix<double> Aeq;
  //Aeq.resize(V.rows(), V.rows());

  // Initially planned using Aeq and Beq below. But using Z_b, Z_bc is much simpler.
  // // for all the handle points keep the diagonal entries in Aeq = 1
  // std::vector<Eigen::Triplet<double> > Aeq_triplets;
  // for (int j=0; j<b.size(); j++) {
  //   Aeq_triplets.push_back(Eigen::Triplet<double>(b(j), b(j), 1.0));
  // }
  // Aeq.setFromTriplets(Aeq_triplets.begin(), Aeq_triplets.end());

  igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);
}

