#include "arap_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/cotmatrix_entries.h>
#include <igl/min_quad_with_fixed.h>
#include <vector>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);
  igl::min_quad_with_fixed_precompute(L, b, Eigen::SparseMatrix<double>(), false, data);

  Eigen::MatrixXd C(F.rows(), 3);
  igl::cotmatrix_entries(V, F, C);

  std::vector<Eigen::Triplet<double>> triplets;

  K.resize(V.rows(), 3 * V.rows());
  for (int f = 0; f < F.rows(); f++) {
    Eigen::RowVector3i face = F.row(f);
    for (int x = 0; x < 3; x++) {
      int i = face((x + 1) % 3);
      int j = face((x + 2) % 3);
      Eigen::RowVector3d edge = V.row(i) - V.row(j);
      Eigen::Vector3d eij = C(f, x) / 3.0 * edge ;

      for (int a = 0; a < 3; a++) {
        int k = face((x + a) % 3);
        for (int b = 0; b < 3; b++) {
          triplets.emplace_back(Eigen::Triplet<double>(i, 3 * k + b, eij(b)));
          triplets.emplace_back(Eigen::Triplet<double>(j, 3 * k + b, -eij(b)));
        }
      }
    }
  }
  K.setFromTriplets(triplets.begin(), triplets.end());

}
