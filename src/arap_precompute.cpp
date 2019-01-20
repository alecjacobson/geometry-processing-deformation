#include "arap_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>
typedef Eigen::Triplet<double> T;

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

  int nV = V.rows();
  std::vector<T> tripletList;

  for(int f = 0; f < F.rows(); f++) {
    for(int a = 0; a < 3; a++) {
      int i = F(f, (a + 1) % 3); 
      int j = F(f, (a + 2) % 3);

      Eigen::Vector3d diff = V.row(i) - V.row(j);
      Eigen::Vector3d eij = L.coeff(i, j) * diff / 6.0;

      for(int ki = 0; ki < 3; ki++) {
        int k = F(f, (a + ki) % 3);
        for(int B = 0; B < 3; B++) {
          tripletList.push_back(T(i, 3 * k + B, eij[B])); // Duplicates will be summed
          tripletList.push_back(T(j, 3 * k + B, -eij[B]));
        }
      }
    }
  }

  K.resize(nV, 3 * nV);
  K.setFromTriplets(tripletList.begin(), tripletList.end());
}
