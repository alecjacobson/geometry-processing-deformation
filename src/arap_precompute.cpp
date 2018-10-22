#include "arap_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix_entries.h>
#include <igl/cotmatrix.h>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  // construct K
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  
  Eigen::MatrixXd C(F.rows(), 3);
  igl::cotmatrix_entries(V, F, C);
  
  for (int f = 0; f < F.rows(); f++) {
    for (int e = 0; e < 3; e++) {
      int i = F(f, e);
      int j = F(f, (e+1)%3);
      // Not sure why I had to use a factor of 1/9 to avoid scaling
      // Theoretically, it should be (1/6)*2 = 1/3
      Eigen::Vector3d eij = (1/9.)*C(f, (e+2)%3) * (V.row(i) - V.row(j));
      for (int ki = 0; ki < 3; ki++) {
        int k = F(f, ki);
        for (int a = 0; a < 3; a++) {
          for (int b = 0; b < 3; b++) {
            tripletList.push_back(T(i, 3*k + b, eij[b]));
            tripletList.push_back(T(j, 3*k + b, -eij[b]));
          }
        }
      }
    }
  }
  
  K.resize(V.rows(), 3*V.rows());
  K.setFromTriplets(tripletList.begin(), tripletList.end());
  
  // prefactorize quadratic form corresponding to global energy
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);
  igl::min_quad_with_fixed_precompute(L, b, Eigen::SparseMatrix<double>(), false, data);
}
