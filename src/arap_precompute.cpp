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
  // REPLACE WITH YOUR CODE
  Eigen::SparseMatrix<double> A, Aeq;
  igl::cotmatrix(V, F, A);
  igl::min_quad_with_fixed_precompute(A, b, Aeq, false, data);

  //scale by 3 result too large
  A = A / 6.0;

  std::vector<T> Tlist;

  Tlist.reserve(F.rows() * 27 * 2);

  for (int idx = 0; idx < F.rows(); idx++)
  {
    for (int pos = 0; pos < 3; pos++)
    {
      int i = F(idx, pos % 3);
      int j = F(idx, (pos + 1) % 3);

      Eigen::RowVector3d u = (V.row(i) - V.row(j)) * A.coeff(i, j);

      for (int x = 0; x < 3; x++)
      {
        int k = F(idx, (pos + x) % 3);

        for (int y = 0; y < 3; y++)
        {
          Tlist.push_back(T(i, 3 * k + y, u(y)));
          Tlist.push_back(T(j, 3 * k + y, -u(y)));
        }
      }
    }
  }

  K.resize(V.rows(), 3 * V.rows());
  K.setFromTriplets(Tlist.begin(), Tlist.end());
}
