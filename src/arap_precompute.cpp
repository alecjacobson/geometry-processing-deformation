#include "arap_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>

#include "minitrace.h"

using namespace Eigen;
using namespace std;

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  // prefactorize L first
  SparseMatrix<double> L, Aeq;
  igl::cotmatrix(V, F, L);
  igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);

  // construct K
  vector<Triplet<double>> triplets;
  triplets.reserve(F.rows() * 3 * 3 * 3 * 2);

  for (int m = 0; m < F.rows(); m++) {
    for (int n = 0; n < 3; n++) {
      int i = F(m, n % 3);
      int j = F(m, (n + 1) % 3);
      RowVector3d eij = (V.row(i) - V.row(j)) * L.coeff(i, j) / 6.0;
      // iterate all three vertices
      for (int a = 0; a < 3; a++) {
        int k = F(m, (n + a) % 3);
        // iterate all three coordinates
        for (int b = 0; b < 3; b++) {
          triplets.push_back(Triplet<double>(i, 3 * k + b, eij(b)));
          triplets.push_back(Triplet<double>(j, 3 * k + b, -eij(b)));
        }
      }
    }
  }

  K.resize(V.rows(), V.rows() * 3);

  // MTR_BEGIN("C++", "setFrom");
  K.setFromTriplets(triplets.begin(), triplets.end());
  // MTR_END("C++", "setFrom");

}
