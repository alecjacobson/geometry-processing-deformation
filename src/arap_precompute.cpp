#include "arap_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <iostream>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  // REPLACE WITH YOUR CODE
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V,F,L);
  igl::min_quad_with_fixed_precompute(L,b,Eigen::SparseMatrix<double>(),false,data);

  K.resize(V.rows(),3*V.rows());
  typedef Eigen::Triplet<double> T;
  std::vector<T> tlist;
  tlist.reserve(3*3*3*2*F.rows());

  for (int i = 0; i < F.rows(); i ++) {
    for (int j = 0; j < 3; j ++) {
      int u,v;
      u = F(i,(j+1)%3);
      v = F(i,(j+2)%3);
      Eigen::RowVector3d e = (V.row(u)-V.row(v))*L.coeffRef(u,v);

      for (int l = 0; l < 3; l++)
      {
        tlist.push_back(T(u, 3 * F(i, 0) + l, e[l]));
        tlist.push_back(T(v, 3 * F(i, 0) + l, -e[l]));
        tlist.push_back(T(u, 3 * F(i, 1) + l, e[l]));
        tlist.push_back(T(v, 3 * F(i, 1) + l, -e[l]));
        tlist.push_back(T(u, 3 * F(i, 2) + l, e[l]));
        tlist.push_back(T(v, 3 * F(i, 2) + l, -e[l]));
      }
    }
  }

  K.setFromTriplets(tlist.begin(),tlist.end());
  K = K / 6.0;
  
}
