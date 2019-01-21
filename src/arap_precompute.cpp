#include "arap_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix_entries.h>
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

  Eigen::SparseMatrix<double> Aeq;
  igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);

  Eigen::MatrixXd cot_mtx;
  igl::cotmatrix_entries(V, F, cot_mtx);

  typedef Eigen::Triplet<double> tuple;
  std::vector<tuple> tuple_list;
  K.resize(V.rows(),3*V.rows());
  K.setZero();

  for (int face_idx = 0; face_idx < F.rows(); face_idx++)
    for (int edge_idx = 0; edge_idx < 3; edge_idx++)
    {
      int i = F(face_idx, (edge_idx + 1) % 3);
      int j = F(face_idx, (edge_idx + 2) % 3);
      Eigen::Vector3d edge = V.row(i) - V.row(j);
      for (int beta = 0; beta < 3; beta++)
      {
        double e_ij = cot_mtx(face_idx,edge_idx) * edge(beta) / 3;
        for (int vert_idx = 0; vert_idx < 3; vert_idx++)
        {
          int k = F(face_idx, vert_idx);
          tuple_list.push_back(tuple(i, 3*k+beta, e_ij));
          tuple_list.push_back(tuple(j, 3*k+beta, -e_ij));
        }
      }
    }
  K.setFromTriplets(tuple_list.begin(), tuple_list.end());
}
