#include "arap_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>
#include <iostream>

void compute_cot_diff(
  const int a, const int b, const Eigen::MatrixXd & V, const double cot, 
  std::vector<Eigen::Triplet<double>> & tripletList)
{
  int num_v = V.rows();
  double diff_x, diff_y, diff_z;
  diff_x = V(a, 0) - V(b, 0);
  diff_y = V(a, 1) - V(b, 1);
  diff_z = V(a, 2) - V(b, 2);
  double val = 0.5 * cot;
  
  tripletList.push_back(Eigen::Triplet<double>(a, b, val * diff_x));
  tripletList.push_back(Eigen::Triplet<double>(b, a, val * -diff_x));
  tripletList.push_back(Eigen::Triplet<double>(a, a, val * diff_x));
  tripletList.push_back(Eigen::Triplet<double>(b, b, val * -diff_x));
  
  tripletList.push_back(Eigen::Triplet<double>(a, num_v + b, val * diff_y));
  tripletList.push_back(Eigen::Triplet<double>(b, num_v + a, val * -diff_y));
  tripletList.push_back(Eigen::Triplet<double>(a, num_v + a, val * diff_y));
  tripletList.push_back(Eigen::Triplet<double>(b, num_v + b, val * -diff_y));
  
  tripletList.push_back(Eigen::Triplet<double>(a, 2*num_v + b, val * diff_z));
  tripletList.push_back(Eigen::Triplet<double>(b, 2*num_v + a, val * -diff_z));
  tripletList.push_back(Eigen::Triplet<double>(a, 2*num_v + a, val * diff_z));
  tripletList.push_back(Eigen::Triplet<double>(b, 2*num_v + b, val * -diff_z));
}

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  int num_v = V.rows();
  Eigen::SparseMatrix<double> L(num_v, num_v);
  igl::cotmatrix(V, F, L);
  
  // Construct K
  K.resize(num_v, 3 * num_v);
  std::vector<Eigen::Triplet<double>> tripletList;
  int v0, v1, v2;
  for (int i = 0; i < F.rows(); i++) {
    v0 = F(i, 0);
    v1 = F(i, 1);
    v2 = F(i, 2);
    compute_cot_diff(v0, v1, V, L.coeffRef(v0,v1), tripletList);
    compute_cot_diff(v1, v2, V, L.coeffRef(v1,v2), tripletList);
    compute_cot_diff(v2, v0, V, L.coeffRef(v2,v0), tripletList);
  }
  K.setFromTriplets(tripletList.begin(), tripletList.end());
  
  /* Setup inputs to min_quad 
   * L = A
   * b = known
   * Aeq = empty matrix
   * pd = false
   * */
  Eigen::SparseMatrix<double> Aeq;
  igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);
}
