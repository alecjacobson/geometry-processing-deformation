#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <iostream>


void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  // REPLACE WITH YOUR CODE
  int vetr_num = V.rows();
  Eigen::SparseMatrix<double> A;
  Eigen::SparseMatrix<double> L;
  Eigen::SparseMatrix<double> M;
  Eigen::SparseMatrix<double> Q;
  Eigen::SparseMatrix<double> Aeq;

  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

  Eigen::SparseMatrix<double> M2(vetr_num, vetr_num); ;
  for (int i = 0; i < vetr_num; i++) {
  	M2.insert(i,i) = 1.0 / M.coeffRef(i,i);
  }
  std::cout << M2 << std::endl;
  std::cout << "succes" << std::endl;
  Q = L.transpose() * M2 * L;
  igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);

}
