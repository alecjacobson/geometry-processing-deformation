#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <iostream>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  // REPLACE WITH YOUR CODE
  Eigen::SparseMatrix<double> A; 
  igl::cotmatrix(V, F, A);
  min_quad_with_fixed_precompute(A, b, Eigen::SparseMatrix<double>(), false, data);
}

