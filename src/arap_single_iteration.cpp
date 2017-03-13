#include "arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  
}
