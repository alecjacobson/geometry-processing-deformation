#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  // REPLACE WITH YOUR CODE
  data.n = V.rows();
}

