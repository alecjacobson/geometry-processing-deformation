#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
    const int size = data.n;
    const int dims = bc.cols();

    min_quad_with_fixed_solve(data, Eigen::MatrixXd::Zero(size, dims),
                                bc, Eigen::MatrixXd::Zero(0, dims), D);
}

