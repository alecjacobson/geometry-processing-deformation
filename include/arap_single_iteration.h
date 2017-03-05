#ifndef ARAP_SINGLE_ITERATION_H
#define ARAP_SINGLE_ITERATION_H
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace igl
{
  template <typename T> struct min_quad_with_fixed_data;
}

// Conduct a single iteration of the local-global solver for minimizing the
// as-rigid-as-possible energy.
//
// Inputs:
//   data  pre-factorized system matrix etc. (see `arap_precompute`
//   K  pre-constructed bi-linear term of energy combining rotations and
//     positions
//   U  #V by dim list of current positions
// Outputs:
//   U  #V by dim list of new positions
void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U);

#endif
