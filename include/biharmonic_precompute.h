#ifndef BIHARMONIC_PRECOMPUTE_H
#define BIHARMONIC_PRECOMPUTE_H
#include <Eigen/Core>

namespace igl
{
  template <typename T> struct min_quad_with_fixed_data;
}

// Precompute data needed to efficiently solve for a biharmonic deformation.
//
// Inputs:
//   V  #V by dim vertex positions
//   F  #F by simplex-size list of element indices
//   b  #b indices into V of handle vertices
// Outputs:
//   data  pre-factorized system matrix etc. (see `igl::min_quad_with_fixed`)
void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data);

#endif
