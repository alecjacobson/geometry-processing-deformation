#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D) {

    // solve the minimization problem with bc constraints using the precomputed data.
    Eigen::MatrixXd Beq;
    Eigen::VectorXd B = Eigen::VectorXd::Zero(data.n);
    D.resize(data.n, 3);
    igl::min_quad_with_fixed_solve(data, B, bc, Beq, D);
}

