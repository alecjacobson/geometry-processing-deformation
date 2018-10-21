#include "arap_single_iteration.h"

#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
    // Local step: find locally optimal rotation matrices R

    const int n = data.n;
    Eigen::MatrixXd R(3*n, 3); // a stack of rotation 3x3 matrix for each vertex
    Eigen::MatrixXd C = K.transpose() * U;
    Eigen::Matrix3d m, r;

    for (int k = 0; k < n; ++k) {
        m = C.block(3*k, 0, 3, 3);
        igl::polar_svd3x3(m, r);
        R.block(3*k, 0, 3, 3) = r;
    }

    // Global step: find distortion minimizing vertex positions U

    Eigen::MatrixXd Beq;
    igl::min_quad_with_fixed_solve(data, K*R, bc, Beq, U);
}