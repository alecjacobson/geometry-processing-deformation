#include "arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>

void arap_single_iteration(
    const igl::min_quad_with_fixed_data<double> & data,
    const Eigen::SparseMatrix<double> & K,
    const Eigen::MatrixXd & bc,
    Eigen::MatrixXd & U)
{
    Eigen::MatrixXd R(3 * U.rows(), 3);
    Eigen::MatrixXd C = U.transpose() * K;

    for (int i = 0; i < U.rows(); i++)
    {
        Eigen::Matrix3d Ri;
        Eigen::Matrix3d Ci = C.block(0, 3 * i, 3, 3).eval();
        Ci /= Ci.maxCoeff();
        igl::polar_svd3x3(Ci, Ri);
        R.block(3 * i, 0, 3, 3) = Ri.transpose();
    }

    Eigen::MatrixXd B = K * R;
    Eigen::VectorXd Beq; // empty constraint
    igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);
}
