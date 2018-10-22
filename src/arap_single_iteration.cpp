#include "arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>

void arap_single_iteration(
        const igl::min_quad_with_fixed_data<double> &data,
        const Eigen::SparseMatrix<double> &K,
        const Eigen::MatrixXd &bc,
        Eigen::MatrixXd &U) {

    //Let's first solve the local step


    // C^T=V^T*K
    // C=(V^T*K)^T so...
    Eigen::MatrixXd C = (U.transpose() * K).transpose();


    Eigen::MatrixXd R(3 * data.n, 3);

    //R
    for (int idx = 0; idx < data.n; idx++) {

        Eigen::Matrix3d Rk;
        Eigen::Matrix3d Ck = C.block(idx * 3, 0, 3, 3);
        igl::polar_svd3x3(Ck, Rk);
        R.block(idx * 3, 0, 3, 3) = Rk;
    }

    //Let's now solve the final step
    Eigen::MatrixXd Beq;
    igl::min_quad_with_fixed_solve(data, K * R, bc, Beq, U);


}
