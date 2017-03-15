#include "arap_single_iteration.h"

#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>

#include <ostream>

void arap_single_iteration(
    const igl::min_quad_with_fixed_data<double> & data,
    const Eigen::SparseMatrix<double> & K,
    const Eigen::MatrixXd & bc,
    Eigen::MatrixXd & U)
{
    const int n = data.n;

    // local - solve for the rotations...
    Eigen::MatrixXd R( 3*n, 3 );
    R = Eigen::MatrixXd::Zero( 3*n, 3 );

    Eigen::MatrixXd C( 3*n, 3 );
    C = (U.transpose() * K).transpose();
    
    // solve min tr( Ct R ), walk each vertexes (k) rotation Rk (3,3)
    // and weighted covariance Ck (3,3)
    for (int k = 0; k < n; ++k)
    {
        Eigen::Matrix3d Rk = Eigen::Matrix3d::Zero();
        Eigen::Matrix3d Ck = C.block( 3*k, 0, 3, 3 );

        
        igl::polar_svd3x3( Ck, Rk );

        // write into R - I love blocks
        R.block( 3*k, 0, 3, 3 ) = Rk;
    }

    // global - use the rotations and solve the quadratic
    Eigen::MatrixXd B = K * R;
    Eigen::VectorXd Beq; // empty equality constraints
    igl::min_quad_with_fixed_solve( data, B, bc, Beq, U );
}
