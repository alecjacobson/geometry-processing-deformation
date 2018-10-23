#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>

/// Whitelisted igl functions
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

void biharmonic_precompute(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const Eigen::VectorXi &b,
        igl::min_quad_with_fixed_data<double> &data) {

    // Firstly, we need to compute Q where
    // Q = L^T M^{-1} L

    //Computing L and M
    Eigen::SparseMatrix<double> M;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);

    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);

    //Let's compute M^{-1} L

    // According to the documentation of eigen
    // https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html
    // this is the recommended solver for very sparse and not too large problems
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(M);
    Eigen::SparseMatrix<double> M_1_L = solver.solve(L).sparseView();

    //Q
    Eigen::SparseMatrix<double> Q = L.transpose() * M_1_L;

    Eigen::SparseMatrix<double> Aeq;
    igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);

}

