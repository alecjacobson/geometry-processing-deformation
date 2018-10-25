#include "biharmonic_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/min_quad_with_fixed.h>

void inverted_massmatrix(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Eigen::SparseMatrix<double> & M)
{
    igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_BARYCENTRIC, M);

    // invert diagonal matrix M
    for (int i = 0; i < M.rows(); i++)
    {
        M.coeffRef(i, i) = 1 / M.coeffRef(i, i);
    }
}

void biharmonic_precompute(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::VectorXi & b,
    igl::min_quad_with_fixed_data<double> & data)
{
    Eigen::SparseMatrix<double> M, L, Q, Aeq;
    igl::cotmatrix(V, F, L);
    inverted_massmatrix(V, F, M);
    Q = L.transpose() * M * L;
    igl::min_quad_with_fixed_precompute(Q, b, Aeq, true, data);
}

