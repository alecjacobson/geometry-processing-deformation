#include "arap_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>

void arap_precompute(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::VectorXi & b,
    igl::min_quad_with_fixed_data<double> & data,
    Eigen::SparseMatrix<double> & K)
{
    Eigen::SparseMatrix<double> L, A, Aeq;
    igl::cotmatrix(V, F, L);
    A = 2 * L; // need to double since igl::min_quad_with_fixed minimizes 0.5Z'*A*Z + Z'*B + C
    igl::min_quad_with_fixed_precompute(A, b, Aeq, false, data);

    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;
    triplets.reserve(18 * F.size());

    for (int i = 0; i < F.rows(); i++)
    {
        Eigen::Vector3i f = F.row(i);
        for (int j = 0; j < 3; j++)
        {
            int v1 = f[j];
            int v2 = f[(j + 1) % 3];
            int v3 = f[(j + 2) % 3];

            Eigen::RowVector3d d = V.row(v1) - V.row(v2);
            Eigen::RowVector3d e = 1 / (double)3 * L.coeffRef(v1, v2) * d;

            for (int k = 0; k < 3; k++)
            {
                triplets.push_back(T(v1, 3 * v1 + k, e[k]));
                triplets.push_back(T(v2, 3 * v1 + k, -e[k]));
                triplets.push_back(T(v1, 3 * v2 + k, e[k]));
                triplets.push_back(T(v2, 3 * v2 + k, -e[k]));
                triplets.push_back(T(v1, 3 * v3 + k, e[k]));
                triplets.push_back(T(v2, 3 * v3 + k, -e[k]));
            }
        }
    }

    K.resize(V.rows(), 3 * V.rows());
    K.setFromTriplets(triplets.begin(), triplets.end());
}
