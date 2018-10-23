#include "arap_precompute.h"
#include "igl/cotmatrix.h"
#include <igl/min_quad_with_fixed.h>


void arap_precompute(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const Eigen::VectorXi &b,
        igl::min_quad_with_fixed_data<double> &data,
        Eigen::SparseMatrix<double> &K) {


    //Similar to biharmonic_precompute

    Eigen::SparseMatrix<double> L;

    igl::cotmatrix(V, F, L);
    Eigen::SparseMatrix<double> Aeq;
    igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);


    //Let's compute K now.
    K.resize(V.rows(), 3 * V.rows());

    typedef Eigen::Triplet<double> E_T;
    std::vector<E_T> Tlist;

    for (int idx = 0; idx < F.rows(); idx++) {
        for (int ij = 0; ij < 3; ij++) {
            int i = F(idx, (ij + 1) % 3);
            int j = F(idx, (ij + 2) % 3);
            Eigen::Vector3d eij = L.coeffRef(i, j) * (V.row(i) - V.row(j));

            for (int kij = 0; kij < 3; kij++) {
                // Let's get the current entry k (opposite to  i,j)
                // This is just shifting based on the ij index being processed.
                int k = F(idx, (ij + kij) % 3);

                for (int b = 0; b < 3; b++) {
                    // pushback seems to be deprecated according to Clang.
                    // let's use emplace_back instead
                    Tlist.emplace_back(E_T(i, 3 * k + b, eij(b) / 6.0));
                    Tlist.emplace_back(E_T(j, 3 * k + b, -eij(b) / 6.0));

                }
            }

        }


    }
    K.setFromTriplets(Tlist.begin(), Tlist.end());


}
