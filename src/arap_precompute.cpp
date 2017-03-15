#include "arap_precompute.h"

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <igl/cotmatrix_entries.h>
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>

#include <vector>

void arap_precompute(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::VectorXi & b,
    igl::min_quad_with_fixed_data<double> & data,
    Eigen::SparseMatrix<double> & K)
{
    // REPLACE WITH YOUR CODE                                                                                         
    const int n = V.rows();
    // columns are edges for tri-row 0,1,2 --> edges in C [1,2],[2,0],[0,1]

    Eigen::MatrixXd C(F.rows(), 3);
    igl::cotmatrix_entries(V, F, C);

    K.resize(n, 3 * n);

    std::vector< Eigen::Triplet<double> > trips;
    const int nFaces = F.rows();
    for (int f = 0; f<nFaces; ++f) // face
    {
        Eigen::RowVector3i E = F.row(f);
        for (int e = 0; e<3; ++e)  // edge e_ij = ( v_i - v_j ) in face -- CCW winding assumed
        {
            for (int v = 0; v<3; ++v)
            {
                // edges are -> [1,2],[2,0],[0,1]                        
                int i = E((e + 1) % 3); // vertex i of current edge
                int j = E((e + 2) % 3); // vertex j of current edge
                int k = F(f, v);

                Eigen::RowVector3d e_ij = V.row(i) - V.row(j);
                double c_ij = C(f, e);
                Eigen::RowVector3d eHat_ij = c_ij*e_ij;

                for (int beta = 0; beta<3; ++beta)
                {
                    trips.push_back(Eigen::Triplet<double>(i, 3 * k + beta, eHat_ij(beta)));
                    trips.push_back(Eigen::Triplet<double>(j, 3 * k + beta, -eHat_ij(beta)));
                }
            }
        }
    }

    K.setFromTriplets( trips.begin(), trips.end() );

    //Eigen::MatrixXd R(3 * n, 3);
    //Eigen::MatrixXd Ct(3 * n, 3);
    //Ct = V.transpose() * K;

    // now prep for the global step
    Eigen::SparseMatrix<double> L(n, n);
    igl::cotmatrix(V, F, L);

    Eigen::SparseMatrix<double> Aeq;
    igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);
}

