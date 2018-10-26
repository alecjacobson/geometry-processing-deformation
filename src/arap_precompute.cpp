#include "arap_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/cotmatrix_entries.h>
#include <igl/min_quad_with_fixed.h>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
    typedef Eigen::Triplet<double> Tr;
    std::vector<Tr> trips;
    Eigen::SparseMatrix<double> L, X;
    Eigen::MatrixXd C;
    
    //L, C, and precompute data
    igl::cotmatrix(V, F, L);
    igl::cotmatrix_entries(V, F, C);
    igl::min_quad_with_fixed_precompute(L, b, X, false, data);
    
    //compute K by looping over all faces
    for (int i = 0; i < F.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            Eigen::Vector3d e = (C(i, j) * (V.row(F(i,(j + 2) % 3)) - V.row(F(i,(j + 1) % 3)))) / 3.0;
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    trips.push_back(Tr(F(i, (j + 2) % 3), l + F(i, k) * 3, e(l)));
                    trips.push_back(Tr(F(i, (j + 1) % 3), l + F(i, k) * 3, -e(l)));
                }
            }
        }
    }
    
    K.resize(V.rows(), V.rows() * 3);
    K.setFromTriplets(trips.begin(), trips.end());
}
