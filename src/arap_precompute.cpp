#include "arap_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/arap.h>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
    // Construct L
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V,F,L);

    igl::min_quad_with_fixed_precompute(L, b, Eigen::SparseMatrix<double>(), false, data);

    /*
     * Note: When constructing K it’s easiest to iterate over all half-edges in the mesh
     * (by iterating over all faces and then each of the three edges).
     * Each half-edge ij contributes terms tying vi,vj to each of the (three)
     * rotations Rk that apply against their difference.
     */

    std::vector<Eigen::Triplet<double>> tripletList;

    // Construct K
    K.resize(V.rows(), 3 * V.rows());

    for (int f = 0; f < F.rows(); f++) {
        for (int k = 0; k < 3; k++) {
            for (int b = 0; b < 3; b++) {
                // edge 01
                tripletList.push_back(Eigen::Triplet<double>(F(f,0), 3*F(f,k)+b , L.coeff(F(f,0),F(f,1)) * ( V(F(f,0),b)-V(F(f,1),b))/6));
                tripletList.push_back(Eigen::Triplet<double>(F(f,1), 3*F(f,k)+b , L.coeff(F(f,0),F(f,1)) * (-V(F(f,0),b)+V(F(f,1),b))/6));
                // edge 12
                tripletList.push_back(Eigen::Triplet<double>(F(f,1), 3*F(f,k)+b , L.coeff(F(f,1),F(f,2)) * ( V(F(f,1),b)-V(F(f,2),b))/6));
                tripletList.push_back(Eigen::Triplet<double>(F(f,2), 3*F(f,k)+b , L.coeff(F(f,1),F(f,2)) * (-V(F(f,1),b)+V(F(f,2),b))/6));
                // edge 20
                tripletList.push_back(Eigen::Triplet<double>(F(f,2), 3*F(f,k)+b , L.coeff(F(f,2),F(f,0)) * ( V(F(f,2),b)-V(F(f,0),b))/6));
                tripletList.push_back(Eigen::Triplet<double>(F(f,0), 3*F(f,k)+b , L.coeff(F(f,2),F(f,0)) * (-V(F(f,2),b)+V(F(f,0),b))/6));
            }
        }
    }

    K.setFromTriplets(tripletList.begin(), tripletList.end());
}
