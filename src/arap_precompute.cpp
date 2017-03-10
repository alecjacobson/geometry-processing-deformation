#include "arap_precompute.h"
#include <iostream>
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

    std::cout << "Numerical issues start?" << std::endl;
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V,F,L);
    L = -L;
    std::cout << "Precompute?" << std::endl;
    igl::min_quad_with_fixed_precompute(L,b,Eigen::SparseMatrix<double>(), true,data);
    std::cout << "Numerical issues?" << std::endl;

    Eigen::MatrixXd C;
    igl::cotmatrix_entries(V,F,C);
    std::vector<Eigen::Triplet<double>> trips;

    std::cout << "C: " << std::endl << C << std::endl;
    K = Eigen::SparseMatrix<double>(V.rows(), 3 * V.rows());
    for(int k = 0; k < F.rows(); ++k) {
        auto&& f = F.row(k);
        for(int i = 0; i < 3; ++i) {
            auto&& a = V.row((i+0)%3);
            auto&& b = V.row((i+1)%3);
            auto&& c = V.row((i+2)%3);

            Eigen::RowVector3d ba = b-a;
            Eigen::RowVector3d ca = c-a;
            Eigen::RowVector3d cb = c-b;

            std::cout << "CB: " << cb << std::endl;

            double cot = C(k,i);

            for(int l = 0; l < 3; ++l) {
                for(int j = 0; j < 3; ++j) {
                    trips.emplace_back(f(l),3*f(l)+j,cot * cb(j));
                }
            }
        }
    }

    K.setFromTriplets(trips.begin(),trips.end());




}
