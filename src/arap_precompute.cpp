#include "arap_precompute.h"
#include <iostream>
#include <igl/cotmatrix.h>
#include <igl/cotmatrix_entries.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/edges.h>
#include <igl/polar_svd3x3.h>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{

    //std::cout << "Numerical issues start?" << std::endl;
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V,F,L);
    L = -L;
    igl::min_quad_with_fixed_precompute(L,b,Eigen::SparseMatrix<double>(), true,data);


    Eigen::MatrixXd C;
    igl::cotmatrix_entries(V,F,C);
    C.noalias() = C / 3.0;



    const int n = V.rows();

    K = Eigen::SparseMatrix<double>(n, 3 * n);


    std::vector<Eigen::Triplet<double>> trips;



    for(int f_ = 0; f_ < F.rows(); ++f_) {
        auto&& f = F.row(f_);
        for(int a = 0; a < 3; ++a) {
                auto&& k = f(a);
                for(int b = 0; b < 3; ++b) {
                    auto&& c = C(f_,b);
                    auto&& i = f((b+1)%3);
                    auto&& j = f((b+2)%3);
                    for(int dim = 0; dim < 3; ++dim) {
                        trips.emplace_back(i,3*k+dim,c * (V(i,dim) - V(j,dim)));
                        trips.emplace_back(j,3*k+dim,c * (V(j,dim) - V(i,dim)));
                    }
                }
        }
    }
    K.setFromTriplets(trips.begin(),trips.end());

}
