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
    igl::min_quad_with_fixed_data<double> herdata;

    //std::cout << "Numerical issues start?" << std::endl;
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V,F,L);
    L = -L;
    igl::min_quad_with_fixed_precompute(L,b,Eigen::SparseMatrix<double>(), true,data);


    std::cout << "Laplacian term: " << std::endl;
    std::cout << V.transpose() * L * V << std::endl;

    Eigen::MatrixXd C;
    igl::cotmatrix_entries(V,F,C);
    C.noalias() = C;

    Eigen::MatrixXi E;
    igl::edges(F,E);


    const int n = V.rows();

    K = Eigen::SparseMatrix<double>(n, 3 * n);
    //Eigen::SparseMatrix<double> D(E.rows(),n);


    std::map<std::array<int,2>,int> Einv;
    std::vector<Eigen::Triplet<double>> trips;
    for(int i = 0; i < E.rows(); ++i) {
        if(E(i,0) < E(i,1)) {
            Einv[std::array<int,2>{{E(i,0),E(i,1)}}] = i;
        } else {
            Einv[std::array<int,2>{{E(i,1),E(i,0)}}] = i;
            std::swap(E(i,0), E(i,1));
        }
        trips.emplace_back(i,E(i,0),-1);
        trips.emplace_back(i,E(i,1),1);

    }


    std::vector<Eigen::Triplet<double>> htrips;



    for(int f_ = 0; f_ < F.rows(); ++f_) {
        auto&& f = F.row(f_);

        for(int a = 0; a < 3; ++a) {
            auto&& k = f((a));
            auto&& i = f((a+1)%3);
            auto&& j = f((a+2)%3);

            int ij;
            if(i < j) {
                ij = Einv[std::array<int,2>{{i,j}}];
            } else {
                ij = Einv[std::array<int,2>{{j,i}}];
            }

            double cot = C(f_,a);

            htrips.emplace_back(ij,ij,cot);


        }


    }

    Eigen::SparseMatrix<double> D(E.rows(),n);
    Eigen::SparseMatrix<double> H(E.rows(),E.rows());
    D.setFromTriplets(trips.begin(),trips.end());
    H.setFromTriplets(htrips.begin(),htrips.end());
    Eigen::SparseMatrix<double> mL(n,n);
    Eigen::SparseMatrix<double> DV3(E.rows(),3*n);
    Eigen::SparseMatrix<double> R(3*n,3);


    trips.clear();
    htrips.clear();

    Eigen::SparseMatrix<double> adp(n, 3*n);
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < 3; ++j) {
            adp.
        }
    }

    trips.clear();
    for (int k=0; k<D.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(D,k); it; ++it) {
            for(int i = 0; i < 3; ++i) {
                trips.emplace_back(it.row(), 3*it.col() + i, it.value() * V(it.col(),i));
            }
        }
    }
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < 3; ++j) {
            htrips.emplace_back(3*i+j,j,1);
        }
    }
    DV3.setFromTriplets(trips.begin(),trips.end());
    R.setFromTriplets(htrips.begin(),htrips.end());
    K = D.transpose() * H * DV3;

    std::cout << "Boring DV: " << std::endl;
    std::cout << (D * V).eval() << std::endl;
    std::cout << "special DV: " << std::endl;
    std::cout << (DV3 * R).eval() << std::endl;

    std::cout << "M<yK: " << std::endl << K << std::endl;
    std::cout << "mylap: " << std::endl << (V.transpose() * K * R).eval() << std::endl;


}
