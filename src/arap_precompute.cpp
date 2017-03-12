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

    std::cout << "Laplacian term: " << std::endl;
    std::cout << V.transpose() * L * V << std::endl;

    Eigen::MatrixXd C;
    igl::cotmatrix_entries(V,F,C);
    C.noalias() = C;

    std::cout << "Cot entries: " << std::endl << C << std::endl;

    Eigen::MatrixXi E;
    igl::edges(F,E);

    std::vector<Eigen::Triplet<double>> trips;

    const int n = V.rows();

    K = Eigen::SparseMatrix<double>(n, 3 * n);
    //Eigen::SparseMatrix<double> D(E.rows(),n);


    Eigen::MatrixXd EV(E.rows(),3);
    std::map<std::array<int,2>,int> Einv;
    for(int i = 0; i < EV.rows(); ++i) {
        if(E(i,0) < E(i,1)) {
            Einv[std::array<int,2>{{E(i,0),E(i,1)}}] = i;
        } else {
            Einv[std::array<int,2>{{E(i,1),E(i,0)}}] = i;
            std::swap(E(i,0), E(i,1));
        }

        auto&& a = V.row(E(i,0));
        auto&& b = V.row(E(i,1));
        EV.row(i) = b-a;
        //EV.row(i) /= EV.row(i).squaredNorm();
    }

    std::cout << "EV: " << std::endl;
    std::cout << EV << std::endl;
    std::cout << std::endl;

    Eigen::MatrixXd R(3*n,3);

    Eigen::Vector3d Le;
    Eigen::Vector3d Ke;


    Eigen::MatrixXd mL(n,n);
    mL.setZero();
    Eigen::MatrixXd mD(EV.rows(),n);
    Eigen::VectorXd mH(EV.rows());
    mD.setZero();
    mH.setZero();
    Eigen::MatrixXd D(3*n,n);
    Eigen::MatrixXd DD(3*n,n);
    Eigen::VectorXd H(3*n);
    D.setZero();
    H.setZero();
    for(int f_ = 0; f_ < F.rows(); ++f_) {
        auto&& f = F.row(f_);
        Le.setZero();
        Ke.setZero();
        std::cout << "Face: " << f << std::endl;

        Eigen::Vector3d a = (V.row(f(1))-V.row(f(0))).transpose();
        Eigen::Vector3d b = (V.row(f(2))-V.row(f(0))).transpose();
        Eigen::Matrix3d A;
        double area = .5 * b.cross(a).norm();
        for(int a = 0; a < 3; ++a) {
            auto&& k = f((a));
            auto&& i = f((a+1)%3);
            auto&& j = f((a+2)%3);

            int ij;
            double sign = 1;
            if(i < j) {
                ij = Einv[std::array<int,2>{{i,j}}];
            } else {
                sign = -1;
                ij = Einv[std::array<int,2>{{j,i}}];
            }
            mD(ij,i) =-1;
            mD(ij,j) = 1;

            Eigen::RowVector3d e = sign * EV.row(ij);

            double cot = C(f_,a);

            Le(a) = .5 * cot * EV.row(ij).norm();
            mH(ij) += cot;


            for(int a = 0; a < 3; ++a) {
                H(3*i + a) += cot;
                double t = EV.row(ij)(a);
                //double t = EV.row(ij)(a);
                D(3*i + a,i) += t;
                D(3*i + a,j) -= t;
            }

            Eigen::RowVector3d v = .5 * cot * EV.row(ij);

            mL(i,j) -= cot;
            mL(j,i) -= cot;
            mL(i,i) += cot;
            mL(j,j) += cot;
            for(int b = 0; b < 3; ++b) {
                int ib = i + n * b;
                int jb = j + n * b;
                trips.emplace_back(i,ib,v(b));
                trips.emplace_back(j,jb,-v(b));
            }
            A.col(a) = e;
        }

        Eigen::JacobiSVD<Eigen::Matrix3d> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
        std::cout << "Singular values: " << svd.singularValues().transpose() << std::endl;
        Eigen::Vector3d uinv = 1.0 / svd.singularValues().array();
        for(int i = 0 ; i < uinv.rows(); ++i) {
            if(std::abs(uinv(i)) > 1e5) {
                uinv(i) = 0;
            }
        }

        A = svd.matrixU() * uinv.asDiagonal() * svd.matrixV().transpose();
        for(int a = 0; a < 3; ++a) {
            auto&& k = f((a));
            auto&& i = f((a+1)%3);
            auto&& j = f((a+2)%3);

            for(int a = 0; a < 3; ++a) {
                DD(3*i + a,i) += A(a,i);
                DD(3*i + a,j) -= A(a,i);
            }

        }

        
    }
    std::cout << "D: " << std::endl << D << std::endl;
    std::cout << "DD: " << std::endl << DD << std::endl;
    std::cout << "DD: " << std::endl << (DD.transpose() * D ).eval() << std::endl;
    std::cout << "H: " << std::endl << H << std::endl;
    std::cout << "mD: " << std::endl << mD << std::endl;
    std::cout << "mH: " << std::endl << mH << std::endl;
    std::cout << "Lap: \n" << L << std::endl;
    std::cout << "MYLAP: " << std::endl << mL << std::endl;
    mL = mD.transpose() * mH.asDiagonal() * mD;
    std::cout << "MYLAP: " << std::endl << mL << std::endl;
    mL = D.transpose() * H.asDiagonal() * D;
    std::cout << "MYLAP: " << std::endl << mL << std::endl;


    K.setFromTriplets(trips.begin(),trips.end());
    return;


    {


        Eigen::MatrixXd U = V;
        Eigen::MatrixXd R(3*U.rows(),3);
        Eigen::MatrixXd  C = K.transpose() * U;
        std::cout << "K: " << std::endl << K << std::endl;
        std::cout << "U: " << std::endl << U << std::endl;
        std::cout << "C: " << std::endl << C << std::endl;

        for(int i = 0; i < U.rows(); ++i) {
            Eigen::Matrix3d c;
            for(int j = 0; j < 3; ++j) {
                c.col(j) = C.row(i + j * U.rows());
            }
            std::cout << "Baby c: " << std::endl << c << std::endl;
            Eigen::Matrix3d r;
            igl::polar_svd3x3<Eigen::Matrix3d>(c,r);
            std::cout << "r: " << std::endl << r << std::endl;
            for(int j = 0; j < 3; ++j) {
                R.row(i + U.rows() * j) = r.col(j);
            }
        }
        std::cout << "VTKR:"<<std::endl;
        std::cout << V.transpose() * K * R << std::endl;
    }
}
