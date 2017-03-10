#include "arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>
#include <iostream>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
    
    std::cout << "Input size: " << U.rows() << " " << U.cols() << std::endl;

    Eigen::MatrixXd R(3*U.rows(),3);
    Eigen::VectorXd V(U.size());
    for(int i = 0; i < U.rows(); ++i) {
        V.block(3*i,0,3,1) = U.row(i).transpose();
    }

    std::cout << "K: " << std::endl << K << std::endl;
    std::cout << "U: " << std::endl << U << std::endl;

    Eigen::MatrixXd  C = K.transpose() * U;
    std::cout << "C: " << std::endl << C << std::endl;

    for(int i = 0; i < U.rows(); ++i) {
        auto a = C.block(3*i,0,3,3);
        Eigen::Matrix3d r;
        igl::polar_svd3x3<Eigen::Matrix3d>(a,r);
        R.block(3*i,0,3,3) = r;
        std::cout << a << std::endl << std::endl;
    }
    return;

    Eigen::MatrixXd KR = (K * R);

  for(int i = 0; i < 3; ++i) {
      Eigen::VectorXd Y = bc.col(i);
      Eigen::VectorXd B = KR.col(i);
      Eigen::VectorXd Beq = Eigen::VectorXd::Zero(data.n);
      Eigen::VectorXd Z = U.col(i);
      igl::min_quad_with_fixed_solve(data,B,Y,Beq,Z);
      U.col(i) = Z;
  }
}
