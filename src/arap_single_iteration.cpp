#include "arap_single_iteration.h"
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

    std::cout << "K shape: " << K.rows() << " " << K.cols() << std::endl;
    std::cout << "U shape: " << U.rows() << " " << U.cols() << std::endl;

    Eigen::MatrixXd  C = K.transpose() * U;

    for(int i = 0; i < U.rows(); ++i) {
        auto a = C.block(3*i,i,3,3);
        Eigen::Matrix3d r;
        igl::polar_svd3x3<Eigen::Matrix3d>(a,r);
        R.block(3*i,i,3,3) = r;;
    }

  Eigen::SparseMatrix<double> KR = K * R;
  //igl::min_quad_with_fixed_solve(data,B,bc,Beq,U);

  for(int i = 0; i < 3; ++i) {
      auto Y = bc.col(i);
      auto B = KR.col(i);
      auto Beq = Eigen::SparseMatrix<double>();
      Eigen::VectorXd Z = U.col(i);
      igl::min_quad_with_fixed_solve(data,B,Y,Beq,Z);
      U.col(i) = Z;
  }
}
