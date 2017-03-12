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


    Eigen::MatrixXd  C = K.transpose() * U;

    for(int i = 0; i < U.rows(); ++i) {
        Eigen::Matrix3d a;
        for(int j = 0; j < 3; ++j) {
            a.col(j) = C.row(i + j * U.rows());
        }
        Eigen::Matrix3d r;
        igl::polar_svd3x3<Eigen::Matrix3d>(a,r);
        for(int j = 0; j < 3; ++j) {
            R.row(i + U.rows() * j) = r.col(j);
        }
    }

    Eigen::MatrixXd KR =   (K * R);
  for(int i = 0; i < 3; ++i) {
      Eigen::VectorXd Y = bc.col(i);
      Eigen::VectorXd B = KR.col(i);
      Eigen::VectorXd Beq;
      Eigen::VectorXd Z ;
      igl::min_quad_with_fixed_solve(data,B,Y,Beq,Z);
      U.col(i) = Z;
  }
  /*
    for(int i = 0; i < U.rows(); ++i) {
       U.row(i).transpose()  = V.block(3*i,0,3,1) ;
    }
    */

}
