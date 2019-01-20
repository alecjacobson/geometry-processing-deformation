#include "arap_single_iteration.h"
#include <Eigen/LU>
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <fstream>

void closest_rotation(
  const Eigen::Matrix3d & M,
  Eigen::Matrix3d & R)
{
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::Matrix3d U=svd.matrixU();
  Eigen::Matrix3d V=svd.matrixV();
  Eigen::Matrix3d omega=Eigen::Matrix3d::Identity();
  R=U*V.transpose();
  if (R.determinant()<0) omega(2,2)=-1;
  R=(U*omega*V.transpose()).transpose();
}

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  printf("start itr\n");
  Eigen::MatrixXd CT=U.transpose()*K,C,id=Eigen::MatrixXd::Identity(3,3),R(3*data.n,3);
  Eigen::Matrix3d Cinv,RR;
  for (int i=0; i<data.n; i++){
  	C=CT.block(0,3*i,3,3);
  	closest_rotation(C.transpose(),RR);
  	R.block(3*i,0,3,3)=RR;
  }
  printf("check\n");
  const Eigen::MatrixXd B=K*R/3;
  printf("%u %u\n",B.cols(),bc.cols());
  igl::min_quad_with_fixed_solve(data,B,bc,Eigen::VectorXd(),U);
  // REPLACE WITH YOUR CODE
  printf("end itr\n");
}
