#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h> 
#include <fstream>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  printf("start bipre\n");
  Eigen::SparseMatrix<double> L,M;
  igl::cotmatrix(V,F,L);
  igl::massmatrix(V,F,igl::MassMatrixType(1),M);
  Eigen::VectorXd diag=M.diagonal();
  diag=(1.0/diag.array()).matrix();
  typedef Eigen::Triplet<double> T;
  std::vector<T> list;
  Eigen::SparseMatrix<double> M1;
  for (int i=0; i<diag.size(); i++){
    list.push_back(T(i,i,diag(i)));
  }
  M1.resize(diag.size(),diag.size());
  M1.setFromTriplets(list.begin(), list.end()); 
  Eigen::SparseMatrix<double> A=L.transpose()*M1*L;
  igl::min_quad_with_fixed_precompute(A,b,Eigen::SparseMatrix<double>(),false,data);
  // REPLACE WITH YOUR CODE
  //data.n = V.rows();
  printf("end bipre\n");
}

