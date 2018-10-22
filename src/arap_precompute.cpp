#include "arap_precompute.h"
#include <igl/cotmatrix_entries.h>
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>
#include <fstream>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  printf("start pre\n");
  Eigen::MatrixXd L;
  igl::cotmatrix_entries(V,F,L);
  Eigen::SparseMatrix<double> LL;
  igl::cotmatrix(V,F,LL);
  igl::min_quad_with_fixed_precompute(LL,b,Eigen::SparseMatrix<double>(),false,data);
  Eigen::RowVector3d tp;
  K.resize(V.rows(),3*V.rows());
  for (int i=0; i<F.rows(); i++){
  	for (int j=0; j<3; j++){
  	  tp=L(i,j)*(V.row(F(i,(j+1)%3))-V.row(F(i,(j+2)%3)));
  	  for (int k=0; k<3; k++)
  	  	for (int dim=0; dim<3; dim++){
          K.coeffRef(F(i,(j+1)%3),3*F(i,k)+dim)+=tp(dim);
          K.coeffRef(F(i,(j+2)%3),3*F(i,k)+dim)-=tp(dim);
		}
	}
  }
  // REPLACE WITH YOUR CODE
  printf("end pre\n");
}
