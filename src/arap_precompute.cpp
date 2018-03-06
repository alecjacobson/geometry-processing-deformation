#include "arap_precompute.h"
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
  // REPLACE WITH YOUR CODE
  int n = V.rows();
  K.resize(n, 3 * n);
  
  Eigen::MatrixXd e;
  igl::cotmatrix_entries(V,F,e);
  //[1, 2] [2, 0] [0, 1]
  
  for (int row = 0; row < F.rows(); row++){
  	//iterate through each edge
  	for (int col = 0; col = F.cols(); col++){
  		int i = F(row, col);
  		int j = F(row, (col+1)%F.cols());
		int k_p[3] = {i, j , (col+2)%3};
		
		
		Eigen::MatrixXd e_ij = e(row, (i+2)%F.cols()) * (V.row(i) - V.row(j)).transpose();
		
		for (int k_i = 0; k_i < 3; k_i++)
		{
			int k = k_p[k_i];
			for (int beta = 1; beta < 4; beta++){
		  		K.coeffRef(i, 3 * k + beta) = e_ij(beta - 1, 0);
		  		K.coeffRef(j, 3 * k + beta) = - e_ij(beta - 1, 0);
	  		}
	  	}
  	}
  	
  }
  
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);
  Eigen::SparseMatrix<double> dummy(0,0);
  igl::min_quad_with_fixed_precompute(L, b, dummy ,false, data);
}
