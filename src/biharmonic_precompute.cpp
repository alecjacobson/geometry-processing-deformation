#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

#include <igl/doublearea.h>
#include <igl/edge_lengths.h>

//Calculate the invserse of a mass matrix
void inverse_massmatrix(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & M)
{
  // Add your code here
    Eigen::MatrixXd areas;
    Eigen::MatrixXd l;
    igl::edge_lengths(V, F, l);
    igl::doublearea(l, areas);
    M.resize(F.maxCoeff() + 1, F.maxCoeff() + 1);
    for (int row = 0; row < F.rows(); ++row){
        for (int col = 0; col < F.cols(); ++col){
            int i = F(row, col);
            M.coeffRef(i,i) += 1/3.0 * 0.5 * areas(row);
        }
    }
    for (int i = 0; i < (F.maxCoeff() + 1); ++i){
    	M.coeffRef(i,i) = 1/M.coeffRef(i,i);
    }

}

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  // REPLACE WITH YOUR CODE
  data.n = V.rows();
  Eigen::SparseMatrix<double> M;
  Eigen::SparseMatrix<double> L;
  inverse_massmatrix(V, F, M);
  igl::cotmatrix(V, F, L);
  Eigen::SparseMatrix<double> Ltrans;
  Eigen::SparseMatrix<double> Q = L.transpose() * M * L;
  
  
  Eigen::SparseMatrix<double> dummy(0,0);
  igl::min_quad_with_fixed_precompute(Q, b, dummy ,true, data);
}





