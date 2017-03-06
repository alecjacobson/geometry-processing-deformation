#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  // REPLACE WITH YOUR CODE





    Eigen::SparseMatrix<double> L,M,Q;

    igl::cotmatrix(V,F,L);

    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT, M);

    for(int i = 0; i < M.rows(); ++i) {
        auto&& v = M.coeffRef(i,i);
        v = 1.0/v;
    }

    Q = L.transpose() * M * L;



    igl::min_quad_with_fixed_precompute(Q,b,Eigen::SparseMatrix<double>(), true,data);


}

