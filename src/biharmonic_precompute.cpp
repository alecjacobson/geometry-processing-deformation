#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V,F,L);

    Eigen::SparseMatrix<double> M;
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT, M);

    Eigen::SparseMatrix<double> M_inverse;
    igl::invert_diag(M,M_inverse);

    Eigen::SparseMatrix<double> Q;
    Q = L.transpose()*M_inverse*L;

    igl::min_quad_with_fixed_precompute(Q, b, Eigen::SparseMatrix<double>(), true, data);
}

