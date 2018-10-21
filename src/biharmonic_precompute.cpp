#include "biharmonic_precompute.h"

#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/min_quad_with_fixed.h>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{    
    int nv = V.rows();
    Eigen::SparseMatrix<double> L, M, I, Q, Aeq;
    igl::cotmatrix(V, F, L);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

    I.resize(nv, nv);
    I.setIdentity();
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
    M = solver.compute(M).solve(I);

    Q.resize(nv, nv);
    Q = L.transpose() * M * L;

    igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);
}