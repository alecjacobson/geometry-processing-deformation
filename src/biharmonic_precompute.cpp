#include "biharmonic_precompute.h"
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
    Eigen::SparseMatrix<double> M, L, Q, X, M_L;
    Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > S;
    
    //compute bi-Laplacian Q = L_transpose * M_inverse * L
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
    igl::cotmatrix(V, F, L);
    S.compute(M);
    M_L = S.solve(L).sparseView();
    Q = L.transpose() * M_L;
    
    //precompute data
    igl::min_quad_with_fixed_precompute(Q, b, X, false, data);
}
