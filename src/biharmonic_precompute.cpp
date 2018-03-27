#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include "igl/massmatrix.h"
#include "igl/cotmatrix.h"
#include "igl/invert_diag.h"


void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{

    //get number of vertices
    int numV = F.maxCoeff() + 1;
    
    //Q defined as in the README
    //L cotangent matrix
    //M mass matrix
    Eigen::SparseMatrix<double> Q, L, M, M_inv;
    
    Eigen::SparseMatrix<double> Aeq;
    
    L.resize(numV,numV);
    M.resize(numV,numV);
    M_inv.resize(numV,numV);
    Q.resize(numV,numV);
    
    //Compute cotangent, mass and inverted-mass matrices
    igl::cotmatrix(V,F,L);
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
    igl::invert_diag(M,M_inv);
    
    Q = L.transpose() * M_inv * L;
    
    //Precompute the data
    min_quad_with_fixed_precompute(Q,b,Aeq,1,data);

}

