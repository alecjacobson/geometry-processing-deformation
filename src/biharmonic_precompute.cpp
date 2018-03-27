#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include "igl/massmatrix.h"
#include "igl/cotmatrix.h"
#include "igl/invert_diag.h"
#include <iostream>
using namespace std;

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{

    
    int numV = F.maxCoeff() + 1;
    int numH = b.size();
    Eigen::VectorXi newB;
    Eigen::SparseMatrix<double> Q, L, M, M_inv;
    
    Eigen::SparseMatrix<double> Aeq;
    
    L.resize(numV,numV);
    M.resize(numV,numV);
    M_inv.resize(numV,numV);
    Q.resize(numV,numV);
    
    igl::cotmatrix(V,F,L);
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
    igl::invert_diag(M,M_inv);
    
    Q = L.transpose() * M_inv * L;
    //Quh = Q.block(0,numV-numH, numV-numH,numH);
     //Qhh = Q.block(numV-numH,numV-numH, numH,numH);
    
    //B.resize(numV-numH, 1);
    //B = Quh * b;
    
    //Aeq.resize(1,numV);
    //cout << "hi" << endl;
    min_quad_with_fixed_precompute(Q,b,Aeq,1,data);
    //cout << "hi" << endl;
}

