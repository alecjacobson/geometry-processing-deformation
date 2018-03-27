#include "arap_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
    //Compute number of vertices and faces
    int numV = F.maxCoeff() + 1;
    int numF = F.rows();
    
    //Triplet list for computing K
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    //tripletList.reserve(F.rows()*54);
    
    Eigen::MatrixXd vals(3,2);
    
    //Makes it easier to reference other vertices
    vals(0,0) = 0;
    vals(0,1) = 1;
    vals(1,0) = 1;
    vals(1,1) = 2;
    vals(2,0) = 2;
    vals(2,1) = 0;

    //Cotangent matrix
    Eigen::SparseMatrix<double> L;
    
    //Placeholder matrix
    Eigen::SparseMatrix<double> Aeq;
    
    L.resize(numV,numV);
    K.resize(numV,3*numV);
    
    //Compute cotangent matrix
    igl::cotmatrix(V,F,L);
    L = L /6.0;

    //Precompute the data, L is not positive-semidefinite
    min_quad_with_fixed_precompute(L,b,Aeq,false,data);
    
    //Need to compute K now
    for (int fNo = 0; fNo < numF; fNo ++) {
        for (int eNo = 0; eNo < 3; eNo ++) {
            Eigen::Vector3d edgeDiff;
            //compute e_ij as in README
            edgeDiff.array() = V.row(F(fNo,vals(eNo,1))).array() - V.row(F(fNo,vals(eNo,0))).array();
            edgeDiff.array() = edgeDiff.array() * L.coeffRef(F(fNo,vals(eNo,1)),F(fNo,vals(eNo,0)));
            
            //Adds to triplet list
            for (int k = 0; k < 3; k ++) {
                for (int beta = 0; beta < 3; beta ++) {
                    tripletList.push_back(T(F(fNo,vals(eNo,1)),3*F(fNo,k) + beta, edgeDiff(beta) ));
                    
                    tripletList.push_back(T(F(fNo,vals(eNo,0)),3*F(fNo,k) + beta, -edgeDiff(beta) ));
                }
            }
        }
    }
    
    //Create sparse matrix
    K.setFromTriplets(tripletList.begin(), tripletList.end());
}
