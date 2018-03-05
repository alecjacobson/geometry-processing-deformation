#include "arap_precompute.h"

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
    Eigen::SparseMatrix<double> L, Aeq;
    igl::cotmatrix(V,F,L);
    igl::min_quad_with_fixed_precompute(L,b,Aeq,false,data);
    
    Eigen::MatrixXd C;
    igl::cotmatrix_entries(V,F,C);
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(F.rows()*3*6);
    
    K.resize(V.rows(),3*V.rows());
    for(int f = 0; f < F.rows(); f++){
        for(int k = 0; k < 3; k++){
            int j = (k+1)%3;
            int i = (k+2)%3;
            
            Eigen::RowVector3d e = C(f,k)*(V.row(F(f,i)) - V.row(F(f,j)));
            for(int kk = 0; kk < 3; kk++){
                tripletList.push_back(T(F(f,i),3*F(f,kk),e(0)));
                tripletList.push_back(T(F(f,i),3*F(f,kk)+1,e(1)));
                tripletList.push_back(T(F(f,i),3*F(f,kk)+2,e(2)));
                
                tripletList.push_back(T(F(f,j),3*F(f,kk),-e(0)));
                tripletList.push_back(T(F(f,j),3*F(f,kk)+1,-e(1)));
                tripletList.push_back(T(F(f,j),3*F(f,kk)+2,-e(2)));
            }
        }
    }
    K.setFromTriplets(tripletList.begin(), tripletList.end());
    K/=3.0;
}
