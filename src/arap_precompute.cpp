#include "arap_precompute.h"

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  // cot matrix
  Eigen::SparseMatrix<double> L;  
  igl::cotmatrix(V,F,L);

  // data
  // from min_quad_with_fixed_precompute code
  Eigen::SparseMatrix<double> Aeq;
  igl::min_quad_with_fixed_precompute(L,b,Aeq,false,data);

  // cotmatrix entries
  // Directly accessible by face idx here
  Eigen::MatrixXd CE;
  igl::cotmatrix_entries(V, F, CE);

  // Build K
  K.resize(V.rows(),3*V.rows());
  // God, this took so much time!

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  tripletList.reserve(F.rows()*3*6);

  for (int f=0; f<F.rows(); f++){
    for (int v=0; v<3; v++){
      // for triangle meshes
      int i = F(f, (v + 1)%3);
      int j = F(f, (v + 2)%3);
      // cot_v corresponds to the opposite half edge

      Eigen::Vector3d eij = CE(f,v) * (V.row(i) - V.row(j));
      eij = eij/3.0;

      for(int k=0; k<3; k++){
        for(int t=0; t<3; t++){
          tripletList.push_back(T(i, 3*F(f,k) + t, eij(t)));
          tripletList.push_back(T(j, 3*F(f,k) + t, -eij(t)));
        }
      }      
    }
  }
  
  K.setFromTriplets(tripletList.begin(), tripletList.end());
}
