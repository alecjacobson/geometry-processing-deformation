#include "arap_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>

using namespace std;
void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  // REPLACE WITH YOUR CODE
	Eigen::SparseMatrix<double> L; 
  igl::cotmatrix(V,F,L);

  std::vector< Eigen::Triplet<double> > tripletList;
  tripletList.reserve(F.rows() * 3 * 3 * 2);
  for (int FIdx = 0; FIdx < F.rows(); FIdx++){
  	for (int ii = 0; ii < 3; ii++){
  		int i = F(FIdx, (ii+1) % 3);
  		int j = F(FIdx, (ii+2) % 3);
  		int k = F(FIdx, ii);
  		for (int dim = 0; dim < 3; dim++){
  			tripletList.emplace_back(i, 3*k+dim, L.coeffRef(i,j) * (V(i,dim) - V(j,dim)));
  			tripletList.emplace_back(j, 3*k+dim, L.coeffRef(i,j) * (V(j,dim) - V(i,dim)));
  		}
  	}
  }

  K.resize(V.rows(), V.rows()*3);
  K.setFromTriplets(tripletList.begin(),tripletList.end());

  Eigen::SparseMatrix<double> Aeq; // empty matrix
  igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);
}
