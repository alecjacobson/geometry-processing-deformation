#include "arap_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/cotmatrix_entries.h>
#include <iostream>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  // REPLACE WITH YOUR CODE
  int num_vect = V.rows();;
  std::cout << num_vect << std::endl;
  typedef Eigen::Triplet<double> T;
  Eigen::SparseMatrix<double> Aeq;
  Eigen::SparseMatrix<double> L;
  Eigen::MatrixXd C;
  std::vector<T> tripletList;

  igl::cotmatrix_entries(V, F, C);
  igl::cotmatrix(V, F, L);
  igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);
  std::cout << C.rows() << std::endl;
  std::cout << L.rows() << std::endl;
  for (int i = 0; i < F.rows(); i++){
	   for (int j = 0; j < 3; j++){ 
		    int ka = (j + 1) % 3;
		    int kb = (j + 2) % 3;
	 	    int a = F(i,ka);
	 	    int b = F(i,kb);
        std::cout << a << std::endl;
		    Eigen::RowVector3d e = C(i,j) * (V.row(a)-V.row(b)); 

		    for (int k = 0; k < 3; k++) {
		      tripletList.push_back(T(a, 3 * F(i,0) + k, e(k)));
          tripletList.push_back(T(a, 3 * F(i,1) + k, e(k)));
          tripletList.push_back(T(a, 3 * F(i,2) + k, e(k)));
          tripletList.push_back(T(b, 3 * F(i,0) + k, e(k)));
          tripletList.push_back(T(b, 3 * F(i,1) + k, e(k)));
          tripletList.push_back(T(b, 3 * F(i,2) + k, e(k)));
          }
      }
    }
  K.resize(num_vect, 3 * num_vect);
  K.setFromTriplets(tripletList.begin(), tripletList.end());
}
