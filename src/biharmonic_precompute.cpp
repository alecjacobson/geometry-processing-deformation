#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include "igl/massmatrix.h"
#include "igl/cotmatrix.h"
#include <iostream>
#include <Eigen/SparseCholesky>

using namespace std;

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  int num_points = V.rows();

  Eigen::SparseMatrix<double> Laplacian;
  igl::cotmatrix(V, F, Laplacian);

  Eigen::SparseMatrix<double> Mass;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, Mass);

  // Invert diagonal Mass matrix
  Eigen::VectorXd diagonal = Mass.diagonal();
  if (not(Mass.nonZeros() == Mass.rows()))
  {
    cout << "Bad assumption that the Mass Matrix is diagonal" << endl;
  }
  Eigen::SparseMatrix<double> Mass_inverse(num_points, num_points); 
  Mass_inverse.reserve(num_points);
  for (int i = 0; i < num_points; i++) 
  {
    Mass_inverse.insert(i, i) = 1.0 / diagonal[i];
  }

  // Aeq is a is a sparse matrix of linear equality constraint coefficients (none)
  Eigen::SparseMatrix<double> Aeq;

  // Laplacian always has at least one 0 eigen value, so it's not PD
  Eigen::SparseMatrix<double> Q = Laplacian.transpose() * Mass_inverse * Laplacian;
  auto positive_definite = false;

  // b is the boundary conditions (given)
  // We shall save the system matrix factorization into "data"
  igl::min_quad_with_fixed_precompute(Q, b, Aeq, positive_definite, data);
}

