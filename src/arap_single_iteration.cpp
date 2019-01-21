#include "arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  // REPLACE WITH YOUR CODE

	// First, compute the big C matrix. By distributing matrix transposes, we know that this
	// is just the tranpose of the K matrix times the current U. This product is probably not sparse.
	Eigen::MatrixXd C_big = K.transpose()*U;

	// Initialize the big R matrix in which we'll later stack each individual rotation
	Eigen::MatrixXd R_big(3*data.n, 3);

	// Loop over each vertex to compute the individual rotation matrices
	for (int ii = 0; ii < data.n; ii++)
	{
		// First, extract the relevant part of the big C matrix. Remember to skip 3x times the
		// rows of the big matrix at each iteration.
		Eigen::Matrix3d C_local = C_big.block(3*ii, 0, 3, 3);

		// Put this into the igl SVD function to get the closest rotation
		Eigen::Matrix3d R_local;
		igl::polar_svd3x3(C_local, R_local);

		// Finally, insert this into the same slot in the big R matrix
		R_big.block(3*ii, 0, 3, 3) = R_local;
	}

	// Now do the matrix-matrix multiplication
	Eigen::MatrixXd B = K*R_big;

	// Solve, similar to the biharmonic case
	Eigen::MatrixXd Beq;

	igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);

	return;

}
