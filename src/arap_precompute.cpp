#include "arap_precompute.h"
#include <igl/cotmatrix_entries.h>
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  // REPLACE WITH YOUR CODE

	// First, do the factorization on the laplacian matrix, which corresponds to the first term in
	// the global minimization function.

	// Get cotangent laplacian
	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);

	// Do precomputation, same as in the biharmonic case
	Eigen::SparseMatrix<double> Aeq;

	// Factorize
	igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);


	// We know we're going to need the cotangents of angles that appear in the cotengent laplacian.
	// Looks like igl has a function for this...

	Eigen::MatrixXd Cij (F.rows(), 3);
	igl::cotmatrix_entries(V, F, Cij);

	// Initialize triplets that will be used to fill the K matrix later
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;

	std::vector<int> V_set_1 = {1, 2, 0}, V_set_2 = {2, 0, 1};

	// As suggested in the readme, first loop though all faces
	for (int ii = 0; ii < F.rows(); ii++)
	{
		// Extract the vertices on this face
		Eigen::VectorXi V_ii = F.row(ii);

		// Now loop over each edge on this face
		for (int jj = 0; jj < 3; jj++)
		{
			// Extract the vertices of this edge, but in the order that agrees with the order
			// that cotmatrix_entries returns them in. So, first time round we need vertices (1,2),
			// second time we need (2,0), and then finally we need (0,1). To make this easier I made
			// the vectors V_set_1 and V_set_2 above to extract the right vertices easily.
			int Vi = V_ii(V_set_1[jj]);
			int Vj = V_ii(V_set_2[jj]);

			// Compute the difference between these
			Eigen::VectorXd V_diff = V.row(Vi) - V.row(Vj);

			// Extract this extry in the cotangent entries
			double cij = Cij(ii, jj);

			// Compute the product of these two to enter into K. Since the cotangent entries returned
			// are scaled by 1/2 already, we need to further divide by 3 to get the 1/6 factor needed.
			Eigen::VectorXd e = cij*V_diff/3.0;

			// The alpha loop
			for (int alpha = 0; alpha < 3; alpha++)
			{
				// Extract this vertex
				int V_alpha = F(ii, alpha);

				// The beta loop
				for (int beta = 0; beta < 3; beta++)
				{
					// Insert these entries as per readme
					triplets.push_back(T(Vi, 3*V_alpha + beta, e(beta)));
					triplets.push_back(T(Vj, 3*V_alpha + beta, -e(beta)));

					// Note: I'm hoping that if multiple values are inserted into the same K entry
					// via triplets, they'll just add up...potential place to debug
				}
			} 
		}
	}

	// Create K matrix
	K.resize(V.rows(), 3*V.rows());
	K.setFromTriplets(triplets.begin(), triplets.end());

	return;

}
