#include "arap_precompute.h"
#include <igl/cotmatrix_entries.h>
#include <igl/cotmatrix.h>
#include <igl/min_quad_with_fixed.h>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K) {
    
    // K represents the bilinear form combining unknown vertex positions and
    // unknown rotations
    int n = V.rows();
    K.resize(n, 3*n);
    
    // the factor of 1/3 comes from the fact that the 
    // cotmatrix_entries function populates C with 1/2 cotangent entries but
    // we need 1/6 entries.
    Eigen::MatrixXd cot(F.rows(), 3);
    igl::cotmatrix_entries(V,F,cot);
    cot *= 1.0/3.0;
    
    // initialize a triplet list to be used in populating K.
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(K.rows()*2);
    
    // the below nested for loops are structured to correspond to the summations in the 
    // readme at (27). They loop through each edge on each face and compute the corresponding K
    // entries based on the weighted difference vectors between the incident vertices on the
    // given edge. Each such entry is added to the tripletList along with its corresponding index.
    for (int f = 0; f < F.rows(); f++) {
        for (int e = 0; e < 3; e++) {
            
            // get the indices for the vertices incident on the current
            // edge.
            int i = F(f, (e + 1) % 3);
            int j = F(f, (e + 2) % 3);
            
            for (int v = 0; v < 3; v++) {
                
                // get the index for the current vertex on the current face.
                int k = F(f, v);
                
                // compute the weighted difference vectors.
                double c = cot(f,e);
                Eigen::Vector3d eij = c*(V.row(i) - V.row(j));
                
                // add entries to the triplet list.
                for (int beta = 0; beta < 3; beta++) {
                    tripletList.push_back(T(i, 3*k + beta, eij(beta)));
                    tripletList.push_back(T(j, 3*k + beta, -eij(beta)));
                }
            }
        }
    }
    
    // use tripletList to populate K.
    K.setFromTriplets(tripletList.begin(), tripletList.end());
    
    // populate the data struct to solve the global step by constructing the
    // cotangent Laplacian and passing it to min_quad_with_fixed_precompute.
    Eigen::SparseMatrix<double> L;
    Eigen::SparseMatrix<double> Aeq;
    igl::cotmatrix(V,F,L);
    igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);
}
