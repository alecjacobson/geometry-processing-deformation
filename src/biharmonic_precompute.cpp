#include "biharmonic_precompute.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/min_quad_with_fixed.h>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{

  /*
    Can measure distortion of the deformation d = (x - x̑) with linear approximation if d is small.
    Gradient methods give a poor approximation, so instead do a Laplacian based energy measure
    
    min_x ∫ ‖∆d‖² dA

    Our Laplace operator computes the locally integrated Laplacian of the given function specified by per-vertex values f
    Need to integrate the square of the point-wise Laplacian, which can be approximated with M⁻¹Lf

    So our optimization:
    min_x ∫ ‖∆d‖² dA ≈ tr(DᵀLᵀM⁻ᵀMM⁻¹LD) = tr(DᵀLᵀM⁻¹LD) = tr(DᵀQD)

    For which can be split into separate parts, where the handle vertices are separate:

    D = / D_u \
        \ D_h / 
	
    Where after shuffling around the Q and D's arrive at:
    D_u = Q_uu⁻¹ Q_uh D_h

    Which the Q_uu⁻¹ Q_uh do not change unless the handle vertices are changed.
    We will use the min_quad_with_fixed() method to precompute this, and it takes a parameter of known indices of the solution (D_h)
    so the D does not need to be explicitly split up by us
   */

  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
  // need the inverse of the massmatrix, but is sparse and only has entries along the diagonal.
  // can just scale each entry mᵢⱼ = 1/mᵢⱼ which gives an inverse. Eigen has this built in as cwiseInverse() which is nice
  Eigen::SparseMatrix<double> Q = L.transpose() * M.cwiseInverse() * L;
  
  // b are our handle vertices as indices into V
  Eigen::SparseMatrix<double> Aeq(0,0);
  if(!igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data))
    printf("ERROR mqwf precompute\n");
}

