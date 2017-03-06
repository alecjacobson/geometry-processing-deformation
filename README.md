# Geometry Processing – Deformation

> **To get started:** Fork this repository then issue
> 
>     git clone --recursive http://github.com/[username]/geometry-processing-deformation.git
>

## Installation, Layout, and Compilation

See
[introduction](http://github.com/alecjacobson/geometry-processing-introduction).

## Execution

Once built, you can execute the assignment from inside the `build/` by running
on a given mesh:

    ./deformation [path to mesh.obj]

When the mesh is blue, the system is in "place handles" mode. Click on the mesh
to select vertex locations for control point handles. After pressing space to switch to
deformation mode, drag the handles. Pressing `m` will switch between different
deformation methods.

![](images/knight-deformation.gif)

## Background

In this assignment we explore smooth deformation of an existing shape. [Shape
deformation](https://en.wikipedia.org/wiki/Deformation_(mechanics)) has many
applications in geometry process; we will focus on interactive **_handle-based
deformation_**. In this setup, the user repositions a sparse set of points and the
goal is to propagate the transformation at these "handles" to the rest of the
shape. To be interactive we should aim for computing the new deformation of the
shape at around 30 [frames per
second](https://en.wikipedia.org/wiki/Frame_rate).

Shape deformation is the transformation from its _rest shape_ to a
new/current/deformed shape. If the position of a point on some 3D rest shape
 given by $\hat{\x} ∈ \R³$ then we will write that the unknown position
on the deformed shape is given by $\x ∈ \R³$. We can write this point's
displacement vector as $\x - \hat{\x} =: \d ∈ \R³$.

The propagation of the handles' deformation can be thought of in two
complementary ways:

 1. as a [scattered data
 interpolation](https://en.wikipedia.org/wiki/Multivariate_interpolation#Irregular_grid_.28scattered_data.29)
 problem, where handles provide sparse samples of an unknown [displacement
 field](https://en.wikipedia.org/wiki/Displacement_field_(mechanics)), or
 2. as a shape optimization problem, where we try to define a _new_ shape that
 retains the details of the old shape but fulfills the handle constraints.

In the following discussion we will take advantage of the ability to switch
between thinking of the unknowns as the positions of the deformed shape ($\x$)
and displacements ($\d$). These views are equivalent, but often one or the
other provides a better intuitive understanding.

#### Continuity

We will limit ourselves to _continuous_ deformations of shapes. That is, the
shape will not tear, crack or change its topological features. 

If we represent our shape discretely as a [triangle
mesh](https://en.wikipedia.org/wiki/Triangle_mesh) (e.g., with _rest_ vertices
in $\hat{\V}
∈ \R^{n × 3}$ and faces in $F ∈ \{1,…,n\}^{m × 3}$, then we can _trivially_
ensure a continuous deformation by determining new vertex positions $\V$. The
topology (connectivity) of the mesh ($F$) will not change.

#### Generic Distortion Minimization

A rest surface $\hat{S}$
[immersed](https://en.wikipedia.org/wiki/Immersion_(mathematics)) in $\R³$ can
be described as a mapping $\hat{\x}$ from _some_ 2D parameteric domain $Ω$. For any
parameters $u$ and $v$, $\hat{\x}$ describes the 3D position:

\\[
\hat{\x}(u,v) ∈ \R³.
\\]

Similarly the deformed surface can be represented as a position function $\x:
Ω → \R³$. The displacement vector field is thus a function taking the
difference: $\d(u,v) = \x(u,v) - \hat{\x}(u,v)$.

For the handle-based deformation problem we would like to find a new surface
(defined by $\x$) that:

 1. adds as little _distortion_ as possible to the shape, and 
 2. satisfies the users constraints at selected handle positions.

We can cast this as an energy optimization problem. Suppose we have
energy [functional](https://en.wikipedia.org/wiki/Functional_(mathematics))
$E(\x)$ that 
measures the amount of distortion between the new shape ($\x$) and the rest
shape $\hat{x}$, then we could optimize for the best possible shape $\x$ by
minimizing $E$:

\\[
\min_\x E(\x). 
\\]

To ensure that the user's $k$ handle points are interpolated we add the
constraints:

\\[
\text{ subject to } \x(u_i, v_i) = \g_i \ ∀ i = \{1, … , k\},
\\]
where $\g_i$ is the position of the $i$th control point handle.

While the constraints are straightforward, we have many choices for how to
formulate the energy function $E$. A natural choice is to measure distortion in
an egalitarian way by integrating a _local_ measure of distortion at all points
on the surface:

\\[
\min_\x ∫_Ω ‖ e(\x) ‖² \ dA \quad \text{ subject to } \x(u_i, v_i) = \g_i \ ∀ i = \{1, … , k\},
\\]

where $e$ is a vector- or scalar- valued function measuring local (unsquared)
distortion. We will now consider different choices for $e$.

### Linear Methods

If we assume that the deformation between the rest shape given by $\hat{\x}$
and the new shape given by $\x$ is _small_ then we can measure the distortion
of the deformation in terms of the smoothness of the displacement field. This
simplest methods will integrate the magnitude of derivatives of the
displacement field ($\d$): if the displacement field has large variations or
sudden changes then it is inducing a lot of distortion.

#### Gradient-based energy

Let us first consider minimizing the integral of squared variation of the
displacement field:

\\[
\min_\d ∫_Ω ‖ ∇\d ‖_F^2 \ dA \quad \text{ subject to } \d_i =
\g_i-\hat{\x}_i \ ∀ i = \{1, … , k\},
\\]
where 
$∇\d = \left(\begin{array}{ccc}
\frac{∂d^x}{∂u} & \frac{∂d^y}{∂u} & \frac{∂d^z}{∂u} \\
\frac{∂d^x}{∂v} & \frac{∂d^y}{∂v} & \frac{∂d^z}{∂v}
\end{array}\right)$
the [Jacobian](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant)
matrix of the displacement field $\d$ 

> ##### Deformation Gradient
>
> If $\I ∈ \R^{3 × 3}$ is the identity matrix, then the quantity $\F := \I +
> ∇\d$ is referred to as the [deformation
> gradient](https://en.wikipedia.org/wiki/Finite_strain_theory#Deformation_gradient_tensor)
> in the [mechanics](https://en.wikipedia.org/wiki/Continuum_mechanics)
> community.

This is simply the familiar [Dirichlet
energy](https://en.wikipedia.org/wiki/Dirichlet's_energy) applied to each
coordinate function of the displacement field _independently_.

We can discretize this over our triangle mesh surface the same way we have in
smoothing and parameterization assignments:

\\[
\min_{\D} \tr{\D^\transpose \L \D} \quad \text{subject to } \D_\text{handles} =
\g_\text{handles} - \hat{\V}_\text{handles},
\\]
where the rows of $\g_\text{handles} ∈ \R^{k × 3}$ contains the new positions
of the $k$ control point handles.

While easy to implement, this method suffers from a couple immediate problems:

 1. it is not smooth at constraints, and 
 2. then _influence_ of handles dies off too quickly.

![Two regions of selected points (purple) are constrained. The gradient-based
energy produces a _continuous_ displacement but is not _smooth_ at the inner
constrained points (sharp crease). The Laplacian-based energy produces a
displacement with continuous positions and
derivatives.](images/bump-k-harmonic.jpg)

By minimizing the Dirichlet energy, each coordinate of "displacement field" is
a [harmonic function](https://en.wikipedia.org/wiki/Harmonic_function).
Intuitively (however abstractly) we can think of each function as _diffusing_
the user's constraints as if they were
[heat](https://en.wikipedia.org/wiki/Heat_equation) values. As such, each
function diffuses quickly to an average, slowly varying value over most of the
domain. As a displacement field a constant value for each coordinate function
would mean a translation: most of the shape is simply translated.

![](images/bunny-harmonic.gif)

The gradient operator ($∇$) is a [_linear_
operator](https://en.wikipedia.org/wiki/Linear_map). We can alternatively view
our minimization above in terms of the unknown positions $\x$:

\\[
\min_\d ∫_Ω ‖ ∇\d ‖_F^2 \ dA ⇒ 
\min_\x ∫_Ω ‖ ∇(\x - \hat{\x}) ‖_F^2 \ dA ⇒ 
\min_\x ∫_Ω ‖ \underbrace{∇\x}_\text{after} -
\underbrace{∇\hat{\x}}_\text{before} ‖_F^2 \ dA.
\\]
If we think of the gradient of the position function $∇\x$ (with respect to the
underlying parameterization $u,v$) as a local geometric [feature
descriptor](https://en.wikipedia.org/wiki/Feature_(computer_vision)) then this
energy can be re-understood as measuring the difference in this feature before
and after the deformation. This is very sensible as we are trying to measure
distortion. We would expect that a low-distortion deformation would correspond
with a small change to local features.

Unfortunately, the gradient of the position function $\x$ is a _poor_,
first-order local feature.

#### Laplacian-based energy

If we model distortion as the change in a local feature descriptor, then a
natural local and _relative_ descriptor would be one that compared the position
of some point on the shape to the average of its local neighborhood. We have
studied an operator that computes this in the smoothing assignment. The
[Laplace(-Beltrami) operator](https://en.wikipedia.org/wiki/Laplace_operator)
can be derived as taking exactly the difference of a functions value at a point
and the average (i.e., [centroid](https://en.wikipedia.org/wiki/Centroid)) of
an infinitesimal region around that point:

\\[
∆ f(\x) = \lim_{|B(\x)| → 0} \frac{1}{|B(\x))|} ∫_{B(\x)} f(\z) \;d\z - f(\x)
\\]

(see, e.g., "Differential coordinates for local mesh morphing and deformation"
[Alexa et al. 2003] and expanded upon in "Laplacian Surface Editing"
[Sorkine et al. 2004]).

When applied to the embedding function $\x$ the Laplace operator computes the
difference in _position_ between a point and its local neighborhood. This
_vector_ points in the
[normal](https://en.wikipedia.org/wiki/Normal_(geometry)) direction and its
magnitude corresponds inversely with [how flat the surface is
locally](https://en.wikipedia.org/wiki/Mean_curvature).

Let's replace the gradient feature above with this _second-order_ feature
descriptor and massage our optimization problem back in terms of
displacements:

\\[
\min_\x ∫_Ω ‖ \underbrace{∆\x}_\text{after} -
\underbrace{∆\hat{\x}}_\text{before} ‖^2 \ dA ⇒
\min_\x ∫_Ω ‖ ∆(\x - \hat{\x}) ‖^2 \ dA ⇒ 
\min_\d ∫_Ω ‖ ∆\d ‖^2 \ dA.
\\]

Just as we can show that harmonic functions ($∆\d = 0$) minimize the Dirichlet
energy, we can use [calculus of
variations](https://en.wikipedia.org/wiki/Calculus_of_variations) apply
[Green's identity](https://en.wikipedia.org/wiki/Green's_identities) _twice_ to
show that minimizers of the squared-Laplacian energy are _bi-_harmonic
functions ($∆∆\d = 0$ or $∆²\d = 0$). Obviously all harmonic functions are also
biharmonic functions. This implies that the space of biharmonic functions is
strictly _larger_. In particular, this will allow use to ensure continuity of
first derivatives across constrained values at handles.

![](images/bunny-biharmonic.gif)

The fact that each coordinate of the displacement field $\d$ will be a
bi-harmonic function is not so elucidating, but by minimizing the
squared-Laplacian energy integrated over the domain, we can say that it is
as-_harmonic_-as-possible. Harmonic functions include constant functions (i.e.,
translations) but also any displacement that _bends_ in one direction by an
equal and opposite amount as it bends in the other direction: $∆\d = 0 → ∂\d/∂u
+ ∂\d/∂v = 0 → ∂\d/∂u = -∂\d/∂v$.

To discretize this energy we can make use of our discrete Laplacian operator
$\L$. This matrix computes the locally _integrated_ Laplacian of a given
function specified by per-vertex values $\f$. Now we would like to integrate
the square of the _point-wise_ Laplacian. We can approximate the point-wise
Laplacian by the local [integral
average](https://en.wikipedia.org/wiki/Mean_of_a_function) of the Laplacian
$\M^{-1}\L \f$. Integrating this over the mesh we have our complete approximate
of the energy:

\\[
\begin{align}
∫_Ω ‖∆\d‖² \;dA &≈ \tr{ \D^\transpose \L^\transpose \M^{-\transpose} \M \M^{-1} \L
\D}\\
&= \tr{ \D^\transpose \underbrace{\L^\transpose \M^{-1} \L}_{\Q} \D },
\end{align}
\\]
where $\M ∈ \R^{n × n}$ is the mass-matrix for the given mesh and $\Q ∈
\R^{n×n}$ can be thought of as the bi-Laplacian matrix.

> #### $k$-harmonic
> 
> The logical continuation of harmonic and biharmonic deformation is to consider
> _triharmonic_ ($∆³\d = 0$) and _tetraharmonic_ ($∆⁴\d = 0$) and so on. It's
> straightforward to implement these, though there are diminishing returns and
> increasing costs.

##### Precomputation

With out loss of generality, assume that the rows of the unknown displacements
$\D$ have been sorted so that displacements corresponding to handle vertices
are in the bottom part:

\\[
\D = \left(\begin{array}{c}
\D_\text{u} \\
\D_\text{h}
\end{array}
\right)
\\]

Since the displacements at handles are _known_ before the optimization, we can
separate the knowns and unknowns in the energy:

\\[
\min_{\D_\text{u}}
\tr{(\D_\text{u}^\transpose \ \D_\text{h}^\transpose)
\left(\begin{array}{cc}
\Q_\text{u,u} & \Q_\text{u,h} \\
\Q_\text{h,u} & \Q_\text{h,h} 
\end{array}\right)
\left(\begin{array}{c}
  \D_\text{u} \\
  \D_\text{h}
\end{array}
\right)} \\
\min_{\D_\text{u}}
\tr{\D_\text{u}^\transpose \Q_\text{u,u} \D_\text{u} +
2 \D_\text{u}^\transpose \Q_\text{u,h} \D_\text{h} + 
\underbrace{\D_\text{h}^\transpose \Q_\text{h,h}
\D_\text{h}}_\text{constant}} \\
\min_{\D_\text{u}} 
\tr{
\D_\text{u}^\transpose \Q_\text{u,u} \D_\text{u} +
2 \D_\text{u}^\transpose \Q_\text{u,h} \D_\text{h}}
\\]
where $\Q_\text{u,u} ∈ \R^{(n-k) × (n-k)}$ is the quadratic coefficients matrix
corresponding to the unknown displacements.

This quadratic optimization problem may solved by setting all partial
derivatives with respect to degrees of freedom in $\D_\text{u}$ to zero:

\\[
2 \Q_\text{u,u} \D_\text{u} + 2 \Q_\text{u,h} \D_\text{h}
= 0 → \D_\text{u} = \Q_\text{u,u}^{-1} \Q_\text{u,h} \D_\text{h}
\\]

If we don't change _which_ vertices are handles, but only change the positions
of the selected handles, then only $\D_\text{h}$ changes above. In particular,
the matrix $\Q_\text{u,u}$ is unchanged. Therefore, we can
[prefactorize](https://en.wikipedia.org/wiki/Cholesky_decomposition) it so that
_applying its inverse_ is fast (`igl::min_quad_with_fixed` does this for you).

> Actually the entire term $\Q_\text{u,u}^{-1} \Q_\text{u,h} =: \W ∈ \R^{(n-k)
> × k}$ does not change. The columns of $\W$ reveal how unknown displacements
> respond to each handle point's displacement (see also "An intuitive framework
> for real-time freeform modeling" [Botsch & Kobbelt 2004]). This provides a
> gateway to the relationship with linear blend skinning and automatic
> (biharmonic) weighting functions (see "Bounded Biharmonic Weights for
> Real-Time Deformation" [Jacobson et al. 2011]).

##### Trouble in paradise 

Biharmonic displacements work well for small deformations and deformations that
do not imply a large rotation. However, the Laplacian of the position function
$∆\x$ as a feature descriptor is _not_ rotation invariant. This problem is true
for all linear differential features including the gradient of the embedding
function $∇\x$ considered above (see "On linear variational surface deformation
methods" [Botsch & Sorkine 2008]).

![](images/knight-biharmonic-large-rotation.gif)

This means that if the user transforms all of the handle locations by a rigid
transformation $\T$ these energies will _not_ measure zero for a displacement
equivalent to applying the rigid transformation $\T$ to the entire shape. We
would like _global_ rotation invariance, but we would also like this property
to apply _locally_ so that parts of the shape can easily rotate.

### As-rigid-as-possible

In the scenario where each handle
$i$ are perfectly transformed by a [rigid
transformation](https://en.wikipedia.org/wiki/Rigid_transformation) $\x_i =
\Rot \hat{\x}_i + \t$, where $\Rot ∈ SO(3) ⊂ \R^{3×3}$ is a rotation matrix and
$\t∈\R^3$ is a translation vector. If an
[oracle](https://en.wikipedia.org/wiki/Oracle) could only tell us this
particular rigid transformation then we could repair 
the gradient-based energy above by pre-rotating the rest shape by this
transformation:

\\[
\begin{align}
∫_Ω ‖ ∇ \x - ∇(\Rot \hat{\x} + \t) ‖² \;dA 
  &= ∫_Ω ‖ ∇ \x - ∇(\Rot \hat{\x}) - ∇\t‖² \;dA \\
  &= ∫_Ω ‖ ∇ \x - \Rot ∇\hat{\x} ‖² \;dA,
\end{align}
\\]

where the translation vector $\t$ falls out because a translation has constant
gradient.

We do not know the rotation $\Rot$ ahead of time, but we could be as generous
as possible and use the "best" rotation $\Rot ← \argmin_\Rot ∫_Ω ‖ ∇\x -  \Rot
∇\hat{\x} ‖² \;dA$:

\\[
∫_Ω \left\|∇\x - \left( \argmin_\Rot ∫_Ω ‖ ∇\x -  \Rot ∇\hat{\x} ‖² \;dA \right)∇\hat{\x}\right\|² \;dA.
\\]

If we treat $\Rot$ as a degree of freedom along with the unknown positions
$\x$, we can unify this into an optimization over $\x$ and $\Rot$:

\\[
\min_{\x,\Rot∈SO(3)} ∫_Ω \left\|∇\x - \Rot ∇\hat{\x}\right\|² \;dA.
\\]

Optimizing this energy will ensure _global_ rotation invariance. To ensure
_local_ rotation invariance, we can replace $\Rot ∈ SO(3)$ with a spatially
varying _function_ $\Rot : Ω → SO(3)$ that outputs a "best" rotation for any
point on the shape (see "A simple geometric model for elastic deformations"
[Chao et al. 2010]). In this way, the optimal rotation will be locally rigid
everywhere, or _as-rigid-as-possible_ (ARAP).

![](images/knight-arap-large-rotation.gif)

> For embedded solid shapes, we can take the rest shape given by $\hat{\x}$ as
> the parameterization $Ω$, so that $∇\hat{\x} = \I$. This allows us to rewrite
> the as-rigid-as-possible energy as the square of the difference between the
> [deformation gradient](#deformationgradient) and the closest rotation:
>
> \\[
> ∫_Ω \left\|∇\x - \Rot ∇\hat{\x}\right\|² \;dA \\
> ∫_Ω \left\|(∇\x + \I - \I) - \Rot \I \right\|² \;dA \\
> ∫_Ω \left\|(\I + ∇\x - ∇\hat{x}) - \Rot \right\|² \;dA \\
> ∫_Ω \left\|(\I + ∇\d) - \Rot \right\|² \;dA \\
> ∫_Ω \left\|\F - \Rot \right\|² \;dA \\
> \\]
>
> This form provides a bridge between the as-rigid-as-possible energy common in
> geometry processing to _corotated linear elasticity_ used in
> mechanics/physically-based simulation (made explicit in "A simple geometric
> model for elastic deformations" [Chao et al. 2010]). See Section 3.4 of "FEM
> Simulation of 3D Deformable Solids" [Sifakis 2012] for a graphics-mechanics
> perspective.
>

#### Discrete as-rigid-as-possible energy

For a triangle mesh with displacing vertices, the gradient of the embedding
function is constant inside each triangle. In this way we can write the raw
gradient energy above as a double sum over all half-edges $ij$ of all faces $f$
in the mesh:

\\[
½ ∫_Ω ‖ ∇ \x - ∇\hat{\x}‖² \;dA = ½ ∑\limits_{f ∈ F} ∑\limits_{ ij ∈ f} c_{ij} ‖
(\v_i-\v_j) - (\hat{\v}_i-\hat{\v}_j)‖²,
\\]
where $c_{ij}$ is cotangent of the angle opposite half-edge $ij$.

To inject localized best fit rotations, we will assign an unknown rotation
matrix $\Rot_k$ to each vertex $k$ of the mesh and accounts for a third of the
energy integrated over incident triangles:
\\[
½ ∫_Ω ‖ ∇ \x - \Rot ∇\hat{\x}‖² \;dA = 
⅙ ∑\limits_{k=1}^n ∑\limits_{ ij ∈ F(k)} 
c_{ij} ‖ (\v_i-\v_j) - \Rot_k (\hat{\v}_i-\hat{\v}_j)‖²,
\\]
where $F(k)$ is the set of all faces incident on the $k$-th vertex.

#### Optimization

The simplest method for optimizing the ARAP energy is by alternating between

 1. finding the optimal rotations $\Rot_k$ assuming the vertex positions $\V$ are
 fixed, and
 2. finding the optimal vertex positions $\V$ assuming all rotations $\Rot_k$ are
 fixed.

Each rotation $\Rot_k$ only affects the local energy and doesn't interact with
the _other_ rotations. So each can be optimized _locally_. In contrast, the
mesh vertex positions $\V$ depend on each other requiring a _global_ solve. In
the geometry processing literature, this is known as a local-global
optimization (see "As-rigid-as-possible surface modeling" [Sorkine & Alexa
2007]). It is also known as "alternating block coordinate descent" because we
have separated the variables into disjoint sets ($\V,\Rot_1,…,\Rot_n$) and taking
the optimal descent direction for each independently.

Observing the discrete energy above we can see that the energy is quadratic in
$\V$ and quadratic in each $\Rot_k$. Let's start by separating the terms that are
quadratic and linear in $\V$:

\\[
⅙ \underbrace{∑\limits_{k=1}^n ∑\limits_{ ij ∈ F(k)}  c_{ij} (\v_i-\v_j)^\transpose(\v_i-\v_j)}_\text{quadratic}
+
⅙ \underbrace{∑\limits_{k=1}^n ∑\limits_{ ij ∈ F(k)}  c_{ij} (\v_i-\v_j)^\transpose \Rot_k (\hat{\v}_i-\hat{\v}_j)}_\text{linear}
\\]

if we stack the rotation matrices $\Rot_k$ into large matrix $\Rot ∈ \R^{3n ×
3}$ then we can write this energy in matrix form as:

\\[
\tr{ \V^\transpose \L \V } + \tr{ \V^\transpose \K \Rot },
\\]
where $\L ∈ \R^{n × n}$ is the familiar cotangent discrete Laplacian matrix and
$\K ∈ \R^{n × 3n}$ sparse matrix containing cotangents multiplied against
differences across edges in the rest mesh (e.g., $\hat{\v}_i - \hat{\v}_j$).

##### Local step

Minimizing this energy with respect $\R$ corresponds to minimizing:

\\[
\tr{ \underbrace{\V^\transpose \K}_{\C^\transpose} \Rot },
\\]
where $\C ∈ \R^{3n × 3}$ stacks weighted covariance matrices $\C_k ∈ \R^{3 ×
3}$ for each region _covered_ by the corresponding rotation $\Rot_k$. We have
seen this problem before in the registration assignment. For each $\C_k$,
$\Rot_k$ will be the closest rotation matrix solved via [singular value
decomposition](https://en.wikipedia.org/wiki/Singular_value_decomposition).

##### Global step

Minimizing the energy above with respect to $\V$ corresponds to solving a
Dirichlet energy-like minimization problem:

\\[
\min_\V \tr{ \V^\transpose \L \V } + \tr{ \V^\transpose \B }
\\]

where $\K \Rot =: \B ∈ \R^{n × 3}$ is a matrix of rotated vertex gradients.
Adding the handle constraints to the corresponding rows of $\V$ this is easily
minimized by setting all partial derivatives with respect to the unknowns in
$\V$ equal to zero (as in the linear methods above).

##### Implementation

In order to facilitate interactive deformation we would like our local and
global iterations to be computed as quickly as possible. Since the quadratic
form in the global step is the _same_ regardless of the current rotations or
current handle positions, we can
[prefactorize](https://en.wikipedia.org/wiki/Cholesky_decomposition) it (again,
as above). The matrix $\K$ also does not depend on the rotations, current
positions or handle positions, so we can pre-build this matrix. This way
computing $\C$ and $\B$ is a simple matrix-matrix multiplication.

> **Note:** When constructing $\K$ it's easiest to iterate over _all_
> half-edges in the mesh (by iterating over all faces and then each of the
> three edges). Each half-edge $ij$ _contributes_ terms tying $\v_i,\v_j$ to
> _each_ of the (three) rotations $\Rot_k$ that apply against their difference
> (see "Fast Automatic Skinning Transformations" [Jacobson et al. 2012]).

##### Still trouble in paradise

![](images/knight-arap-high-vs-low-resolution.gif)

The as-rigid-as-possible deformation method for surfaces described above has a
number of remaining problems:

 1. like its gradient-based energy ancestor the deformation is not smooth at
 constraints;
 2. the energy punishes _bending_ of the surface--which is good--but does so in
 a way that diminishes as the mesh becomes higher and higher resolution, in
 otherwords, the discrete energy is mesh-resolution dependent; and
 3. the energy is biased by the original combinatorics of the mesh (even in
 flat regions).

## Tasks

### Blacklist

 - `igl::arap`
 - `igl::arap_linear_block`
 - `igl::covariance_scatter_matrix`
 - `igl::harmonic`

### Whitelist

 - `igl::polar_svd3x3` (or your previous assignment's `closest_rotation`)
 - `igl::min_quad_with_fixed` 
 - `igl::cotmatrix_entries` 
 - `igl::cotmatrix`  (or your previous implementation)
 - `igl::massmatrix`  (or your previous implementation)

### `src/biharmonic_precompute.cpp`

Precompute data needed to efficiently solve for a biharmonic deformation given
a mesh with vertices `V` and faces `F` and a list of selected vertices as
indices `b` into `V`. The output should be a prefacorized system using the
`data` struct employed by `igl::min_quad_with_fixed`.

### `src/biharmonic_solve.cpp`

Given precomputation `data` and a list of handle _displacements_ determine
_displacements_ for all vertices in the mesh.

### `src/arap_precompute.cpp`

Precompute data needed to efficiently conduct local-global iterations for an
arap deformation. This includes the `data` struct employed by
`igl::min_quad_with_fixed` to solve the global step" and constructing the
bi-linear form `K` that mixes rotation degrees of freedom with unknown
positions for preparing the covariance matrices of the local step and the
linear term of the global step.

### `src/arap_single_iteration.cpp`

Given precomputed data (`data` and `K`), handle _positions_ `bc` and current
positions of all vertices `U`, conduct a _single_ iteration of the local-global
solver for minimizing the as-rigid-as-possible energy. Output the _positions_
of all vertices of the mesh (by overwriting `U`).
