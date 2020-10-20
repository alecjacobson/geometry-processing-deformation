# Geometry Processing - Deformation

> **To get started:** Clone this repository by issuing
> 
>     git clone --recursive http://github.com/alecjacobson/geometry-processing-deformation.git
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
 given by <img src="./tex/64bd5626ade5f8fe061404694b17f7bf.svg?invert_in_darkmode" align=middle width=48.49298024999999pt height=26.76175259999998pt/> then we will write that the unknown position
on the deformed shape is given by <img src="./tex/6b02f3286d916a27c6682d25c52b7558.svg?invert_in_darkmode" align=middle width=48.49298024999999pt height=26.76175259999998pt/>. We can write this point's
displacement vector as <img src="./tex/e4745e8d5ce6cf4769e509f167a46cf2.svg?invert_in_darkmode" align=middle width=115.54736819999998pt height=26.76175259999998pt/>.

The propagation of the handles' deformation can be thought of in two
complementary ways:

 1. as a [scattered data
 interpolation](https://en.wikipedia.org/wiki/Multivariate_interpolation#Irregular_grid_.28scattered_data.29)
 problem, where handles provide sparse samples of an unknown [displacement
 field](https://en.wikipedia.org/wiki/Displacement_field_(mechanics)), or
 2. as a shape optimization problem, where we try to define a _new_ shape that
 retains the details of the old shape but fulfills the handle constraints.

In the following discussion we will take advantage of the ability to switch
between thinking of the unknowns as the positions of the deformed shape (<img src="./tex/b0ea07dc5c00127344a1cad40467b8de.svg?invert_in_darkmode" align=middle width=9.97711604999999pt height=14.611878600000017pt/>)
and displacements (<img src="./tex/a1a0237b674867e56447da2abd68e758.svg?invert_in_darkmode" align=middle width=10.502226899999991pt height=22.831056599999986pt/>). These views are equivalent, but often one or the
other provides a better intuitive understanding.

#### Continuity

We will limit ourselves to _continuous_ deformations of shapes. That is, the
shape will not tear, crack or change its topological features. 

If we represent our shape discretely as a [triangle
mesh](https://en.wikipedia.org/wiki/Triangle_mesh) (e.g., with _rest_ vertices
in <img src="./tex/d6e5fa2ddb5d55a7af4013b6c29ca3ef.svg?invert_in_darkmode" align=middle width=71.47063604999998pt height=32.14621080000002pt/> and faces in <img src="./tex/96e281ef91b68940c4b79f1ea930a478.svg?invert_in_darkmode" align=middle width=132.49039814999998pt height=26.76175259999998pt/>, then we can _trivially_
ensure a continuous deformation by determining new vertex positions <img src="./tex/26eb59da31fb48cb17abfe4c6dc80375.svg?invert_in_darkmode" align=middle width=14.554737449999989pt height=22.55708729999998pt/>. The
topology (connectivity) of the mesh (<img src="./tex/b8bc815b5e9d5177af01fd4d3d3c2f10.svg?invert_in_darkmode" align=middle width=12.85392569999999pt height=22.465723500000017pt/>) will not change.

#### Generic Distortion Minimization

A rest surface <img src="./tex/72a5310fd675ff6d9c8adc711d650863.svg?invert_in_darkmode" align=middle width=11.44981034999999pt height=32.054807400000016pt/>
[immersed](https://en.wikipedia.org/wiki/Immersion_(mathematics)) in <img src="./tex/fabbfebc4049d77e28eefb36851e7538.svg?invert_in_darkmode" align=middle width=18.424726649999986pt height=26.76175259999998pt/> can
be described as a mapping <img src="./tex/f72190ded391e03cf6edf05c4ad12939.svg?invert_in_darkmode" align=middle width=9.97711604999999pt height=24.200985600000003pt/> from _some_ 2D parametric domain <img src="./tex/9f531c9f3f1ebeef802ced46eabb0336.svg?invert_in_darkmode" align=middle width=11.87217899999999pt height=22.465723500000017pt/>. For any
parameters <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/> and <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/>, <img src="./tex/f72190ded391e03cf6edf05c4ad12939.svg?invert_in_darkmode" align=middle width=9.97711604999999pt height=24.200985600000003pt/> describes the 3D position:

<p align="center"><img src="./tex/3b604eb5b64ec224dc5956482b70e436.svg?invert_in_darkmode" align=middle width=91.9405476pt height=18.312383099999998pt/></p>


Similarly the deformed surface can be represented as a position function <img src="./tex/5bc85ecba8cd06a77697da9524374358.svg?invert_in_darkmode" align=middle width=79.54304325pt height=26.76175259999998pt/>. The displacement vector field is thus a function taking the
difference: <img src="./tex/4c919b4505ccf53ee821116684896780.svg?invert_in_darkmode" align=middle width=186.6435747pt height=24.65753399999998pt/>.

For the handle-based deformation problem we would like to find a new surface
(defined by <img src="./tex/b0ea07dc5c00127344a1cad40467b8de.svg?invert_in_darkmode" align=middle width=9.97711604999999pt height=14.611878600000017pt/>) that:

 1. adds as little _distortion_ as possible to the shape, and 
 2. satisfies the users constraints at selected handle positions.

We can cast this as an energy optimization problem. Suppose we have
energy [functional](https://en.wikipedia.org/wiki/Functional_(mathematics))
<img src="./tex/4a31f2cc615655655447fc20a0124abe.svg?invert_in_darkmode" align=middle width=35.84472704999999pt height=24.65753399999998pt/> that 
measures the amount of distortion between the new shape (<img src="./tex/b0ea07dc5c00127344a1cad40467b8de.svg?invert_in_darkmode" align=middle width=9.97711604999999pt height=14.611878600000017pt/>) and the rest
shape <img src="./tex/dbe3dc5c9d8feca52f6522e1ae39393f.svg?invert_in_darkmode" align=middle width=9.720341399999992pt height=23.744328300000017pt/>, then we could optimize for the best possible shape <img src="./tex/b0ea07dc5c00127344a1cad40467b8de.svg?invert_in_darkmode" align=middle width=9.97711604999999pt height=14.611878600000017pt/> by
minimizing <img src="./tex/84df98c65d88c6adf15d4645ffa25e47.svg?invert_in_darkmode" align=middle width=13.08219659999999pt height=22.465723500000017pt/>:

<p align="center"><img src="./tex/2c6c06e3f37c9dab7f581c9263f8f9de.svg?invert_in_darkmode" align=middle width=70.54795275pt height=22.1917806pt/></p>


To ensure that the user's <img src="./tex/63bb9849783d01d91403bc9a5fea12a2.svg?invert_in_darkmode" align=middle width=9.075367949999992pt height=22.831056599999986pt/> handle points are interpolated we add the
constraints:

<p align="center"><img src="./tex/bc532968cf0518d3dad40ef4d5193e61.svg?invert_in_darkmode" align=middle width=289.51655039999997pt height=16.438356pt/></p>

where <img src="./tex/1f4a081144b46d01e7e016680d6ecb4b.svg?invert_in_darkmode" align=middle width=14.10290474999999pt height=14.611878600000017pt/> is the position of the <img src="./tex/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode" align=middle width=5.663225699999989pt height=21.68300969999999pt/>th control point handle.

While the constraints are straightforward, we have many choices for how to
formulate the energy function <img src="./tex/84df98c65d88c6adf15d4645ffa25e47.svg?invert_in_darkmode" align=middle width=13.08219659999999pt height=22.465723500000017pt/>. A natural choice is to measure distortion in
an egalitarian way by integrating a _local_ measure of distortion at all points
on the surface:

<p align="center"><img src="./tex/1d844cb3b790419cd582211584015405.svg?invert_in_darkmode" align=middle width=444.2427792pt height=37.3519608pt/></p>


where <img src="./tex/8cd34385ed61aca950a6b06d09fb50ac.svg?invert_in_darkmode" align=middle width=7.654137149999991pt height=14.15524440000002pt/> is a vector- or scalar- valued function measuring local (unsquared)
distortion. We will now consider different choices for <img src="./tex/8cd34385ed61aca950a6b06d09fb50ac.svg?invert_in_darkmode" align=middle width=7.654137149999991pt height=14.15524440000002pt/>.

### Linear Methods

If we assume that the deformation between the rest shape given by <img src="./tex/f72190ded391e03cf6edf05c4ad12939.svg?invert_in_darkmode" align=middle width=9.97711604999999pt height=24.200985600000003pt/>
and the new shape given by <img src="./tex/b0ea07dc5c00127344a1cad40467b8de.svg?invert_in_darkmode" align=middle width=9.97711604999999pt height=14.611878600000017pt/> is _small_ then we can measure the distortion
of the deformation in terms of the smoothness of the displacement field. This
simplest methods will integrate the magnitude of derivatives of the
displacement field (<img src="./tex/a1a0237b674867e56447da2abd68e758.svg?invert_in_darkmode" align=middle width=10.502226899999991pt height=22.831056599999986pt/>): if the displacement field has large variations or
sudden changes then it is inducing a lot of distortion.

#### Gradient-based energy

Let us first consider minimizing the integral of squared variation of the
displacement field:

<p align="center"><img src="./tex/70614799c2b109c5fc4ec7338d423e05.svg?invert_in_darkmode" align=middle width=434.70442949999995pt height=37.3519608pt/></p>

where 
<img src="./tex/936d2de13b1622e90e11c12793af5bcf.svg?invert_in_darkmode" align=middle width=196.41939899999997pt height=49.99863329999999pt/>
the [Jacobian](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant)
matrix of the displacement field <img src="./tex/a1a0237b674867e56447da2abd68e758.svg?invert_in_darkmode" align=middle width=10.502226899999991pt height=22.831056599999986pt/> 

> ##### Deformation Gradient
>
> If <img src="./tex/2cea0ab85371ca3fac2a6a872ffcc01e.svg?invert_in_darkmode" align=middle width=62.51135549999999pt height=26.76175259999998pt/> is the identity matrix, then the quantity <img src="./tex/6d3d74299b738dc4a5970460702d795c.svg?invert_in_darkmode" align=middle width=89.83977914999998pt height=22.831056599999986pt/>
> is referred to as the [deformation
> gradient](https://en.wikipedia.org/wiki/Finite_strain_theory#Deformation_gradient_tensor)
> in the [mechanics](https://en.wikipedia.org/wiki/Continuum_mechanics)
> community.

This is simply the familiar [Dirichlet
energy](https://en.wikipedia.org/wiki/Dirichlet's_energy) applied to each
coordinate function of the displacement field _independently_.

We can discretize this over our triangle mesh surface the same way we have in
smoothing and parameterization assignments:

<p align="center"><img src="./tex/a8dde6bf89d57624cd15508bd0cde91d.svg?invert_in_darkmode" align=middle width=420.6178086pt height=26.707794299999996pt/></p>

where the rows of <img src="./tex/eab13ffc1e0b2895a5199d248cd85d44.svg?invert_in_darkmode" align=middle width=109.43735384999998pt height=27.91243950000002pt/> contains the new positions
of the <img src="./tex/63bb9849783d01d91403bc9a5fea12a2.svg?invert_in_darkmode" align=middle width=9.075367949999992pt height=22.831056599999986pt/> control point handles.

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

The gradient operator (<img src="./tex/289bc3c68d369342d98da6ff43ec5b9f.svg?invert_in_darkmode" align=middle width=13.69867124999999pt height=22.465723500000017pt/>) is a [_linear_
operator](https://en.wikipedia.org/wiki/Linear_map). We can alternatively view
our minimization above in terms of the unknown positions <img src="./tex/b0ea07dc5c00127344a1cad40467b8de.svg?invert_in_darkmode" align=middle width=9.97711604999999pt height=14.611878600000017pt/>:

<p align="center"><img src="./tex/339fe5271130cfc4ee90508bd46fbf0c.svg?invert_in_darkmode" align=middle width=552.6760998pt height=44.94072825pt/></p>

If we think of the gradient of the position function <img src="./tex/313c2d91a02a7a93ef69bf70775408a8.svg?invert_in_darkmode" align=middle width=23.67578729999999pt height=22.465723500000017pt/> (with respect to the
underlying parameterization <img src="./tex/cfecde842a36413fb233cf4913fbcb8f.svg?invert_in_darkmode" align=middle width=25.27401689999999pt height=14.15524440000002pt/>) as a local geometric [feature
descriptor](https://en.wikipedia.org/wiki/Feature_(computer_vision)) then this
energy can be re-understood as measuring the difference in this feature before
and after the deformation. This is very sensible as we are trying to measure
distortion. We would expect that a low-distortion deformation would correspond
with a small change to local features.

Unfortunately, the gradient of the position function <img src="./tex/b0ea07dc5c00127344a1cad40467b8de.svg?invert_in_darkmode" align=middle width=9.97711604999999pt height=14.611878600000017pt/> is a _poor_,
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

<p align="center"><img src="./tex/3446b76fb53ce9bc1c4a197b8c95320e.svg?invert_in_darkmode" align=middle width=331.80687704999997pt height=40.2286731pt/></p>


(see, e.g., "Differential coordinates for local mesh morphing and deformation"
[Alexa et al. 2003] and expanded upon in "Laplacian Surface Editing"
[Sorkine et al. 2004]).

When applied to the embedding function <img src="./tex/b0ea07dc5c00127344a1cad40467b8de.svg?invert_in_darkmode" align=middle width=9.97711604999999pt height=14.611878600000017pt/> the Laplace operator computes the
difference in _position_ between a point and its local neighborhood. This
_vector_ points in the
[normal](https://en.wikipedia.org/wiki/Normal_(geometry)) direction and its
magnitude corresponds inversely with [how flat the surface is
locally](https://en.wikipedia.org/wiki/Mean_curvature).

Let's replace the gradient feature above with this _second-order_ feature
descriptor and massage our optimization problem back in terms of
displacements:

<p align="center"><img src="./tex/fb6285e09a41b3c869b994aa16252b0b.svg?invert_in_darkmode" align=middle width=542.01507195pt height=44.94072825pt/></p>


Just as we can show that harmonic functions (<img src="./tex/e1658554850fd2353953938449bf4076.svg?invert_in_darkmode" align=middle width=54.33773894999999pt height=22.831056599999986pt/>) minimize the Dirichlet
energy, we can use [calculus of
variations](https://en.wikipedia.org/wiki/Calculus_of_variations) apply
[Green's identity](https://en.wikipedia.org/wiki/Green's_identities) _twice_ to
show that minimizers of the squared-Laplacian energy are _bi-_harmonic
functions (<img src="./tex/5d67c468fd74d4fe59de53e17f2c5798.svg?invert_in_darkmode" align=middle width=68.03641184999998pt height=22.831056599999986pt/> or <img src="./tex/8846588dfe14c3c5bc1d829dc2394df8.svg?invert_in_darkmode" align=middle width=61.712199449999986pt height=26.76175259999998pt/>). Obviously all harmonic functions are also
biharmonic functions. This implies that the space of biharmonic functions is
strictly _larger_. In particular, this will allow use to ensure continuity of
first derivatives across constrained values at handles.

![](images/bunny-biharmonic.gif)

The fact that each coordinate of the displacement field <img src="./tex/a1a0237b674867e56447da2abd68e758.svg?invert_in_darkmode" align=middle width=10.502226899999991pt height=22.831056599999986pt/> will be a
bi-harmonic function is not so elucidating, but by minimizing the
squared-Laplacian energy integrated over the domain, we can say that it is
as-_harmonic_-as-possible. Harmonic functions include constant functions (i.e.,
translations) but also any displacement that _bends_ in one direction by an
equal and opposite amount as it bends in the other direction: <img src="./tex/dbc21813d9ca34f91bbaedb3ac9e74df.svg?invert_in_darkmode" align=middle width=378.35527619999993pt height=24.65753399999998pt/>.

To discretize this energy we can make use of our discrete Laplacian operator
<img src="./tex/80637df1ca7533740cc7b3fdd1ab540b.svg?invert_in_darkmode" align=middle width=11.36979854999999pt height=22.55708729999998pt/>. This matrix computes the locally _integrated_ Laplacian of a given
function specified by per-vertex values <img src="./tex/47b0192f8f0819d64bce3612c46d15ea.svg?invert_in_darkmode" align=middle width=7.56844769999999pt height=22.831056599999986pt/>. Now we would like to integrate
the square of the _point-wise_ Laplacian. We can approximate the point-wise
Laplacian by the local [integral
average](https://en.wikipedia.org/wiki/Mean_of_a_function) of the Laplacian
<img src="./tex/e36ca5ba76e286a7816020307d4de83d.svg?invert_in_darkmode" align=middle width=54.53182514999998pt height=26.76175259999998pt/>. Integrating this over the mesh we have our complete approximate
of the energy:

<p align="center"><img src="./tex/833473d09a2476a4260b20f306e6bcee.svg?invert_in_darkmode" align=middle width=309.10938465pt height=103.1059425pt/></p>


where <img src="./tex/e81a1a11e934fd54074805470b4f9ad1.svg?invert_in_darkmode" align=middle width=76.43450429999999pt height=26.17730939999998pt/> is the mass-matrix for the given mesh and <img src="./tex/235813b75c42f1089ad90b8356b44c60.svg?invert_in_darkmode" align=middle width=72.69021704999999pt height=26.17730939999998pt/> can be thought of as the bi-Laplacian matrix.

> #### <img src="./tex/63bb9849783d01d91403bc9a5fea12a2.svg?invert_in_darkmode" align=middle width=9.075367949999992pt height=22.831056599999986pt/>-harmonic
> 
> The logical continuation of harmonic and biharmonic deformation is to consider
> _triharmonic_ (<img src="./tex/cf16cb899bde3ad37f07b6760174aa93.svg?invert_in_darkmode" align=middle width=61.712199449999986pt height=26.76175259999998pt/>) and _tetraharmonic_ (<img src="./tex/c5ee46d8fa24a32ddd4d58079f5fac85.svg?invert_in_darkmode" align=middle width=61.712199449999986pt height=26.76175259999998pt/>) and so on. It's
> straightforward to implement these, though there are diminishing returns and
> increasing costs.

##### Precomputation

With out loss of generality, assume that the rows of the unknown displacements
<img src="./tex/17104becada06c6cda0447c33ec6c846.svg?invert_in_darkmode" align=middle width=14.49764249999999pt height=22.55708729999998pt/> have been sorted so that displacements corresponding to handle vertices
are in the bottom part:

<p align="center"><img src="./tex/171e558dac589aea526c441cc3bd5d21.svg?invert_in_darkmode" align=middle width=99.63451740000001pt height=39.452455349999994pt/></p>


Since the displacements at handles are _known_ before the optimization, we can
separate the knowns and unknowns in the energy:

<p align="center"><img src="./tex/ff66603fb510911121e4d20bd0b6a0ae.svg?invert_in_darkmode" align=middle width=327.81992925pt height=39.452455349999994pt/></p>

<p align="center"><img src="./tex/489388061205e7f18d25ea7f6a10b2eb.svg?invert_in_darkmode" align=middle width=361.4211018pt height=60.10065104999999pt/></p>

<p align="center"><img src="./tex/15034b5c7aa725b109bdfdf112f433db.svg?invert_in_darkmode" align=middle width=248.11643835pt height=27.056663699999994pt/></p>

where <img src="./tex/e00f2f4fea4951cb0548cf6343680709.svg?invert_in_darkmode" align=middle width=147.38239559999997pt height=29.190975000000005pt/> is the quadratic coefficients matrix
corresponding to the unknown displacements.

This quadratic optimization problem may solved by setting all partial
derivatives with respect to degrees of freedom in <img src="./tex/bf05505d70cec031bbcd89b8ec055cca.svg?invert_in_darkmode" align=middle width=21.75795434999999pt height=22.55708729999998pt/> to zero:

<p align="center"><img src="./tex/9b59d73088e022e117547943a5dca140.svg?invert_in_darkmode" align=middle width=336.7116192pt height=20.50407645pt/></p>


If we don't change _which_ vertices are handles, but only change the positions
of the selected handles, then only <img src="./tex/7ed12ef03df047be423d80cfadd1a328.svg?invert_in_darkmode" align=middle width=21.75795434999999pt height=22.55708729999998pt/> changes above. In particular,
the matrix <img src="./tex/7af39e007598c4b27ab79ad5cfc84f3d.svg?invert_in_darkmode" align=middle width=32.44294184999999pt height=22.55708729999998pt/> is unchanged. Therefore, we can
[prefactorize](https://en.wikipedia.org/wiki/Cholesky_decomposition) it so that
_applying its inverse_ is fast (`igl::min_quad_with_fixed` does this for you).

> Actually the entire term <img src="./tex/b0363fba948f10fb0aaa23b65a3d9332.svg?invert_in_darkmode" align=middle width=198.2628681pt height=29.190975000000005pt/> does not change. The columns of <img src="./tex/380c103b60c66d6420ec8923cdc6e6e8.svg?invert_in_darkmode" align=middle width=19.80585089999999pt height=22.55708729999998pt/> reveal how unknown displacements
> respond to each handle point's displacement (see also "An intuitive framework
> for real-time freeform modeling" [Botsch & Kobbelt 2004]). This provides a
> gateway to the relationship with linear blend skinning and automatic
> (biharmonic) weighting functions (see "Bounded Biharmonic Weights for
> Real-Time Deformation" [Jacobson et al. 2011]).

##### Trouble in paradise 

Biharmonic displacements work well for small deformations and deformations that
do not imply a large rotation. However, the Laplacian of the position function
<img src="./tex/52cdbba1eb4fddc47169867a1330fa0d.svg?invert_in_darkmode" align=middle width=23.67578729999999pt height=22.465723500000017pt/> as a feature descriptor is _not_ rotation invariant. This problem is true
for all linear differential features including the gradient of the embedding
function <img src="./tex/313c2d91a02a7a93ef69bf70775408a8.svg?invert_in_darkmode" align=middle width=23.67578729999999pt height=22.465723500000017pt/> considered above (see "On linear variational surface deformation
methods" [Botsch & Sorkine 2008]).

![](images/knight-biharmonic-large-rotation.gif)

This means that if the user transforms all of the handle locations by a rigid
transformation <img src="./tex/02f380174e367c8935a57f86907fc7da.svg?invert_in_darkmode" align=middle width=13.15061054999999pt height=22.55708729999998pt/> these energies will _not_ measure zero for a displacement
equivalent to applying the rigid transformation <img src="./tex/02f380174e367c8935a57f86907fc7da.svg?invert_in_darkmode" align=middle width=13.15061054999999pt height=22.55708729999998pt/> to the entire shape. We
would like _global_ rotation invariance, but we would also like this property
to apply _locally_ so that parts of the shape can easily rotate.

### As-rigid-as-possible

In the scenario where each handle
<img src="./tex/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode" align=middle width=5.663225699999989pt height=21.68300969999999pt/> are perfectly transformed by a [rigid
transformation](https://en.wikipedia.org/wiki/Rigid_transformation) <img src="./tex/fd86595e93fdb45c43400aff59c1fe09.svg?invert_in_darkmode" align=middle width=94.43821364999998pt height=24.200985600000003pt/>, where <img src="./tex/43b9b5d6dee02a08443aa04710f133ad.svg?invert_in_darkmode" align=middle width=136.46551874999997pt height=26.76175259999998pt/> is a rotation matrix and
<img src="./tex/c74b6d60aba5c61cba32777b4a37a171.svg?invert_in_darkmode" align=middle width=45.86742269999999pt height=26.76175259999998pt/> is a translation vector. If an
[oracle](https://en.wikipedia.org/wiki/Oracle) could only tell us this
particular rigid transformation then we could repair 
the gradient-based energy above by pre-rotating the rest shape by this
transformation:

<p align="center"><img src="./tex/e01b09619c0a59d53b386a0c42bd096e.svg?invert_in_darkmode" align=middle width=421.99689345pt height=81.279264pt/></p>


where the translation vector <img src="./tex/f40598ec49a99f9a93c399f7dacc6d3e.svg?invert_in_darkmode" align=middle width=7.35155849999999pt height=20.87411699999998pt/> falls out because a translation has constant
gradient.

We do not know the rotation <img src="./tex/6423e0d54c2545769ad013e5f6a4cf94.svg?invert_in_darkmode" align=middle width=14.17800779999999pt height=22.55708729999998pt/> ahead of time, but we could be as generous
as possible and use the "best" rotation <img src="./tex/91f23d9c6124d195c0b48f522fdfed62.svg?invert_in_darkmode" align=middle width=256.22221019999995pt height=26.76175259999998pt/>:

<p align="center"><img src="./tex/79f4938c10172efd2fc2caa43c07d4a4.svg?invert_in_darkmode" align=middle width=380.8330416pt height=42.80410035pt/></p>


If we treat <img src="./tex/6423e0d54c2545769ad013e5f6a4cf94.svg?invert_in_darkmode" align=middle width=14.17800779999999pt height=22.55708729999998pt/> as a degree of freedom along with the unknown positions
<img src="./tex/b0ea07dc5c00127344a1cad40467b8de.svg?invert_in_darkmode" align=middle width=9.97711604999999pt height=14.611878600000017pt/>, we can unify this into an optimization over <img src="./tex/b0ea07dc5c00127344a1cad40467b8de.svg?invert_in_darkmode" align=middle width=9.97711604999999pt height=14.611878600000017pt/> and <img src="./tex/6423e0d54c2545769ad013e5f6a4cf94.svg?invert_in_darkmode" align=middle width=14.17800779999999pt height=22.55708729999998pt/>:

<p align="center"><img src="./tex/aa60ec43b66b731035d119c86b06f4ce.svg?invert_in_darkmode" align=middle width=230.55203655pt height=37.3519608pt/></p>


Optimizing this energy will ensure _global_ rotation invariance. To ensure
_local_ rotation invariance, we can replace <img src="./tex/9ef0dc2c2b529f2bad33e83fb198c711.svg?invert_in_darkmode" align=middle width=79.2965943pt height=24.65753399999998pt/> with a spatially
varying _function_ <img src="./tex/eaff36d2409c381d58217e2c1c5daaf3.svg?invert_in_darkmode" align=middle width=110.34665729999999pt height=24.65753399999998pt/> that outputs a "best" rotation for any
point on the shape (see "A simple geometric model for elastic deformations"
[Chao et al. 2010]). In this way, the optimal rotation will be locally rigid
everywhere, or _as-rigid-as-possible_ (ARAP).

![](images/knight-arap-large-rotation.gif)

> For embedded solid shapes, we can take the rest shape given by <img src="./tex/f72190ded391e03cf6edf05c4ad12939.svg?invert_in_darkmode" align=middle width=9.97711604999999pt height=24.200985600000003pt/> as
> the parameterization <img src="./tex/9f531c9f3f1ebeef802ced46eabb0336.svg?invert_in_darkmode" align=middle width=11.87217899999999pt height=22.465723500000017pt/>, so that <img src="./tex/348265fccc6882ffad052563b1b4d5f7.svg?invert_in_darkmode" align=middle width=52.762342049999994pt height=24.200985600000003pt/>. This allows us to rewrite
> the as-rigid-as-possible energy as the square of the difference between the
> [deformation gradient](#deformationgradient) and the closest rotation:
>
<p align="center"><img src="./tex/d193c12c52f9e359a987f95f7d706bc0.svg?invert_in_darkmode" align=middle width=155.7017814pt height=37.3519608pt/></p>

<p align="center"><img src="./tex/e0c54b5d77de41f88771d5e13b71ecb9.svg?invert_in_darkmode" align=middle width=206.50050629999998pt height=37.3519608pt/></p>

<p align="center"><img src="./tex/dfb06d4e02cd236a3bd3de5f7f9914d1.svg?invert_in_darkmode" align=middle width=215.25633359999998pt height=37.3519608pt/></p>

<p align="center"><img src="./tex/6d62d4838d11477e5b92dd613227424c.svg?invert_in_darkmode" align=middle width=172.5966231pt height=37.3519608pt/></p>

<p align="center"><img src="./tex/31a4b157b288a413347c42fcb9d7abdd.svg?invert_in_darkmode" align=middle width=120.2451129pt height=37.3519608pt/></p>

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
gradient energy above as a double sum over all half-edges <img src="./tex/e5a8bc7bac1dd7d337c9e609a4ae3f99.svg?invert_in_darkmode" align=middle width=13.373644349999989pt height=21.68300969999999pt/> of all faces <img src="./tex/190083ef7a1625fbc75f243cffb9c96d.svg?invert_in_darkmode" align=middle width=9.81741584999999pt height=22.831056599999986pt/>
in the mesh:

<p align="center"><img src="./tex/aa9affd62a894e0e4af8cb56257f3bcb.svg?invert_in_darkmode" align=middle width=445.61633325pt height=44.383808699999996pt/></p>

where <img src="./tex/536cf759f040345e1526f87290d35fdb.svg?invert_in_darkmode" align=middle width=17.86921289999999pt height=14.15524440000002pt/> is cotangent of the angle opposite half-edge <img src="./tex/e5a8bc7bac1dd7d337c9e609a4ae3f99.svg?invert_in_darkmode" align=middle width=13.373644349999989pt height=21.68300969999999pt/>.

To inject localized best fit rotations, we will assign an unknown rotation
matrix <img src="./tex/b45c84b1466bdbc8a65b014b93aeac34.svg?invert_in_darkmode" align=middle width=21.44403689999999pt height=22.55708729999998pt/> to each vertex <img src="./tex/63bb9849783d01d91403bc9a5fea12a2.svg?invert_in_darkmode" align=middle width=9.075367949999992pt height=22.831056599999986pt/> of the mesh and accounts for a third of the
energy integrated over incident triangles:
<p align="center"><img src="./tex/ecfd8636a9708578ca1f9717093e099b.svg?invert_in_darkmode" align=middle width=499.25208795pt height=48.7901106pt/></p>

where <img src="./tex/f9648d8a82a1a6fefe1152ab19eb75e2.svg?invert_in_darkmode" align=middle width=34.714717949999994pt height=24.65753399999998pt/> is the set of all faces incident on the <img src="./tex/63bb9849783d01d91403bc9a5fea12a2.svg?invert_in_darkmode" align=middle width=9.075367949999992pt height=22.831056599999986pt/>-th vertex.

> ##### Where do rotations _live_?
> We have assigned a rotation for each vertex: there are <img src="./tex/55a049b8f161ae7cfeb0197d75aff967.svg?invert_in_darkmode" align=middle width=9.86687624999999pt height=14.15524440000002pt/> rotation matrices
> as auxiliary degrees of freedom. This is in contrast to assigning
> rotations per face (as in "A Local/Global Approach to Mesh Parameterization"
> [Liu et al. 2008]). Per-face--or more generally per-element--rotations work
> well for [codimension](https://en.wikipedia.org/wiki/Codimension) zero
> objects (triangle meshes in <img src="./tex/3177e934cf575c08431076a1a5479ba5.svg?invert_in_darkmode" align=middle width=18.424726649999986pt height=26.76175259999998pt/> or tetrahedral meshes in <img src="./tex/fabbfebc4049d77e28eefb36851e7538.svg?invert_in_darkmode" align=middle width=18.424726649999986pt height=26.76175259999998pt/>). But for
> triangle mesh surfaces in <img src="./tex/fabbfebc4049d77e28eefb36851e7538.svg?invert_in_darkmode" align=middle width=18.424726649999986pt height=26.76175259999998pt/> (i.e., codimension 1), per-face rotations
> would lead to "crumpling") because _bending_ along edges would not be
> measured.

#### Optimization

The simplest method for optimizing the ARAP energy is by alternating between

 1. finding the optimal rotations <img src="./tex/b45c84b1466bdbc8a65b014b93aeac34.svg?invert_in_darkmode" align=middle width=21.44403689999999pt height=22.55708729999998pt/> assuming the vertex positions <img src="./tex/26eb59da31fb48cb17abfe4c6dc80375.svg?invert_in_darkmode" align=middle width=14.554737449999989pt height=22.55708729999998pt/> are
 fixed, and
 2. finding the optimal vertex positions <img src="./tex/26eb59da31fb48cb17abfe4c6dc80375.svg?invert_in_darkmode" align=middle width=14.554737449999989pt height=22.55708729999998pt/> assuming all rotations <img src="./tex/b45c84b1466bdbc8a65b014b93aeac34.svg?invert_in_darkmode" align=middle width=21.44403689999999pt height=22.55708729999998pt/> are
 fixed.

Each rotation <img src="./tex/b45c84b1466bdbc8a65b014b93aeac34.svg?invert_in_darkmode" align=middle width=21.44403689999999pt height=22.55708729999998pt/> only affects the local energy and doesn't interact with
the _other_ rotations. So each can be optimized _locally_. In contrast, the
mesh vertex positions <img src="./tex/26eb59da31fb48cb17abfe4c6dc80375.svg?invert_in_darkmode" align=middle width=14.554737449999989pt height=22.55708729999998pt/> depend on each other requiring a _global_ solve. In
the geometry processing literature, this is known as a local-global
optimization (see "As-rigid-as-possible surface modeling" [Sorkine & Alexa
2007]). It is also known as "alternating block coordinate descent" because we
have separated the variables into disjoint sets (<img src="./tex/4f90daf6c058a218597ba681f5a2ad5d.svg?invert_in_darkmode" align=middle width=102.24653009999999pt height=22.55708729999998pt/>) and taking
the optimal descent direction for each independently.

Observing the discrete energy above we can see that the energy is quadratic in
<img src="./tex/26eb59da31fb48cb17abfe4c6dc80375.svg?invert_in_darkmode" align=middle width=14.554737449999989pt height=22.55708729999998pt/> and quadratic in each <img src="./tex/b45c84b1466bdbc8a65b014b93aeac34.svg?invert_in_darkmode" align=middle width=21.44403689999999pt height=22.55708729999998pt/>. Let's start by separating the terms that are
quadratic and linear in <img src="./tex/26eb59da31fb48cb17abfe4c6dc80375.svg?invert_in_darkmode" align=middle width=14.554737449999989pt height=22.55708729999998pt/>:

<p align="center"><img src="./tex/0c777522918da74b454b740b803648ee.svg?invert_in_darkmode" align=middle width=535.5726469499999pt height=75.23751015pt/></p>


if we stack the rotation matrices <img src="./tex/b45c84b1466bdbc8a65b014b93aeac34.svg?invert_in_darkmode" align=middle width=21.44403689999999pt height=22.55708729999998pt/> into large matrix <img src="./tex/0ebe293c2c02704d6822722036907ba3.svg?invert_in_darkmode" align=middle width=77.64646229999998pt height=26.76175259999998pt/> then we can write this energy in matrix form as:

<p align="center"><img src="./tex/d35f8116cbe53cf1b33de4e1e8b07cc3.svg?invert_in_darkmode" align=middle width=186.20986559999997pt height=20.5316694pt/></p>

where <img src="./tex/e9a4432c2ea64a0aa58762f307040216.svg?invert_in_darkmode" align=middle width=69.85918335pt height=26.17730939999998pt/> is the familiar cotangent discrete Laplacian matrix and
<img src="./tex/7b41f6063ebf3ef4b4158a9d6b70e423.svg?invert_in_darkmode" align=middle width=79.85920799999998pt height=26.76175259999998pt/> sparse matrix containing cotangents multiplied against
differences across edges in the rest mesh (e.g., <img src="./tex/507c868c2c2348c985e69eb4db954e99.svg?invert_in_darkmode" align=middle width=51.622728299999984pt height=24.200985600000003pt/>).

> ###### I'm so confused. What's in the <img src="./tex/558e1b6b0d61666c16dd87622253a301.svg?invert_in_darkmode" align=middle width=14.817277199999989pt height=22.55708729999998pt/> matrix?
> 
> Let's take it slow. The <img src="./tex/558e1b6b0d61666c16dd87622253a301.svg?invert_in_darkmode" align=middle width=14.817277199999989pt height=22.55708729999998pt/> matrix is represents the
> [bilinear form](https://en.wikipedia.org/wiki/Bilinear_form) that combines unknown 
> vertex positions and unknown rotations. We have identified above that we can
> write this in summation or matrix form:
<p align="center"><img src="./tex/07573d67f1dee5f5173bf41452345533.svg?invert_in_darkmode" align=middle width=380.01614969999997pt height=48.7901106pt/></p>

> but how did we get here?
> 
> Let's start with the summation form. The _constants_ of this formula are the
> <img src="./tex/536cf759f040345e1526f87290d35fdb.svg?invert_in_darkmode" align=middle width=17.86921289999999pt height=14.15524440000002pt/> terms and the <img src="./tex/d6a7afdc954cfc35dc55794b86741573.svg?invert_in_darkmode" align=middle width=65.23005719999999pt height=24.65753399999998pt/> terms. Since these always
> appear together, let us merge them into weighted edge difference vectors
> <img src="./tex/c003b2ce58f183fd6cb75b5d5a4d7dbf.svg?invert_in_darkmode" align=middle width=169.16249295pt height=26.76175259999998pt/>:
> 
<p align="center"><img src="./tex/592a0f0ed8bc9a5657f6571ecfe36c04.svg?invert_in_darkmode" align=middle width=211.05717105pt height=53.57849639999999pt/></p>

> the inner term in the summation is an [inner
> product](https://en.wikipedia.org/wiki/Inner_product_space); that is, a
> [scalar](https://en.wikipedia.org/wiki/Scalar_(mathematics)). Let's expose this
> by expanding the matrix-vector products of the inner-product:
<p align="center"><img src="./tex/1dd58d1b9dc13e98774999860ecb15e9.svg?invert_in_darkmode" align=middle width=264.74640555pt height=51.251300099999995pt/></p>

> 
> If our mesh is stored as a vertex list and face list, it's not easy/efficient
> to loop over per-vertex rotations (outer sum) and then over all half-edges of
> incident faces (second sum). Instead, let's rearrange these sums to loop over
> all faces first, then the half-edges of that face, and then over all
> per-vertex rotations that involve this half-edge:
<p align="center"><img src="./tex/f6fefd3e59025fc478024025a416cd49.svg?invert_in_darkmode" align=middle width=359.83764134999996pt height=51.251300099999995pt/></p>

> where the third sum is over all rotations <img src="./tex/63bb9849783d01d91403bc9a5fea12a2.svg?invert_in_darkmode" align=middle width=9.075367949999992pt height=22.831056599999986pt/> such that the half-edge <img src="./tex/e5a8bc7bac1dd7d337c9e609a4ae3f99.svg?invert_in_darkmode" align=middle width=13.373644349999989pt height=21.68300969999999pt/>
> belongs to the half-edges of the faces incident on the <img src="./tex/63bb9849783d01d91403bc9a5fea12a2.svg?invert_in_darkmode" align=middle width=9.075367949999992pt height=22.831056599999986pt/>-th vertex: <img src="./tex/1e243979a1bae020c8e691cf9127b0b0.svg?invert_in_darkmode" align=middle width=94.60652024999999pt height=24.65753399999998pt/>. Well, this means <img src="./tex/63bb9849783d01d91403bc9a5fea12a2.svg?invert_in_darkmode" align=middle width=9.075367949999992pt height=22.831056599999986pt/> can either be <img src="./tex/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode" align=middle width=5.663225699999989pt height=21.68300969999999pt/> or <img src="./tex/36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode" align=middle width=7.710416999999989pt height=21.68300969999999pt/> or the third vertex of
> the <img src="./tex/190083ef7a1625fbc75f243cffb9c96d.svg?invert_in_darkmode" align=middle width=9.81741584999999pt height=22.831056599999986pt/>-th face.
> 
> Now let's turn our attention back to the
> [summand](https://en.wiktionary.org/wiki/summand). The terms indexed by <img src="./tex/dbbd12c1d7f968c7fca71ae001318ee6.svg?invert_in_darkmode" align=middle width=10.57650494999999pt height=14.15524440000002pt/>
> never _mix_. That is, we never add/multiply <img src="./tex/e0f64800fe0ae613446bd42442c44adb.svg?invert_in_darkmode" align=middle width=17.10375479999999pt height=21.839370299999988pt/>, <img src="./tex/670480dfa7e7478306a3a9ed75ddd6e7.svg?invert_in_darkmode" align=middle width=16.170587399999988pt height=25.70766330000001pt/>, and <img src="./tex/1ea4537baf1fe60d2302b062739421c0.svg?invert_in_darkmode" align=middle width=27.04576544999999pt height=31.780732499999996pt/>
> unless <img src="./tex/73c1b1cd55ca71bd8dfa618f76b20f92.svg?invert_in_darkmode" align=middle width=71.76369584999999pt height=22.831056599999986pt/>. This implies that we can write this summation in matrix form
> as:
<p align="center"><img src="./tex/1ed2a049354a064ee9bb1537eae9e54a.svg?invert_in_darkmode" align=middle width=247.16293469999997pt height=18.84197535pt/></p>

> where <img src="./tex/162f945ac5aaf149d6350368558f12e6.svg?invert_in_darkmode" align=middle width=63.74933234999998pt height=22.648391699999998pt/> is <img src="./tex/dbbd12c1d7f968c7fca71ae001318ee6.svg?invert_in_darkmode" align=middle width=10.57650494999999pt height=14.15524440000002pt/>-th column of <img src="./tex/26eb59da31fb48cb17abfe4c6dc80375.svg?invert_in_darkmode" align=middle width=14.554737449999989pt height=22.55708729999998pt/>, <img src="./tex/c314ab46bad6349ca204152ede0c71de.svg?invert_in_darkmode" align=middle width=70.18771979999998pt height=26.76175259999998pt/> is the <img src="./tex/dbbd12c1d7f968c7fca71ae001318ee6.svg?invert_in_darkmode" align=middle width=10.57650494999999pt height=14.15524440000002pt/>-th column of
> <img src="./tex/6423e0d54c2545769ad013e5f6a4cf94.svg?invert_in_darkmode" align=middle width=14.17800779999999pt height=22.55708729999998pt/> and <img src="./tex/9592f8c60998c87a195f1c6740394179.svg?invert_in_darkmode" align=middle width=146.22890865pt height=26.76175259999998pt/> are sparse matrices.
> 
> _Further_, the constant term <img src="./tex/4d0f907770e1712d010d3eee63ba1edd.svg?invert_in_darkmode" align=middle width=18.40954664999999pt height=31.780732499999996pt/> in the summation acts the same on
> <img src="./tex/ea54d53814c759cf48a9e892e9cf6689.svg?invert_in_darkmode" align=middle width=97.92750164999998pt height=31.780732499999996pt/> for any <img src="./tex/dbbd12c1d7f968c7fca71ae001318ee6.svg?invert_in_darkmode" align=middle width=10.57650494999999pt height=14.15524440000002pt/> value. This implies that <img src="./tex/cf1582cc520ea6b124848c981fecddd7.svg?invert_in_darkmode" align=middle width=109.58855819999998pt height=22.55708729999998pt/>,
> so we can reduce the matrix form to:
<p align="center"><img src="./tex/9b77a40d95b72320c406a1b6c870b56e.svg?invert_in_darkmode" align=middle width=76.20982545pt height=16.438356pt/></p>

> 
> Finally, we can answer what is in each entry <img src="./tex/869feffc92ee93a54a3d464c23e77556.svg?invert_in_darkmode" align=middle width=30.76879409999999pt height=22.465723500000017pt/>, where <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/> chooses the
> row and <img src="./tex/83550feba9d733e0c2be45729d81afaa.svg?invert_in_darkmode" align=middle width=81.67979159999999pt height=22.831056599999986pt/> chooses the column of <img src="./tex/558e1b6b0d61666c16dd87622253a301.svg?invert_in_darkmode" align=middle width=14.817277199999989pt height=22.55708729999998pt/>. Treating our nested summation as
> nested for-loops, we can increment entries in <img src="./tex/558e1b6b0d61666c16dd87622253a301.svg?invert_in_darkmode" align=middle width=14.817277199999989pt height=22.55708729999998pt/> 
<p align="center"><img src="./tex/bff89ee5a15e9095c48a95fb4cd53fc4.svg?invert_in_darkmode" align=middle width=129.608259pt height=51.9950244pt/></p>

> for each half-edge encountered. 




##### Local step

Minimizing this energy with respect <img src="./tex/6423e0d54c2545769ad013e5f6a4cf94.svg?invert_in_darkmode" align=middle width=14.17800779999999pt height=22.55708729999998pt/> corresponds to minimizing:

<p align="center"><img src="./tex/4b842f6f84b269bced34145c54e1c9c1.svg?invert_in_darkmode" align=middle width=104.36630489999999pt height=59.1786591pt/></p>

where <img src="./tex/c16082bf5205d2ee536939eb65f66ad6.svg?invert_in_darkmode" align=middle width=77.12134979999999pt height=26.76175259999998pt/> stacks weighted covariance matrices <img src="./tex/f96969c1d8bd033a8b11bed326d758a3.svg?invert_in_darkmode" align=middle width=77.08327109999999pt height=26.76175259999998pt/> for each region _covered_ by the corresponding rotation <img src="./tex/b45c84b1466bdbc8a65b014b93aeac34.svg?invert_in_darkmode" align=middle width=21.44403689999999pt height=22.55708729999998pt/>. We have
seen this problem before in the registration assignment. For each <img src="./tex/386ffa8656300cc181104143b61e47c5.svg?invert_in_darkmode" align=middle width=20.91892439999999pt height=22.55708729999998pt/>,
<img src="./tex/b45c84b1466bdbc8a65b014b93aeac34.svg?invert_in_darkmode" align=middle width=21.44403689999999pt height=22.55708729999998pt/> will be the closest rotation matrix solved via [singular value
decomposition](https://en.wikipedia.org/wiki/Singular_value_decomposition).

##### Global step

Minimizing the energy above with respect to <img src="./tex/26eb59da31fb48cb17abfe4c6dc80375.svg?invert_in_darkmode" align=middle width=14.554737449999989pt height=22.55708729999998pt/> corresponds to solving a
Dirichlet energy-like minimization problem:

<p align="center"><img src="./tex/3b7e5fdbf9086f610c3e542d9bceb3a6.svg?invert_in_darkmode" align=middle width=196.23279225pt height=25.4128281pt/></p>


where <img src="./tex/a0940c9ca54a8591066f8419b39e04e3.svg?invert_in_darkmode" align=middle width=125.84246564999998pt height=26.76175259999998pt/> is a matrix of rotated vertex gradients.
Adding the handle constraints to the corresponding rows of <img src="./tex/26eb59da31fb48cb17abfe4c6dc80375.svg?invert_in_darkmode" align=middle width=14.554737449999989pt height=22.55708729999998pt/> this is easily
minimized by setting all partial derivatives with respect to the unknowns in
<img src="./tex/26eb59da31fb48cb17abfe4c6dc80375.svg?invert_in_darkmode" align=middle width=14.554737449999989pt height=22.55708729999998pt/> equal to zero (as in the linear methods above).

##### Implementation

In order to facilitate interactive deformation we would like our local and
global iterations to be computed as quickly as possible. Since the quadratic
form in the global step is the _same_ regardless of the current rotations or
current handle positions, we can
[prefactorize](https://en.wikipedia.org/wiki/Cholesky_decomposition) it (again,
as above). The matrix <img src="./tex/558e1b6b0d61666c16dd87622253a301.svg?invert_in_darkmode" align=middle width=14.817277199999989pt height=22.55708729999998pt/> also does not depend on the rotations, current
positions or handle positions, so we can pre-build this matrix. This way
computing <img src="./tex/12d3ebda1a212bd89197298f60cf3ce1.svg?invert_in_darkmode" align=middle width=13.652895299999988pt height=22.55708729999998pt/> and <img src="./tex/ff44d867a998c08241beb49b30148782.svg?invert_in_darkmode" align=middle width=13.44741914999999pt height=22.55708729999998pt/> is a simple matrix-matrix multiplication.

> **Note:** When constructing <img src="./tex/558e1b6b0d61666c16dd87622253a301.svg?invert_in_darkmode" align=middle width=14.817277199999989pt height=22.55708729999998pt/> it's easiest to iterate over _all_
> half-edges in the mesh (by iterating over all faces and then each of the
> three edges). Each half-edge <img src="./tex/e5a8bc7bac1dd7d337c9e609a4ae3f99.svg?invert_in_darkmode" align=middle width=13.373644349999989pt height=21.68300969999999pt/> _contributes_ terms tying <img src="./tex/26299117fe9555933d0e92948be82590.svg?invert_in_darkmode" align=middle width=38.83742114999999pt height=14.611878600000017pt/> to
> _each_ of the (three) rotations <img src="./tex/b45c84b1466bdbc8a65b014b93aeac34.svg?invert_in_darkmode" align=middle width=21.44403689999999pt height=22.55708729999998pt/> that apply against their difference
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
bilinear form `K` that mixes rotation degrees of freedom with unknown
positions for preparing the covariance matrices of the local step and the
linear term of the global step.

### `src/arap_single_iteration.cpp`

Given precomputed data (`data` and `K`), handle _positions_ `bc` and current
positions of all vertices `U`, conduct a _single_ iteration of the local-global
solver for minimizing the as-rigid-as-possible energy. Output the _positions_
of all vertices of the mesh (by overwriting `U`).
