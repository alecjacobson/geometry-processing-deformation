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

## Background

### Linear Methods

Measure distortion with respect to small displacements

#### Gradient-based energy

- Dependent on global coordinate system

- Not smooth at constraints 

- Each coordinate of "displacement field" is a harmonic function (picture of
  arrows and picture of x-coordinate)
  - Each coordinate is a heat diffusion toward some average 
  - As-_constant_-as-possible

#### Laplacian-based energy

- Each coordinate of "displacement field" is a biharmonic function 
  - larger space of functions than harmonic
  - includes functions that are smooth at the constraints
  - As-_harmonic_-as-possible (not simply harmonic since we have to satisfy more
    constraints): locally harmonic roughly corresponds to locally
    "fair"/non-ripply (technically, harmonic means bending in one direction
    equals opposite amount of bending in the other direction. [Saddle]()
    shaped)

- Corresponds to higher-order interpolation à la Hermite splines 

Dependent on local coordinate system, but still confused by large rotations.

#### $k$-harmonic

The logical continuation of harmonic and biharmonic deformation is to consider
triharmonic and tetraharmonic and so on. It's straightforward to extend these,
though there are diminishing returns and increasing costs.

### Non-Linear Methods

Take a step back. Why didn't we like the gradient based energy? Global
coordinate frame, non-smoothness at constraints. Put aside non-smoothness at
constraints. Let's try to tackle large rotations.

**Deformation gradient** for 2D shapes

#### As-rigid-as-possible energy

Warm-up: For 2D shapes

> ##### Volumes

##### Surfaces

Overlapping integration regions

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

### k-harmonic precomputation

### k-harmonic deformation

### covariance scatter matrix

### arap precompuation

### arap deformation
