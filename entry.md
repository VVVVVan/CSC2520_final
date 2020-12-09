# Cubic Stylization - final project

This is the final project of course CSC2520: Geometry Processing in Fall 2020. The chosen paper is [Cubic Stylization](https://www.dgp.toronto.edu/projects/cubic-stylization/). Cubic stylization is a 3D stylization tool. The input is a manifold triangle mesh and the output is a cubic triangle mesh. 

## References

> **Paper:** See ["Cubic Stylization"
[Liu, H.T.D. and Jacobson, A., 2019]](https://www.dgp.toronto.edu/projects/cubic-stylization/)
>
> **Source Code:** See [Cubic Stylization](https://github.com/HTDerekLiu/CubicStylization_Cpp)
>
> **Deformation guide:** See [Geometry Processing - Deformation](https://github.com/alecjacobson/geometry-processing-deformation).

## Contents

* `data` contains avaiable meshes
* `images` contains example outputs
* `include` contains all the `.h` files
* `paper` contains the Cubic Stylization paper
* `source` contains all the `.cpp` files
* `entry.md` is for tutorial
* `main.cpp` for tutorial demo
* Final results in [Video](https://youtu.be/xv5lCofoemo)

## How to run the code

1. Cloning the repository with the RECURSIVE option

```
git clone --recursive https://github.com/VVVVVan/CSC2520_final.git
```

(or download the `final.zip` in markus.)

2. compile the application, please type these commands in the terminal

```
mkdir build
cd build
cmake ..
make
```

(if it doesn't compile in maxOS, try change first line of the `CMakeLists.txt` to  

```
cmake_minimum_required(VERSION 3.0)
set(CMAKE_C_COMPILER "/usr/bin/clang")
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
```

)

3. execution by given mesh

```
./cubicstylization [meshName]
```

4. Pressing [space] would give the cubic stylization. Pressing [>]or[.] will increase <img src="https://render.githubusercontent.com/render/math?math=\lambda"> and [<] or[,] will decrease <img src="https://render.githubusercontent.com/render/math?math=\lambda">. A small window also pops up which contains the mesh file name and <img src="https://render.githubusercontent.com/render/math?math=\lambda"> value.
![Example outputs with no mesh inputs](./images/spot.png "Example outputs with no mesh inputs")

## Implement

### Summary

The paper introduces a 3D stylization algorithm that can turn an input shape into the style of a cube. The cubic style sculptures are realized by minimizing the *as-rigid-as-possible* energy with an <img src="https://render.githubusercontent.com/render/math?math=l^1">-regularization or rotated surface normals (Algorithm 1). To accelerate ARAP deformation, they include a hierarchical approach by deforming a low-resolution model and recovering the details back afterward (Algorithm 2). They also expose many variants to incorporate artistic controls.

### What is implemented

For this project, Algorithm 1 is implemented to have the basic function, turning the shape of the triangle mesh to cubic.

**Pseudocode** of Algorithm 1 from [paper](https://www.dgp.toronto.edu/projects/cubic-stylization/):

```c++
// Algorithm 1: 
Cubic Stylization(V,F) {
    U <- V; // U denote tilde(V)
    while not converge do {
        R <- local-step(V, U, lambda);
        U <- global-step(R);
    }
}
```

The implementation can not only turn mesh with and without boundaries into cubic but also turn mesh with boundaries and non-orientable surfaces into cubic.

![Example outputs for ogre which has boundaries](./images/ogre_ori.png "Example outputs for ogre which has boundaries")
![Example outputs for ogre which has boundaries](./images/ogre_cubic.png "Example outputs for ogre which has boundaries")

![Example outputs for non-orientable surfaces](./images/kb_ori.png "Example outputs for non-orientable surfaces")
![Example outputs for non-orientable surfaces](./images/kb_cubic.png "Example outputs for non-orientable surfaces")

### Explanation

We need two steps to solve a *as-rigid-as-possible* problem, local step and global step (details in [Geometry Processing - Deformation](https://github.com/alecjacobson/geometry-processing-deformation)). Thus, the following is the explanation of each step.

#### Pre-process

The pro-process step for the algorithm is to compute all the necessary variable for the equation (equation 1 in [paper](https://www.dgp.toronto.edu/projects/cubic-stylization/)):

![](./images/eq1.png)

(The first term is ASAP term and the second term is CUBENESS. The ASAP term could write as <img src="https://render.githubusercontent.com/render/math?math=tr(V^TLV)"> + <img src="https://render.githubusercontent.com/render/math?math=tr(V^TKR)"> as discussed in [class](https://github.com/alecjacobson/geometry-processing-deformation).)

where  

* <img src="https://render.githubusercontent.com/render/math?math=R_i"> is a 3-by-3 rotation matrix, update in local step below
* <img src="https://render.githubusercontent.com/render/math?math=N_i"> is the "spokes and rims" edges of the <img src="https://render.githubusercontent.com/render/math?math=i">th vertex, compute by getting the face first (by `igl::vertex_triangle_adjacency`) and then vertex in each face
* <img src="https://render.githubusercontent.com/render/math?math=L"> is cotangent discrete Laplacian matrix, compute directly by function `igl::cotmatrix`
* <img src="https://render.githubusercontent.com/render/math?math=K"> is cotangents multiplied against differences across edges in the rest mesh, compute directly by function `arap_rhs` for `igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS` type
* <img src="https://render.githubusercontent.com/render/math?math=w_{ij}"> is the cotangent weight, compute by finding corresponding cotangent discrete for each pair of adjacent edges
* ![](./images/dij.png) is the edge vectors between vertices <img src="https://render.githubusercontent.com/render/math?math=i,j"> , compute by getting vectors between each pair of adjacent edges
* ![](./images/dijt.png) is the edge vectors between vertices <img src="https://render.githubusercontent.com/render/math?math=i,j">  in deformed states, update in local step below
* <img src="https://render.githubusercontent.com/render/math?math=\lambda"> is the parameter that used to control cubeness and initialized as 0.0
* <img src="https://render.githubusercontent.com/render/math?math=a_i"> is the barycentric area of vertex <img src="https://render.githubusercontent.com/render/math?math=i"> , compute by diagonal of `igl::massmatrix`
* <img src="https://render.githubusercontent.com/render/math?math=\hat{n}_i"> is the unit area-weighted normal vector of a vertex <img src="https://render.githubusercontent.com/render/math?math=i"> , compute directly by `igl::per_vertex_normals`

#### Local Step

The local step is to find the optimal <img src="https://render.githubusercontent.com/render/math?math=R_i"> that minimize the energy above, such that (equation 2 in [paper](https://www.dgp.toronto.edu/projects/cubic-stylization/)):

![](./images/eq2.png)

which could be further written as (equation 3 in [paper](https://www.dgp.toronto.edu/projects/cubic-stylization/)):

![](./images/eq3.png)

This equation is a standard ADMM formulation, and could be solved using [scaled-form ADMM updates](https://stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf) (equation 4,5,6,7 in [paper](https://www.dgp.toronto.edu/projects/cubic-stylization/)):

![](./images/eq4.png)

![](./images/eq5.png)

![](./images/eq6.png)

![](./images/eq7.png)

where <img src="https://render.githubusercontent.com/render/math?math=\rho"> is the penality and $u$ is the scaled dual variable.

Before go into optimizing the local step, there are some local-step parameters need to be initilized according to sec 3.1 of the [paper](https://www.dgp.toronto.edu/projects/cubic-stylization/). Initialize <img src="https://render.githubusercontent.com/render/math?math=z">,<img src="https://render.githubusercontent.com/render/math?math=u"> with all zeros by `setZero()`, and ![](./images/pa.png).  

##### Update <img src="https://render.githubusercontent.com/render/math?math=R_i">

![](./images/eq4.png)

The equation could wirte as:

![](./images/eq1-1.png)

And the optimal <img src="https://render.githubusercontent.com/render/math?math=R_i"> = <img src="https://render.githubusercontent.com/render/math?math=V_i"><img src="https://render.githubusercontent.com/render/math?math=U_i^T"> where <img src="https://render.githubusercontent.com/render/math?math=V_i">, <img src="https://render.githubusercontent.com/render/math?math=U_i"> is the singluar valur decomposition of <img src="https://render.githubusercontent.com/render/math?math=M_i"> compute by `Eigen::JacobiSVD`. Changing the sign of the column of <img src="https://render.githubusercontent.com/render/math?math=U_i"> to ensure that det(<img src="https://render.githubusercontent.com/render/math?math=R_i">) is positive.

##### Update <img src="https://render.githubusercontent.com/render/math?math=z">

![](./images/eq5.png)

This equation is an instance of the *lasso* problem and could solve with a *shrinkage* step:

![](./images/eq5-1.png)

Based on [sec 4.4.3 of ADMM paper](https://stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf), the shrinkage could also write as:  

![](./images/eq5-2.png)

which is easy to implement.

##### update <img src="https://render.githubusercontent.com/render/math?math=\tilde{u}">

![](./images/eq6.png)

This is very straightforward and no special methods included.

##### update <img src="https://render.githubusercontent.com/render/math?math=\rho, u">

![](./images/eq7.png)

This `update()` function is based on [sec 3.4.1 of ADMM paper](https://stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf), such that:

![](./images/eq7-1.png)

after that <img src="https://render.githubusercontent.com/render/math?math=u"> is updated accordingly, such that if <img src="https://render.githubusercontent.com/render/math?math=\rho"> increase, <img src="https://render.githubusercontent.com/render/math?math=u"> decrease by <img src="https://render.githubusercontent.com/render/math?math=\tau^{decr}">; if <img src="https://render.githubusercontent.com/render/math?math=\rho"> decrease, <img src="https://render.githubusercontent.com/render/math?math=u"> increase by <img src="https://render.githubusercontent.com/render/math?math=\tau^{incr}">.

****
Hints and Ref:

> * Store the half edges/ pair of neighbor vertexes could help compute the <img src="https://render.githubusercontent.com/render/math?math=\tilde{D_i}"> in local step
> * <img src="https://render.githubusercontent.com/render/math?math=d_{ij}"> is <img src="https://render.githubusercontent.com/render/math?math=vj"> - <img src="https://render.githubusercontent.com/render/math?math=vi">!
> * Need to notice the orientation/size of each vertex during implementing
> * Could get first term of <img src="https://render.githubusercontent.com/render/math?math=W_i"> before going into ADMM loops to save some time
> * Need copy <img src="https://render.githubusercontent.com/render/math?math=z^k"> first before updating <img src="https://render.githubusercontent.com/render/math?math=z">, <img src="https://render.githubusercontent.com/render/math?math=z^k"> needed for <img src="https://render.githubusercontent.com/render/math?math=\rho"> updates  
> * Read [ADMM paper](https://stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf) or [ADMM slide](https://web.stanford.edu/~boyd/papers/pdf/admm_slides.pdf) for more details
> * [Example ADMM lasso code](https://web.stanford.edu/~boyd/papers/admm/lasso/lasso.html)

#### Global Step

The global step is to find the optimal <img src="https://render.githubusercontent.com/render/math?math=V"> that minimize the energy, such that:

![](./images/gs.png)

The global step is easier than local step. Thanks to `igl::min_quad_with_fixed`, it helps prefactorize and solve the Quadratic Energy Minimization problems.

****

> * If <img src="https://render.githubusercontent.com/render/math?math=n"> denotes number of vertex in the mesh: ![](./images/k.png) and ![](./images/R.png). Need resize <img src="https://render.githubusercontent.com/render/math?math=R"> to time correctly and resize <img src="https://render.githubusercontent.com/render/math?math=B"> for updating <img src="https://render.githubusercontent.com/render/math?math=V">.
> * Read [Quadratic Energy Minimization](https://libigl.github.io/tutorial/#laplacian) or review [Geometry Processing - Deformation](https://github.com/alecjacobson/geometry-processing-deformation) for more details.

### In my code

#### `include/cubic_style_data.h`  

Contains all the parameters and variables that is needed to pass through functions like ADMM parameters and matrixes for ASAP.

#### `source/cubic_style_precomputation.cpp`  

Given `V`, `F`, `data`, compute all the pre-processing discussed above, including matrices implementation for local step (i.e `$D_i, N, L, K$`)and the `data` struct employed by igl::min_quad_with_fixed to solve the global step.  

#### `source/cubic_style_single_iteration.cpp`  

Given `U`, `data`, compute the optimization step, including local (scaled-form ADMM updates) and global steps. Minimizing the as-rigid-as-possible plus $l^1$ -regularization energy. Output position of deformed vertexes in `U`.  

#### `main.cpp`

Contains short demo that could change the mesh to cubic by pressing [space] and change <img src="https://render.githubusercontent.com/render/math?math=\lambda"> value by pressing [>,.] or [<,,]. It also pops up a window that contains mesh file name and <img src="https://render.githubusercontent.com/render/math?math=\lambda"> value.

## Future steps

* Try to implement Affine Progressive Meshes
* Adding mouse control for more artistic controls
