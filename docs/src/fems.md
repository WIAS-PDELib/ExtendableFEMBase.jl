# Implemented Finite Elements

This page provides an overview of the finite element type hierarchy and lists all finite elements currently implemented in ExtendableFEMBase.

## The Finite Element Type Hierarchy

Finite elements in ExtendableFEMBase are organized as leaves in an abstract type hierarchy. The complete type tree is as follows:

```
AbstractFiniteElement
├─ AbstractH1FiniteElement
│  ├─ AbstractH1FiniteElementWithCoefficients
│  │  ├─ H1P1TEB
│  │  └─ H1BR
│  ├─ H1CR
│  ├─ H1MINI
│  ├─ L2P0
│  ├─ L2P1
│  ├─ H1P1
│  ├─ H1P2
│  ├─ H1P2B
│  ├─ H1P3
│  ├─ H1Pk
│  ├─ H1Q1
│  └─ H1Q2
├─ AbstractHcurlFiniteElement
│  ├─ HCURLN0
│  └─ HCURLN1
└─ AbstractHdivFiniteElement
   ├─ HDIVBDM1
   ├─ HDIVBDM2
   ├─ HDIVRT0
   ├─ HDIVRT1
   ├─ HDIVRTk
   └─ HDIVRTkENRICH
```


#### Remarks
- Each finite element type depends on one, two, or three parameters. The first parameter is always the number of components (`ncomponents`), which determines whether the element is scalar- or vector-valued. Some elements also require the parameter `edim <: Int` if their structure differs by spatial dimension. Arbitrary order elements require a third parameter specifying the polynomial order.
- Each finite element provides a set of basis functions in reference coordinates for each applicable `AbstractElementGeometry`, as well as degree-of-freedom (dof) maps for each mesh entity.
- Discontinuous (broken) finite elements can be created using the `broken` switch in the [`FESpace`](@ref) constructor.
- The element type determines how basis functions are transformed from local to global coordinates and how `FunctionOperators` are evaluated.
- Additional continuity properties of element types lead to more specialized basis function sets:
    - `AbstractH1FiniteElement` types provide evaluations of nonzero basis functions on faces/boundary faces.
    - `AbstractHdivFiniteElement` types provide evaluations of nonzero normal fluxes of basis functions on faces/boundary faces.
    - `AbstractHcurlFiniteElement` types provide evaluations of nonzero tangential fluxes of basis functions on edges/boundary edges.
- Each finite element has its own standard interpolation routine `interpolate!` (see [Finite Element Interpolations](@ref)), which can be applied to a function with the signature `function(result, qpinfo)`. The specific interpolation behavior is described for each element below.


## List of Implemented Finite Elements

The following table summarizes all finite elements currently implemented in ExtendableFEMBase and indicates the reference geometries on which they are available. For each entry, the dofmap pattern for cell degrees of freedom is shown in brackets, along with the number of local degrees of freedom for a vector-valued realization. Click on an FEType to view more details.

| FEType | Triangle2D | Parallelogram2D | Tetrahedron3D | Parallelepiped3D |
| :----------------: | :----------------: |  :----------------: |  :----------------: |  :----------------: | 
| AbstractH1FiniteElementWithCoefficients |   |   |   |   |
| [`H1BR`](@ref) | ✓ (N1f1, 9) | ✓ (N1f1, 12) | ✓ (N1f1, 16) |   |
| [`H1P1TEB`](@ref) | ✓ (N1f1, 9) |   | ✓ (N1e1, 18) |   |
| AbstractH1FiniteElement |   |   |   |   |
| [`H1BUBBLE`](@ref) | ✓ (I1, 2) | ✓ (I1, 2) | ✓ (I1, 3) |   |
| [`H1CR`](@ref) | ✓ (F1, 6) | ✓ (F1, 8) | ✓ (F1, 12) |   |
| [`H1MINI`](@ref) | ✓ (N1I1, 8) | ✓ (N1I1, 10) | ✓ (N1I1, 15) |   |
| [`L2P0`](@ref) | ✓ (I1, 2) | ✓ (I1, 2) | ✓ (I1, 3) | ✓ (I1, 3) |
| [`L2P1`](@ref) | ✓ (I3, 6) | ✓ (I3, 6) | ✓ (I4, 12) | ✓ (I4, 12) |
| [`H1P1`](@ref) | ✓ (N1, 6) |  | ✓ (N1, 12) |  |
| [`H1P2`](@ref) | ✓ (N1F1, 12) |  | ✓ (N1E1, 30) |   |
| [`H1P2B`](@ref) | ✓ (N1F1I1, 14) |   |   |   |
| [`H1P3`](@ref) | ✓ (N1F2I1, 20) |   | ✓ (N1E2F1, 60)  |   |
| [`H1Pk`](@ref) | ✓ (order-dep) |   |   |   |
| [`H1Q1`](@ref) | ✓ (N1, 6) | ✓ (N1, 8) | ✓ (N1, 12) | ✓ (N1, 24) |
| [`H1Q2`](@ref) | ✓ (N1F1, 12) | ✓ (N1F1I1, 18) | ✓ (N1E1, 30) |   |
| AbstractHcurlFiniteElement |   |   |   |   |
| [`HCURLN0`](@ref) | ✓ (f1, 3) | ✓ (f1, 4) | ✓ (e1, 6) |   |
| [`HCURLN1`](@ref) | ✓ (f1, 6) |  |  |   |
| AbstractHdivFiniteElement |   |   |   |   |
| [`HDIVBDM1`](@ref) | ✓ (f2, 6) | ✓ (f2, 8) | ✓ (f3, 12) |   |
| [`HDIVBDM2`](@ref) | ✓ (f3i3, 12) |   |   |   |
| [`HDIVRT0`](@ref) | ✓ (f1, 3) | ✓ (f1, 4) | ✓ (f1, 4) | ✓ (f1, 6) |
| [`HDIVRT1`](@ref) | ✓ (f2i2, 8) |   | ✓ (f3i3, 15) |   |
| [`HDIVRTk`](@ref) | ✓ (order-dep) |   |  |   |
| [`HDIVRTkENRICH`](@ref) | ✓ (order-dep) |  | ✓ (order-dep) | |


Note: The dofmap pattern describes how local degrees of freedom are associated with grid entities and provides insight into the continuity properties of the element. Here, "N" or "n" denotes nodes, "F" or "f" denotes faces, "E" or "e" denotes edges, and "I" denotes interior degrees of freedom (i.e., those without continuity across elements). Capital letters indicate that each component has its own degree of freedom, while lowercase letters mean only one degree of freedom is associated with the entity. For example, "N1f1" (as in the Bernardi-Raugel element) means that each node has one dof per component and each face has a single dof. Typically, finite elements involving lowercase letters are only defined for vector-valued cases (i.e., the number of components must match the element dimension), while those with only capital letters are available for any number of components.


## H1-conforming finite elements

### P0 finite element

Piecewise constant finite element that has one degree of freedom on each cell of the grid. (It is masked as a H1-conforming finite element, because it uses the same operator evaluations.)

The interpolation of a given function into this space preserves the cell integrals.

```@docs
L2P0
```

### P1 finite element

The lowest-order Courant finite element that has a degree of freedom on each vertex of the grid. On simplices the
basis functions coincide with the linear barycentric coordinates. Only the L2P1 element is also defined on quads.

The interpolation of a given function into this space performs point evaluations at the nodes.

```@docs
L2P1
H1P1
```

### Q1 finite element

The lowest-order finite element that has a degree of freedom on each vertex of the grid. On simplices the
basis functions coincide with the linear barycentric coordinates. This element is also defined on quads.

The interpolation of a given function into this space performs point evaluations at the nodes.

```@docs
H1Q1
```


### MINI finite element

The mini finite element adds cell bubles to the P1 element that are e.g. beneficial to define inf-sup stable finite element pairs for the Stokes problem.

The interpolation of a given function into this space performs point evaluations at the nodes and preserves its
cell integral.

```@docs
H1MINI
```


### P1TEB finite element

This element adds tangent-weighted edge bubbles to the P1 finite element and therefore is only available as a vector-valued element.

The interpolation of a given function into this space performs point evaluations at the nodes and preserves face integrals of its tangential flux.

```@docs
H1P1TEB
```


### Bernardi-Raugel (BR) finite element

The Bernardi-Raugel adds normal-weighted face bubbles to the P1 finite element and therefore is only available
as a vector-valued element.

The interpolation of a given function into this space performs point evaluations at the nodes and preserves face integrals of its normal flux.

```@docs
H1BR
```


### P2 finite element

The P2 finite element method on simplices equals quadratic polynomials. On the Triangle2D shape the degrees of freedom
are associated with the three vertices and the three faces of the triangle. On the Tetrahedron3D shape the degrees of freedom are associated with the four verties and the six edges.

The interpolation of a given function into this space performs point evaluations at the nodes and preserves its face/edge integrals in 2D/3D.

```@docs
H1P2
```


### Q2 finite element

A second order finite element. On simplices it equals the P2 finite element, and on Quadrilateral2D it has 9 degrees of freedom (vertices, faces and one cell bubble).

The interpolation of a given function into this space performs point evaluations at the nodes and preserves lowest order face moments and (only on quads) also the cell integreal mean.

```@docs
H1Q2
```

### P2B finite element

The P2B finite element adds additional cell bubles (in 2D and 3D) and face bubbles (only in 3D) that are e.g. used to define inf-sup stable finite element pairs for the Stokes problem.

The interpolation of a given function into this space performs point evaluations at the nodes and preserves its cell and face integrals in 2D and also edge integrals in 3D.

```@docs
H1P2B
```

### P3 finite element

The P3 finite element method on simplices equals cubic polynomials. On the Triangle2D shape the degrees of freedom
are associated with the three vertices, the three faces (double dof) of the triangle and the cell itself (one cell bubble).

The interpolation of a given function into this space performs point evaluations at the nodes and preserves cell and face integrals in 2D.

```@docs
H1P3
```

### Pk finite element (experimental)

The Pk finite element method generically generates polynomials of arbitrary order k on simplices (Edge1D, Triangle2D so far).

The interpolation of a given function into this space performs point evaluations at the nodes and preserves cell and face integrals in 2D (moment order depends on the order and the element dimension).

```@docs
H1Pk
```

### Crouzeix-Raviart (CR) finite element

The Crouzeix-Raviart element associates one lowest-order function with each face. On the Triangle2D shape, the basis function of a face is one minus two times the nodal basis function of the opposite node. 

The interpolation of a given function into this space preserves its face integrals.

```@docs
H1CR
```



## Hdiv-conforming finite elements

These Raviart-Thomas and Brezzi-Douglas-Marini finite elements of lower order and their standard interpolations are available:

```@docs
HDIVRT0
HDIVBDM1
HDIVRT1
HDIVBDM2
HDIVRTk
HDIVRTkENRICH
```

## Hcurl-conforming finite elements

So far only the lowest order Nedelec element is available in 2D and 3D. On Triangle2D it has one degree of freedom for each face (i.e. the rotated RT0 element), on Tetrahedron3D it has one degree of freedom associated to each of the six edges.

Its standard interpolation of a given functions preserves its tangential face/edge integrals.

```@docs
HCURLN0
HCURLN1
```


## Incomplete finite elements without approximation power

```@docs
H1BUBBLE
```
