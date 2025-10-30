# Finite Element Interpolations

## Source functions and QPInfo

The functions that can be interpolated with the methods below are expected to have a certain interface, i.e.:
```julia
function f!(result, qpinfo) end
```
The qpinfo argument communicates vast information of the current quadrature point:

| qpinfo child       | Type               | Description         |
| :----------------  | :----------------  |  :---------------- |
| qpinfo.x           | Vector{Real}       | space coordinates of quadrature point |
| qpinfo.time        | Real               | current time |
| qpinfo.item        | Integer            | current item that contains qpinfo.x |
| qpinfo.cell        | Integer            | cell number (when reasonable) |
| qpinfo.region      | Integer            | region number of item |
| qpinfo.xref        | Vector{Real}       | reference coordinates within item of qpinfo.x |
| qpinfo.volume      | Real               | volume of item |
| qpinfo.normal      | Vector{Real}       | normal vector (when reasonable) |
| qpinfo.params      | Vector{Any}        | parameters that can be transferred via keyword arguments |
| qpinfo.grid        | ExtendableGrid     | full access to grid |


## Standard Interpolations

Each finite element type provides a standard interpolation routine that can be applied to user-defined source functions. By default, interpolation is performed over all cells, but it can also be restricted to faces or edges using an appropriate `AssemblyType`.

```@docs
interpolate!
```

Additionally, you can transfer finite element functions from one grid to another using the `lazy_interpolate!` routine, which interpolates between different meshes.


```@docs
lazy_interpolate!
```

```@docs
compute_interpolation_jacobian
```

The following function continuously interpolates finite element function into a H1Pk space by
point evaluations at the Lagrange nodes of the H1Pk element (averaged over all neighbours).

```@docs
continuify
```

## Nodal Evaluations

Plotting routines require nodal values, i.e., the values of a finite element function at the mesh nodes. The generic `nodevalues!` function evaluates any finite element function at the grid nodes, averaging values if the function is discontinuous. For H1-conforming finite elements and identity evaluations, the `nodevalues_view` function can provide a direct view into the coefficient field, avoiding unnecessary memory allocations.


```@docs
nodevalues!
nodevalues
nodevalues_view
nodevalues_subset!
```



## Displace Mesh

Nodal values (e.g. of a FEVector that discretizes a displacement) can be used to displace the mesh.

```@docs
displace_mesh!
displace_mesh
```
