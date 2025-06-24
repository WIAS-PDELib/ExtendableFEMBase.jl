
# FESpace

A **finite element space** (FESpace) represents the set of all functions that can be expressed using a particular finite element type on a given computational grid. In ExtendableFEMBase, constructing a finite element space requires only specifying the finite element type and the grid; all necessary degree-of-freedom (dof) mappings are generated automatically on first access. 

```@autodocs
Modules = [ExtendableFEMBase]
Pages = ["finiteelements.jl"]
Order   = [:type, :function]
```

## DofMaps

```@autodocs
Modules = [ExtendableFEMBase]
Pages = ["dofmaps.jl"]
Order   = [:type, :function]
```


The following DofMap subtypes are available and are used as keys to access the dofmap via ```FESpace[DofMap]``` (which is equivalent to ```FESpace.dofmaps[DofMap]```).

| DofMap             | Explanation                                       |
| :----------------: | :------------------------------------------------ | 
| CellDofs           | degrees of freedom for on each cell               | 
| FaceDofs           | degrees of freedom for each face                  | 
| EdgeDofs           | degrees of freedom for each edge (in 3D)          | 
| BFaceDofs          | degrees of freedom for each boundary face         |
| BEdgeDofs          | degrees of freedom for each boundary edge (in 3D) |
