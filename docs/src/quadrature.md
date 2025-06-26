# Quadrature

A quadrature rule consists of a set of points (coordinates of evaluation points in the reference geometry) and associated weights. Constructors are provided for various `AbstractElementGeometries` (from ExtendableGrids) and for different orders; some element types support generic formulas for arbitrary order. See below for a detailed list.

```@autodocs
Modules = [ExtendableFEMBase]
Pages = ["quadrature.jl"]
Order   = [:type, :function]
```

#### Accumulating Vector (internal, for completeness)

Internally, global integration uses an accumulating vector and calls cell-wise integration routines.

```@docs
AccumulatingVector
```
