# FEVector

A `FEVector` represents a block-structured vector used in finite element computations. It consists of one or more `FEVectorBlock`s, each associated with a specific `FESpace`. All blocks share a common one-dimensional array that stores the global degrees of freedom (DoFs). Each block can only write to a region of the array specified by offsets, ensuring that each `FESpace` manages its own DoFs within the global vector.


```@autodocs
Modules = [ExtendableFEMBase]
Pages = ["fevector.jl"]
Order   = [:type, :function]
```
