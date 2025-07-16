# FEMatrix

A `FEMatrix` represents a block-structured finite element matrix, where each block (a `FEMatrixBlock`) corresponds to a pair of finite element spaces and operates on a submatrix of a shared `ExtendableSparseMatrix`.

```@autodocs
Modules = [ExtendableFEMBase]
Pages = ["fematrix.jl"]
Order   = [:type, :function]
```
