
# FEEvaluator

The `FEEvaluator` provides a unified interface for evaluating finite element basis functions, their derivatives, and related quantities for a given function operator, quadrature rule, and mesh entity. It manages the storage and reuse of basis evaluations both on the reference element (where derivatives are computed via automatic differentiation) and on the current mesh item. The mesh item context can be updated dynamically using the `update!` function.


```@autodocs
Modules = [ExtendableFEMBase]
Pages = ["feevaluator.jl"]
Order   = [:type, :function]
```
