# PointEvaluator

Point evaluators provide a convenient interface to evaluate finite element functions (`FEVector`) at arbitrary spatial points, not restricted to mesh nodes or quadrature points. This is useful for post-processing, visualization, or extracting solution values at specific locations.

```@autodocs
Modules = [ExtendableFEMBase]
Pages = ["point_evaluator.jl"]
Order   = [:type, :function]
```
