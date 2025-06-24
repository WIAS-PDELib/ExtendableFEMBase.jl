# Function Operators

## StandardFunctionOperators

`StandardFunctionOperators` are abstract types that represent fundamental (linear) operators, such as the identity, gradient, divergence, and others. These operators provide a unified interface for evaluating finite element basis functions in various ways.

### List of Primitive Operators

| StandardFunctionOperator                                    | Description                                              | Mathematically                                                            |
| :--------------------------------------------------- | :------------------------------------------------------- | :------------------------------------------------------------------------ |
| Identity                                             | identity                                                 | ``v \rightarrow v``                                                       |
| IdentityComponent{c}                                 | identity of c-th component                               | ``v \rightarrow v_c``                                                     |
| NormalFlux                                           | normal flux (function times normal)                      | ``v \rightarrow v \cdot \vec{n}`` (only ON_FACES)                         |
| TangentFlux                                          | tangent flux (function times tangent)                    | ``v \rightarrow v \cdot \vec{t}`` (only ON_EDGES)                         |
| Gradient                                             | gradient/Jacobian (as a vector)                          | ``v \rightarrow \nabla v``                                                |
| SymmetricGradient                                    | symmetric part of the gradient                           | ``v \rightarrow Voigt(\mathrm{sym}(\nabla v))``                           |
| Divergence                                           | divergence                                               | ``v \rightarrow \mathrm{div}(v) = \nabla \cdot v``                        |
| CurlScalar                                           | curl operator 1D to 2D (rotated gradient)                | ``v \rightarrow [-dv/dx_2,dv/dx_1]``                                      |
| Curl2D                                               | curl operator 2D to 1D                                   | ``v \rightarrow dv_1/dx_2 - dv_2/dx_1``                                   |
| Curl3D                                               | curl operator 3D to 3D                                   | ``v \rightarrow \nabla \times v``                                         |
| Hessian                                              | Hesse matrix = all 2nd order derivatives (as a vector)   | ``v \rightarrow D^2 v``      (e.g. in 2D: xx,xy,yx,yy for each component) |
| SymmetricHessian{a}                                  | symmetric part of Hesse matrix, offdiagonals scaled by a | ``v \rightarrow sym(D^2 v)`` (e.g. in 2D: xx,yy,a*xy for each component)  |
| Laplacian                                            | Laplace Operator (diagonal of Hessian)                   | ``v \rightarrow \Delta v``   (e.g. in 2D: xx,yy for each component)       |



!!! note

    The transformation from the reference domain to the physical domain differs for each finite element class. As a result, the evaluation of each function operator must be implemented specifically for every finite element class. Not all function operators are currently available for every dimension or element type, but new implementations are added as needed or upon request.
    
    Additionally, function operators can be combined with user-defined kernels to postprocess/construct more advanced operators from the available primitives (for example, the deviatoric part of a tensor).


```@autodocs
Modules = [ExtendableFEMBase]
Pages = ["functionoperators.jl"]
Order   = [:type, :function]
```


## ReconstructionOperators

Special operators are provided to evaluate a primitive operator on a reconstructed version of a test function. These are useful for advanced discretizations and post-processing.

```@autodocs
Modules = [ExtendableFEMBase]
Pages = ["reconstructionoperators.jl"]
Order   = [:type, :function]
```

### Divergence-Free Reconstruction Operators

For gradient-robust discretizations of certain classical non-divergence-conforming ansatz spaces, reconstruction operators are available that map a discretely divergence-free H1 function to a pointwise divergence-free H(div) function. Currently, such operators are implemented for the vector-valued Crouzeix-Raviart (H1CR) and Bernardi–Raugel (H1BR) finite element types, as well as for the P2-bubble (H1P2B) element in two dimensions.

**Example:** `Reconst{HDIVRT0{d}, Identity}` reconstructs the identity operator into HDIVRT0, and is available for `H1BR{d}` and `H1CR{d}` for `d = 1, 2`.



## Operator Pairs (Experimental)

Two function operators can be combined into an `OperatorPair`, allowing you to provide two operators in each argument of an assembly pattern. Both operators must be well-defined on the relevant element geometries and finite element spaces, and their actions must be compatible with the input and result fields. This feature is experimental and may have limitations in some cases. An `OperatorTriple` is also available for combining three operators.
