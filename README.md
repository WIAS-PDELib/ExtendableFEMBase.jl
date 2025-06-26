[![Build status](https://github.com/WIAS-PDELib/ExtendableFEMBase.jl/workflows/linux-macos-windows/badge.svg)](https://github.com/WIAS-PDELib/ExtendableFEMBase.jl/actions)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://wias-pdelib.github.io/ExtendableFEMBase.jl/stable/index.html)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://wias-pdelib.github.io/ExtendableFEMBase.jl/dev/index.html)
[![DOI](https://zenodo.org/badge/667751152.svg)](https://zenodo.org/doi/10.5281/zenodo.10563410)
[![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl)

# ExtendableFEMBase

ExtendableFEMBase.jl provides foundational data structures and tools for assembling custom finite element solvers in Julia. It is designed for flexibility, efficiency, and extensibility, and serves as the low-level backend for [ExtendableFEM.jl](https://github.com/WIAS-PDELib/ExtendableFEM.jl). All functionality is built on top of [ExtendableGrids.jl](https://github.com/WIAS-PDELib/ExtendableGrids.jl).

## Features

- Wide range of finite element types (H1, Hdiv, Hcurl, etc.)
- Flexible finite element spaces (`FESpace`)
- Block-structured matrices and vectors (`FEMatrix`, `FEVector`)
- Primitive and composite function operators (e.g., gradient, divergence)
- Efficient basis evaluation and quadrature routines
- Interpolation and reconstruction operators
