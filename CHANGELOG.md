# CHANGES

## v1.1.1 June 16, 2025
  - docstring improvements
    
## v1.1.0 May 15, 2025
  - overhaul of interpolators: each finite element
    interpolation is composed of the three new
    structures NodalInterpolator, MomentInterpolator, FunctionalInterpolator
  - improved show functions
  - coffset function for FESpace
  - first, last functions for FEMatrixBlocks and FEVectorBlocks
  - improved submatrix function, submatrix function for FEMatrixBlock
  
## v1.0.0 April 7, 2025
  - doc improvements
  - several bugfixes regarding finite elements defined on subgrids
  - integrate now reacts on regions parameter and has proper docstrings

## v0.8.1 November 5, 2024
  - fixed dimension in nodevals output (xdim was determined wrongly and caused huge arrays)
  - examples in the documentation reactivated

## v0.8 November 1, 2024
  - started changelog
  - fist registered version since repository move (fixed some URL links in the doc and readme)
  - extension for UnicodePlots for the function unicode_gridplot and unicode_scalarplot (without use of Term for now)
  - fixed some qpinfo.item and qpinfo.cell updates in some interpolation routine (thanks to @Da-Be-Ru)

## October 28, 2024

Moved repository from https://github.com/chmerdon/ExtendableFEMBase.jl to https://github.com/WIAS-PDELib/ExtendableFEMBase.jl.
[WIAS-PDELib](https://github.com/WIAS-PDELib/) is a github organization created to collectively manage the Julia packages developed under
the lead of the [WIAS Numerical Mathematics and Scientific Computing](https://wias-berlin.de/research/rgs/fg3)  research group.
According to the [github docs on repository transfer](https://docs.github.com/en/repositories/creating-and-managing-repositories/transferring-a-repository#whats-transferred-with-a-repository),
all links to the previous repository location are automatically redirected to the new location, and all relationships with forks stay intact.
