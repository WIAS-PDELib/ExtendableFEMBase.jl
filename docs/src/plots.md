# Plots

## GridVisualize and PlutoVista

Visualization of finite element solutions is possible, for example, via [Nodal Evaluations](@ref) and the plotting routines from
[ExtendableGrids.jl](https://github.com/WIAS-PDELib/ExtendableGrids.jl).
For interactive Pluto notebooks, it is recommended to use [PlutoVista.jl](https://github.com/j-fu/PlutoVista.jl) as the backend for high-quality plots.

## UnicodePlots

For quick, in-terminal visualization, several UnicodePlots-based plotters are available via the `ExtendableFEMBaseUnicodePlotsExt` extension.
This extension is loaded automatically when UnicodePlots is available.

```@docs
unicode_gridplot
unicode_scalarplot
```
