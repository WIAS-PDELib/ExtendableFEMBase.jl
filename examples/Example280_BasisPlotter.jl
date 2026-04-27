#=

# 280 : Basis-Plotter
([source code](@__SOURCE_URL__))

This example plots all the basis functions of a H1 finite element on Edge1D or Triangle2D
as unicode plots. This is the result with the default parameters (dim = 1, order = 3):

![](https://github.com/chmerdon/ExtendableFEMBase.jl/blob/master/docs/src/assets/example280.png?raw=true")

=#

module Example280_BasisPlotter

using ExtendableFEMBase
using ExtendableGrids
using GridVisualize
using UnicodePlots
using Term

## everything is wrapped in a main function
function main(; dim = 1, order = 3, Plotter = UnicodePlots)

    ## generate two grids
    @assert dim in [1, 2] "dim must be 1 or 2"
    refgeom = dim == 1 ? Edge1D : Triangle2D
    xgrid = reference_domain(refgeom)

    ## set finite element type and get some information
    FEType = H1Pk{1, dim, order}
    ndofs = get_ndofs(ON_CELLS, FEType, refgeom)
    FEType = H1Pk{ndofs, dim, order}

    ## generate FEVector with ncomponents = ndofs
    ## that will carry one basis function in each component
    FEFunc = FEVector(FESpace{FEType}(xgrid))
    coffsets = ExtendableFEMBase.get_local_coffsets(FEType, ON_CELLS, refgeom)
    for j in 1:ndofs
        FEFunc[1][j + coffsets[j]] = 1
    end

    ## interpolate on finer grid
    xgrid_plot = simplexgrid(0:0.01:1)
    I = FEVector(FESpace{H1P1{ndofs}}(xgrid_plot))
    lazy_interpolate!(I[1], FEFunc, [(1, Identity)])

    ## plot
    nodevals = nodevalues_view(I[1])
    plt = GridVisualizer(; Plotter = Plotter, layout = (1, 1), size = (800, 600))
    colors = [:red, :green, :blue, :white, :yellow, :cyan, :magenta]
    for j in 1:ndofs
        GridVisualize.scalarplot!(plt[1, 1], xgrid_plot, nodevals[j]; Plotter = Plotter, clear = false, limits = (-1, 1.5), label = "dof $j", color = colors[j])
    end
    return reveal(plt)
end
end
