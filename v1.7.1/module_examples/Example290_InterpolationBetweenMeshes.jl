#=

# 290 : Interpolation Between Meshes
([source code](@__SOURCE_URL__))

This example demonstrates the interpolation between meshes feature. Here, we interpolate a function with the P2 element of a coarse triangulation and then interpolate
this P2 function on two uniform refinements into some P1 function. Then, both finite element functions are plotted.

The computed solution for the default parameters looks like this:

![](example290.png)
=#

module Example290_InterpolationBetweenMeshes

using ExtendableFEMBase
using ExtendableGrids
using GridVisualize
using UnicodePlots, Term
using Test #

## function to interpolate
function u!(result, qpinfo)
    x = qpinfo.x
    result[1] = sin(4 * pi * x[1]) * sin(4 * pi * x[2])
    return result[2] = cos(4 * pi * x[1]) * cos(4 * pi * x[2])
end

## everything is wrapped in a main function
function main(; ν = 1.0e-3, nrefs = 3, Plotter = UnicodePlots)

    ## generate two grids
    xgrid1 = uniform_refine(grid_unitsquare(Triangle2D), nrefs)
    xgrid2 = uniform_refine(xgrid1, 3; store_parents = true)

    @show xgrid1 xgrid2

    ## set finite element types for the two grids
    FEType1 = H1Pk{2, 2, 2}
    FEType2 = H1Pk{2, 2, 1}

    ## generate coressponding finite element spaces and FEVectors
    FES1 = FESpace{FEType1}(xgrid1)
    FES2 = FESpace{FEType2}(xgrid2)
    FEFunction1 = FEVector(FES1)
    FEFunction2 = FEVector(FES2)

    ## interpolate function onto first grid
    @time interpolate!(FEFunction1[1], u!)
    @time interpolate!(FEFunction2[1], u!)

    ## interpolate onto other grid
    @time lazy_interpolate!(FEFunction2[1], FEFunction1)
    @time lazy_interpolate!(FEFunction2[1], FEFunction1; use_cellparents = true)

    ## plot
    plt = GridVisualizer(; Plotter = Plotter, layout = (2, 2), clear = true, resolution = (1000, 1000))
    scalarplot!(plt[1, 1], FEFunction1[1], levels = 11, title = "u_h ($FEType1, coarse grid)")
    scalarplot!(plt[1, 2], FEFunction2[1], levels = 11, title = "u_h ($FEType2, fine grid)")
    gridplot!(plt[2, 1], xgrid1, title = "coarse grid", markersize = 0)
    gridplot!(plt[2, 2], xgrid2, title = "fine grid", markersize = 0)
    reveal(plt)
    return plt
end

function generateplots(dir = pwd(); Plotter = nothing, kwargs...)
    plt = main(; Plotter = Plotter, kwargs...)
    scene = GridVisualize.reveal(plt)
    return GridVisualize.save(joinpath(dir, "example290.png"), scene; Plotter = Plotter)
end

## check that the lazy interpolation with and without parent cell search
## gives the same result
function runtests(;
        nrefs = 2,
        nsub = 2
    )
    xgrid1 = uniform_refine(grid_unitsquare(Triangle2D), nrefs)
    xgrid2 = uniform_refine(xgrid1, nsub; store_parents = true)

    FEType1 = H1Pk{2, 2, 2}
    FEType2 = H1Pk{2, 2, 1}

    FES1 = FESpace{FEType1}(xgrid1)
    FES2 = FESpace{FEType2}(xgrid2)

    FEFunction1 = FEVector(FES1)
    interpolate!(FEFunction1[1], u!)

    FEFunction2 = FEVector(FES2)
    FEFunction2_pc = FEVector(FES2)

    lazy_interpolate!(FEFunction2[1], FEFunction1)
    lazy_interpolate!(FEFunction2_pc[1], FEFunction1; use_cellparents = true)

    @info "ndofs = coarse: $(FES1.ndofs), fine: $(FES2.ndofs)"
    @test norm(FEFunction2.entries - FEFunction2_pc.entries, Inf) < 1.0e-12

    ## check that the restriction from the higher to the lower degree element
    ## on the same grid computed via lazy_interpolate! coincides with the
    ## direct nodal interpolation into the lower degree element
    ## (the lower degree nodes are a subset of the higher degree nodes)
    function u1!(result, qpinfo)
        x = qpinfo.x
        return result[1] = sin(4 * pi * x[1]) * sin(4 * pi * x[2])
    end
    FES_p3 = FESpace{H1Pk{1, 2, 3}}(xgrid1)
    FES_p2 = FESpace{H1Pk{1, 2, 2}}(xgrid1)

    V_p3 = FEVector(FES_p3)
    interpolate!(V_p3[1], u1!)

    V_p2_lazy = FEVector(FES_p2)
    lazy_interpolate!(V_p2_lazy[1], V_p3)

    V_p2_direct = FEVector(FES_p2)
    interpolate!(V_p2_direct[1], u1!)

    return @test norm(V_p2_lazy.entries - V_p2_direct.entries, Inf) < 1.0e-12
end #hide

end
