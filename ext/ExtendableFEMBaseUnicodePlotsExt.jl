module ExtendableFEMBaseUnicodePlotsExt

using ExtendableGrids: ExtendableGrid, simplexgrid, ON_CELLS, ON_FACES, ON_EDGES, CellNodes, CellGeometries, BFaceNodes, BFaceRegions, Coordinates
using ExtendableFEMBase: FEVector, FEVectorBlock, FESpace, lazy_interpolate!, num_cells, local_celledgenodes, get_ncomponents, H1P1, IdentityComponent, standard_kernel
using UnicodePlots: UnicodePlots, BrailleCanvas, Plot, heatmap, lineplot,
    lineplot!, lines!

import ExtendableFEMBase: unicode_gridplot, unicode_scalarplot

function unicode_gridplot(
        xgrid::ExtendableGrid;
        title = "gridplot",
        resolution = (40, 20),
        autoscale = true,
        color = (200, 200, 200),
        bface_color = (255, 0, 0),
        CanvasType = BrailleCanvas,
        plot_based = ON_CELLS,   # or ON_FACES/ON_EDGES
        kwargs...,
    )
    Base.depwarn(
        "The unicode extension of ExtendableFEMBase will be deprecated in a future release. A unicode gridplot can now be achieved
         via GridVisualize and UnicodePlots as an backend, e.g., use gridplot(xgrid; Plotter = UnicodePlots) instead of unicode_gridplot(xgrid)",
        :unicode_gridplot; force = true
    )
    coords = xgrid[Coordinates]
    ex = extrema(view(coords, 1, :))
    ey = extrema(view(coords, 2, :))

    if autoscale
        wx = ex[2] - ex[1]
        wy = ey[2] - ey[1]
        rescale = wx / wy * (resolution[1] / (2 * resolution[2]))
        if rescale > 1
            resolution = (resolution[1], Int(ceil(resolution[2] / rescale)))
        else
            resolution = (Int(ceil(resolution[1] / rescale)), resolution[2])
        end
    end

    canvas = CanvasType(
        resolution[2], resolution[1],            # number of rows and columns (characters)
        origin_y = ey[1], origin_x = ex[1],              # position in virtual space
        height = ey[2] - ey[1], width = ex[2] - ex[1]; blend = false, kwargs...
    )    # size of the virtual space

    ## plot all edges
    if plot_based in [ON_FACES, ON_EDGES]
        ## plot all edges via FaceNodes
        facenodes = xgrid[FaceNodes]
        nfaces = size(facenodes, 2)
        for j in 1:nfaces
            lines!(canvas, coords[1, facenodes[1, j]], coords[2, facenodes[1, j]], coords[1, facenodes[2, j]], coords[2, facenodes[2, j]]; color = color)
        end
    elseif plot_based == ON_CELLS
        ## plot all edges via CellNodes and local_celledgenodes
        cellnodes = xgrid[CellNodes]
        cellgeoms = xgrid[CellGeometries]
        ncells = num_cells(xgrid)
        for j in 1:ncells
            cen = local_celledgenodes(cellgeoms[j])
            for k in 1:size(cen, 2)
                lines!(canvas, coords[1, cellnodes[cen[1, k], j]], coords[2, cellnodes[cen[1, k], j]], coords[1, cellnodes[cen[2, k], j]], coords[2, cellnodes[cen[2, k], j]]; color = color)
            end
        end
    end

    ## plot BFaces again and color BFaceRegions
    bfacenodes = xgrid[BFaceNodes]
    bfaceregions = xgrid[BFaceRegions]
    ebfr = extrema(bfaceregions)
    nbfaces = size(bfacenodes, 2)
    for j in 1:nbfaces
        cscale = (bfaceregions[j] - ebfr[1]) / (ebfr[2] - ebfr[1] + 1)
        c = Int.(round.(bface_color .* cscale))
        lines!(canvas, coords[1, bfacenodes[1, j]], coords[2, bfacenodes[1, j]], coords[1, bfacenodes[2, j]], coords[2, bfacenodes[2, j]]; color = c)                 # pixel space
    end
    plot = Plot(canvas; title = title, kwargs...)
    return plot
end


function unicode_scalarplot(
        u::FEVectorBlock;
        components = 1:get_ncomponents(u),
        abs = false,
        nrows = 1,                  # only used for 2D plots
        resolution = (40, 40),
        colormap = :viridis,        # only used for 2D plots
        title = u.name,
        kwargs...
    )
    Base.depwarn(
        "The unicode extension of ExtendableFEMBase will be deprecated in a future release. A unicode scalarplot can now be achieved
         via GridVisualize and UnicodePlots as an backend, e.g., use scalarplot(u; Plotter = UnicodePlots)
         instead of unicode_scalarplot(u)",
        :unicode_scalarplot; force = true
    )

    xgrid = u.FES.xgrid
    coords = xgrid[Coordinates]
    dim = size(coords, 1)

    if dim == 1
        ex = extrema(view(coords, 1, :))
        X = LinRange(ex[1], ex[2], resolution[1])
        xgrid_plot = simplexgrid(X)
    elseif dim == 2
        ex = extrema(view(coords, 1, :))
        ey = extrema(view(coords, 2, :))

        X = LinRange(ex[1], ex[2], resolution[1])
        Y = LinRange(ey[1], ey[2], resolution[2])
        xgrid_plot = simplexgrid(X, Y)
    else
        @warn "sorry, no unicode_scalarplat available for dimension $dim"
        return nothing
    end

    I = [FEVector(FESpace{H1P1{1}}(xgrid_plot)) for c in components]
    if abs
        lazy_interpolate!(I[1][1], [u], [(1, Identity)]; postprocess = (result, input, qpinfo) -> (result[1] = sqrt(sum(view(input, components) .^ 2))), not_in_domain_value = 0)
    else
        for c in 1:length(components)
            lazy_interpolate!(I[c][1], [u], [(1, IdentityComponent{components[c]})]; postprocess = standard_kernel, not_in_domain_value = 0)
        end
    end

    if dim == 1
        plt = nothing
        ylim = extrema(I[1].entries)
        for c in 2:length(components)
            e = extrema(I[c].entries)
            ylim = (min(ylim[1], e[1]), max(ylim[2], e[2]))
        end
        for c in 1:length(components)
            if c == 1
                plt = lineplot(X, view(I[c][1]), ylim = ylim, xlabel = "x", name = title * (length(components) == 1 ? "" : "[$(components[c])]"), height = resolution[2], width = resolution[1])
            else
                lineplot!(plt, X, view(I[c][1]), name = title * "[$(components[c])]")
            end
        end
        return plt
    elseif dim == 2
        plts = [
            heatmap(
                    reshape(view(I[c][1]), (resolution[1], resolution[2]))',
                    xlabel = "x",
                    ylabel = "y",
                    xfact = (ex[2] - ex[1]) / (resolution[1] - 1),
                    yfact = (ey[2] - ey[1]) / (resolution[2] - 1),
                    xoffset = ex[1],
                    yoffset = ey[1],
                    title = title * (length(components) == 1 ? "" : "[$(components[c])]"),
                    colormap = colormap,
                ) for c in 1:length(components)
        ]
        return length(plts) == 1 ? plts[1] : plts
        ## the following command looks nicer, but requires the Term-based extension of UnicodePlots
        ## return UnicodePlots.gridplot(map(i -> plts[i], 1:length(components)); layout = (nrows, nothing))
    end
end

end
