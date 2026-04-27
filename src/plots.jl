"""
    scalarplot!(vis, feVectorBlock::FEVectorBlock, operator = Identity; kwargs...)

A standard scalarplot of (the operator evaluation of) a finite element vector.

All kwargs of the calling method are transferred to the scalarplot in this method.
"""
function GridVisualize.scalarplot!(vis::Union{Nothing, Dict{Symbol, Any}}, feVectorBlock::FEVectorBlock, operator = Identity; abs = false, component = 1, kwargs...)
    return GridVisualize.scalarplot!(vis, feVectorBlock.FES.dofgrid, view(nodevalues(feVectorBlock, operator; abs = abs), component, :); kwargs...)
end
function GridVisualize.scalarplot(feVectorBlock::FEVectorBlock, operator = Identity; abs = false, component = 1, kwargs...)
    return GridVisualize.scalarplot(feVectorBlock.FES.dofgrid, view(nodevalues(feVectorBlock, operator; abs = abs), component, :); kwargs...)
end


"""
    broken_scalarplot!(vis, feVectorBlock::FEVectorBlock, operator = Identity; kwargs...)

A "broken" scalarplot (the operator evaluation of) a broken finite element vector.
Instead of averaging the discontinuous values on the grid nodes, each grid cell is plotted
independently. Thus, a discontinuous plot is generated.

All kwargs of the calling method are transferred to the scalarplot in this method.
"""
function broken_scalarplot!(vis, feVectorBlock::FEVectorBlock, operator = Identity; kwargs...)

    dofgrid = feVectorBlock.FES.dofgrid
    cell_nodes = dofgrid[CellNodes]
    coords = dofgrid[Coordinates]

    all_values = nodevalues(feVectorBlock, operator; cellwise = true) # cellwise evaluation of the FE
    all_coords = @views coords[:, cell_nodes[:]]
    all_cells = reshape(1:length(all_values), size(all_values))

    return GridVisualize.scalarplot!(vis, simplexgrid(all_coords, all_cells, dofgrid[CellRegions]), view(all_values, :); kwargs...)
end
function broken_scalarplot(feVectorBlock::FEVectorBlock, operator = Identity; kwargs...)
    vis = GridVisualizer(; kwargs...)
    broken_scalarplot!(vis, feVectorBlock, operator = operator; kwargs...)
    return reveal(vis)
end

"""
    vectorplot!(vis, feVectorBlock::FEVectorBlock, operator = Identity; kwargs...)

A standard vectorplot of (the operator evaluation of) a finite element vector.

All kwargs of the calling method are transferred to the vectorplot in this method.
"""
function GridVisualize.vectorplot!(p, feVectorBlock::FEVectorBlock, operator = Identity; title = feVectorBlock.name, kwargs...)
    return GridVisualize.vectorplot!(p, feVectorBlock.FES.dofgrid, eval_func_bary(PointEvaluator([(1, operator)], [feVectorBlock])); title = title, kwargs...)
end
function GridVisualize.vectorplot(feVectorBlock::FEVectorBlock, operator = Identity; title = feVectorBlock.name, kwargs...)
    return GridVisualize.vectorplot(feVectorBlock.FES.dofgrid, eval_func_bary(PointEvaluator([(1, operator)], [feVectorBlock])); title = title, kwargs...)
end
