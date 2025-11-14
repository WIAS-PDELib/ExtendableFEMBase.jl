"""
````
function lazy_interpolate!(
        target::FEVectorBlock{T1, Tv, Ti},
        source,
        operators = [(1, Identity)];
        postprocess! = standard_kernel,
        xtrafo! = nothing,
        items = [],
        resultdim = get_ncomponents(eltype(target.FES)),
        not_in_domain_value = 1.0e30,
        start_cell = 1,
        only_localsearch = false,
        use_cellparents::Bool = false,
        eps = 1.0e-13,
        kwargs...) where {T1, Tv, Ti}
````

Interpolates (operator-evaluations of) the given FEVector source (or an array of FEVectorBlocks)
into the finite element space assigned to the `target` FEVectorBlock.

The interpolation is performed using a point evaluation pattern and cell search. If `CellParents` information
is available in the target grid, enabling `use_cellparents=true` can improve the efficiency of the search.

A custom postprocessing function can be provided via the `postprocess!` argument, which should have the interface:
    postprocess!(result, input, qpinfo)
where `result` is the output buffer, `input` is the operator evaluation, and `qpinfo` provides quadrature point information.

If the source and target grids have different coordinate dimensions, a coordinate transformation function
`xtrafo!` must be provided, with the interface:
    xtrafo!(x_source, x)
which maps coordinates `x` from the target grid to coordinates in the source grid.

If a point cannot be found in the source grid, the value `not_in_domain_value` is used as the function value.
The `items` argument can be used to restrict the interpolation to specific target cells.

# Arguments
- `target::FEVectorBlock`: The target finite element vector block.
- `source`: The source array of FEVectorBlocks.
- `operators`: Array of operator argument tuples (source block tag, operator type)  (default: `[(1, Identity)]`).

# Keyword Arguments
- `postprocess!`: Function to postprocess operator evaluations (default: `standard_kernel`).
- `xtrafo!`: Optional coordinate transformation function (default: `nothing`).
- `items`: List of target cells to interpolate (default: `[]`).
- `resultdim`: Result dimension (default: `get_ncomponents(eltype(target.FES))`).
- `not_in_domain_value`: Value assigned if a point is not found in the source domain (default: `1.0e30`).
- `start_cell`: Starting cell index for cell search (default: `1`).
- `only_localsearch`: Restrict cell search to local neighborhood (default: `false`).
- `use_cellparents`: Use parent cell information for search (default: `false`).
- `eps`: Tolerance for cell search (default: `1.0e-13`).
- `kwargs...`: Additional keyword arguments passed to `interpolate!`.

# Notes
- Discontinuous quantities at target grid vertices are evaluated in the first found cell of the source grid; no averaging is performed.
- The function is not the most efficient for large-scale problems due to its reliance on pointwise cell search.

"""
function lazy_interpolate!(
        target::FEVectorBlock{T1, Tv, Ti},
        source,
        operators = [(1, Identity)];
        postprocess = standard_kernel,
        xtrafo = nothing,
        items = [],
        resultdim = get_ncomponents(eltype(target.FES)),
        not_in_domain_value = 1.0e30,
        start_cell = 1,
        only_localsearch = false,
        use_cellparents::Bool = false,
        eps = 1.0e-13,
        kwargs...
    ) where {T1, Tv, Ti}

    # wrap point evaluation into function that is put into normal interpolate!
    xgrid = source[1].FES.xgrid
    xdim_source::Int = size(xgrid[Coordinates], 1)
    xdim_target::Int = size(target.FES.xgrid[Coordinates], 1)
    if xdim_source != xdim_target
        @assert xtrafo !== nothing "grids have different coordinate dimensions, need xtrafo!"
    end
    PE = PointEvaluator(postprocess, operators, source; TCoeff = T1)
    xref = zeros(Tv, xdim_source)
    x_source = zeros(Tv, xdim_source)
    cell::Int = start_cell
    lastnonzerocell::Int = start_cell
    same_cells::Bool = xgrid == target.FES.xgrid
    CF::CellFinder{Tv, Ti} = CellFinder(xgrid)

    if same_cells || use_cellparents
        xCellParents::Array{Ti, 1} = same_cells ? (1:num_cells(target.FES.xgrid)) : target.FES.xgrid[CellParents]
        #@show xCellParents
        function point_evaluation_parentgrid!(result, qpinfo)
            x = xtrafo !== nothing ? xtrafo(x_source, qpinfo.x) : qpinfo.x
            #@show qpinfo.x qpinfo.cell
            cell = gFindLocal!(xref, CF, x; icellstart = xCellParents[qpinfo.cell], eps = eps, trybrute = !only_localsearch)
            return evaluate_bary!(result, PE, xref, cell)
        end
        fe_function = point_evaluation_parentgrid!
    else
        function point_evaluation_arbitrarygrids!(result, qpinfo)
            x = xtrafo !== nothing ? xtrafo(x_source, qpinfo.x) : qpinfo.x
            cell = gFindLocal!(xref, CF, x; icellstart = lastnonzerocell, eps = eps, trybrute = !only_localsearch)
            if cell == 0
                return fill!(result, not_in_domain_value)
            else
                evaluate_bary!(result, PE, xref, cell)
                lastnonzerocell = cell
                return cell
            end
        end
        fe_function = point_evaluation_arbitrarygrids!
    end

    return interpolate!(target, ON_CELLS, fe_function; resultdim = resultdim, items = items, kwargs...)
end
