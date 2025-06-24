#########################
# PUBLIC INTERPOLATIONS #
#########################

# remap boundary face interpolation to faces by using BFaceFaces (if there is no special function by the finite element defined)
function ExtendableGrids.interpolate!(target::FEVectorBlock, FES::FESpace, ::Type{ON_BFACES}, source; items = items, kwargs...)
    if length(items) == 0
        items = FES.dofgrid[BFaceFaces]
    else
        items = FES.dofgrid[BFaceFaces][items]
    end
    return interpolate!(target, FES, ON_FACES, source; items = items, kwargs...)
end

"""
````
function ExtendableGrids.interpolate!(
    target::FEVectorBlock{T, Tv, Ti},
    AT::Type{<:AssemblyType},
    source;
    items = [],
    kwargs...
) where {T, Tv, Ti}
````

Interpolate a function or data into the finite element space (i.e. computes the coefficients) associated with `target`, using the specified assembly type.

# Arguments
- `target::FEVectorBlock`: The block of the FE vector to store the interpolated coefficients.
- `AT::Type{<:AssemblyType}`: The assembly type specifying where interpolation is performed (e.g., `ON_CELLS`).
- `source`: The function or callable to interpolate. Should have the signature `source!(result, qpinfo)`.

# Keyword Arguments
- `items`: List of mesh entities (cells, faces, etc.) to interpolate on. If empty, all entities of the specified type are used.
- `kwargs...`: Additional keyword arguments passed to lower-level routines (e.g., `bonus_quadorder`, `time`).

# Notes
- For "broken" FE spaces, interpolation is performed in a continuous auxiliary space and then mapped to the broken space.
- The `source!` function is called at each quadrature point and should fill `result` with the function values at that point.
- The `qpinfo` argument provides information about the current quadrature point, including coordinates, weights, and possibly time.

"""
function ExtendableGrids.interpolate!(
        target::FEVectorBlock{T, Tv, Ti},
        AT::Type{<:AssemblyType},
        source;
        items = [],
        kwargs...
    ) where {T, Tv, Ti}

    FEType = eltype(target.FES)
    return if target.FES.broken == true
        ## interpolate continuously
        FESc = FESpace{FEType}(target.FES.dofgrid)
        Targetc = FEVector{T}(FESc)
        interpolate!(Targetc[1], FESc, AT, source; items = items, kwargs...)
        celldofs = target.FES[CellDofs]
        if items == []
            items = 1:num_sources(celldofs)
        end

        ## copy continuous dofs to broken dofs
        function barrier(xCellDofs, xCellDofsc, items)
            for cell in items
                for k in 1:num_targets(xCellDofs, cell)
                    dof = xCellDofs[k, cell]
                    dofc = xCellDofsc[k, cell]
                    target[dof] = Targetc.entries[dofc]
                end
            end
            return
        end

        barrier(celldofs, FESc[CellDofs], items)

    else
        interpolate!(target, target.FES, AT, source; items = items, kwargs...)
    end
end

"""
````
function ExtendableGrids.interpolate!(target::FEVectorBlock,
	 source::Function;
	 kwargs...)
````

see interpolate!(target, ON_CELLS, source; kwargs...)
"""
function ExtendableGrids.interpolate!(target::FEVectorBlock, source; kwargs...)
    return interpolate!(target, ON_CELLS, source; kwargs...)
end


"""
````
function nodevalues_subset!(
	target::AbstractArray{<:Real,2},
	source::AbstractArray{T,1},
	FE::FESpace{Tv,Ti,FEType,AT},
	operator::Type{<:AbstractFunctionOperator} = Identity;
	regions::Array{Int,1} = [0],
	abs::Bool = false,
	factor = 1,
	nodes = [],				  
	target_offset::Int = 0,   # start to write into target after offset
	zero_target::Bool = true, # target vector is zeroed
	continuous::Bool = false)
````

Evaluate (an operator of) a finite element function (given by the coefficient vector `source` and FE space `FE`) at a specified subset of nodes, and write the results into `target`.

For each node in `nodes`, the function is evaluated (optionally with `operator`) and the result is written to the corresponding column of `target`. If `continuous` is `false`, values are averaged over all neighboring cells; if `true`, only one cell is used per node. If `abs` is `true`, the Euclidean norm is computed instead of the raw values.

# Arguments
- `target`: Output array to store the evaluated values (size: result dimension × number of nodes).
- `source`: Coefficient vector for the FE function.
- `FE`: The finite element space.
- `operator`: Function operator to apply (default: `Identity`).
- `abs`: If `true`, store the Euclidean norm of the result at each node (default: `false`).
- `factor`: Scaling factor applied to the result (default: `1`).
- `nodes`: List of node indices to evaluate (default: all nodes).
- `regions`: List of region indices to restrict evaluation (default: all regions).
- `target_offset`: Offset for writing into `target` (default: `0`).
- `source_offset`: Offset for reading from `source` (default: `0`).
- `zero_target`: If `true`, zero out `target` before writing (default: `true`).
- `continuous`: If `true`, evaluate only once per node; otherwise, average over all neighboring cells (default: `false`).

# Notes
- The result dimension is determined by the FE space and operator and the `abs` argument.
"""
function nodevalues_subset!(
        target::AbstractArray{T, 2},
        source::AbstractArray{T, 1},
        FE::FESpace{Tv, Ti, FEType, AT},
        operator::Type{<:AbstractFunctionOperator} = Identity;
        abs::Bool = false,
        factor = 1,
        nodes = [],
        regions::Array{Int, 1} = [0],
        target_offset::Int = 0,
        source_offset::Int = 0,
        zero_target::Bool = true,
        continuous::Bool = false
    ) where {T, Tv, Ti, FEType, AT}

    xgrid = FE.dofgrid
    xItemGeometries = xgrid[CellGeometries]
    xItemRegions::GridRegionTypes{Ti} = xgrid[CellRegions]
    xItemDofs::DofMapTypes{Ti} = FE[CellDofs]
    xItemNodes::Adjacency{Ti} = xgrid[CellNodes]

    EG = xgrid[UniqueCellGeometries]
    ndofs4EG::Array{Int, 1} = Array{Int, 1}(undef, length(EG))
    qf = Array{QuadratureRule, 1}(undef, length(EG))
    basisevaler::Array{FEEvaluator{T, Tv, Ti}, 1} = Array{FEEvaluator{T, Tv, Ti}, 1}(undef, length(EG))
    for j in 1:length(EG)
        qf[j] = VertexRule(EG[j])
        basisevaler[j] = FEEvaluator(FE, operator, qf[j]; T = T)
        ndofs4EG[j] = size(basisevaler[j].cvals, 2)
    end
    cvals_resultdim::Int = size(basisevaler[1].cvals, 1)
    target_resultdim::Int = abs ? 1 : cvals_resultdim
    @assert size(target, 1) >= target_resultdim "too small target dimension"

    # setup basisevaler for each unique cell geometries
    EG = xgrid[UniqueCellGeometries]

    if zero_target
        fill!(target, 0)
    end

    if regions == [0]
        try
            regions = Array{Int, 1}(Base.unique(xItemRegions[:]))
        catch
            regions = [xItemRegions[1]]
        end
    end
    nregions = length(regions)

    nnodes::Int = num_sources(xgrid[Coordinates])
    nneighbours::Int = 0

    nnodes = length(nodes)
    xNodeCells = xgrid[NodeCells]
    countedneighbours::Int = 0
    i::Int = 0
    iEG::Int = 1
    temp::Array{T, 1} = zeros(T, cvals_resultdim)
    localT::Array{T, 1} = zeros(T, cvals_resultdim)
    for n in 1:nnodes
        node = nodes[n]

        # get number of neighbours
        nneighbours = continuous ? 1 : num_targets(xNodeCells, node)
        countedneighbours = 0

        for b in 1:nneighbours
            item = xNodeCells[b, node]
            if xItemRegions[item] in regions
                countedneighbours += 1

                ## find local node index
                i = 1
                while xItemNodes[i, item] != node
                    i += 1
                end

                # find index for CellType
                if length(EG) > 1
                    itemET = xItemGeometries[item]
                    for j in 1:length(EG)
                        if itemET == EG[j]
                            iEG = j
                            break
                        end
                    end
                end

                # update FEbasisevaler
                update_basis!(basisevaler[iEG], item)

                fill!(localT, 0)
                for dof_i in 1:ndofs4EG[iEG]
                    dof = xItemDofs[dof_i, item]
                    eval_febe!(temp, basisevaler[iEG], dof_i, i)
                    for k in 1:cvals_resultdim
                        localT[k] += source[source_offset + dof] * temp[k]
                    end
                end
                localT .*= factor

                if abs
                    for k in 1:cvals_resultdim
                        target[1 + target_offset, n] += localT[k]^2
                    end
                else
                    for k in 1:cvals_resultdim
                        target[k + target_offset, n] += localT[k]
                    end
                end
            end
        end
        if countedneighbours > 0
            for k in 1:cvals_resultdim
                target[k + target_offset, n] /= countedneighbours
            end
        end

        if abs
            target[1 + target_offset, n] = sqrt(target[1 + target_offset, n])
        end

    end
    return
end


"""
````
function nodevalues!(
	target::AbstractArray{<:Real,2},
	source::AbstractArray{T,1},
	FE::FESpace{Tv,Ti,FEType,AT},
	operator::Type{<:AbstractFunctionOperator} = Identity;
	kwargs...)
````

calls the nodevalues_subset! function with all nodes of the dofgrid, see nodevalues_subset!(target, source, FE, operator; nodes = 1:num_nodes(FE.dofgrid), kwargs...)
"""
function nodevalues!(
        target::AbstractArray{T, 2},
        source::AbstractArray{T, 1},
        FE::FESpace{Tv, Ti, FEType, AT},
        operator::Type{<:AbstractFunctionOperator} = Identity;
        kwargs...
    ) where {T, Tv, Ti, FEType, AT}

    nodevalues_subset!(target, source, FE, operator; nodes = 1:num_nodes(FE.dofgrid), kwargs...)

    return nothing
end

"""
````
piecewise_nodevalues!(
        target::AbstractArray{T, 2},
        source::AbstractArray{T, 1},
        FE::FESpace{Tv, Ti, FEType, AT},
        operator::Type{<:AbstractFunctionOperator} = Identity;
        abs::Bool = false,
        factor = 1,
        regions::Array{Int, 1} = [0],
        target_offset::Int = 0,
        source_offset::Int = 0,
        zero_target::Bool = true,
        continuous::Bool = false
    )
````

Evaluate a finite element function (given by the coefficient vector `source` and FE space `FE`) at all nodes, but store the results in a piecewise (cellwise) fashion, i.e., for each cell, the values at its local nodes are written into `target`. The result is organized so that each column of `target` corresponds to a cell, and each row corresponds to the values at the cell's nodes.

# Arguments
- `target`: Output array to store the evaluated values (size: result dimension × number of nodes per cell, number of cells).
- `source`: Coefficient vector for the FE function.
- `FE`: The finite element space.
- `operator`: Function operator to apply (default: `Identity`).
- `abs`: If `true`, store the Euclidean norm of the result at each node (default: `false`).
- `factor`: Scaling factor applied to the result (default: `1`).
- `regions`: List of region indices to restrict evaluation (default: all regions).
- `target_offset`: Offset for writing into `target` (default: `0`).
- `source_offset`: Offset for reading from `source` (default: `0`).
- `zero_target`: If `true`, zero out `target` before writing (default: `true`).
- `continuous`: If `true`, evaluate only once per node; otherwise, average over all neighboring cells (default: `false`).

# Notes
- The result dimension is determined by the FE space, the operator, and the `abs` argument.
"""
function piecewise_nodevalues!(
        target::AbstractArray{T, 2},
        source::AbstractArray{T, 1},
        FE::FESpace{Tv, Ti, FEType, AT},
        operator::Type{<:AbstractFunctionOperator} = Identity;
        abs::Bool = false,
        factor = 1,
        regions::Array{Int, 1} = [0],
        target_offset::Int = 0,
        source_offset::Int = 0,
        zero_target::Bool = true
    ) where {T, Tv, Ti, FEType, AT}

    xgrid = FE.dofgrid
    xItemGeometries = xgrid[CellGeometries]
    xItemRegions::GridRegionTypes{Ti} = xgrid[CellRegions]
    xItemDofs::DofMapTypes{Ti} = FE[CellDofs]
    xItemNodes::Adjacency{Ti} = xgrid[CellNodes]
    nitems = num_sources(xItemNodes)
    target_resultdim::Int = 0

    if regions == [0]
        try
            regions = Array{Int, 1}(Base.unique(xItemRegions[:]))
        catch
            regions = [xItemRegions[1]]
        end
    end
    nregions = length(regions)

    # setup basisevaler for each unique cell geometries
    EG = xgrid[UniqueCellGeometries]

    if zero_target
        fill!(target, 0)
    end

    function barrier(EG, qf, BE)
        node::Int = 0
        dof::Ti = 0
        weights::Array{T, 1} = qf.w
        nweights = length(weights)
        basisvals = BE.cvals
        ndofs = size(basisvals, 2)
        cvals_resultdim::Int = size(basisvals, 1)
        temp::Array{T, 1} = zeros(T, cvals_resultdim)
        localT::Array{T, 1} = zeros(T, cvals_resultdim)
        target_resultdim = abs ? 1 : cvals_resultdim
        @assert size(target, 1) >= target_resultdim "too small target dimension"

        for item in 1:nitems
            for r in 1:nregions
                # check if item region is in regions
                if xItemRegions[item] == regions[r] && xItemGeometries[item] == EG

                    # update FEbasisevaler
                    update_basis!(BE, item)

                    for i in eachindex(weights) # vertices
                        node = xItemNodes[i, item]
                        fill!(localT, 0)

                        for dof_i in 1:ndofs
                            dof = xItemDofs[dof_i, item]
                            eval_febe!(temp, BE, dof_i, i)
                            for k in 1:cvals_resultdim
                                localT[k] += source[source_offset + dof] * temp[k]
                                #target[k+target_offset,node] += temp[k] * source[source_offset + dof]
                            end
                        end
                        localT .*= factor
                        if abs
                            for k in 1:cvals_resultdim
                                target[target_offset + i, item] += localT[k]^2
                            end
                        else
                            for k in 1:cvals_resultdim
                                target[target_offset + i + (k - 1) * nweights, item] += localT[k]
                            end
                        end
                    end

                    if abs
                        for i in 1:nweights
                            target[i + target_offset, item] = sqrt(target[i + target_offset, item])
                        end
                    end
                    break # region for loop
                end # if in region
            end # region for loop
        end # item for loop
        return
    end # barrier

    for j in 1:length(EG)
        qf = VertexRule(EG[j])
        BE = FEEvaluator(FE, operator, qf; T = T)
        barrier(EG[j], qf, BE)
    end

    return nothing
end


"""
````
function nodevalues!(
        target::AbstractArray{<:Real,2},
        source::AbstractArray{T,1},
        FE::FESpace{Tv,Ti,FEType,AT},
        operator::Type{<:AbstractFunctionOperator} = Identity;
        abs::Bool = false,
        factor = 1,
        regions::Array{Int,1} = [0],
        target_offset::Int = 0,
        source_offset::Int = 0,
        zero_target::Bool = true,
        continuous::Bool = false
    )
````

Evaluate a finite element function (given by the coefficient vector `source` and FE space `FE`) at all nodes of the grid, applying an optional function operator, and write the results into `target`.

By default, the function is evaluated at every node in the mesh. If `continuous` is `false`, the value at each node is averaged over all neighboring cells (suitable for discontinuous quantities). If `continuous` is `true`, the value is taken from a single cell per node (suitable for continuous quantities). The result can optionally be scaled, offset, or restricted to specific regions.

# Arguments
- `target`: Output array to store the evaluated values (size: result dimension × number of nodes).
- `source`: Coefficient vector for the FE function.
- `FE`: The finite element space.
- `operator`: Function operator to apply at each node (default: `Identity`).

# Keyword Arguments
- `abs`: If `true`, store the Euclidean norm of the result at each node (default: `false`).
- `factor`: Scaling factor applied to the result (default: `1`).
- `regions`: List of region indices to restrict evaluation (default: all regions).
- `target_offset`: Offset for writing into `target` (default: `0`).
- `source_offset`: Offset for reading from `source` (default: `0`).
- `zero_target`: If `true`, zero out `target` before writing (default: `true`).
- `continuous`: If `true`, evaluate only once per node; otherwise, average over all neighboring cells (default: `false`).

# Notes
- The result dimension is determined by the FE space, the operator, and the `abs` argument.
- The function modifies `target` in-place.
- For vector-valued or higher-dimensional results, the first dimension of `target` corresponds to the result dimension.

"""
function nodevalues!(
        target::AbstractArray{T, 2},
        source::AbstractArray{T, 1},
        FE::FESpace{Tv, Ti, FEType, AT},
        operator::Type{<:AbstractFunctionOperator} = Identity;
        abs::Bool = false,
        factor = 1,
        regions::Array{Int, 1} = [0],
        target_offset::Int = 0,
        source_offset::Int = 0,
        zero_target::Bool = true,
        continuous::Bool = false
    ) where {T, Tv, Ti, FEType, AT}

    xgrid = FE.dofgrid
    xItemGeometries = xgrid[CellGeometries]
    xItemRegions::GridRegionTypes{Ti} = xgrid[CellRegions]
    xItemDofs::DofMapTypes{Ti} = FE[CellDofs]
    xItemNodes::Adjacency{Ti} = xgrid[CellNodes]
    nitems = num_sources(xItemNodes)
    target_resultdim::Int = 0

    if regions == [0]
        try
            regions = Array{Int, 1}(Base.unique(xItemRegions[:]))
        catch
            regions = [xItemRegions[1]]
        end
    end
    nregions = length(regions)

    # setup basisevaler for each unique cell geometries
    EG = xgrid[UniqueCellGeometries]

    if zero_target
        fill!(target, 0)
    end

    nnodes::Int = num_sources(xgrid[Coordinates])
    nneighbours::Array{Int, 1} = zeros(Int, nnodes)
    flag4node::Array{Bool, 1} = zeros(Bool, nnodes)

    function barrier(EG, qf, BE)
        node::Int = 0
        dof::Ti = 0
        weights::Array{T, 1} = qf.w
        basisvals = BE.cvals
        ndofs = size(basisvals, 2)
        cvals_resultdim::Int = size(basisvals, 1)
        temp::Array{T, 1} = zeros(T, cvals_resultdim)
        localT::Array{T, 1} = zeros(T, cvals_resultdim)
        target_resultdim = abs ? 1 : cvals_resultdim
        @assert size(target, 1) >= target_resultdim "too small target dimension"

        for item in 1:nitems
            for r in 1:nregions
                # check if item region is in regions
                if xItemRegions[item] == regions[r] && xItemGeometries[item] == EG

                    # update FEbasisevaler
                    update_basis!(BE, item)

                    for i in eachindex(weights) # vertices
                        node = xItemNodes[i, item]
                        fill!(localT, 0)
                        if continuous == false || flag4node[node] == false
                            nneighbours[node] += 1
                            flag4node[node] = true
                            for dof_i in 1:ndofs
                                dof = xItemDofs[dof_i, item]
                                eval_febe!(temp, BE, dof_i, i)
                                for k in 1:cvals_resultdim
                                    localT[k] += source[source_offset + dof] * temp[k]
                                    #target[k+target_offset,node] += temp[k] * source[source_offset + dof]
                                end
                            end
                            localT .*= factor
                            if abs
                                for k in 1:cvals_resultdim
                                    target[1 + target_offset, node] += localT[k]^2
                                end
                            else
                                for k in 1:cvals_resultdim
                                    target[k + target_offset, node] += localT[k]
                                end
                            end
                        end
                    end
                    break # region for loop
                end # if in region
            end # region for loop
        end # item for loop
        return
    end # barrier

    for j in 1:length(EG)
        qf = VertexRule(EG[j])
        BE = FEEvaluator(FE, operator, qf; T = T)
        barrier(EG[j], qf, BE)
    end


    if continuous == false
        for n in 1:nnodes, k in 1:target_resultdim
            if nneighbours[n] > 0
                target[k + target_offset, n] /= nneighbours[n]
            end
        end
    end

    if abs
        for n in 1:nnodes
            target[1 + target_offset, n] = sqrt(target[1 + target_offset, n])
        end
    end

    return nothing
end

"""
````
    nodevalues(
        source::FEVectorBlock,
        operator::Type{<:AbstractFunctionOperator} = Identity;
        continuous = "auto",
        nodes = [],
        cellwise = false,
        abs = false,
        kwargs...
    )
````

Evaluate a finite element function (given by the coefficient vector `source`) at nodes of the grid, applying an optional function operator, and return the result as a newly allocated array of the appropriate size.

This function provides a flexible interface for extracting nodal or cellwise values from a finite element solution. By default, it evaluates at all nodes, but a subset of nodes can be specified. The result can be returned in a cellwise (piecewise) layout if desired.

# Arguments
- `source`: The `FEVectorBlock` containing the coefficients of the FE function.
- `operator`: The function operator to apply at each node (default: `Identity`).

# Keyword Arguments
- `continuous`: If `"auto"`, automatically choose continuous/discontinuous evaluation based on FE type and operator. If `true`, evaluate only once per node; if `false`, average over all neighboring cells. Is ignored when `cellwise` is `true`.
- `nodes`: List of node indices to evaluate (default: all nodes).
- `cellwise`: If `true`, return values in a cellwise (piecewise) layout (default: `false`).
- `abs`: If `true`, return the Euclidean norm at each node (default: `false`).
- `kwargs...`: Additional keyword arguments passed to lower-level routines (e.g., `regions`, `factor`, `target_offset`, `zero_target`, etc.).

# Returns
- A newly allocated array containing the evaluated values, with shape depending on the options chosen.

# Notes
- If `nodes` is empty, all nodes are evaluated.
- If `cellwise` is `true`, the result is organized per cell (suitable for discontinuous or element-wise quantities).
- The result dimension is determined by the FE space, the operator, and the `abs` argument.
"""
function nodevalues(
        source::FEVectorBlock{T, Tv, Ti, FEType, APT},
        operator::Type{<:AbstractFunctionOperator} = Identity;
        continuous = "auto",
        nodes = [],
        cellwise = false,
        abs = false,
        kwargs...
    ) where {T, Tv, Ti, APT, FEType}
    if continuous == "auto"
        if FEType <: AbstractH1FiniteElement && operator == Identity && !source.FES.broken && !(FEType <: H1CR)
            continuous = true
        else
            continuous = false
        end
    end
    if abs
        nvals = 1
    else
        xdim = size(source.FES.dofgrid[Coordinates], 1)
        ncomponents = get_ncomponents(eltype(source.FES))
        nvals = Length4Operator(operator, xdim, ncomponents)
    end
    if nodes == []
        if cellwise
            target = zeros(T, nvals * max_num_targets_per_source(source.FES.dofgrid[CellNodes]), num_cells(source.FES.dofgrid))
            piecewise_nodevalues!(target, source.entries, source.FES, operator; abs = abs, kwargs...)
        else
            target = zeros(T, nvals, num_nodes(source.FES.dofgrid))
            nodevalues!(target, source.entries, source.FES, operator; continuous = continuous, abs = abs, kwargs...)
        end
    else
        @assert cellwise == false
        target = zeros(T, nvals, length(nodes))
        nodevalues_subset!(target, source.entries, source.FES, operator; continuous = continuous, abs = abs, nodes = nodes, kwargs...)
    end
    return target
end

"""
````
function nodevalues_view(
	source::FEVectorBlock,
	operator::Type{<:AbstractFunctionOperator} = Identity)
````

Return a vector of views into the nodal values of the given finite element function, allowing direct access to the underlying coefficient storage for each component.

This function provides efficient, zero-copy access to the nodal values of an `FEVectorBlock` for unbroken H1-conforming finite element spaces with the identity operator. Each entry in the returned vector is a view into the coefficients corresponding to one component of the FE function, optionally restricted to a subset of nodes.

# Arguments
- `source`: The `FEVectorBlock` containing the coefficients of the FE function.
- `operator`: The function operator to apply (must be `Identity` for direct views; default: `Identity`).

# Keyword Arguments
- `nodes`: List of node indices to view (default: all nodes).

# Returns
- A vector of `SubArray` views, one for each component, directly referencing the coefficients for the specified nodes.

# Notes
- Only available for unbroken H1-conforming elements and the Identity operator.

"""
function nodevalues_view(source::FEVectorBlock{T, Tv, Ti, FEType, APT}, operator::Type{<:AbstractFunctionOperator} = Identity; nodes = [0]) where {T, Tv, Ti, APT, FEType}

    if (FEType <: AbstractH1FiniteElement) && (operator == Identity) && (source.FES.broken == false)
        # give a direct view without computing anything
        ncomponents = get_ncomponents(FEType)
        array_of_views = []
        offset::Int = source.offset
        coffset::Int = source.FES.coffset
        if nodes == [0]
            nodes = 1:num_nodes(source.FES.dofgrid)
        end
        for k in 1:ncomponents
            push!(array_of_views, view(source.entries, offset .+ nodes))
            offset += coffset
        end
        return array_of_views
    else
        @error "nodevalues_view not available for FEType = $FEType and operator = $operator"
    end
end


function nodevalues(xgrid::ExtendableGrid{Tv, Ti}, f::Function; T = Float64, resultdim = 1, time = 0) where {Tv, Ti}
    xCoordinates::Array{Tv, 2} = xgrid[Coordinates]
    nnodes::Int = size(xCoordinates, 2)
    QP = QPInfos(xgrid; time = time)
    nodevals = zeros(T, resultdim, nnodes)
    for j in 1:nnodes
        QP.x = view(xCoordinates, :, j)
        f(view(nodevals, :, j), QP)
    end
    return nodevals
end

"""
````
function continuify(
	source::FEVectorBlock,
	operator = Identity;
	abs::Bool = false,
	broken = false,
	order = "auto",
	factor = 1,
	regions::Array{Int,1} = [0]) where {T,Tv,Ti,FEType,APT}
````

Interpolate the evaluation of an operator applied to a finite element function onto a continuous Lagrange finite element space (`H1Pk`), returning a new `FEVector` with the interpolated values.

This function performs nodal interpolation of the (possibly vector-valued) result of applying `operator` to the FE function represented by `source`. The result is a new FE function in a continuous Lagrange space of the specified order and dimension. If `broken = true`, the interpolation is performed in a piecewise (discontinuous) fashion.

# Arguments
- `source`: The `FEVectorBlock` containing the coefficients of the original FE function.
- `operator`: The function operator to apply before interpolation (default: `Identity`).

# Keyword Arguments
- `abs`: If `true`, interpolate the Euclidean norm of the result (default: `false`).
- `broken`: If `true`, generate a piecewise (discontinuous) interpolation (default: `false`).
- `order`: Polynomial order of the target Lagrange space (default: `"auto"`, which chooses an appropriate order).
- `factor`: Scaling factor applied to the result before interpolation (default: `1`).
- `regions`: List of region indices to restrict interpolation (default: all regions).

# Returns
- A new `FEVector` in a continuous (or broken) Lagrange space, containing the interpolated values.

"""
function continuify(
        source::FEVectorBlock{T, Tv, Ti, FEType, APT},
        operator::Type{<:AbstractFunctionOperator} = Identity;
        abs::Bool = false,
        broken = false,
        order = "auto",
        factor = 1,
        regions::Array{Int, 1} = [0]
    ) where {T, Tv, Ti, FEType, APT}

    FE = source.FES
    xgrid = FE.dofgrid
    xItemGeometries = xgrid[CellGeometries]
    xItemRegions::GridRegionTypes{Ti} = xgrid[CellRegions]
    xItemDofs::DofMapTypes{Ti} = FE[CellDofs]
    ncomponents = get_ncomponents(FEType)
    xdim = size(xgrid[Coordinates], 1)
    nitems::Int = num_sources(xItemDofs)

    if regions == [0]
        try
            regions = Array{Int, 1}(Base.unique(xItemRegions[:]))
        catch
            regions = [xItemRegions[1]]
        end
    end

    # setup basisevaler for each unique cell geometries
    EG = xgrid[UniqueCellGeometries]
    if order == "auto"
        order = max(get_polynomialorder(FEType, EG[1]) + QuadratureOrderShift4Operator(operator), 1)
    end
    cvals_resultdim::Int = Length4Operator(operator, xdim, ncomponents)
    target_resultdim::Int = abs ? 1 : cvals_resultdim


    name = "$operator(" * source.name * ")"
    edim = dim_element(EG[1])
    if edim in 1:2
        FETypeC = H1Pk{cvals_resultdim, edim, order}
    else
        if order == 1
            FETypeC = H1P1{cvals_resultdim}
        elseif order == 2
            FETypeC = H1P2{cvals_resultdim, 3}
        elseif order == 2
            FETypeC = H1P3{cvals_resultdim, 3}
        else
            @error "continuify target order > 3 currently not available in 3D"
        end
    end
    FEScont = FESpace{FETypeC}(xgrid; broken = broken)
    target = FEVector(FEScont)
    xItemDofsC::DofMapTypes{Ti} = target[1].FES[CellDofs]
    target_offset::Int = broken ? Int(get_ndofs(ON_CELLS, FETypeC, EG[1]) / cvals_resultdim) : target[1].FES.coffset

    if order > 2
        @warn "continuify may not work correctly if target order is larger than 2 currently"
    end
    @debug "Interpolating $(source.name) ($FEType) >> $(target[1].name) ($FETypeC)"

    #subset_handler! = get_basissubset(ON_CELLS, target[1].FES, EG[1])
    #subset_ids = Array{Int,1}(1 : get_ndofs(ON_CELLS, FETypeC, EG[1]))

    nneighbours::Array{Int, 1} = broken ? [] : zeros(Int, target_offset)

    function barrier(EG, qf, BE)

        basisvals = BE.cvals
        ndofs::Int = size(basisvals, 2)
        nregions::Int = length(regions)
        dof::Ti = 0
        dofc::Ti = 0
        temp::Array{T, 1} = zeros(T, cvals_resultdim)
        localT::Array{T, 1} = zeros(T, cvals_resultdim)
        weights::Array{T, 1} = qf.w

        for item in 1:nitems
            for r in 1:nregions
                # check if item region is in regions
                if xItemRegions[item] == regions[r] && xItemGeometries[item] == EG

                    # update FEbasisevaler
                    update_basis!(BE, item)
                    #subset_handler!(subset_ids, item)

                    for i in eachindex(weights) # dofs
                        fill!(localT, 0)
                        for dof_i in 1:ndofs
                            dof = xItemDofs[dof_i, item]
                            eval_febe!(temp, BE, dof_i, i)
                            for k in 1:cvals_resultdim
                                localT[k] += source[dof] * temp[k]
                            end
                        end
                        localT .*= factor
                        dofc = xItemDofsC[i, item]
                        if !broken
                            nneighbours[dofc] += 1
                        end
                        if abs
                            for k in 1:cvals_resultdim
                                target.entries[dofc + (k - 1) * target_offset] += localT[k]^2
                            end
                        else
                            for k in 1:cvals_resultdim
                                target.entries[dofc + (k - 1) * target_offset] += localT[k]
                            end
                        end
                    end
                    break # region for loop
                end # if in region
            end # region for loop
        end # item for loop
        return
    end

    for j in 1:length(EG)
        qf = VertexRule(EG[j], order)
        BE = FEEvaluator(FE, operator, qf; T = T)
        cvals_resultdim = size(BE.cvals, 1)
        barrier(EG[j], qf, BE)
    end

    if !broken
        for dofc in 1:target_offset, k in 1:target_resultdim
            target.entries[dofc + (k - 1) * target_offset] /= nneighbours[dofc]
        end
    end

    return target
end


"""
````
function displace_mesh!(xgrid::ExtendableGrid, source::FEVectorBlock; magnify = 1)
````
Displace all nodes of the given grid by adding a vector-valued finite element field as a displacement, scaled by an optional magnification factor.

This function modifies the coordinates of `xgrid` in-place by adding the nodal values of the FE function represented by `source` (typically a displacement field) to each node. The displacement can be scaled by the `magnify` parameter. After the update, all cached geometric quantities in the grid are invalidated and recomputed.

# Arguments
- `xgrid`: The `ExtendableGrid` whose node coordinates will be updated.
- `source`: An `FEVectorBlock` representing the displacement field (must be vector-valued and defined on the grid).

# Keyword Arguments
- `magnify`: Scaling factor for the displacement field (default: `1`).

"""
function displace_mesh!(xgrid::ExtendableGrid, source::FEVectorBlock; magnify = 1)
    nnodes = size(xgrid[Coordinates], 2)
    nodevals = zeros(eltype(xgrid[Coordinates]), get_ncomponents(Base.eltype(source.FES)), nnodes)
    nodevalues!(nodevals, source, Identity)
    xgrid[Coordinates] .+= magnify * nodevals

    # remove all keys from grid components that might have changed and need a reinstantiation
    ExtendableGrids.update!(xgrid, CellVolumes)
    ExtendableGrids.update!(xgrid, FaceVolumes)
    ExtendableGrids.update!(xgrid, EdgeVolumes)
    ExtendableGrids.update!(xgrid, FaceNormals)
    ExtendableGrids.update!(xgrid, EdgeTangents)
    ExtendableGrids.update!(xgrid, BFaceVolumes)
    return ExtendableGrids.update!(xgrid, BEdgeVolumes)
end


"""
````
function displace_mesh(xgrid::ExtendableGrid, source::FEVectorBlock; magnify = 1)
````
Return a new grid with node coordinates displaced by a vector-valued finite element field, optionally scaled by a magnification factor.

# Arguments
- `xgrid`: The `ExtendableGrid` to be copied and displaced.
- `source`: An `FEVectorBlock` representing the displacement field (must be vector-valued and defined on the grid).

# Keyword Arguments
- `magnify`: Scaling factor for the displacement field (default: `1`).

# Returns
- A new `ExtendableGrid` with displaced node coordinates.

"""
function displace_mesh(xgrid::ExtendableGrid, source::FEVectorBlock; kwargs...)
    xgrid_displaced = deepcopy(xgrid)
    displace_mesh!(xgrid_displaced, source; kwargs...)
    return xgrid_displaced
end
