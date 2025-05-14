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
function ExtendableGrids.interpolate!(target::FEVectorBlock,
	 AT::Type{<:AssemblyType},
	 source!::Function;
	 items = [],
	 bonus_quadorder = 0,
	 time = 0,
	 kwargs...)
````

Interpolates the given source into the finite elements space assigned to the target FEVectorBlock with the specified AssemblyType
(usually ON_CELLS). 

The source functions should adhere to the interface
```julia
	source!(result, qpinfo)
```
The qpinfo argument communicates vast information of the current quadrature/evaluation point.

The bonus_quadorder argument can be used to steer the quadrature order of integrals that needs to be computed
for the interpolation (the default quadrature order corresponds to the polynomial order of the finite element).
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
	 items = [],
	 bonus_quadorder = 0,
	 time = 0,
	 kwargs...)
````

Interpolates the given source function into the finite element space assigned to the target FEVectorBlock. 
	
The source functions should adhere to the interface
```julia
	source!(result, qpinfo)
```
The qpinfo argument communicates vast information of the current quadrature/evaluation point.

The bonus_quadorder argument can be used to steer the quadrature order of integrals that needs to be computed
for the interpolation (the default quadrature order corresponds to the polynomial order of the finite element).
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

Evaluates the finite element function with the coefficient vector source (interpreted as a coefficient vector for the FESpace FE)
and the specified FunctionOperator at the specified list of nodes of the grid and writes the values in that order into target.
Node values for nodes that are not part of the specified regions (default = all regions) are set to zero.
Discontinuous (continuous = false) quantities are evaluated in all neighbouring cells (in the specified regions)
of each node and then averaged. Continuous (continuous = true) quantities are only evaluated once at each node.
"""
function nodevalues_subset!(
        target::AbstractArray{T, 2},
        source::AbstractArray{T, 1},
        FE::FESpace{Tv, Ti, FEType, AT},
        operator::Type{<:AbstractFunctionOperator} = Identity;
        abs::Bool = false,
        factor = 1,
        nodes = [],
        regions::Array{Int, 1} = [0], # are ignored at the moment
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
	regions::Array{Int,1} = [0],
	abs::Bool = false,
	factor = 1,
	target_offset::Int = 0,   # start to write into target after offset
	zero_target::Bool = true, # target vector is zeroed
	continuous::Bool = false)
````

Evaluates the finite element function with the coefficient vector source (interpreted as a coefficient vector for the FESpace FE)
and the specified FunctionOperator at all the nodes of the (specified regions of the) grid and writes the values into target.
Discontinuous (continuous = false) quantities are evaluated in all neighbouring cells of each node and then averaged. Continuous
(continuous = true) quantities are only evaluated once at each node.
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
	source::FEVectorBlock,
	operator::Type{<:AbstractFunctionOperator} = Identity;
	regions::Array{Int,1} = [0],
	abs::Bool = false,
	factor = 1,
	cellwise = false,		  # return cellwise nodevalues ncells x nnodes_on_cell
	target_offset::Int = 0,   # start to write into target after offset
	zero_target::Bool = true, # target vector is zeroed
	continuous::Bool = false)
````

Evaluates the finite element function with the coefficient vector source
and the specified FunctionOperator at all the nodes of the (specified regions of the) grid and writes the values into target.
Discontinuous (continuous = false) quantities are evaluated in all neighbouring cells of each node and then averaged. Continuous
(continuous = true) quantities are only evaluated once at each node.
"""
function nodevalues!(target, source::FEVectorBlock, operator::Type{<:AbstractFunctionOperator} = Identity; cellwise = false, kwargs...)
    return if cellwise
        piecewise_nodevalues!(target, source.entries, source.FES, operator; source_offset = source.offset, kwargs...)
    else
        nodevalues!(target, source.entries, source.FES, operator; source_offset = source.offset, kwargs...)
    end
end


"""
````
function nodevalues(
	source::FEVectorBlock,
	operator::Type{<:AbstractFunctionOperator} = Identity;
	regions::Array{Int,1} = [0],
	abs::Bool = false,
	factor = 1,
	nodes = [],				  
	cellwise = false,		  # return cellwise nodevalues ncells x nnodes_on_cell (only if nodes == [])
	target_offset::Int = 0,   # start to write into target after offset
	zero_target::Bool = true, # target vector is zeroed
	continuous::Bool = false)
````

Evaluates the finite element function with the coefficient vector source
and the specified FunctionOperator at the specified list of nodes of the grid (default = all nodes)
and writes the values in that order into target. Nodes that are not part of the specified regions (default = all regions)
are set to zero.
Discontinuous (continuous = false) quantities are evaluated in all neighbouring cells of each node and then averaged. Continuous
(continuous = true) quantities are only evaluated once at each node.
"""
function nodevalues(source::FEVectorBlock{T, Tv, Ti, FEType, APT}, operator::Type{<:AbstractFunctionOperator} = Identity; continuous = "auto", nodes = [], cellwise = false, abs = false, kwargs...) where {T, Tv, Ti, APT, FEType}
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
            piecewise_nodevalues!(target, source.entries, source.FES, operator; continuous = continuous, source_offset = source.offset, abs = abs, kwargs...)
        else
            target = zeros(T, nvals, num_nodes(source.FES.dofgrid))
            nodevalues!(target, source.entries, source.FES, operator; continuous = continuous, source_offset = source.offset, abs = abs, kwargs...)
        end
    else
        @assert cellwise == false
        target = zeros(T, nvals, length(nodes))
        nodevalues_subset!(target, source.entries, source.FES, operator; continuous = continuous, source_offset = source.offset, abs = abs, nodes = nodes, kwargs...)
    end
    return target
end

"""
````
function nodevalues_view(
	source::FEVectorBlock,
	operator::Type{<:AbstractFunctionOperator} = Identity)
````

Returns a vector of views of the nodal values of the source block (currently works for unbroken H1-conforming elements) that directly accesses the coefficients.
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

interpolates operator evaluation of source into a FE function of FEType H1Pk, i.e., Lagrange interpolation of arbitrary 
operator evaluations of the source finite element type, broken = true generates a piecewise interpolation
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
Moves all nodes of the grid by adding the displacement field in source (expects a vector-valued finite element)
times a magnify value.
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
Returns a new grid by adding the displacement field in source (expects a vector-valued finite element)
to the coordinates of the provided xgrid times a magnify value.
"""
function displace_mesh(xgrid::ExtendableGrid, source::FEVectorBlock; kwargs...)
    xgrid_displaced = deepcopy(xgrid)
    displace_mesh!(xgrid_displaced, source; kwargs...)
    return xgrid_displaced
end
