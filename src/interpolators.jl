Base.get!(FES::FESpace, entity::Type{<:AssemblyType}) = get!(() -> init_interpolator!(FES, entity), FES.interpolators, entity)
Base.getindex(FES::FESpace, entity::Type{<:AssemblyType}) = get!(FES, entity)
Base.setindex!(FES::FESpace, v, E::Type{<:AssemblyType}) = FES.interpolators[E] = v

struct NodalInterpolator{EFT} <: AbstractInterpolationOperator
    evaluator::EFT
end

function NodalInterpolator(FES::FESpace{T}, grid = FES.dofgrid; component_offset = FES.coffset, time = 0, kwargs...) where {T}
    FEType = eltype(FES)
    ncomponents = get_ncomponents(FEType)
    nnodes = num_nodes(grid)
    offset4component = 0:component_offset:(ncomponents * component_offset)
    xNodeCells = atranspose(grid[CellNodes])
    xCoordinates = grid[Coordinates]
    xCellRegions = grid[CellRegions]
    result = zeros(T, ncomponents)
    QP = QPInfos(grid; time, kwargs...)

    function evaluate!(target, exact_function!, items)
        if isempty(items)
            items = 1 : nnodes
        end
        for j in items
            cell = xNodeCells[1, j]
            QP.item = cell
            QP.cell = cell
            QP.region = xCellRegions[cell]
            QP.x .= view(xCoordinates, :, j)
            exact_function!(result, QP)
            for k in 1:ncomponents
                target[j + offset4component[k]] = result[k]
            end
        end
    end

    return NodalInterpolator(evaluate!)
end

struct MomentInterpolator{EFT} <: AbstractInterpolationOperator
    evaluator::EFT
end

function MomentInterpolator(FE::FESpace{Tv, Ti, FEType, APT}, AT::Type{<:AssemblyType}, xgrid = FE.dofgrid; FEType_ref = "auto", bestapprox = false, order = 0, items = [], kwargs...) where {T, Tv, Ti, FEType <: AbstractH1FiniteElement, APT}

    xItemVolumes::Array{Tv, 1} = xgrid[GridComponentVolumes4AssemblyType(AT)]
    xItemNodes::Adjacency{Ti} = xgrid[GridComponentNodes4AssemblyType(AT)]
    xItemDofs::DofMapTypes{Ti} = Dofmap4AssemblyType(FE, AT)
    EGs = xgrid[GridComponentUniqueGeometries4AssemblyType(AT)]

    @assert length(EGs) == 1 "ensure_moments! currently only works with grids with a single element geometry"
    EG = EGs[1]

    nitems::Int = num_sources(xItemNodes)
    if items == []
        items = 1:nitems
    end
    ncomponents::Int = get_ncomponents(FEType)
    edim::Int = dim_element(EG)
    order_FE = get_polynomialorder(FEType, EG)
    coffset::Int = get_ndofs(AT, FEType, EG) / ncomponents
    interior_offset = interior_dofs_offset(AT, FEType, EG)

    ## get basis for moments
    ## here we assume that the FEType looks like a H1Pk element on EG (which is true for all H1Pk elements)
    if order == 0
        FEType_moments = L2P0{ncomponents}
    elseif order == 1
        FEType_moments = H1P1{ncomponents}
    elseif order == 2
        FEType_moments = H1P2{ncomponents, edim}
    else
        FEType_moments = H1Pk{ncomponents, edim, order}
    end

    if FEType_ref == "auto"
        if AT == ON_CELLS
            FEType_ref = FEType
        else
            if edim == 2 && order == 0
                FEType_ref = H1P3{ncomponents, edim} # order + 3
            elseif edim == 1
                FEType_ref = H1Pk{ncomponents, edim, order + 2}
            else
                @error "not yet supported"
            end
        end
    end

    moments_basis! = get_basis(ON_CELLS, FEType_moments, EG)
    nmoments::Int = get_ndofs_all(ON_CELLS, FEType_moments, EG)
    xgrid_ref = reference_domain(EG)
    nmoments4c::Int = nmoments / ncomponents
    idofs = zeros(Int, 0)
    for c in 1:ncomponents, m in 1:nmoments4c
        push!(idofs, (c - 1) * coffset + interior_offset + m)
    end

    MOMxBASIS::Array{Float64, 2} = zeros(Float64, 0, 0)
    FE_onref = FESpace{FEType_ref}(xgrid_ref)
    moments_eval::Matrix{Float64} = zeros(Float64, 0, 0)
    if (bestapprox) # interior dofs are set by best-approximation
        refbasis = get_basis(ON_CELLS, FEType_ref, EG)
        ndofs_ref = get_ndofs(ON_CELLS, FEType_ref, EG)
        refbasis_vals = zeros(Tv, ndofs_ref, ncomponents)

        function refbasis_times_refbasis(result, qpinfo)
            x = qpinfo.x
            refbasis(refbasis_vals, x)
            for j in 1:ndofs_ref, k in 1:ndofs_ref
                result[(k - 1) * ndofs_ref + j] = dot(view(refbasis_vals, j, :), view(refbasis_vals, k, :))
            end
            return nothing
        end

        MOMxBASIS = reshape(integrate(xgrid_ref, ON_CELLS, refbasis_times_refbasis, ndofs_ref^2; quadorder = 2 * order_FE, kwargs...), (ndofs_ref, ndofs_ref))
        MOMxBASIS ./= xgrid_ref[CellVolumes][1]

        ## extract quadratic matrix for interior dofs
        MOMxINTERIOR = zeros(length(idofs), length(idofs))
        for j in 1:length(idofs), k in 1:length(idofs)
            MOMxINTERIOR[j, k] = MOMxBASIS[idofs[j], idofs[k]]
        end
        moments_eval = zeros(Float64, size(MOMxBASIS, 1), ncomponents)
        moments_basis! = get_basis(ON_CELLS, FEType_ref, EG)
        MOMxBASIS = MOMxBASIS[:, idofs]
    else # interior dofs are set by moments
        ## calculate moments times basis functions

        refbasis = get_basis(ON_CELLS, FEType_ref, EG)
        momentbasis = get_basis(ON_CELLS, FEType_moments, EG)
        ndofs_ref = get_ndofs(ON_CELLS, FEType_ref, EG)
        ndofs_moment = get_ndofs(ON_CELLS, FEType_moments, EG)
        refbasis_vals = zeros(Tv, ndofs_ref, ncomponents)
        momentbasis_vals = zeros(Tv, ndofs_moment, ncomponents)

        function momentbasis_times_refbasis(result, qpinfo)
            x = qpinfo.x
            refbasis(refbasis_vals, x)
            momentbasis(momentbasis_vals, x)
            for j in 1:ndofs_ref, k in 1:ndofs_moment
                result[(k - 1) * ndofs_ref + j] = dot(view(refbasis_vals, j, :), view(momentbasis_vals, k, :))
            end
            return nothing
        end

        MOMxBASIS_temp = integrate(xgrid_ref, ON_CELLS, momentbasis_times_refbasis, ndofs_ref * ndofs_moment; quadorder = 2 * order_FE)
        MOMxBASIS = reshape(MOMxBASIS_temp, (ndofs_ref, ndofs_moment))
        MOMxBASIS ./= xgrid_ref[CellVolumes][1]

        ## extract quadratic matrix for interior dofs
        #MOMxINTERIOR = view(MOMxBASIS,idofs,1:nmoments)  # inverting this in line 228 does not work in julia 1.6 !!!
        MOMxINTERIOR = zeros(length(idofs), size(MOMxBASIS, 2))
        for j in 1:length(idofs), k in 1:size(MOMxBASIS, 2)
            MOMxINTERIOR[j, k] = MOMxBASIS[idofs[j], k]
        end

        moments_eval = zeros(Float64, nmoments, ncomponents)
    end

    ### get permutation of dofs on reference EG and real cells
    subset_handler = get_basissubset(AT, FE, EG)
    current_subset = Array{Int, 1}(1:size(MOMxBASIS, 1))
    #doforder_ref::Array{Int,1} = FE_onref[CellDofs][:,1]
    invA::Array{Float64, 2} = inv(MOMxINTERIOR)

    ## evaluator for moments of exact_function
    result_f::Vector{Tv} = zeros(Tv, ncomponents)
    function f_times_moments(exact_function!)
        function closure(result, qpinfo)
            fill!(moments_eval, 0)
            moments_basis!(moments_eval, qpinfo.xref)
            exact_function!(result_f, qpinfo)
            fill!(result, 0)
            if (bestapprox)
                for m in 1:nmoments, k in 1:ncomponents
                    result[m] += result_f[k] * moments_eval[idofs[m], k]
                end
            else
                for m in 1:nmoments, k in 1:ncomponents
                    result[m] += result_f[k] * moments_eval[m, k]
                end
            end
        end
    end

    # integrate moments of exact_function over edges
    edgemoments::Array{Tv, 2} = zeros(Tv, nmoments, nitems)

    function evaluate!(target, exact_function!, items)
        if isempty(items)
            items = 1 : nitems
        end
        integrate!(edgemoments, xgrid, AT, f_times_moments(exact_function!); quadorder = 2 * order_FE, items = items, kwargs...)
        for item::Int in items
            if subset_handler != NothingFunction
                subset_handler(current_subset, item)
            end
            for m::Int in 1:nmoments, exdof in 1:interior_offset, c in 1:ncomponents
                localdof = coffset * (c - 1) + exdof
                edgemoments[m, item] -= target[xItemDofs[localdof, item]] * MOMxBASIS[current_subset[localdof], m] * xItemVolumes[item]
            end
            for m::Int in 1:nmoments
                localdof = idofs[m]
                target[xItemDofs[localdof, item]] = 0
                for n::Int in 1:nmoments
                    target[xItemDofs[localdof, item]] += invA[n, m] * edgemoments[n, item]
                end
                target[xItemDofs[localdof, item]] /= xItemVolumes[item]
            end
        end
    end

    return MomentInterpolator(evaluate!)
end
