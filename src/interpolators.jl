get_interpolator(FES::FESpace, entity::Type{<:AssemblyType}) = get!(() -> init_interpolator!(FES, entity), FES.interpolators, entity)

struct NodalInterpolator{EFT} <: AbstractInterpolationOperator
    evaluate!::EFT
end

function NodalInterpolator(FES::FESpace{T}, grid = FES.dofgrid; broken = FES.broken, component_offset = FES.coffset, time = 0, kwargs...) where {T}
    FEType = eltype(FES)
    ncomponents = get_ncomponents(FEType)
    offset4component = 0:component_offset:(ncomponents * component_offset)
    xCoordinates = grid[Coordinates]
    xCellRegions = grid[CellRegions]
    result = zeros(T, ncomponents)
    QP = QPInfos(grid; time, kwargs...)

    if broken
        ## FE space is broken
        xCellNodes = grid[CellNodes]
        xCellDofs = FES[CellDofs]
        ncells = num_cells(grid)
        function evaluate_broken!(target, exact_function!, items; kwargs...)
            if isempty(items)
                items = 1 : ncells
            end
            for cell in items
                nnodes_on_cell = num_targets(xCellNodes, cell)
                QP.item = cell
                QP.cell = cell
                QP.region = xCellRegions[cell]
                for n in 1:nnodes_on_cell
                    j = xCellNodes[n, cell]
                    QP.x .= view(xCoordinates, :, j)
                    exact_function(result, QP)
                    for k in 1:ncomponents
                        target[xCellDofs[1, cell] + n - 1 + (k - 1) * nnodes_on_cell] = result[k]
                    end
                end
            end
        end
        return NodalInterpolator(evaluate_broken!)
    else
        ## FE space is continuous in node, so only one evaluation is required
        nnodes = num_nodes(grid)
        xNodeCells = atranspose(grid[CellNodes])
        function evaluate!(target, exact_function!, items; kwargs...)
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
end

struct MomentInterpolator{EFT} <: AbstractInterpolationOperator
    evaluate!::EFT
end

function MomentInterpolator(FE::FESpace{Tv, Ti, FEType, APT}, AT::Type{<:AssemblyType}, xgrid = FE.dofgrid; FEType_ref = "auto", bestapprox = false, order = 0, items = [], kwargs...) where {Tv, Ti, FEType <: AbstractH1FiniteElement, APT}

    itemvolumes = xgrid[GridComponentVolumes4AssemblyType(AT)]
    itemnodes = xgrid[GridComponentNodes4AssemblyType(AT)]
    itemregions = xgrid[GridComponentRegions4AssemblyType(AT)]
    itemdofs = Dofmap4AssemblyType(FE, AT)
    has_normals = true
    if AT <: ON_FACES
        itemnormals = xgrid[FaceNormals]
    elseif AT <: ON_BFACES
        itemnormals = xgrid[FaceNormals][:, xgrid[BFaceFaces]]
    else
        has_normals = false
    end
    EGs = xgrid[GridComponentUniqueGeometries4AssemblyType(AT)]

    @assert length(EGs) == 1 "MomentInterpolator currently only works with grids with a single element geometry"
    EG = EGs[1]

    nitems::Int = num_sources(itemnodes)
    if items == []
        items = 1:nitems
    end
    ncomponents::Int = get_ncomponents(FEType)
    edim::Int = dim_element(EG)
    order_FE = get_polynomialorder(FEType, EG)
    coffset::Int = get_ndofs(AT, FEType, EG) / ncomponents
    interior_offset = interior_dofs_offset(AT, FEType, EG)
    @assert interior_offset >= 0 "This FEType seems to missing a definition for interior_dofs_offset!"

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
                @error "MomentInterpolator with order $order for edim $edim not yet supported"
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
    moments_eval::Matrix{Float64} = zeros(Float64, 0, 0)
    refbasis! = get_basis(ON_CELLS, FEType_ref, EG)
    ndofs_ref = get_ndofs(ON_CELLS, FEType_ref, EG)
    refbasis_vals = zeros(Tv, ndofs_ref, ncomponents)

    if (bestapprox) # interior dofs are set by best-approximation

        ## get mass matrix of reference basis
        function refbasis_times_refbasis(result, qpinfo)
            x = qpinfo.x
            fill!(refbasis_vals, 0)
            refbasis!(refbasis_vals, x)
            for j in 1:ndofs_ref, k in 1:ndofs_ref
                result[(k - 1) * ndofs_ref + j] = dot(view(refbasis_vals, j, :), view(refbasis_vals, k, :))
            end
            return nothing
        end

        MOMxBASIS_temp = integrate(xgrid_ref, ON_CELLS, refbasis_times_refbasis, ndofs_ref * ndofs_ref; quadorder = 2 * order_FE)
        if ndofs_ref == 1
            MOMxBASIS_temp = [MOMxBASIS_temp]
        end
        MOMxBASIS = reshape(MOMxBASIS_temp, (ndofs_ref, ndofs_ref))
        MOMxBASIS ./= xgrid_ref[CellVolumes][1]

        ## extract quadratic matrix for interior dofs
        MOMxINTERIOR = zeros(length(idofs), length(idofs))
        for j in 1:length(idofs), k in 1:length(idofs)
            MOMxINTERIOR[j, k] = MOMxBASIS[idofs[j], idofs[k]]
        end
        moments_eval = zeros(Float64, size(MOMxBASIS, 1), ncomponents)
        moments_basis! = get_basis(ON_CELLS, FEType_ref, EG)
        MOMxBASIS = MOMxBASIS[:, idofs]
    else # interior dofs are set by moments#

        ## calculate moments times basis functions
        ndofs_moment = get_ndofs(ON_CELLS, FEType_moments, EG)
        momentbasis_vals = zeros(Tv, ndofs_moment, ncomponents)

        function momentbasis_times_refbasis(result, qpinfo)
            x = qpinfo.x
            refbasis!(refbasis_vals, x)
            moments_basis!(momentbasis_vals, x)
            for j in 1:ndofs_ref, k in 1:ndofs_moment
                result[(k - 1) * ndofs_ref + j] = dot(view(refbasis_vals, j, :), view(momentbasis_vals, k, :))
            end
            return nothing
        end

        MOMxBASIS_temp = integrate(xgrid_ref, ON_CELLS, momentbasis_times_refbasis, ndofs_ref * ndofs_moment; quadorder = 2 * order_FE)
        if ndofs_ref == 1
            MOMxBASIS_temp = [MOMxBASIS_temp]
        end
        MOMxBASIS = reshape(MOMxBASIS_temp, (ndofs_ref, ndofs_moment))
        MOMxBASIS ./= xgrid_ref[CellVolumes][1]

        ## extract quadratic matrix for interior dofs
        MOMxINTERIOR = zeros(length(idofs), size(MOMxBASIS, 2))
        for j in 1:length(idofs), k in 1:size(MOMxBASIS, 2)
            MOMxINTERIOR[j, k] = MOMxBASIS[idofs[j], k]
        end

        moments_eval = zeros(Float64, nmoments, ncomponents)
    end

    ### get permutation of dofs on reference EG and real cells
    subset_handler = get_basissubset(AT, FE, EG)
    current_subset = Array{Int, 1}(1:size(MOMxBASIS, 1))
    invA::Array{Float64, 2} = inv(MOMxINTERIOR)

    # prepare integration of moments
    L2G = L2GTransformer(EG, xgrid, AT)
    current_quadorder = order_FE + order
    QF = QuadratureRule{Tv, EG}(current_quadorder)
    edgemoments = zeros(Tv, nmoments, nitems)
    result_f = zeros(Tv, ncomponents)
    QP = QPInfos(xgrid; time = 0, kwargs...)

    function assembly_loop!(target, edgemoments, items, exact_function!, QF, L2G)
        weights, xref = QF.w, QF.xref
        nweights = length(weights)
        fill!(edgemoments, 0)
        for item in items
            QP.region = itemregions[item]
            QP.item = item
            if has_normals
                QP.normal .= view(itemnormals, :, item)
            end
            QP.volume = itemvolumes[item]
            update_trafo!(L2G, item)
            
            ## integrate moments of function
            for qp in 1:nweights
                QP.xref = xref[qp]
                eval_trafo!(QP.x, L2G, xref[qp])

                exact_function!(result_f, QP)
                if (bestapprox)
                    fill!(refbasis_vals, 0)
                    refbasis!(refbasis_vals, xref[qp])
                    for m in 1:nmoments, k in 1:ncomponents
                        edgemoments[m, item] += result_f[k] * refbasis_vals[idofs[m], k] * weights[qp]
                    end
                else
                    fill!(moments_eval, 0)
                    moments_basis!(moments_eval, xref[qp])
                    for m in 1:nmoments, k in 1:ncomponents
                        edgemoments[m, item] += result_f[k] * moments_eval[m, k] * weights[qp]
                    end
                end
            end
            for m = 1:nmoments
                edgemoments[m, item] *= itemvolumes[item]
            end
        end
        for item in items
            if subset_handler != NothingFunction
                subset_handler(current_subset, item)
            end
            for m in 1:nmoments, exdof in 1:interior_offset, c in 1:ncomponents
                localdof = coffset * (c - 1) + exdof
                edgemoments[m, item] -= target[itemdofs[localdof, item]] * MOMxBASIS[current_subset[localdof], m] * itemvolumes[item]
            end
            for m in 1:nmoments
                localdof = idofs[m]
                target[itemdofs[localdof, item]] = 0
                for n in 1:nmoments
                    target[itemdofs[localdof, item]] += invA[n, m] * edgemoments[n, item]
                end
                target[itemdofs[localdof, item]] /= itemvolumes[item]
            end
        end
        return nothing
    end

    function evaluate!(target, exact_function!, items; time = 0, quadorder = current_quadorder, params = [], bonus_quadorder = 0, kwargs...)
        new_quadorder = quadorder + bonus_quadorder
        if new_quadorder !== current_quadorder
            QF = QuadratureRule{Tv, EG}(new_quadorder)
            quadorder = new_quadorder
        end
        QP.params = params
        QP.time = time
        if isempty(items)
            items = 1:nitems
        end
        assembly_loop!(target, edgemoments, items, exact_function!, QF, L2G)
        return nothing
    end

    return MomentInterpolator(evaluate!)
end
