#################
# INTERPOLATORS #
#################
#
# this file contains interpolating structures to compose full finite element interpolations
# and to be carried in the interpolators backpack of the finite element space
#
# there are three interpolators:
# - NodalInterpolator (for point evaluating degrees of freedom)
# - FunctionalInterpolator (for functional degrees of freedom and dual basis functions -> cheap evaluations)
# - MomentInterpolator (for functional degrees of freedom and non-dual basis functions -> solves local problems)

get_interpolator(FES::FESpace, entity::Type{<:AssemblyType}) = get!(() -> init_interpolator!(FES, entity), FES.interpolators, entity)

struct NodalInterpolator{EFT} <: AbstractInterpolationOperator
    evaluate!::EFT
end


"""
````
function NodalInterpolator(FES::FESpace{T}, grid = FES.dofgrid; broken = FES.broken, component_offset = FES.coffset, kwargs...)
````

Constructs an interpolator structure that has an evaluate! function that evaluates some
function at the requested nodes (items) of the grid. In broken mode only ON_CELL the items
refer to cells and all nodes of these cells are evaluated.
"""
function NodalInterpolator(
    FES::FESpace{T},
    grid = FES.dofgrid;
    broken = FES.broken,
    component_offset = FES.coffset,
    kwargs...) where {T}

    FEType = eltype(FES)
    ncomponents = get_ncomponents(FEType)
    xCoordinates = grid[Coordinates]
    xCellRegions = grid[CellRegions]
    result = zeros(T, ncomponents)
    QP = QPInfos(grid; time = 0.0)

    if broken
        ## FE space is broken
        xCellNodes = grid[CellNodes]
        xCellDofs = FES[CellDofs]
        ncells = num_cells(grid)
        function evaluate_broken!(target, exact_function!, items; time = 0, params = [], kwargs...)
            QP.time = time
            QP.params = params === nothing ? [] : params
            if isempty(items)
                items = 1 : ncells
            end
            for cell in items
                if cell < 1
                    continue
                end
                nnodes_on_cell = num_targets(xCellNodes, cell)
                QP.item = cell
                QP.cell = cell
                QP.region = xCellRegions[cell]
                for n in 1:nnodes_on_cell
                    j = xCellNodes[n, cell]
                    QP.x .= view(xCoordinates, :, j)
                    exact_function!(result, QP)
                    for k in 1:ncomponents
                        target[xCellDofs[1, cell] + n - 1 + (k - 1) * nnodes_on_cell] = result[k]
                    end
                end
            end
        end
        return NodalInterpolator(evaluate_broken!)
    else
        ## FE space is continuous in node, so only one evaluation is required
        offset4component = 0:component_offset:(ncomponents * component_offset)
        nnodes = num_nodes(grid)
        xNodeCells = atranspose(grid[CellNodes])
        function evaluate!(target, exact_function!, items; time = 0, params = [], kwargs...)
            QP.time = time
            QP.params = params === nothing ? [] : params
            if isempty(items)
                items = 1 : nnodes
            end
            for j in items
                if j < 1
                    continue
                end
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


"""
````
function MomentInterpolator(FE::FESpace{Tv, Ti, FEType, APT}, AT::Type{<:AssemblyType}, xgrid = FE.dofgrid; operator = Identity, FEType_ref = :auto, FEType_moments = :auto, moments_operator = operator, moments_dofs = Int[], bestapprox = false, order = 0, coffset::Int = -1, componentwise = true, kwargs...) where {Tv, Ti, FEType <: AbstractFiniteElement, APT}
````

Constructs an interpolation structure that has an evaluate! function that sets the interior degrees of freedom
such that the moments of the given function are preserved up to the given order. For this small local problems
are solved with a mass matrix of interior basis functions (evaluated with operator) and (by moments_dofs selected)
basis functions of the moments of type FEType_ref and with moments_operator. In the bestapprox mode
the mass matrix instead is formed from the scalar product of the interior basis functions.
"""
function MomentInterpolator(
    FE::FESpace{Tv, Ti, FEType, APT},
    AT::Type{<:AssemblyType},
    xgrid = FE.dofgrid;
    operator = Identity,
    FEType_ref = :auto,
    FEType_moments = :auto,
    moments_operator = operator,
    moments_dofs = Int[],
    bestapprox = false,
    order = 0,
    coffset::Int = -1,
    componentwise = true,
    kwargs...) where {Tv, Ti, FEType <: AbstractFiniteElement, APT}

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
    ncomponents::Int = get_ncomponents(FEType)
    edim::Int = dim_element(EG)
    order_FE = get_polynomialorder(FEType, EG)
    if coffset == -1
        coffset = FEType <: AbstractH1FiniteElement ? Int(get_ndofs(AT, FEType, EG) / ncomponents) : Int(0)
    end
    interior_offset = interior_dofs_offset(AT, FEType, EG)
    @assert interior_offset >= 0 "This FEType seems to missing a definition for interior_dofs_offset!"

    ## reference basis for FE on EG
    ## here we assume that the FEType looks like a H1Pk element on EG (which is true for all H1Pk elements)
    if FEType_ref == :auto
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

    ## get basis for moments
    if FEType_moments == :auto
        if order == 0
            FEType_moments = L2P0{ncomponents}
        elseif order == 1
            FEType_moments = H1P1{ncomponents}
        elseif order == 2
            FEType_moments = H1P2{ncomponents, edim}
        else
            FEType_moments = H1Pk{ncomponents, edim, order}
        end
    end

    ## check if mass matrix can be computed once or needs to be recomputed on every mesh
    if FEType_ref <: AbstractH1FiniteElement &&
        !(FEType_ref <: AbstractH1FiniteElementWithCoefficients) &&
        FEType_moments <: AbstractH1FiniteElement &&
         moments_operator == Identity
        fixed_mass_matrix = true
    else
        fixed_mass_matrix = false
    end

    moments_basis! = get_basis(AT, FEType_moments, EG)
    nmoments::Int = get_ndofs_all(AT, FEType_moments, EG)
    if isempty(moments_dofs)
        moments_dofs = Array{Int,1}(1:nmoments)
    else
        nmoments = length(moments_dofs)
    end
    xgrid_ref = reference_domain(EG)
    idofs = zeros(Int, 0)
    if bestapprox 
        FEType_moments = FEType_ref
        append!(idofs, interior_offset+1:get_ndofs(AT, FEType_ref, EG))
        nmoments = length(idofs)
    else
        if FEType_ref <: AbstractH1FiniteElement && componentwise
            nmoments4c::Int = nmoments / ncomponents
            for c in 1:ncomponents, m in 1:nmoments4c
                push!(idofs, (c - 1) * coffset + interior_offset + m)
            end
        else
            ## FE is assumed to be vector-valued
            for m in 1:nmoments
                push!(idofs, interior_offset + m)
            end
        end
    end

    MOMxBASIS::Array{Float64, 2} = zeros(Float64, 0, 0)
    moments_eval::Matrix{Float64} = zeros(Float64, 0, 0)
    refbasis! = get_basis(AT, FEType_ref, EG)
    ndofs_ref = get_ndofs_all(AT, FEType_ref, EG)
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

        ## integrate ON_CELLS here since we are on the refenrece domain of EG!
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
        moments_basis! = get_basis(AT, FEType_ref, EG)
        MOMxBASIS = MOMxBASIS[:, idofs]
    else # interior dofs are set by preserving moments

        ## calculate moments times basis functions
        ndofs_moment = get_ndofs_all(AT, FEType_moments, EG)
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
        MOMxINTERIOR = zeros(length(idofs), nmoments)
        for j in 1:length(idofs), k in 1:nmoments
            MOMxINTERIOR[j, k] = MOMxBASIS[idofs[j], moments_dofs[k]]
        end

        moments_eval = zeros(Float64, nmoments, ncomponents)
    end

    ### get permutation of dofs on reference EG and real cells
    subset_handler = get_basissubset(AT, FE, EG)
    current_subset = Array{Int, 1}(1:size(MOMxBASIS, 1))

    ## inverse of MOMxINTERIOR
    ## (here we assume that it is the same on every cell, which is true for H1Pk elements)
    if fixed_mass_matrix
        invA::Array{Float64, 2} = inv(MOMxINTERIOR)
    end

    # prepare integration of moments
    L2G = L2GTransformer(EG, xgrid, AT)
    current_quadorder = 2*order_FE
    QF = QuadratureRule{Tv, EG}(current_quadorder)
    f_moments = zeros(Tv, nmoments)
    result_f = zeros(Tv, ncomponents)
    QP = QPInfos(xgrid; time = 0.0, kwargs...)

    # prepare mass matrix integration
    FEB = FEEvaluator(FE, operator, QF; AT = AT, T = Tv)
    if bestapprox
        FEB_moments = FEB
    else
        FE_moments = FESpace{FEType_moments}(xgrid)
        FEB_moments = FEEvaluator(FE_moments, moments_operator, QF; AT = AT, T = Tv)
    end
    basisval = zeros(Tv, ncomponents)
    interiordofs = zeros(Int, length(idofs))

    function assembly_loop!(target, f_moments, items, exact_function!, QF, L2G, FEB, FEB_moments)
        weights, xref = QF.w, QF.xref
        nweights = length(weights)
        for item::Int in items
            if item < 1
                continue
            end
            QP.region = itemregions[item]
            QP.item = item
            if has_normals
                QP.normal .= view(itemnormals, :, item)
            end
            QP.volume = itemvolumes[item]
            update_trafo!(L2G, item)
            if moments_operator !== Identity
                update_basis!(FEB_moments, item)
            end
            
            ## integrate moments of function
            fill!(f_moments, 0)
            for qp in 1:nweights
                QP.xref = xref[qp]
                eval_trafo!(QP.x, L2G, xref[qp])

                exact_function!(result_f, QP)
                if (bestapprox)
                    for m in 1:nmoments, k in 1:ncomponents
                        f_moments[m] += result_f[k] * FEB.cvals[k, idofs[m], qp] * weights[qp]
                    end
                else
                    for m in 1:nmoments, k in 1:ncomponents
                        f_moments[m] += result_f[k] * FEB_moments.cvals[k, moments_dofs[m], qp] * weights[qp]
                    end
                end
            end

            ## solve linear system for free/interior dofs
            if fixed_mass_matrix
                if subset_handler != NothingFunction
                    subset_handler(current_subset, item)
                end
                ## subtract moments of fixed dofs
                for m in 1:nmoments, exdof in 1:interior_offset, c in 1:ncomponents
                    localdof = coffset * (c - 1) + exdof
                    f_moments[m] -= target[itemdofs[localdof, item]] * MOMxBASIS[current_subset[localdof], m]
                end
                ## solve for free dofs
                for m in 1:nmoments
                    localdof = idofs[m]
                    target[itemdofs[localdof, item]] = 0
                    for n in 1:nmoments
                        target[itemdofs[localdof, item]] += invA[n, m] * f_moments[n]
                    end
                end
            else
                update_basis!(FEB, item)
                ## subtract moments of fixed dofs
                for dof in 1:interior_offset
                    for i in 1:nweights
                        eval_febe!(basisval, FEB, dof, i)
                        for m in 1:nmoments, k = 1:ncomponents
                            f_moments[m] -= basisval[k] * FEB_moments.cvals[k, moments_dofs[m], i] * target[itemdofs[dof, item]] * weights[i]
                        end
                    end
                end
                ## recompute mass matrix of interior dofs
                fill!(MOMxINTERIOR, 0)
                if bestapprox
                    for dof in 1:nmoments
                        for i in 1:nweights
                            for m in 1:nmoments, k in 1:ncomponents
                                MOMxINTERIOR[m, dof] += FEB.cvals[k, idofs[dof], i] * FEB.cvals[k, idofs[m], i] * weights[i]
                            end
                        end
                        interiordofs[dof] = itemdofs[idofs[dof], item]
                    end
                else
                    for dof in 1:nmoments
                        for i in 1:nweights
                            for m in 1:nmoments, k in 1:ncomponents
                                MOMxINTERIOR[m, dof] += FEB.cvals[k, idofs[dof], i] * FEB_moments.cvals[k, moments_dofs[m], i] * weights[i]
                            end
                        end
                        interiordofs[dof] = itemdofs[idofs[dof], item]
                    end
                end
                ## solve for free dofs
                target[interiordofs] = MOMxINTERIOR \ f_moments
            end
        end

        return nothing
    end

    function evaluate!(target, exact_function!, items; time = 0, quadorder = current_quadorder, params = [], bonus_quadorder = 0, kwargs...)
        new_quadorder = quadorder + bonus_quadorder
        if new_quadorder !== current_quadorder
            QF = QuadratureRule{Tv, EG}(new_quadorder)
            quadorder = new_quadorder
            FEB = FEEvaluator(FE, operator, QF; T = Tv)
            FEB_moments = FEEvaluator(FE_moments, moments_operator, QF; T = Tv)
        end
        QP.params = params === nothing ? [] : params
        QP.time = time
        if isempty(items)
            items = 1:nitems
        end
        assembly_loop!(target, f_moments, items, exact_function!, QF, L2G, FEB, FEB_moments)
        return nothing
    end

    return MomentInterpolator(evaluate!)
end

struct FunctionalInterpolator{EFT} <: AbstractInterpolationOperator
    evaluate!::EFT
end


"""
````
function FunctionalInterpolator(
    functionals!::Function,
    FE::FESpace{Tv, Ti, FEType, APT},
    AT::Type{<:AssemblyType} = ON_FACES,
    xgrid = FE.dofgrid;
    operator = NormalFlux, nfluxes = 0, dofs = [], kwargs...) where {Tv, Ti, FEType <: AbstractFiniteElement, APT}
````

Constructs an interpolation structure that has an evaluate! function that determines the interior degrees of freedom
(or the specified local dofs) by evaluating functionals! (the result dimension nfluxes should correspond to the number of dofs,
both are set to the number of interior dofs as defaults).
The functionals are corrected by the operator evaluations of the fixed dofs.
The mean parameter decides if the functional is divided by the area in the end (if mean = true).
"""
function FunctionalInterpolator(
    functionals!::Function,
    FE::FESpace{Tv, Ti, FEType, APT},
    AT::Type{<:AssemblyType} = ON_FACES,
    xgrid = FE.dofgrid;
    bonus_quadorder = 0,
    operator = NormalFlux,
    nfluxes = 0,
    dofs = [],
    mean = false,
    kwargs...) where {Tv, Ti, FEType <: AbstractFiniteElement, APT}
    
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

    # prepare integration of moments
    if nfluxes == 0
        if isempty(dofs)
            nfluxes = max_num_targets_per_source(itemdofs)
        else
            nfluxes = length(dofs)
        end
    end
    if isempty(dofs)
        dofs = 1 : nfluxes
    end
    ncomponents::Int = get_ncomponents(FEType)
    order_FE = get_polynomialorder(FEType, EG)
    L2G = L2GTransformer(EG, xgrid, AT)
    current_quadorder = order_FE + bonus_quadorder
    QF = QuadratureRule{Tv, EG}(current_quadorder)
    f_fluxes = zeros(Tv, nfluxes)
    result_f = zeros(Tv, ncomponents)
    QP = QPInfos(xgrid; time = 0.0, kwargs...)
    nitems = size(itemnodes, 2)

    ## prepare evaluation of fixed dofs
    interior_offset = interior_dofs_offset(AT, FEType, EG)
    FEB = FEEvaluator(FE, operator, QF; AT = AT, T = Tv, L2G = L2G)

    function assembly_loop!(target, f_fluxes, items, exact_function!, QF, L2G, FEB)
        weights, xref = QF.w, QF.xref
        nweights = length(weights)
        for item::Int in items
            if item < 1
                continue
            end
            for m = 1:nfluxes
                target[itemdofs[dofs[m], item]] = 0
            end
            QP.region = itemregions[item]
            QP.item = item
            if has_normals
                QP.normal .= view(itemnormals, :, item)
            end
            QP.volume = itemvolumes[item]
            update_trafo!(L2G, item)
            if interior_offset > 0
                update_basis!(FEB, item)
            end
            
            ## compute fluxes of function
            for qp in 1:nweights
                fill!(f_fluxes, 0)
                QP.xref = xref[qp]
                eval_trafo!(QP.x, L2G, xref[qp])
                exact_function!(result_f, QP)
                functionals!(f_fluxes, result_f, QP)

                ## subtract flux of fixed dofs
                if interior_offset > 0
                    for m = 1 : nfluxes, dof = 1 : interior_offset
                        f_fluxes[m] -= target[itemdofs[dof, item]] * FEB.cvals[m, dof, qp] 
                    end
                end

                weight = weights[qp]
                if !mean
                    weight *= itemvolumes[item]
                end

                ## set fluxes to dofs
                for m = 1:nfluxes
                    target[itemdofs[dofs[m], item]] += f_fluxes[m] * weight
                end
            end
        end
        return nothing
    end

    function evaluate!(target, exact_function!, items; time = 0, quadorder = current_quadorder, params = [], bonus_quadorder = 0, kwargs...)
        new_quadorder = quadorder + bonus_quadorder
        if new_quadorder !== current_quadorder
            QF = QuadratureRule{Tv, EG}(new_quadorder)
            FEB = FEEvaluator(FE, operator, QF; AT = AT, T = Tv)
            quadorder = new_quadorder
        end
        QP.params = params === nothing ? [] : params
        QP.time = time
        if isempty(items)
            items = 1:nitems
        end
        assembly_loop!(target, f_fluxes, items, exact_function!, QF, L2G, FEB)
        return nothing
    end

    return FunctionalInterpolator(evaluate!)
end



function slice(VTA::VariableTargetAdjacency, items = [], only_unique::Bool = true)
    subitems = zeros(Int, 0)
    if items == []
        items = 1:num_sources(VTA)
        subitems = VTA.colentries
    else
        for item in items
            append!(subitems, VTA[:, item])
        end
        if only_unique
            subitems = unique(subitems)
        end
    end
    return subitems
end

function slice(VTA::Array{<:Signed, 2}, items = [], only_unique::Bool = true)
    #=
	subitems = zeros(Int,0)
	if items == []
		items = 1 : size(VTA,2)
	end
	for item in items
		append!(subitems, VTA[:,item])
	end
	if only_unique
		subitems = unique(subitems)
	end
	=#
    return unique(view(VTA, :, items))
end
