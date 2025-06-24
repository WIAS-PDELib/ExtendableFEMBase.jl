#####################
# SegmentIntegrator #
#####################


mutable struct SegmentIntegrator{Tv <: Real, UT, KFT <: Function, EG}
    u_args::Array{UT, 1}
    ops_args::Array{DataType, 1}
    kernel::KFT
    integrator::Any
    parameters::Dict{Symbol, Any}
end

segment_geometry(::SegmentIntegrator{Tv, UT, KFT, EG}) where {Tv, UT, KFT, EG} = EG


default_segint_kwargs() = Dict{Symbol, Tuple{Any, String}}(
    :name => ("SegmentIntegrator", "name for operator used in printouts"),
    :resultdim => (0, "dimension of result field (default = length of arguments)"),
    :matrix_mode => (false, "integrator integrates basis functions of FEspace separately to assembly a matrix that maps solution to segment integrations (requires that kernel is linear)"),
    :entry_tolerance => (0, "threshold to add entry to sparse matrix (only in matrix_mode)"),
    :params => (nothing, "array of parameters that should be made available in qpinfo argument of kernel function"),
    :factor => (1, "factor that should be multiplied during assembly"),
    :quadorder => ("auto", "quadrature order"),
    :bonus_quadorder => (0, "additional quadrature order added to quadorder"),
    :verbosity => (0, "verbosity level"),
)


function SegmentIntegrator(EG, kernel, u_args, ops_args; Tv = Float64, kwargs...)
    parameters = Dict{Symbol, Any}(k => v[1] for (k, v) in default_segint_kwargs())
    _update_params!(parameters, kwargs)
    @assert length(u_args) == length(ops_args)
    return SegmentIntegrator{Tv, typeof(u_args[1]), typeof(kernel), EG}(u_args, ops_args, kernel, nothing, parameters)
end


"""
````
function SegmentIntegrator(
	EG::ElementGeometry,
	[kernel!::Function],
	oa_args::Array{<:Tuple{<:Any, DataType},1};
	kwargs...)
````

Construct a `SegmentIntegrator` for integrating over segments of the given element geometry `EG`.

At each quadrature point, the specified operator evaluations are performed, optionally postprocessed by the provided `kernel!` function, and accumulated using the quadrature weights. If no kernel is given, the operator arguments are integrated directly.

# Arguments
- `EG::ElementGeometry`: The segment geometry over which to integrate.
- `kernel!::Function` (optional): A function of the form `kernel!(result, eval_args, qpinfo)` that postprocesses operator evaluations at each quadrature point. If omitted, a default kernel is used.
- `oa_args::Array{<:Tuple{Any, DataType},1}`: Array of tuples `(tag, FunctionOperator)`, where `tag` identifies the unknown or vector position, and `FunctionOperator` specifies the operator to evaluate.

# Keyword arguments:
$(_myprint(default_segint_kwargs()))

# Kernel Function Interface

    kernel!(result, eval_args, qpinfo)

- `result`: Preallocated array to store the kernel output.
- `eval_args`: Array of operator evaluations at the current quadrature point.
- `qpinfo`: `QPInfos` struct with information about the current quadrature point.

# Usage

After construction, call `initialize!` to prepare the integrator for a given solution, then use `integrate_segment!` to perform integration over a segment.


"""
function SegmentIntegrator(EG, kernel, oa_args::Array{<:Tuple{Union{<:Any, Int}, DataType}, 1}; kwargs...)
    u_args = [oa[1] for oa in oa_args]
    ops_args = [oa[2] for oa in oa_args]
    return SegmentIntegrator(EG, kernel, u_args, ops_args; kwargs...)
end

function SegmentIntegrator(EG, oa_args::Array{<:Tuple{Union{<:Any, Int}, DataType}, 1}; kwargs...)
    u_args = [oa[1] for oa in oa_args]
    ops_args = [oa[2] for oa in oa_args]
    return SegmentIntegrator(EG, standard_kernel, u_args, ops_args; kwargs...)
end


"""
````
function initialize!(
	O::SegmentIntegrator{T, UT},
	sol;
	time = 0,
	kwargs...)
````

Initializes the SegmentIntegrator for the specified solution.
"""
function initialize!(O::SegmentIntegrator{T, UT}, sol; time = 0, kwargs...) where {T, UT}
    _update_params!(O.parameters, kwargs)
    if UT <: Integer
        ind_args = O.u_args
    else
        ind_args = [findfirst(==(u), sol.tags) for u in O.u_args]
    end
    FES_args = [sol[j].FES for j in ind_args]
    FETypes_args = [eltype(F) for F in FES_args]

    AT = ON_CELLS
    gridAT = EffAT4AssemblyType(get_AT(FES_args[1]), AT)
    xgrid = FES_args[1].xgrid
    itemregions = xgrid[CellRegions]

    ## prepare quadrature formulae
    SG = segment_geometry(O)
    EG = xgrid[UniqueCellGeometries][1]
    dimfill = dim_element(EG) - dim_element(SG)
    @assert dimfill >= 0
    bonus_quadorder = O.parameters[:bonus_quadorder]
    if O.parameters[:quadorder] == "auto"
        polyorder = maximum([get_polynomialorder(FE, EG) for FE in FETypes_args])
        minderiv = minimum([NeededDerivative4Operator(op) for op in O.ops_args])
        quadorder = polyorder - minderiv + bonus_quadorder
    else
        quadorder = O.parameters[:quadorder] + bonus_quadorder
    end
    qf_SG = QuadratureRule{T, SG}(quadorder)
    if O.parameters[:verbosity] > 1
        @info "...... integrating on $SG with quadrature order $quadorder"
    end
    if dimfill > 0
        new_xref = Array{Array{T, 1}, 1}(undef, length(qf_SG.xref))
        for i in 1:length(qf_SG.xref)
            new_xref[i] = zeros(T, dim_element(EG))
        end
        QF = QuadratureRule{Float64, EG}(new_xref, qf_SG.w)
    else
        QF = qf_SG
    end

    ## FE basis evaluator for EG
    nargs = length(FES_args)
    BE_args = [FEEvaluator(FES_args[j], O.ops_args[j], QF; AT = AT) for j in 1:nargs]

    ## L2G map for EG
    L2G = L2GTransformer(EG, xgrid, gridAT)

    ## parameter structure
    QPinfo = QPInfos(xgrid; time = time, params = O.parameters[:params])

    FEATs_args = [EffAT4AssemblyType(get_AT(FES_args[j]), AT) for j in 1:nargs]
    itemdofs_args::Array{Union{Adjacency{Int32}, SerialVariableTargetAdjacency{Int32}}, 1} = [FES_args[j][Dofmap4AssemblyType(FEATs_args[j])] for j in 1:nargs]

    ## prepare operator infos
    op_lengths_args = [size(BE_args[j].cvals, 1) for j in 1:nargs]

    op_offsets_args = [0]
    append!(op_offsets_args, cumsum(op_lengths_args))
    resultdim::Int = O.parameters[:resultdim]
    if resultdim == 0
        resultdim = op_offsets_args[end]
        O.parameters[:resultdim] = resultdim
    end
    input_args = zeros(T, op_offsets_args[end])
    kernel_result::Array{T, 1} = zeros(T, O.parameters[:resultdim])
    ndofs_args::Array{Int, 1} = [size(BE.cvals, 2) for BE in BE_args]
    weights::Array{T, 1} = QF.w
    xrefSG = qf_SG.xref
    kernel = O.kernel
    entry_tol = O.parameters[:entry_tolerance]

    if O.parameters[:matrix_mode]
        if O.parameters[:verbosity] > 1
            @info "$(O.parameters[:name]) configured for matrix assembly"
        end
        Aloc = Vector{Matrix{T}}(undef, nargs)
        for j in 1:nargs
            Aloc[j] = zeros(T, resultdim, ndofs_args[j])
        end
        function integrator_matrix!(
                A::AbstractSparseArray{Tv},
                w::Array{Array{Tv, 1}, 1},    # world coordinates
                b::Array{Array{Tv, 1}, 1},    # barycentric coordinates (w.r.t. item geometry)
                item,                       # cell in which the segment lies (completely)
                segment_id,                 # segment number
            ) where {Tv}

            # calculate new quadrature points
            xref = BE_args[1].xref
            for i in 1:length(xref)
                fill!(xref[i], 0)
                for k in 1:length(xref[1])
                    for j in 1:(length(b) - 1)
                        xref[i][k] += xrefSG[i][j] * b[j][k]
                    end
                    xref[i][k] += (1 - sum(xrefSG[i])) * b[end][k]
                end
            end

            # update basis evaluations on new quadrature points
            for id in 1:nargs
                relocate_xref!(BE_args[id], xref)
                update_basis!(BE_args[id], item)
            end

            # compute volume of segment
            if SG <: AbstractElementGeometry1D
                QPinfo.volume = sqrt((w[1][1] - w[2][1])^2 + (w[1][2] - w[2][2])^2)
            else
                @error "This segment geometry is not implemented!"
            end

            QPinfo.region = itemregions[item]
            QPinfo.item = item

            for qp in eachindex(weights)

                update_trafo!(L2G, item)
                eval_trafo!(QPinfo.x, L2G, xref[qp])

                # evaluate operator
                fill!(input_args, 0)
                for id in 1:nargs
                    for j in 1:ndofs_args[id]
                        for d in 1:op_lengths_args[id]
                            input_args[d + op_offsets_args[id]] = BE_args[id].cvals[d, j, qp]
                        end

                        kernel(kernel_result, input_args, QPinfo)

                        # accumulate
                        for d in 1:resultdim
                            Aloc[id][d, j] += kernel_result[d] * weights[qp]
                        end
                    end
                end
            end
            ## add local matrices to global matrix
            for id in 1:nargs
                Aloc[id] .*= QPinfo.volume
                for j in 1:ndofs_args[id]
                    dof_j = itemdofs_args[id][j, item]
                    for k in 1:resultdim
                        dof_k = (segment_id - 1) * resultdim + k
                        if abs(Aloc[id][k, j]) > entry_tol
                            rawupdateindex!(A, +, Aloc[id][k, j], dof_k, dof_j)
                        end
                    end
                end
                fill!(Aloc[id], 0)
            end
            return nothing
        end
        O.integrator = integrator_matrix!
    else
        if O.parameters[:verbosity] > 1
            @info "$(O.parameters[:name]) configured for default segment evaluation"
        end
        function integrator!(
                result::AbstractArray{Tv, 1},
                w, #::Array{Array{Tv, 1}, 1},    # world coordinates
                b, #::Array{Array{Tv, 1}, 1},    # barycentric coordinates (w.r.t. item geometry)
                item, # cell in which the segment lies (completely)
            ) where {Tv}

            # calculate new quadrature points
            xref = BE_args[1].xref
            for i in 1:length(xref)
                fill!(xref[i], 0)
                for k in 1:length(xref[1])
                    for j in 1:(length(b) - 1)
                        xref[i][k] += xrefSG[i][j] * b[j][k]
                    end
                    xref[i][k] += (1 - sum(xrefSG[i])) * b[end][k]
                end
            end

            # update basis evaluations on new quadrature points
            for id in 1:nargs
                relocate_xref!(BE_args[id], xref)
                update_basis!(BE_args[id], item)
            end

            # compute volume of segment
            if SG <: AbstractElementGeometry1D
                QPinfo.volume = sqrt((w[1][1] - w[2][1])^2 + (w[1][2] - w[2][2])^2)
            else
                @error "This segment geometry is not implemented!"
            end

            QPinfo.region = itemregions[item]
            QPinfo.item = item

            fill!(result, 0)
            for qp in eachindex(weights)

                update_trafo!(L2G, item)
                eval_trafo!(QPinfo.x, L2G, xref[qp])

                # evaluate operator
                fill!(input_args, 0)
                for id in 1:nargs
                    for j in 1:ndofs_args[id]
                        dof_j = itemdofs_args[id][j, item]
                        for d in 1:op_lengths_args[id]
                            input_args[d + op_offsets_args[id]] += sol[id][dof_j] * BE_args[id].cvals[d, j, qp]
                        end
                    end
                end

                kernel(kernel_result, input_args, QPinfo)

                # accumulate
                for j in 1:resultdim
                    result[j] += kernel_result[j] * weights[qp]
                end
            end

            # multiply volume
            result .*= QPinfo.volume
            return nothing
        end
        O.integrator = integrator!
    end


    return nothing
end


"""
````
function integrate_segment!(
	result::Array{T,1},
	SI::SegmentIntegrator,
	w::Array{Array{T,1},1},
	b::Array{Array{T,1},1},
	item
	) where {T}
````

Integrate a segment with world coordinates w and barycentric coordinates b in the cell with the given item number.
"""
function integrate_segment!(
        result::Array{T, 1},
        SI::SegmentIntegrator,
        w::Array{Array{T, 1}, 1},    # world coordinates
        b::Array{Array{T, 1}, 1},    # barycentric coordinates (w.r.t. item geometry)
        item, # cell in which the segment lies (completely)
    ) where {T}

    return SI.integrator(result, w, b, item)
end
