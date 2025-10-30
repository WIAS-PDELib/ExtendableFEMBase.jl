mutable struct PointEvaluator{Tv <: Real, TCoeff <: Real, UT, KFT <: Function}
    u_args::Array{UT, 1}
    ops_args::Array{DataType, 1}
    kernel::KFT
    BE_args::Any
    L2G::Any
    CF::Any
    lastitem::Int
    eval_selector::Any
    evaluator_bary::Any
    evaluator::Any
    xref::Vector{Tv}
    parameters::Dict{Symbol, Any}
end

function Base.show(io::IO, PE::PointEvaluator)
    println(io, "PointEvaluator")
    println(io, "--------------")
    println(io, "Unknowns: ", PE.u_args)
    println(io, "Operators: ", PE.ops_args)
    println(io, "Kernel function: ", typeof(PE.kernel))
    println(io, "Parameters: ", PE.parameters)
    return
end

default_peval_kwargs() = Dict{Symbol, Tuple{Any, String}}(
    :name => ("PointEvaluator", "name for operator used in printouts"),
    :resultdim => (0, "dimension of result field (default = length of operators)"),
    :params => (nothing, "array of parameters that should be made available in qpinfo argument of kernel function"),
    :verbosity => (0, "verbosity level"),
)


"""
$(TYPEDSIGNATURES)

Construct a `PointEvaluator` object for evaluating operator expressions at arbitrary points in a finite element space.

A `PointEvaluator` can be used to evaluate one or more operator evaluations (e.g., function values, gradients) at arbitrary points, optionally postprocessed by a user-supplied kernel function. The evaluation is performed with respect to a given solution and its finite element basis.

After construction, the `PointEvaluator` must be initialized with a solution using `initialize!`. Evaluation at a point is then performed using `evaluate!` or `evaluate_bary!`.

# Arguments
- `kernel!` (optional): Postprocessing function for operator evaluations. Should have the form `kernel!(result, eval_args, qpinfo)`.
- `oa_args`: Array of operator argument tuples `(source block tag, operator type)`.
- `sol` (optional): Solution object for immediate initialization.

# Keyword arguments:
$(_myprint(default_peval_kwargs()))

# Kernel Function Interface

    kernel!(result, eval_args, qpinfo)

- `result`: Preallocated array to store the kernel output.
- `eval_args`: Array of operator evaluations at the current evaluation point.
- `qpinfo`: `QPInfos` struct with information about the current evaluation point.

# Usage

After construction, call `initialize!` to prepare the evaluator for a given solution, then use `evaluate!` or `evaluate_bary!` to perform point evaluations.

"""
function PointEvaluator(kernel, u_args, ops_args, sol = nothing; Tv = Float64, TCoeff = Float64, kwargs...)
    parameters = Dict{Symbol, Any}(k => v[1] for (k, v) in default_peval_kwargs())
    _update_params!(parameters, kwargs)
    @assert length(u_args) == length(ops_args)
    PE = PointEvaluator{Tv, TCoeff, typeof(u_args[1]), typeof(kernel)}(u_args, ops_args, kernel, nothing, nothing, nothing, 1, nothing, nothing, nothing, zeros(Tv, 2), parameters)
    if sol !== nothing
        initialize!(PE, sol)
    end
    return PE
end

function PointEvaluator(kernel::Function, oa_args::Array{<:Tuple{<:Any, DataType}, 1}, sol = nothing; kwargs...)
    u_args = [oa[1] for oa in oa_args]
    ops_args = [oa[2] for oa in oa_args]
    return PointEvaluator(kernel, u_args, ops_args, sol; kwargs...)
end

function PointEvaluator(oa_args::Array{<:Tuple{<:Any, DataType}, 1}, sol = nothing; kwargs...)
    return PointEvaluator(standard_kernel, oa_args, sol; kwargs...)
end


"""
````
function initialize!(
	O::PointEvaluator,
	sol;
	time = 0,
	kwargs...)
````

Initializes the given `PointEvaluator` for a specified solution (FEVector or vector of FEVectorBlocks).

This function prepares the `PointEvaluator` for evaluation by associating it with the provided solution vector. It sets up the necessary finite element basis evaluators, local-to-global transformations, and cell finder structures for the underlying grid.

# Arguments
- `O::PointEvaluator`: The `PointEvaluator` instance to initialize.
- `sol`: The solution object (e.g., array of FEVectorBlocks) to be used for evaluations.

# Keyword Arguments
- `time`: (default: `0`) Time value to be passed to the quadrature point info structure.
$(_myprint(default_peval_kwargs()))

# Notes
- This function must be called before using `evaluate!` or `evaluate_bary!` with the `PointEvaluator`.
Initializes the given `PointEvaluator` for a specified solution (FEVector or vector of FEVectorBlocks).

This function prepares the `PointEvaluator` for evaluation by associating it with the provided solution vector. It sets up the necessary finite element basis evaluators, local-to-global transformations, and cell finder structures for the underlying grid.

# Arguments
- `O::PointEvaluator`: The `PointEvaluator` instance to initialize.
- `sol`: The solution object (e.g., array of FEVectorBlocks) to be used for evaluations.

# Keyword Arguments
- `time`: (default: `0`) Time value to be passed to the quadrature point info structure.
$(_myprint(default_peval_kwargs()))

# Notes
- This function must be called before using `evaluate!` or `evaluate_bary!` with the `PointEvaluator`.
"""
function initialize!(O::PointEvaluator{T, TCoeff, UT}, sol; time = 0, kwargs...) where {T, TCoeff, UT}
    _update_params!(O.parameters, kwargs)
    if UT <: Integer
        ind_args = O.u_args
    else
        ind_args = [findfirst(==(u), sol.tags) for u in O.u_args]
    end
    FES_args = [sol[j].FES for j in ind_args]
    nargs = length(FES_args)
    xgrid = FES_args[1].xgrid
    Ti = eltype(xgrid[CellNodes])
    EGs = xgrid[UniqueCellGeometries]
    AT = ON_CELLS
    gridAT = EffAT4AssemblyType(get_AT(FES_args[1]), AT)
    xgrid = FES_args[1].xgrid
    itemregions = xgrid[CellRegions]
    itemgeometries = xgrid[CellGeometries]

    O.CF = CellFinder(xgrid)
    O.xref = zeros(T, size(xgrid[Coordinates], 1))

    O.BE_args = Array{Array{<:FEEvaluator, 1}, 1}([])
    O.L2G = []
    for EG in EGs
        ## FE basis evaluator for EG
        push!(O.BE_args, [FEEvaluator(FES_args[j], O.ops_args[j], QuadratureRule{T, EG}(0); AT = AT) for j in 1:nargs])

        ## L2G map for EG
        push!(O.L2G, L2GTransformer(EG, xgrid, gridAT))
    end

    ## parameter structure
    QPinfo = QPInfos(xgrid; time = time, params = O.parameters[:params])

    ## prepare input args
    op_lengths_args = [size(O.BE_args[1][j].cvals, 1) for j in 1:nargs]
    op_offsets_args = [0]
    append!(op_offsets_args, cumsum(op_lengths_args))
    input_args = zeros(TCoeff, op_offsets_args[end])

    FEATs_args = [EffAT4AssemblyType(get_AT(FES_args[j]), AT) for j in 1:nargs]
    itemdofs_args::Array{Union{Adjacency{Ti}, SerialVariableTargetAdjacency{Ti}}, 1} = [FES_args[j][Dofmap4AssemblyType(FEATs_args[j])] for j in 1:nargs]
    kernel = O.kernel

    function eval_selector(item)
        return findfirst(==(itemgeometries[item]), EGs)
    end

    function _evaluate_bary!(
            result,
            BE_args::Array{<:FEEvaluator, 1},
            L2G::L2GTransformer,
            xref,
            item, # cell used to evaluate local coordinates
        )

        for id in 1:nargs
            # update basis evaluations at xref
            relocate_xref!(BE_args[id], xref)

            # update operator evaluation on item
            update_basis!(BE_args[id], item)
        end

        # update QPinfo
        QPinfo.item = item
        QPinfo.region = itemregions[item]
        QPinfo.xref = xref
        update_trafo!(L2G, item)
        eval_trafo!(QPinfo.x, L2G, xref)
        # evaluate operator
        fill!(input_args, 0)
        for id in 1:nargs
            for j in 1:size(BE_args[id].cvals, 2)
                dof_j = itemdofs_args[id][j, item]
                for d in 1:op_lengths_args[id]
                    input_args[d + op_offsets_args[id]] += sol[ind_args[id]][dof_j] * BE_args[id].cvals[d, j, 1]
                end
            end
        end

        ## evaluate kernel
        kernel(result, input_args, QPinfo)

        return nothing
    end

    ## initialize cell finder
    CF = CellFinder(xgrid)
    xref = zeros(T, size(xgrid[Coordinates], 1))
    function _evaluate!(
            result,
            BE_args::Array{<:FEEvaluator, 1},
            L2G::L2GTransformer,
            x,
        )


        # evaluate in barycentric coordinates
        _evaluate_bary!(result, BE_args, L2G, xref, item)

        return nothing
    end
    O.evaluator = _evaluate!
    O.evaluator_bary = _evaluate_bary!
    O.eval_selector = eval_selector

    return nothing
end


"""
````
function evaluate_bary!(
	result,
	PE::PointEvaluator,
	xref, 
	item
	)
````

Evaluates the PointEvaluator at the specified reference coordinates in the cell with the specified item number.
"""
function evaluate_bary!(
        result,
        PE::PointEvaluator,
        xref,
        item,
    )

    ## find cell geometry id
    j = PE.eval_selector(item)

    ## evaluate
    return PE.evaluator_bary(result, PE.BE_args[j], PE.L2G[j], xref, item)
end

"""
````
function evaluate!(
	result,
	PE::PointEvaluator,
	x
	)
````

Evaluates the PointEvaluator at the specified coordinates x.
(To do so it internally calls CellFinder to find the cell and the barycentric
coordinates of x and calls evaluate_bary!.)
"""
function evaluate!(
        result,
        PE::PointEvaluator,
        x;
        kwargs...,
    )
    # find correct cell (start from cell of last evaluation)
    item = gFindLocal!(PE.xref, PE.CF, x; icellstart = PE.lastitem, kwargs...)
    @assert item > 0
    PE.lastitem = item

    ## find cell geometry id
    j = PE.eval_selector(item)

    ## evaluate
    return PE.evaluator_bary(result, PE.BE_args[j], PE.L2G[j], PE.xref, item)
end


"""
````
function eval_func_bary(PE::PointEvaluator)
````

Yields the function (result, xref, item) -> evaluate_bary!(result,PE,xref,item).
"""
function eval_func_bary(PE::PointEvaluator)
    return (result, xref, item) -> evaluate_bary!(result, PE, xref, item)
end

"""
````
function eval_func(PE::PointEvaluator)
````

Yields the function (result, x) -> evaluate!(result,PE,x).
"""
function eval_func(PE::PointEvaluator; kwargs...)
    return (result, x) -> evaluate!(result, PE, x; kwargs...)
end
