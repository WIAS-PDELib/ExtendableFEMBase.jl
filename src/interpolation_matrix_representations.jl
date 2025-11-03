"""
````
function compute_lazy_interpolation_jacobian(
    target_space::FESpace, 
    source_space::FESpace; 
    use_cellparents::Bool = false,
    kwargs...
)
````

Compute the Jacobian of the [`lazy_interpolate!`](@ref) call with respect to the
`source_space` degrees of freedom, i.e. for functions ``v = \\sum_j \\alpha_j \\, \\varphi_j`` of the
`source_space` and the interpolation operator ``I(v) = \\sum_k L_k(v)\\,\\phi_k = \\sum_k L_k\\left(\\sum_j \\alpha_j \\varphi_j\\right) \\, \\phi_k``
into the `target_space`, this function computes the jacobian ``\\left[\\frac{\\partial L_k}{\\partial \\alpha_j}\\right]_{k,\\,j}``
and returns its sparse matrix representation.

# Arguments
- `target_space::FESpace`: Finite element space into which the interpolation ``I(v)`` is directed.
- `source_space::FESpace`: Finite element space from which ``v`` is taken.

# Keyword Arguments
- `use_cellparents`: Use parent cell information if available (can speed up the calculation if the `target_space` is defined on a subgrid of `source_space`).
- `kwargs...`: Additional keyword arguments passed to lower-level `lazy_interpolate!` call.

# Notes
- Since [`lazy_interpolate!`](@ref) is based on evaluating functions from the `source_space`
using a [`ExtendableFEMBase.PointEvaluator`](@ref), this should be used carefully on finer 
grids as this is not the most efficient method, but will work out of the box for any 
`source` and `target` spaces.
- This function can be used for computing prolongation or restriction operators if the `FESpace`s are defined on coarser/finer grids, respectively.

"""
function compute_lazy_interpolation_jacobian(target_space::FESpace, source_space::FESpace; use_cellparents::Bool = false, kwargs...)
    # DifferentiationInterface.jacobian needs a function of signature
    # AbstractVector -> AbstractVector
    function do_interpolation(source_vector::AbstractVector; use_cellparents = use_cellparents)
        T = valtype(source_vector)
        target_vector = FEVector{T}(target_space)[1]
        source_FE_Vector = FEVector{T}(source_space)
        source_FE_Vector.entries .= source_vector

        lazy_interpolate!(target_vector, source_FE_Vector, [(1, Identity)]; use_cellparents, kwargs...)

        return target_vector.entries
    end

    n = ndofs(source_space)

    dense_forward_backend = AutoForwardDiff()
    sparse_forward_backend = AutoSparse(
        dense_forward_backend;
        sparsity_detector = TracerSparsityDetector(),
        coloring_algorithm = GreedyColoringAlgorithm(),
    )

    M = jacobian(x -> do_interpolation(x; use_cellparents), sparse_forward_backend, ones(n))

    return M
end

"""
````
function H1Pk_to_HDIVRT0_interpolator(
    VRT0::FESpace,
    V1::FESpace
)
````

Computes the interpolation matrix from the [`H1Pk`](@ref)-conforming source space `V1`
into the [`HDIVRT0`](@ref)-conforming target space `VRT0` *defined on the same grid* and
returns its sparse matrix representation.

# Arguments
- `VRT0::FESpace`: [`HDIVRT0`](@ref) target space into which the interpolation is directed.
- `V1::FESpace`: [`H1Pk`](@ref) source space from which the interpolant is taken.

"""
function H1Pk_to_HDIVRT0_interpolator(VRT0::FESpace, V1::FESpace)
    FEType = get_FEType(V1)
    facedofs_V1::VariableTargetAdjacency{Int32} = V1[FaceDofs]
    dim = get_ncomponents(V1)
    EGF = dim == 2 ? Edge1D : Triangle2D
    ndofs = get_ndofs(ON_FACES, FEType, EGF)
    order = get_polynomialorder(FEType, EGF)
    face_basis = get_basis(ON_FACES, FEType, EGF)
    result = zeros(ndofs, dim)
    ref_integrate!(result, EGF, order, face_basis)

    F0 = FEMatrix(V1, VRT0)
    FE::ExtendableSparseMatrix{Float64, Int64} = F0.entries
    xgrid = V1.xgrid
    @assert xgrid == VRT0.xgrid
    xFaceVolumes = xgrid[FaceVolumes]
    xFaceNormals = xgrid[FaceNormals]
    nfaces = num_sources(xFaceNormals)
    for face in 1:nfaces
        for dof1 in 1:ndofs
            FE[facedofs_V1[dof1, face], face] = dot(view(result, dof1, :), view(xFaceNormals, :, face)) * xFaceVolumes[face]
        end
    end
    flush!(FE)
    return F0
end
