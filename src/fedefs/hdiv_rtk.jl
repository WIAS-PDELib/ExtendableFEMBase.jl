"""
````
abstract type HDIVRTk{edim, order} <: AbstractHdivFiniteElement where {edim<:Int}
````

Hdiv-conforming vector-valued (ncomponents = edim) Raviart-Thomas space of arbitrary order.

allowed ElementGeometries:
- Triangle2D
"""
abstract type HDIVRTk{edim, order} <: AbstractHdivFiniteElement where {edim <: Int, order <: Int} end
HDIVRTk(edim::Int, order::Int) = HDIVRTk{edim, order}
HDIVRTk(; edim::Int, order::Int) = HDIVRTk{edim, order}

function Base.show(io::Core.IO, FEType::Type{<:HDIVRTk{edim, order}}) where {edim, order}
    return print(io, "HDIVRTk{$edim, $order}")
end

get_ncomponents(FEType::Type{<:HDIVRTk}) = FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:HDIVRTk{edim, order}}, EG::Type{<:AbstractElementGeometry1D}) where {edim, order} = order + 1
get_ndofs(::Type{ON_CELLS}, FEType::Type{<:HDIVRTk{edim, order}}, EG::Type{<:Triangle2D}) where {edim, order} = (order + 1) * num_faces(EG) + order * (order + 1)

get_polynomialorder(::Type{<:HDIVRTk{2, order}}, ::Type{<:AbstractElementGeometry1D}) where {order} = order;
get_polynomialorder(::Type{<:HDIVRTk{2, order}}, ::Type{<:AbstractElementGeometry2D}) where {order} = order + 1;

get_dofmap_pattern(FEType::Type{<:HDIVRTk{2, order}}, ::Type{CellDofs}, EG::Type{<:Triangle2D}) where {order} = "f$(order + 1)" * (order > 0 ? "i$(Int((order + 1) * (order)))" : "")
get_dofmap_pattern(FEType::Type{<:HDIVRTk{2, order}}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry1D}) where {order} = "i$(order + 1)"

isdefined(FEType::Type{<:HDIVRTk}, ::Type{<:Triangle2D}) = true

interior_dofs_offset(::Type{<:ON_CELLS}, ::Type{<:HDIVRTk{2, order}}, ::Type{<:Triangle2D}) where {order} = 3 * (order + 1)


function RTk_normalflux_eval!(order)
    moments_weights = basis.(ShiftedLegendre, (0:order))
    function closure(result, f, qpinfo)
        result[1] = dot(f, qpinfo.normal)
        for j in 1:order
            result[j + 1] = result[1] * moments_weights[j + 1](qpinfo.xref[1])
        end
    end
end
init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}) where {Tv, Ti, FEType <: HDIVRTk, APT} = FunctionalInterpolator(RTk_normalflux_eval!(FEType.parameters[2]), FES, ON_FACES; bonus_quadorder = FEType.parameters[2])
init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}) where {Tv, Ti, FEType <: HDIVRTk, APT} = MomentInterpolator(FES, ON_CELLS; order = FEType.parameters[2] - 1)


function ExtendableGrids.interpolate!(Target::AbstractArray{T, 1}, FE::FESpace{Tv, Ti, HDIVRTk{edim, order}, APT}, ::Type{ON_FACES}, exact_function!; items = [], kwargs...) where {T, Tv, Ti, edim, order, APT}
    get_interpolator(FE, ON_FACES).evaluate!(Target, exact_function!, items; kwargs...)
end

function ExtendableGrids.interpolate!(Target::AbstractArray{T, 1}, FE::FESpace{Tv, Ti, HDIVRTk{edim, order}, APT}, ::Type{ON_CELLS}, exact_function!; items = [], kwargs...) where {T, Tv, Ti, edim, order, APT}
    # delegate cell faces to face interpolation
    subitems = slice(FE.dofgrid[CellFaces], items)
    interpolate!(Target, FE, ON_FACES, exact_function!; items = subitems, kwargs...)

    if order == 0
        return nothing
    end

    # set values of interior functions such that moments of Pk are preserved
    get_interpolator(FE, ON_CELLS).evaluate!(Target, exact_function!, items; kwargs...)
    return
end

# only normalfluxes on faces
function get_basis(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, ::Type{<:HDIVRTk{2, order}}, ::Type{<:AbstractElementGeometry1D}) where {order}
    moments_weights = basis.(ShiftedLegendre, (0:order))
    return function closure(refbasis, xref)
        refbasis[1, 1] = 1
        for j in 1:order
            refbasis[j + 1, 1] = (2 * j + 1) * moments_weights[j + 1](xref[1])
        end
        return
    end
end

function get_basis(::Type{ON_CELLS}, ::Type{HDIVRTk{2, order}}, ::Type{<:Triangle2D}) where {order}
    if order == 1
        interior_basis = get_basis(ON_CELLS, L2P0{1}, Triangle2D)
        ninterior = 1
    elseif order > 0
        interior_basis = get_basis(ON_CELLS, H1Pk{1, 2, order - 1}, Triangle2D)
        ninterior = Int((order + 1) * (order) / 2)
    else
        ninterior = 0
    end
    interior_offset = 3 * (order + 1)
    moments_weights = basis.(ShiftedLegendre, (0:order))
    moment_factors = [convert(Polynomial, m).coeffs[end] for m in moments_weights]
    return function closure(refbasis, xref)
        # RT0 basis
        refbasis[1, 1] = xref[1]
        refbasis[1, 2] = xref[2] - 1
        refbasis[2 + order, 1] = xref[1]
        refbasis[2 + order, 2] = xref[2]
        refbasis[3 + 2 * order, 1] = xref[1] - 1
        refbasis[3 + 2 * order, 2] = xref[2]
        for k in 1:2
            for j in 1:order
                refbasis[1 + j, k] = (-1)^j * (2 * j + 1) * moments_weights[j + 1](1 - xref[1] - xref[2]) * refbasis[1, k]
                refbasis[2 + j + order, k] = (-1)^j * (2 * j + 1) * moments_weights[j + 1](xref[1]) * refbasis[2 + order, k]
                refbasis[3 + j + 2 * order, k] = (-1)^j * (2 * j + 1) * moments_weights[j + 1](xref[2]) * refbasis[3 + 2 * order, k]
            end

            # interior RT2 functions (RT1 interior functions times P1) = 6 dofs
            if order > 0
                interior_basis(view(refbasis, (interior_offset + ninterior + 1):(interior_offset + 2 * ninterior), 2), xref)
                for j in 1:ninterior
                    refbasis[interior_offset + j, k] = 12 * xref[2] * refbasis[1, k] * refbasis[interior_offset + ninterior + j, 2]
                    refbasis[interior_offset + ninterior + j, k] = 12 * xref[1] * refbasis[3 + 2 * order, k] * refbasis[interior_offset + ninterior + j, 2]
                end
            end
        end
        return
    end
end


function get_coefficients(::Type{ON_CELLS}, FE::FESpace{Tv, Ti, <:HDIVRTk{2, order}, APT}, EG::Type{<:Triangle2D}, xgrid) where {Tv, Ti, APT, order}
    xCellFaceSigns::Union{VariableTargetAdjacency{Int32}, Array{Int32, 2}} = xgrid[CellFaceSigns]
    nfaces::Int = num_faces(EG)
    dim::Int = dim_element(EG)
    return function closure(coefficients::Array{<:Real, 2}, cell::Int)
        fill!(coefficients, 1.0)
        # multiplication with normal vector signs (only RT0, RT2, RT4 etc.)
        for j in 1:nfaces, k in 1:dim
            coefficients[k, (order + 1) * j - order] = xCellFaceSigns[j, cell]
            for o in 2:order
                if iseven(o)
                    coefficients[k, (order + 1) * j - (order - o)] = xCellFaceSigns[j, cell]
                end
            end
        end
        return nothing
    end
end
