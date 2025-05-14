"""
````
abstract type HDIVBDM2{edim} <: AbstractHdivFiniteElement where {edim<:Int}
````

Hdiv-conforming vector-valued (ncomponents = edim) Brezzi-Douglas-Marini space of order 2

allowed ElementGeometries:
- Triangle2D
"""
abstract type HDIVBDM2{edim} <: AbstractHdivFiniteElement where {edim <: Int} end
HDIVBDM2(edim::Int) = HDIVBDM2{edim}

function Base.show(io::Core.IO, ::Type{<:HDIVBDM2{edim}}) where {edim}
    return print(io, "HDIVBDM2{$edim}")
end

get_ncomponents(FEType::Type{<:HDIVBDM2}) = FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:HDIVBDM2}, EG::Type{<:AbstractElementGeometry1D}) = 3
get_ndofs(::Type{ON_CELLS}, FEType::Type{<:HDIVBDM2}, EG::Type{<:AbstractElementGeometry2D}) = 3 * num_faces(EG) + 3

get_polynomialorder(::Type{<:HDIVBDM2{2}}, ::Type{<:Edge1D}) = 2;
get_polynomialorder(::Type{<:HDIVBDM2{2}}, ::Type{<:Triangle2D}) = 2;

get_dofmap_pattern(FEType::Type{<:HDIVBDM2{2}}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry}) = "f3i3"
get_dofmap_pattern(FEType::Type{<:HDIVBDM2{2}}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry}) = "i3"

isdefined(FEType::Type{<:HDIVBDM2}, ::Type{<:Triangle2D}) = true

interior_dofs_offset(::Type{<:ON_CELLS}, ::Type{<:HDIVBDM2{2}}, ::Type{<:Triangle2D}) = 9


function BDM2_normalflux_eval!(dim)
    @assert dim == 2 "BDM3 for dim=3 not available yet"
    function closure(result, f, qpinfo)
        result[1] = dot(f, qpinfo.normal)
        result[2] = result[1] * (qpinfo.xref[1] - 1 // dim)
        result[3] = result[1] * (qpinfo.xref[1]^2 - qpinfo.xref[1] + 1 // 6)
    end
end
init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}) where {Tv, Ti, FEType <: HDIVBDM2, APT} = FunctionalInterpolator(BDM2_normalflux_eval!(FEType.parameters[1]), FES, ON_FACES; bonus_quadorder = 2)
init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}) where {Tv, Ti, FEType <: HDIVBDM2{2}, APT} = MomentInterpolator(FES, ON_CELLS; FEType_moments = H1MINI{1,2}, moments_dofs = [1,2,4], moments_operator = CurlScalar)


function ExtendableGrids.interpolate!(Target::AbstractArray{T, 1}, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}, exact_function!; items = [], kwargs...) where {T, Tv, Ti, FEType <: HDIVBDM2, APT}
    get_interpolator(FE, ON_FACES).evaluate!(Target, exact_function!, items; kwargs...)
end

function ExtendableGrids.interpolate!(Target::AbstractArray{T, 1}, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}, exact_function!; items = [], kwargs...) where {T, Tv, Ti, FEType <: HDIVBDM2, APT}
    # delegate cell faces to face interpolation
    subitems = slice(FE.dofgrid[CellFaces], items)
    interpolate!(Target, FE, ON_FACES, exact_function!; items = subitems, kwargs...)

    # set interior dofs such that moments with ∇P1 and curl(b) are preserved
    # (b is the cell bubble)
    get_interpolator(FE, ON_CELLS).evaluate!(Target, exact_function!, items; kwargs...)
end

## only normalfluxes on faces
function get_basis(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, ::Type{<:HDIVBDM2}, ::Type{<:AbstractElementGeometry1D})
    return function closure(refbasis, xref)
        refbasis[1, 1] = 1
        refbasis[2, 1] = 12 * (xref[1] - 1 // 2) # linear normal-flux of BDM2 function
        return refbasis[3, 1] = 180 * (xref[1]^2 - xref[1] + 1 // 6) # quadratic normal-flux of BDM2 function
    end
end

function get_basis(::Type{ON_CELLS}, ::Type{HDIVBDM2{2}}, ::Type{<:Triangle2D})
    return function closure(refbasis, xref)
        refbasis[end] = 1 - xref[1] - xref[2]

        # RT0 basis
        refbasis[1, 1] = xref[1]
        refbasis[1, 2] = xref[2] - 1
        refbasis[4, 1] = xref[1]
        refbasis[4, 2] = xref[2]
        refbasis[7, 1] = xref[1] - 1
        refbasis[7, 2] = xref[2]
        # additional BDM1 functions on faces
        refbasis[2, 1] = 6 * xref[1]
        refbasis[2, 2] = 6 - 12 * xref[1] - 6 * xref[2]
        refbasis[5, 1] = -6 * xref[1]
        refbasis[5, 2] = 6 * xref[2]
        refbasis[8, 1] = 6 * (xref[1] - 1) + 12 * xref[2]
        refbasis[8, 2] = -6 * xref[2]
        for k in 1:2
            # additional BDM2 face functions on faces
            refbasis[3, k] = -15 * ((refbasis[end] - 1 // 2) * refbasis[2, k] + refbasis[1, k])
            refbasis[6, k] = -15 * ((xref[1] - 1 // 2) * refbasis[5, k] + refbasis[4, k])
            refbasis[9, k] = -15 * ((xref[2] - 1 // 2) * refbasis[8, k] + refbasis[7, k])
            # additional BDM2 interior functions
            refbasis[10, k] = xref[2] * refbasis[2, k]
            refbasis[11, k] = refbasis[end] * refbasis[5, k]
            refbasis[12, k] = xref[1] * refbasis[8, k]
        end
        return
    end
end


function get_coefficients(::Type{ON_CELLS}, FE::FESpace{Tv, Ti, <:HDIVBDM2, APT}, EG::Type{<:AbstractElementGeometry2D}, xgrid) where {Tv, Ti, APT}
    xCellFaceSigns::Union{VariableTargetAdjacency{Int32}, Array{Int32, 2}} = xgrid[CellFaceSigns]
    nfaces::Int = num_faces(EG)
    dim::Int = dim_element(EG)
    return function closure(coefficients::Array{<:Real, 2}, cell::Int)
        fill!(coefficients, 1.0)
        # multiplication with normal vector signs (only RT0)
        for j in 1:nfaces, k in 1:dim
            coefficients[k, 3 * j - 2] = xCellFaceSigns[j, cell]
            coefficients[k, 3 * j] = xCellFaceSigns[j, cell]
        end
        return nothing
    end
end
