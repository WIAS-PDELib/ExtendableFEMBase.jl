"""
````
abstract type HCURLN1{edim} <: AbstractHcurlFiniteElement where {edim<:Int}
````

Hcurl-conforming vector-valued (ncomponents = edim) Nedelec space of first kind and order 1.

allowed ElementGeometries:
- Triangle2D
"""
abstract type HCURLN1{edim} <: AbstractHcurlFiniteElement where {edim <: Int} end
HCURLN1(edim::Int) = HCURLN1{edim}

function Base.show(io::Core.IO, ::Type{<:HCURLN1{edim}}) where {edim}
    return print(io, "HCURLN1{$edim}")
end

get_ncomponents(FEType::Type{<:HCURLN1}) = FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_EDGES}, Type{<:ON_BEDGES}, Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:HCURLN1}, EG::Type{<:Edge1D}) = 2
get_ndofs(::Type{ON_CELLS}, FEType::Type{HCURLN1{2}}, EG::Type{<:Triangle2D}) = 2 * num_faces(EG) + 2

get_polynomialorder(::Type{<:HCURLN1{2}}, ::Type{<:AbstractElementGeometry1D}) = 1;
get_polynomialorder(::Type{<:HCURLN1{2}}, ::Type{<:AbstractElementGeometry2D}) = 2;

get_dofmap_pattern(FEType::Type{<:HCURLN1{2}}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry2D}) = "f2i2"
get_dofmap_pattern(FEType::Type{<:HCURLN1{2}}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry1D}) = "i2"

isdefined(FEType::Type{<:HCURLN1}, ::Type{<:Triangle2D}) = true

function N1_tangentflux_eval_2d!(result, f, qpinfo)
    result[1] = -f[1] * qpinfo.normal[2] # rotated normal = tangent
    result[1] += f[2] * qpinfo.normal[1]
    result[2] = result[1] * (qpinfo.xref[1] - 1 // 2)
    return nothing
end
init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}) where {Tv, Ti, FEType <: HCURLN1{2}, APT} = FunctionalInterpolator(N1_tangentflux_eval_2d!, FES, ON_FACES; bonus_quadorder = 1)
init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}) where {Tv, Ti, FEType <: HCURLN1{2}, APT} = MomentInterpolator(FES, ON_CELLS)


function ExtendableGrids.interpolate!(Target::AbstractArray{T, 1}, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_EDGES}, exact_function!; items = [], kwargs...) where {T, Tv, Ti, FEType <: HCURLN1, APT}
    edim = get_ncomponents(FEType)
    return if edim == 3
        # todo
    end
end

function ExtendableGrids.interpolate!(Target::AbstractArray{T, 1}, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}, exact_function!; items = [], kwargs...) where {T, Tv, Ti, FEType <: HCURLN1, APT}
    edim = get_ncomponents(FEType)
    return if edim == 2
        get_interpolator(FE, ON_FACES).evaluate!(Target, exact_function!, items; kwargs...)
    elseif edim == 3
        # delegate face edges to edge interpolation
        subitems = slice(FE.dofgrid[FaceEdges], items)
        interpolate!(Target, FE, ON_EDGES, data; items = subitems, kwargs...)
    end
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}, data; items = [], kwargs...) where {Tv, Ti, FEType <: HCURLN1, APT}
    edim = get_ncomponents(FEType)
    return if edim == 2
        # delegate cell faces to face interpolation
        subitems = slice(FE.dofgrid[CellFaces], items)
        interpolate!(Target, FE, ON_FACES, data; items = subitems, kwargs...)
    elseif edim == 3
        # delegate cell edges to edge interpolation
        subitems = slice(FE.dofgrid[CellEdges], items)
        interpolate!(Target, FE, ON_EDGES, data; items = subitems, kwargs...)
    end

    # set values of interior N1 functions such that P0 moments are preserved
    get_interpolator(FE, ON_CELLS).evaluate!(Target, data, items; kwargs...)

end

# on faces dofs are only tangential fluxes
function get_basis(::Union{Type{<:ON_EDGES}, Type{<:ON_BEDGES}, Type{<:ON_BFACES}, Type{<:ON_FACES}}, ::Type{<:HCURLN1}, ::Type{<:AbstractElementGeometry})
    return function closure(refbasis, xref)
        refbasis[1, 1] = 1                # tangent-flux of N0 function on single face
        return refbasis[2, 1] = 12 * (xref[1] - 1 // 2) # linear tangent-flux of additional N1 edge function
    end
end

function get_basis(::Type{ON_CELLS}, ::Type{HCURLN1{2}}, ::Type{<:Triangle2D})
    return function closure(refbasis, xref)
        ## HCURLN0 basis
        refbasis[1, 1] = 1 - xref[2]
        refbasis[1, 2] = xref[1]
        refbasis[3, 1] = -xref[2]
        refbasis[3, 2] = xref[1]
        refbasis[5, 1] = -xref[2]
        refbasis[5, 2] = xref[1] - 1

        ## additional functions
        for k in 1:2
            # additional N1 edge basis functions
            refbasis[2, k] = -12 * (1 // 2 - xref[1] - xref[2]) * refbasis[1, k]
            refbasis[4, k] = -(12 * (xref[1] - 1 // 2)) * refbasis[3, k]
            refbasis[6, k] = -(12 * (xref[2] - 1 // 2)) * refbasis[5, k]
            # interior functions
            refbasis[7, k] = 12 * xref[2] * refbasis[1, k]
            refbasis[8, k] = 12 * xref[1] * refbasis[5, k]
        end

        return nothing
    end
end

function get_coefficients(::Type{ON_CELLS}, FE::FESpace{Tv, Ti, <:HCURLN1, APT}, EG::Type{<:AbstractElementGeometry2D}, xgrid) where {Tv, Ti, APT}
    xCellFaceSigns = xgrid[CellFaceSigns]
    nfaces = num_faces(EG)
    return function closure(coefficients, cell)
        # multiplication with normal vector signs (only RT0)
        for j in 1:nfaces, k in 1:size(coefficients, 1)
            coefficients[k, 2 * j - 1] = xCellFaceSigns[j, cell]
        end
        return nothing
    end
end
