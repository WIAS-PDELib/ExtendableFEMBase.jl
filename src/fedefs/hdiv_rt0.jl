"""
````
abstract type HDIVRT0{edim} <: AbstractHdivFiniteElement where {edim<:Int}
````

Hdiv-conforming vector-valued (ncomponents = edim) lowest-order Raviart-Thomas space.

allowed ElementGeometries:
- Triangle2D
- Quadrilateral2D
- Tetrahedron3D
- Hexahedron3D
"""
abstract type HDIVRT0{edim} <: AbstractHdivFiniteElement where {edim <: Int} end
HDIVRT0(edim::Int) = HDIVRT0{edim}

function Base.show(io::Core.IO, FEType::Type{<:HDIVRT0{edim}}) where {edim}
    return print(io, "HDIVRT0{$edim}")
end

get_ncomponents(FEType::Type{<:HDIVRT0}) = FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:HDIVRT0}, EG::Type{<:AbstractElementGeometry}) = 1
get_ndofs(::Type{ON_CELLS}, FEType::Type{<:HDIVRT0}, EG::Type{<:AbstractElementGeometry}) = num_faces(EG)

get_polynomialorder(::Type{<:HDIVRT0{2}}, ::Type{<:AbstractElementGeometry1D}) = 0;
get_polynomialorder(::Type{<:HDIVRT0{2}}, ::Type{<:AbstractElementGeometry2D}) = 1;
get_polynomialorder(::Type{<:HDIVRT0{3}}, ::Type{<:AbstractElementGeometry2D}) = 0;
get_polynomialorder(::Type{<:HDIVRT0{3}}, ::Type{<:AbstractElementGeometry3D}) = 1;

get_dofmap_pattern(FEType::Type{<:HDIVRT0}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry}) = "f1"
get_dofmap_pattern(FEType::Type{<:HDIVRT0}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry}) = "i1"

isdefined(FEType::Type{<:HDIVRT0}, ::Type{<:Triangle2D}) = true
isdefined(FEType::Type{<:HDIVRT0}, ::Type{<:Quadrilateral2D}) = true
isdefined(FEType::Type{<:HDIVRT0}, ::Type{<:Tetrahedron3D}) = true
isdefined(FEType::Type{<:HDIVRT0}, ::Type{<:Hexahedron3D}) = true

function RT0_normalflux_eval!(result, f, qpinfo)
    return result[1] = dot(f, qpinfo.normal)
end
init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}) where {Tv, Ti, FEType <: HDIVRT0, APT} = FaceFluxEvaluator(RT0_normalflux_eval!, FES, ON_FACES)

function ExtendableGrids.interpolate!(Target::AbstractArray{T, 1}, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}, exact_function!; items = [], kwargs...) where {T, Tv, Ti, FEType <: HDIVRT0, APT}
    return get_interpolator(FE, ON_FACES).evaluate!(Target, exact_function!, items; kwargs...)
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: HDIVRT0, APT}
    # delegate cell faces to face interpolation
    subitems = slice(FE.dofgrid[CellFaces], items)
    return interpolate!(Target, FE, ON_FACES, exact_function!; items = subitems, kwargs...)
end

# only normalfluxes on faces
function get_basis(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, ::Type{<:HDIVRT0}, ::Type{<:AbstractElementGeometry})
    return function closure(refbasis, xref)
        return refbasis[1, 1] = 1
    end
end

function get_basis(::Type{ON_CELLS}, ::Type{HDIVRT0{2}}, ::Type{<:Triangle2D})
    return function closure(refbasis, xref)
        refbasis[1, 1] = xref[1]
        refbasis[1, 2] = xref[2] - 1
        refbasis[2, 1] = xref[1]
        refbasis[2, 2] = xref[2]
        refbasis[3, 1] = xref[1] - 1
        return refbasis[3, 2] = xref[2]
    end
end

function get_basis(::Type{ON_CELLS}, ::Type{HDIVRT0{2}}, ::Type{<:Quadrilateral2D})
    return function closure(refbasis, xref)
        refbasis[1, 1] = 0
        refbasis[1, 2] = xref[2] - 1
        refbasis[2, 1] = xref[1]
        refbasis[2, 2] = 0
        refbasis[3, 1] = 0
        refbasis[3, 2] = xref[2]
        refbasis[4, 1] = xref[1] - 1
        return refbasis[4, 2] = 0
    end
end

function get_basis(::Type{ON_CELLS}, ::Type{HDIVRT0{3}}, ::Type{<:Tetrahedron3D})
    return function closure(refbasis, xref)
        refbasis[1, 1] = 2 * xref[1]
        refbasis[1, 2] = 2 * xref[2]
        refbasis[1, 3] = 2 * (xref[3] - 1)
        refbasis[2, 1] = 2 * xref[1]
        refbasis[2, 2] = 2 * (xref[2] - 1)
        refbasis[2, 3] = 2 * xref[3]
        refbasis[3, 1] = 2 * xref[1]
        refbasis[3, 2] = 2 * xref[2]
        refbasis[3, 3] = 2 * xref[3]
        refbasis[4, 1] = 2 * (xref[1] - 1)
        refbasis[4, 2] = 2 * xref[2]
        return refbasis[4, 3] = 2 * xref[3]
    end
    # note: factor 2 is chosen, such that normal-flux integrated over faces is 1 again
end

function get_basis(::Type{ON_CELLS}, ::Type{HDIVRT0{3}}, ::Type{<:Hexahedron3D})
    return function closure(refbasis, xref)
        fill!(refbasis, 0)
        refbasis[1, 3] = xref[3] - 1
        refbasis[2, 2] = xref[2] - 1
        refbasis[3, 1] = xref[1]
        refbasis[4, 2] = xref[2]
        refbasis[5, 1] = xref[1] - 1
        return refbasis[6, 3] = xref[3]
    end
end

function get_coefficients(::Type{ON_CELLS}, FE::FESpace{Tv, Ti, <:HDIVRT0, APT}, EG::Type{<:AbstractElementGeometry}, xgrid) where {Tv, Ti, APT}
    xCellFaceSigns = xgrid[CellFaceSigns]
    nfaces = num_faces(EG)
    return function closure(coefficients, cell)
        # multiplication with normal vector signs
        for j in 1:nfaces, k in 1:size(coefficients, 1)
            coefficients[k, j] = xCellFaceSigns[j, cell]
        end
        return nothing
    end
end
