"""
````
abstract type H1CR{ncomponents} <: AbstractH1FiniteElement where {ncomponents<:Int}
````

Crouzeix-Raviart element (only continuous at face centers).

allowed ElementGeometries:
- Triangle2D (piecewise linear, similar to P1)
- Quadrilateral2D (similar to Q1 space)
- Tetrahedron3D (piecewise linear, similar to P1)
"""
abstract type H1CR{ncomponents} <: AbstractH1FiniteElement where {ncomponents <: Int} end
H1CR(ncomponents::Int) = H1CR{ncomponents}

function Base.show(io::Core.IO, ::Type{<:H1CR{ncomponents}}) where {ncomponents}
    return print(io, "H1CR{$ncomponents}")
end

get_ncomponents(FEType::Type{<:H1CR}) = FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:H1CR}, EG::Type{<:AbstractElementGeometry}) = FEType.parameters[1]
get_ndofs(::Type{ON_CELLS}, FEType::Type{<:H1CR}, EG::Type{<:AbstractElementGeometry}) = num_faces(EG) * FEType.parameters[1]

get_polynomialorder(::Type{<:H1CR}, ::Type{<:Edge1D}) = 1; # 0 on continuous edges, but = 1 on edges with jumps
get_polynomialorder(::Type{<:H1CR}, ::Type{<:Triangle2D}) = 1;
get_polynomialorder(::Type{<:H1CR}, ::Type{<:Quadrilateral2D}) = 2;
get_polynomialorder(::Type{<:H1CR}, ::Type{<:Tetrahedron3D}) = 1;

get_dofmap_pattern(FEType::Type{<:H1CR}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry}) = "F1"
get_dofmap_pattern(FEType::Type{<:H1CR}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry}) = "I1"

interior_dofs_offset(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, ::Type{<:H1CR}, ::Type{Edge1D}) = 0
interior_dofs_offset(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, ::Type{<:H1CR}, ::Type{Triangle2D}) = 0

isdefined(FEType::Type{<:H1CR}, ::Type{<:AbstractElementGeometry1D}) = true
isdefined(FEType::Type{<:H1CR}, ::Type{<:Triangle2D}) = true
isdefined(FEType::Type{<:H1CR}, ::Type{<:Quadrilateral2D}) = true
isdefined(FEType::Type{<:H1CR}, ::Type{<:Tetrahedron3D}) = true

init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}) where {Tv, Ti, FEType <: H1CR, APT} = MomentInterpolator(FES, ON_FACES; FEType_ref = L2P0{get_ncomponents(FEType)})

function ExtendableGrids.interpolate!(Target::AbstractArray{T, 1}, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}, exact_function!; items = [], kwargs...) where {T, Tv, Ti, FEType <: H1CR, APT}
    return get_interpolator(FE, ON_FACES).evaluate!(Target, exact_function!, items; kwargs...)
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1CR, APT}
    # delegate cell faces to face interpolation
    subitems = slice(FE.dofgrid[CellFaces], items)
    return interpolate!(Target, FE, ON_FACES, exact_function!; items = subitems, kwargs...)
end

# BEWARE ON FACES
#
# all basis functions on a cell are nonzero on all edges,
# but only the dof associated to the face is evaluated
# when using get_basis on faces
function get_basis(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, ::Type{H1CR{ncomponents}}, ::Type{<:AbstractElementGeometry}) where {ncomponents}
    return function closure(refbasis, xref)
        for k in 1:ncomponents
            refbasis[k, k] = 1
        end
        return
    end
end

function get_basis(::Type{ON_CELLS}, ::Type{<:H1CR{ncomponents}}, ET::Type{<:Triangle2D}) where {ncomponents}
    return function closure(refbasis, xref)
        refbasis[end] = 2 * (xref[1] + xref[2]) - 1
        for k in 1:ncomponents
            refbasis[3 * k - 2, k] = 1 - 2 * xref[2]
            refbasis[3 * k - 1, k] = refbasis[end]
            refbasis[3 * k, k] = 1 - 2 * xref[1]
        end
        return
    end
end

function get_basis(::Type{ON_CELLS}, ::Type{H1CR{ncomponents}}, ET::Type{<:Tetrahedron3D}) where {ncomponents}
    return function closure(refbasis, xref)
        refbasis[end] = 3 * (xref[1] + xref[2] + xref[3]) - 2
        for k in 1:ncomponents
            refbasis[4 * k - 3, k] = 1 - 3 * xref[3]
            refbasis[4 * k - 2, k] = 1 - 3 * xref[2]
            refbasis[4 * k - 1, k] = refbasis[end]
            refbasis[4 * k, k] = 1 - 3 * xref[1]
        end
        return
    end
end

function get_basis(::Type{ON_CELLS}, ::Type{H1CR{ncomponents}}, ET::Type{<:Quadrilateral2D}) where {ncomponents}
    return function closure(refbasis, xref)
        refbasis[1, 1] = xref[1] * (1 - xref[1]) + (1 - xref[2])^2 - 1 // 4
        refbasis[2, 1] = xref[2] * (1 - xref[2]) + xref[1] * xref[1] - 1 // 4
        refbasis[3, 1] = xref[1] * (1 - xref[1]) + xref[2] * xref[2] - 1 // 4
        refbasis[4, 1] = xref[2] * (1 - xref[2]) + (1 - xref[1])^2 - 1 // 4
        for k in 2:ncomponents, j in 1:4
            refbasis[4 * (k - 1) + j, k] = refbasis[j, 1]
        end
        return
    end
end
