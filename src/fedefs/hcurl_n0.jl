"""
````
abstract type HCURLN0{edim} <: AbstractHcurlFiniteElement where {edim<:Int}
````

Hcurl-conforming vector-valued (ncomponents = edim) lowest-order Nedelec space of first kind.

allowed ElementGeometries:
- Triangle2D
- Quadrilateral2D
- Tetrahedron3D
"""
abstract type HCURLN0{edim} <: AbstractHcurlFiniteElement where {edim <: Int} end
HCURLN0(edim::Int) = HCURLN0{edim}

function Base.show(io::Core.IO, ::Type{<:HCURLN0{edim}}) where {edim}
    return print(io, "HCURLN0{$edim}")
end

get_ncomponents(FEType::Type{<:HCURLN0}) = FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_EDGES}, Type{<:ON_BEDGES}, Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:HCURLN0}, EG::Type{<:AbstractElementGeometry}) = 1
get_ndofs(::Type{ON_CELLS}, FEType::Type{HCURLN0{2}}, EG::Type{<:AbstractElementGeometry}) = num_faces(EG)
get_ndofs(::Type{ON_CELLS}, FEType::Type{HCURLN0{3}}, EG::Type{<:AbstractElementGeometry}) = num_edges(EG)

get_polynomialorder(::Type{<:HCURLN0{2}}, ::Type{<:AbstractElementGeometry1D}) = 0;
get_polynomialorder(::Type{<:HCURLN0{2}}, ::Type{<:AbstractElementGeometry2D}) = 1;
get_polynomialorder(::Type{<:HCURLN0{3}}, ::Type{<:AbstractElementGeometry1D}) = 0;
get_polynomialorder(::Type{<:HCURLN0{3}}, ::Type{<:AbstractElementGeometry3D}) = 1;

get_dofmap_pattern(FEType::Type{<:HCURLN0{2}}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry2D}) = "f1"
get_dofmap_pattern(FEType::Type{<:HCURLN0{2}}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry1D}) = "i1"

get_dofmap_pattern(FEType::Type{<:HCURLN0{3}}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry3D}) = "e1"
get_dofmap_pattern(FEType::Type{<:HCURLN0{3}}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry2D}) = "e1"
get_dofmap_pattern(FEType::Type{<:HCURLN0{3}}, ::Union{Type{EdgeDofs}, Type{BEdgeDofs}}, EG::Type{<:AbstractElementGeometry1D}) = "i1"

isdefined(FEType::Type{<:HCURLN0}, ::Type{<:Triangle2D}) = true
isdefined(FEType::Type{<:HCURLN0}, ::Type{<:Quadrilateral2D}) = true
isdefined(FEType::Type{<:HCURLN0}, ::Type{<:Tetrahedron3D}) = true

function ExtendableGrids.interpolate!(Target::AbstractArray{T, 1}, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_EDGES}, data; items = [], kwargs...) where {T, Tv, Ti, FEType <: HCURLN0, APT}
    edim = get_ncomponents(FEType)
    return if edim == 3
        xEdgeTangents = FE.dofgrid[EdgeTangents]
        if items == []
            items = 1:size(xEdgeTangents, 2)
        end
        data_eval = zeros(T, 3)
        function tangentflux_eval3d(result, qpinfo)
            data(data_eval, qpinfo)
            return result[1] = dot(data_eval, view(xEdgeTangents, :, qpinfo.item))
        end
        integrate!(Target, FE.dofgrid, ON_EDGES, tangentflux_eval3d; items = items, kwargs...)
    end
end

function ExtendableGrids.interpolate!(Target::AbstractArray{T, 1}, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}, data; items = [], kwargs...) where {T, Tv, Ti, FEType <: HCURLN0, APT}
    edim = get_ncomponents(FEType)
    return if edim == 2
        xFaceNormals = FE.dofgrid[FaceNormals]
        if items == []
            items = 1:size(xFaceNormals, 2)
        end
        data_eval = zeros(T, 2)
        function tangentflux_eval2d(result, qpinfo)
            data(data_eval, qpinfo)
            result[1] = -data_eval[1] * xFaceNormals[2, qpinfo.item] # rotated normal = tangent
            return result[1] += data_eval[2] * xFaceNormals[1, qpinfo.item]
        end
        integrate!(Target, FE.dofgrid, ON_FACES, tangentflux_eval2d; items = items, kwargs...)
    elseif edim == 3
        # delegate face edges to edge interpolation
        subitems = slice(FE.dofgrid[FaceEdges], items)
        interpolate!(Target, FE, ON_EDGES, data; items = subitems, kwargs...)
    end
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}, data; items = [], kwargs...) where {Tv, Ti, FEType <: HCURLN0, APT}
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
end

# on faces dofs are only tangential fluxes
function get_basis(::Union{Type{<:ON_EDGES}, Type{<:ON_BEDGES}, Type{<:ON_BFACES}, Type{<:ON_FACES}}, ::Type{<:HCURLN0}, ::Type{<:AbstractElementGeometry})
    return function closure(refbasis, xref)
        return refbasis[1, 1] = 1
    end
end

function get_basis(::Type{ON_CELLS}, ::Type{HCURLN0{2}}, ::Type{<:Triangle2D})
    return function closure(refbasis, xref)
        refbasis[1, 1] = 1 - xref[2]
        refbasis[1, 2] = xref[1]
        refbasis[2, 1] = -xref[2]
        refbasis[2, 2] = xref[1]
        refbasis[3, 1] = -xref[2]
        refbasis[3, 2] = xref[1] - 1
        return nothing
    end
end

function get_basis(::Type{ON_CELLS}, ::Type{HCURLN0{2}}, ::Type{<:Quadrilateral2D})
    return function closure(refbasis, xref)
        refbasis[1, 1] = 1 - xref[2]
        refbasis[1, 2] = 0
        refbasis[2, 1] = 0
        refbasis[2, 2] = xref[1]
        refbasis[3, 1] = -xref[2]
        refbasis[3, 2] = 0
        refbasis[4, 1] = 0
        refbasis[4, 2] = xref[1] - 1
        return nothing
    end
end

function get_basis(::Type{ON_CELLS}, ::Type{HCURLN0{3}}, ::Type{<:Tetrahedron3D})
    return function closure(refbasis, xref) # [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]
        refbasis[1, 1] = 1.0 - xref[2] - xref[3]
        refbasis[1, 2] = xref[1]
        refbasis[1, 3] = xref[1]             # edge 1 = [1,2] (y=z=0)        t = (1,0,0)
        refbasis[2, 1] = xref[2]
        refbasis[2, 2] = 1 - xref[3] - xref[1]
        refbasis[2, 3] = xref[2]             # edge 2 = [1,3] (x=z=0)        t = (0,1,0)
        refbasis[3, 1] = xref[3]
        refbasis[3, 2] = xref[3]
        refbasis[3, 3] = 1 - xref[1] - xref[2]   # edge 3 = [1,4] (x=y=0)        t = (0,0,1)
        refbasis[4, 1] = -xref[2]
        refbasis[4, 2] = xref[1]
        refbasis[4, 3] = 0                   # edge 4 = [2,3] (z=0, y=1-x)   t = (-1,1,0)
        refbasis[5, 1] = -xref[3]
        refbasis[5, 2] = 0
        refbasis[5, 3] = xref[1]             # edge 5 = [2,4] (y=0, z=1-x)   t = (-1,0,1)
        refbasis[6, 1] = 0
        refbasis[6, 2] = -xref[3]
        return refbasis[6, 3] = xref[2]             # edge 6 = [3,4] (x=0, z=1-y)   t = (0,-1,1)
    end
end

function get_coefficients(::Type{ON_CELLS}, FE::FESpace{Tv, Ti, <:HCURLN0, APT}, EG::Type{<:AbstractElementGeometry2D}, xgrid) where {Tv, Ti, APT}
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

function get_coefficients(::Type{ON_CELLS}, FE::FESpace{Tv, Ti, <:HCURLN0, APT}, EG::Type{<:AbstractElementGeometry3D}, xgrid) where {Tv, Ti, APT}
    xCellEdgeSigns = xgrid[CellEdgeSigns]
    nedges = num_edges(EG)
    return function closure(coefficients, cell)
        # multiplication with normal vector signs
        for j in 1:nedges, k in 1:size(coefficients, 1)
            coefficients[k, j] = xCellEdgeSigns[j, cell]
        end
        return nothing
    end
end
