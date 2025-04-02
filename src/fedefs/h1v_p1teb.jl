"""
````
abstract type H1P1TEB{edim} <: AbstractH1FiniteElementWithCoefficients where {edim<:Int}
````

vector-valued (ncomponents = edim) element that uses P1 functions + tangential-weighted edge bubbles
as suggested by [Diening, L., Storn, J. & Tscherpel, T., "Fortin operator for the Taylor–Hood element", Num. Math. 150, 671–689 (2022)]

(is inf-sup stable for Stokes if paired with continuous P1 pressure space, less degrees of freedom than MINI)

allowed ElementGeometries:
- Triangle2D
- Tetrahedron3D
"""
abstract type H1P1TEB{edim} <: AbstractH1FiniteElementWithCoefficients where {edim <: Int} end
H1P1TEB(edim::Int) = H1P1TEB{edim}

function Base.show(io::Core.IO, ::Type{<:H1P1TEB{edim}}) where {edim}
    return print(io, "H1P1TEB{$edim}")
end

get_ncomponents(FEType::Type{<:H1P1TEB}) = FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:H1P1TEB{2}}, EG::Type{<:AbstractElementGeometry}) = 1 + num_nodes(EG) * FEType.parameters[1]
get_ndofs(::Type{ON_CELLS}, FEType::Type{<:H1P1TEB{2}}, EG::Type{<:AbstractElementGeometry}) = num_faces(EG) + num_nodes(EG) * FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_EDGES}, Type{<:ON_BEDGES}}, FEType::Type{<:H1P1TEB{3}}, EG::Type{<:AbstractElementGeometry}) = 1 + num_nodes(EG) * FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:H1P1TEB{3}}, EG::Type{<:AbstractElementGeometry}) = num_edges(EG) + num_nodes(EG) * FEType.parameters[1]
get_ndofs(::Type{ON_CELLS}, FEType::Type{<:H1P1TEB{3}}, EG::Type{<:AbstractElementGeometry}) = num_edges(EG) + num_nodes(EG) * FEType.parameters[1]

get_polynomialorder(::Type{<:H1P1TEB{2}}, ::Type{<:Edge1D}) = 2
get_polynomialorder(::Type{<:H1P1TEB{2}}, ::Type{<:Triangle2D}) = 2
get_polynomialorder(::Type{<:H1P1TEB{3}}, ::Type{<:Edge1D}) = 2
get_polynomialorder(::Type{<:H1P1TEB{3}}, ::Type{<:Triangle2D}) = 2
get_polynomialorder(::Type{<:H1P1TEB{3}}, ::Type{<:Tetrahedron3D}) = 2

get_dofmap_pattern(FEType::Type{<:H1P1TEB{2}}, ::Type{CellDofs}, EG::Type{<:Triangle2D}) = "N1f1"
get_dofmap_pattern(FEType::Type{<:H1P1TEB{2}}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:Edge1D}) = "N1i1"
get_dofmap_pattern(FEType::Type{<:H1P1TEB{3}}, ::Type{EdgeDofs}, EG::Type{<:Edge1D}) = "N1i1"
get_dofmap_pattern(FEType::Type{<:H1P1TEB{3}}, ::Type{CellDofs}, EG::Type{<:Tetrahedron3D}) = "N1e1"
get_dofmap_pattern(FEType::Type{<:H1P1TEB{3}}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:Triangle2D}) = "N1e1"

isdefined(FEType::Type{<:H1P1TEB}, ::Type{<:Triangle2D}) = true
isdefined(FEType::Type{<:H1P1TEB}, ::Type{<:Tetrahedron3D}) = true


function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{AT_NODES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1P1TEB, APT}
    nnodes = size(FE.dofgrid[Coordinates], 2)
    return point_evaluation!(Target, FE, AT_NODES, exact_function!; items = items, component_offset = nnodes, kwargs...)
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_EDGES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1P1TEB{2}, APT}
    # delegate edge nodes to node interpolation
    subitems = slice(FE.dofgrid[EdgeNodes], items)
    return interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1P1TEB{3}, APT}
    # delegate edges to edge interpolation
    subitems = slice(FE.dofgrid[FaceEdges], items)
    return interpolate!(Target, FE, ON_EDGES, exact_function!; items = subitems, kwargs...)
end

function ExtendableGrids.interpolate!(Target::AbstractArray{T, 1}, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}, exact_function!; items = [], bonus_quadorder = 0, kwargs...) where {T, Tv, Ti, FEType <: H1P1TEB{2}, APT}
    # delegate face nodes to node interpolation
    subitems = slice(FE.dofgrid[FaceNodes], items)
    interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)

    # preserve face means in tangential direction
    xItemVolumes = FE.dofgrid[FaceVolumes]
    xItemNodes = FE.dofgrid[FaceNodes]
    xItemGeometries = FE.dofgrid[FaceGeometries]
    xFaceNormals = FE.dofgrid[FaceNormals]
    xItemDofs = FE[FaceDofs]
    ncomponents = get_ncomponents(FEType)
    nnodes = size(FE.dofgrid[Coordinates], 2)
    nitems = num_sources(xItemNodes)
    offset = ncomponents * nnodes
    if items == []
        items = 1:nitems
    end

    # compute exact face means
    facemeans = zeros(T, ncomponents, nitems)
    integrate!(facemeans, FE.dofgrid, ON_FACES, exact_function!; quadorder = 2 + bonus_quadorder, items = items, kwargs...)
    P1flux::T = 0
    value::T = 0
    itemEG = Edge1D
    nitemnodes::Int = 0
    for item in items
        itemEG = xItemGeometries[item]
        nitemnodes = num_nodes(itemEG)
        # compute normal flux (minus linear part)
        value = 0
        for c in 1:ncomponents
            P1flux = 0
            for dof in 1:nitemnodes
                P1flux += Target[xItemDofs[(c - 1) * nitemnodes + dof, item]] * xItemVolumes[item] / nitemnodes
            end
            if c == 1
                value -= (facemeans[c, item] - P1flux) * xFaceNormals[2, item]
            else
                value += (facemeans[c, item] - P1flux) * xFaceNormals[1, item]
            end
        end
        # set face bubble value
        Target[offset + item] = value / xItemVolumes[item]
    end
    return
end


function ExtendableGrids.interpolate!(Target::AbstractArray{T, 1}, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_EDGES}, exact_function!; items = [], bonus_quadorder = 0, kwargs...) where {T, Tv, Ti, FEType <: H1P1TEB{3}, APT}
    # delegate face nodes to node interpolation
    subitems = slice(FE.dofgrid[EdgeNodes], items)
    interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)

    # preserve edge means in tangential direction
    xItemVolumes = FE.dofgrid[EdgeVolumes]
    xItemNodes = FE.dofgrid[EdgeNodes]
    xItemGeometries = FE.dofgrid[EdgeGeometries]
    xEdgeTangents = FE.dofgrid[EdgeTangents]
    xItemDofs = FE[EdgeDofs]
    nnodes = size(FE.dofgrid[Coordinates], 2)
    nitems = num_sources(xItemNodes)
    ncomponents = get_ncomponents(FEType)
    offset = ncomponents * nnodes
    if items == []
        items = 1:nitems
    end

    # compute exact face means
    edgemeans = zeros(T, ncomponents, nitems)
    integrate!(edgemeans, FE.dofgrid, ON_EDGES, exact_function!; quadorder = 2 + bonus_quadorder, items = items, kwargs...)
    P1flux::T = 0
    value::T = 0
    itemEG = Edge1D
    nitemnodes::Int = 0
    for item in items
        itemEG = xItemGeometries[item]
        nitemnodes = num_nodes(itemEG)
        # compute normal flux (minus linear part)
        value = 0
        for c in 1:ncomponents
            P1flux = 0
            for dof in 1:nitemnodes
                P1flux += Target[xItemDofs[(c - 1) * nitemnodes + dof, item]] * xItemVolumes[item] / nitemnodes
            end
            value += (edgemeans[c, item] - P1flux) * xEdgeTangents[c, item]
        end
        # set face bubble value
        Target[offset + item] = value / xItemVolumes[item]
    end
    return
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1P1TEB, APT}
    # delegate cell faces to face interpolation
    subitems = slice(FE.dofgrid[CellFaces], items)
    return interpolate!(Target, FE, ON_FACES, exact_function!; items = subitems, kwargs...)
end


############
# 2D basis #
############

function get_basis(AT::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, ::Type{H1P1TEB{2}}, EG::Type{<:Edge1D})
    refbasis_P1 = get_basis(AT, H1P1{2}, EG)
    offset = get_ndofs(AT, H1P1{2}, EG)
    return function closure(refbasis, xref)
        refbasis_P1(refbasis, xref)
        # add face bubble to P1 basis
        refbasis[offset + 1, 1] = 6 * xref[1] * refbasis[1, 1]
        return refbasis[offset + 1, 2] = refbasis[offset + 1, 1]
    end
end

function get_basis(AT::Type{ON_CELLS}, ::Type{H1P1TEB{2}}, EG::Type{<:Triangle2D})
    refbasis_P1 = get_basis(AT, H1P1{2}, EG)
    offset = get_ndofs(AT, H1P1{2}, EG)
    return function closure(refbasis, xref)
        refbasis_P1(refbasis, xref)
        # add face bubbles to P1 basis
        refbasis[offset + 1, 1] = 6 * xref[1] * refbasis[1, 1]
        refbasis[offset + 2, 1] = 6 * xref[2] * xref[1]
        refbasis[offset + 3, 1] = 6 * refbasis[1, 1] * xref[2]
        refbasis[offset + 1, 2] = refbasis[offset + 1, 1]
        refbasis[offset + 2, 2] = refbasis[offset + 2, 1]
        return refbasis[offset + 3, 2] = refbasis[offset + 3, 1]
    end
end

function get_coefficients(::Type{ON_CELLS}, FE::FESpace{Tv, Ti, H1P1TEB{2}, APT}, ::Type{<:Triangle2D}) where {Tv, Ti, APT}
    xFaceNormals::Array{Tv, 2} = FE.dofgrid[FaceNormals]
    xCellFaces = FE.dofgrid[CellFaces]
    return function closure(coefficients::Array{<:Real, 2}, cell)
        fill!(coefficients, 1.0)
        coefficients[1, 7] = -xFaceNormals[2, xCellFaces[1, cell]]
        coefficients[2, 7] = xFaceNormals[1, xCellFaces[1, cell]]
        coefficients[1, 8] = -xFaceNormals[2, xCellFaces[2, cell]]
        coefficients[2, 8] = xFaceNormals[1, xCellFaces[2, cell]]
        coefficients[1, 9] = -xFaceNormals[2, xCellFaces[3, cell]]
        return coefficients[2, 9] = xFaceNormals[1, xCellFaces[3, cell]]
    end
end

function get_coefficients(::Type{<:ON_FACES}, FE::FESpace{Tv, Ti, H1P1TEB{2}, APT}, ::Type{<:Edge1D}) where {Tv, Ti, APT}
    xFaceNormals::Array{Tv, 2} = FE.dofgrid[FaceNormals]
    return function closure(coefficients::Array{<:Real, 2}, face)
        # multiplication of face bubble with normal vector of face
        fill!(coefficients, 1.0)
        coefficients[1, 5] = -xFaceNormals[2, face]
        return coefficients[2, 5] = xFaceNormals[1, face]
    end
end


############
# 3D basis #
############

function get_basis(AT::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, ::Type{H1P1TEB{3}}, EG::Type{<:Triangle2D})
    refbasis_P1 = get_basis(AT, H1P1{3}, EG)
    offset = get_ndofs(AT, H1P1{3}, EG)
    return function closure(refbasis, xref)
        refbasis_P1(refbasis, xref)
        # add edge bubbles to P1 basis
        refbasis[offset + 1, 1] = 6 * xref[1] * refbasis[1, 1]
        refbasis[offset + 2, 1] = 6 * xref[2] * xref[1]
        refbasis[offset + 3, 1] = 6 * refbasis[1, 1] * xref[2]
        for j in 1:3, k in 2:3
            refbasis[offset + j, k] = refbasis[offset + j, 1]
        end
        return
    end
end

function get_basis(AT::Type{ON_CELLS}, ::Type{H1P1TEB{3}}, EG::Type{<:Tetrahedron3D})
    refbasis_P1 = get_basis(AT, H1P1{3}, EG)
    offset = get_ndofs(AT, H1P1{3}, EG)
    return function closure(refbasis, xref)
        refbasis_P1(refbasis, xref)
        # add edge bubbles to P1 basis
        refbasis[offset + 1, 1] = 6 * xref[1] * refbasis[1, 1]
        refbasis[offset + 2, 1] = 6 * xref[2] * refbasis[1, 1]
        refbasis[offset + 3, 1] = 6 * xref[3] * refbasis[1, 1]
        refbasis[offset + 4, 1] = 6 * xref[1] * xref[2]
        refbasis[offset + 5, 1] = 6 * xref[1] * xref[3]
        refbasis[offset + 6, 1] = 6 * xref[2] * xref[3]
        for j in 1:6, k in 2:3
            refbasis[offset + j, k] = refbasis[offset + j, 1]
        end
        return nothing
    end
end

function get_coefficients(::Type{ON_CELLS}, FE::FESpace{Tv, Ti, H1P1TEB{3}, APT}, ::Type{<:Tetrahedron3D}, xgrid) where {Tv, Ti, APT}
    xEdgeTangents::Array{Tv, 2} = xgrid[EdgeTangents]
    xCellEdges = xgrid[CellEdges]
    return function closure(coefficients::Array{<:Real, 2}, cell)
        fill!(coefficients, 1.0)
        for e in 1:6, k in 1:2
            coefficients[k, 8 + e] = xEdgeTangents[k, xCellEdges[e, cell]]
        end
        return nothing
    end
end

function get_coefficients(::Type{<:ON_FACES}, FE::FESpace{Tv, Ti, H1P1TEB{3}, APT}, ::Type{<:Triangle2D}, xgrid) where {Tv, Ti, APT}
    xEdgeTangents::Array{Tv, 2} = xgrid[EdgeTangents]
    xFaceEdges = xgrid[FaceEdges]
    return function closure(coefficients::Array{<:Real, 2}, face)
        fill!(coefficients, 1.0)
        for e in 1:3, k in 1:2
            coefficients[k, 6 + e] = xEdgeTangents[k, xFaceEdges[e, face]]
        end
        return nothing
    end
end
