"""
````
abstract type H1BR{edim} <: AbstractH1FiniteElementWithCoefficients where {edim<:Int}
````

vector-valued (ncomponents = edim) Bernardi--Raugel element
(first-order polynomials + normal-weighted face bubbles)

allowed ElementGeometries:
- Triangle2D (piecewise linear + normal-weighted face bubbles)
- Quadrilateral2D (Q1 space + normal-weighted face bubbles)
- Tetrahedron3D (piecewise linear + normal-weighted face bubbles)
"""
abstract type H1BR{edim} <: AbstractH1FiniteElementWithCoefficients where {edim <: Int} end
H1BR(edim::Int) = H1BR{edim}

function Base.show(io::Core.IO, ::Type{<:H1BR{edim}}) where {edim}
    return print(io, "H1BR{$edim}")
end

get_ncomponents(FEType::Type{<:H1BR}) = FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:H1BR}, EG::Type{<:AbstractElementGeometry}) = 1 + num_nodes(EG) * FEType.parameters[1]
get_ndofs(::Type{ON_CELLS}, FEType::Type{<:H1BR}, EG::Type{<:AbstractElementGeometry}) = num_faces(EG) + num_nodes(EG) * FEType.parameters[1]

get_polynomialorder(::Type{<:H1BR{2}}, ::Type{<:Edge1D}) = 2;
get_polynomialorder(::Type{<:H1BR{2}}, ::Type{<:Triangle2D}) = 2;
get_polynomialorder(::Type{<:H1BR{2}}, ::Type{<:Quadrilateral2D}) = 3;
get_polynomialorder(::Type{<:H1BR{3}}, ::Type{<:Triangle2D}) = 3;
get_polynomialorder(::Type{<:H1BR{3}}, ::Type{<:Tetrahedron3D}) = 3;
get_polynomialorder(::Type{<:H1BR{3}}, ::Type{<:Parallelogram2D}) = 4;
get_polynomialorder(::Type{<:H1BR{3}}, ::Type{<:Hexahedron3D}) = 5;

get_dofmap_pattern(FEType::Type{<:H1BR}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry}) = "N1f1"
get_dofmap_pattern(FEType::Type{<:H1BR}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry}) = "N1i1"

isdefined(FEType::Type{<:H1BR}, ::Type{<:Triangle2D}) = true
isdefined(FEType::Type{<:H1BR}, ::Type{<:Quadrilateral2D}) = true
isdefined(FEType::Type{<:H1BR}, ::Type{<:Tetrahedron3D}) = true

interior_dofs_offset(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, ::Type{H1BR{2}}, ::Type{Edge1D}) = 4
interior_dofs_offset(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, ::Type{H1BR{3}}, ::Type{Triangle2D}) = 9

function BR_normalflux_eval!(result, f, qpinfo)
    return result[1] = dot(f, qpinfo.normal)
end
init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{AT_NODES}) where {Tv, Ti, FEType <: H1BR, APT} = NodalInterpolator(FES)
init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}) where {Tv, Ti, FEType <: H1BR{2}, APT} = FunctionalInterpolator(BR_normalflux_eval!, FES, ON_FACES; operator = NormalFlux, dofs = [5])
init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}) where {Tv, Ti, FEType <: H1BR{3}, APT} = FunctionalInterpolator(BR_normalflux_eval!, FES, ON_FACES; operator = NormalFlux, dofs = [10])


function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{AT_NODES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1BR, APT}
    return get_interpolator(FE, AT_NODES).evaluate!(Target, exact_function!, items; kwargs...)
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_EDGES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1BR, APT}
    # delegate edge nodes to node interpolation
    subitems = slice(FE.dofgrid[EdgeNodes], items)
    return interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)
end

function ExtendableGrids.interpolate!(Target::AbstractVector{T}, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}, exact_function!; items = [], kwargs...) where {T, Tv, Ti, FEType <: H1BR, APT}
    # delegate face nodes to node interpolation
    subitems = slice(FE.dofgrid[FaceNodes], items)
    interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)

    # preserve face means in normal direction
    get_interpolator(FE, ON_FACES).evaluate!(Target, exact_function!, items; kwargs...)
    return
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1BR, APT}
    # delegate cell faces to face interpolation
    subitems = slice(FE.dofgrid[CellFaces], items)
    return interpolate!(Target, FE, ON_FACES, exact_function!; items = subitems, kwargs...)
end


############
# 2D basis #
############

function get_basis(AT::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{H1BR{2}}, EG::Type{<:Edge1D})
    refbasis_P1 = get_basis(AT, H1P1{2}, EG)
    offset = get_ndofs(AT, H1P1{2}, EG)
    return function closure(refbasis, xref)
        refbasis_P1(refbasis, xref)
        # add face bubble to P1 basis
        refbasis[offset + 1, 1] = 6 * xref[1] * refbasis[1, 1]
        return refbasis[offset + 1, 2] = refbasis[offset + 1, 1]
    end
end

function get_basis(AT::Type{ON_CELLS}, FEType::Type{H1BR{2}}, EG::Type{<:Triangle2D})
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


function get_basis(AT::Type{ON_CELLS}, FEType::Type{H1BR{2}}, EG::Type{<:Quadrilateral2D})
    refbasis_P1 = get_basis(AT, H1Q1{2}, EG)
    offset = get_ndofs(AT, H1Q1{2}, EG)
    return function closure(refbasis, xref)
        refbasis_P1(refbasis, xref)
        # add face bubbles to Q1 basis
        refbasis[offset + 1, 1] = 6 * xref[1] * (1 - xref[1]) * (1 - xref[2])
        refbasis[offset + 2, 1] = 6 * xref[2] * xref[1] * (1 - xref[2])
        refbasis[offset + 3, 1] = 6 * xref[1] * xref[2] * (1 - xref[1])
        refbasis[offset + 4, 1] = 6 * xref[2] * (1 - xref[1]) * (1 - xref[2])
        refbasis[offset + 1, 2] = refbasis[offset + 1, 1]
        refbasis[offset + 2, 2] = refbasis[offset + 2, 1]
        refbasis[offset + 3, 2] = refbasis[offset + 3, 1]
        return refbasis[offset + 4, 2] = refbasis[offset + 4, 1]
    end
end

function get_coefficients(::Type{ON_CELLS}, FE::FESpace{Tv, Ti, H1BR{2}, APT}, ::Type{<:Triangle2D}, xgrid) where {Tv, Ti, APT}
    xFaceNormals::Array{Tv, 2} = xgrid[FaceNormals]
    xCellFaces = xgrid[CellFaces]
    return function closure(coefficients::Array{<:Real, 2}, cell)
        fill!(coefficients, 1.0)
        coefficients[1, 7] = xFaceNormals[1, xCellFaces[1, cell]]
        coefficients[2, 7] = xFaceNormals[2, xCellFaces[1, cell]]
        coefficients[1, 8] = xFaceNormals[1, xCellFaces[2, cell]]
        coefficients[2, 8] = xFaceNormals[2, xCellFaces[2, cell]]
        coefficients[1, 9] = xFaceNormals[1, xCellFaces[3, cell]]
        return coefficients[2, 9] = xFaceNormals[2, xCellFaces[3, cell]]
    end
end

function get_coefficients(::Type{ON_CELLS}, FE::FESpace{Tv, Ti, H1BR{2}, APT}, ::Type{<:Quadrilateral2D}, xgrid) where {Tv, Ti, APT}
    xFaceNormals::Array{Tv, 2} = xgrid[FaceNormals]
    xCellFaces = xgrid[CellFaces]
    return function closure(coefficients::Array{<:Real, 2}, cell)
        fill!(coefficients, 1.0)
        coefficients[1, 9] = xFaceNormals[1, xCellFaces[1, cell]]
        coefficients[2, 9] = xFaceNormals[2, xCellFaces[1, cell]]
        coefficients[1, 10] = xFaceNormals[1, xCellFaces[2, cell]]
        coefficients[2, 10] = xFaceNormals[2, xCellFaces[2, cell]]
        coefficients[1, 11] = xFaceNormals[1, xCellFaces[3, cell]]
        coefficients[2, 11] = xFaceNormals[2, xCellFaces[3, cell]]
        coefficients[1, 12] = xFaceNormals[1, xCellFaces[4, cell]]
        return coefficients[2, 12] = xFaceNormals[2, xCellFaces[4, cell]]
    end
end

function get_coefficients(::Type{<:ON_FACES}, FE::FESpace{Tv, Ti, H1BR{2}, APT}, ::Type{<:Edge1D}, xgrid) where {Tv, Ti, APT}
    xFaceNormals::Array{Tv, 2} = xgrid[FaceNormals]
    return function closure(coefficients::Array{<:Real, 2}, face)
        # multiplication of face bubble with normal vector of face
        fill!(coefficients, 1.0)
        coefficients[1, 5] = xFaceNormals[1, face]
        return coefficients[2, 5] = xFaceNormals[2, face]
    end
end

############
# 3D basis #
############

function get_basis(AT::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, ::Type{H1BR{3}}, EG::Type{<:Triangle2D})
    refbasis_P1 = get_basis(AT, H1P1{3}, EG)
    offset = get_ndofs(AT, H1P1{3}, EG)
    return function closure(refbasis, xref)
        refbasis_P1(refbasis, xref)
        # add face bubbles to P1 basis
        refbasis[offset + 1, 1] = 60 * xref[1] * refbasis[1, 1] * xref[2]
        refbasis[offset + 1, 2] = refbasis[offset + 1, 1]
        return refbasis[offset + 1, 3] = refbasis[offset + 1, 1]
    end
end

function get_basis(AT::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, ::Type{H1BR{3}}, EG::Type{<:Quadrilateral2D})
    refbasis_P1 = get_basis(AT, H1Q1{3}, EG)
    offset = get_ndofs(AT, H1Q1{3}, EG)
    return function closure(refbasis, xref)
        refbasis_P1(refbasis, xref)
        # add face bubbles to P1 basis
        refbasis[offset + 1, 1] = 36 * xref[1] * (1 - xref[1]) * (1 - xref[2]) * xref[2]
        refbasis[offset + 1, 2] = refbasis[offset + 1, 1]
        return refbasis[offset + 1, 3] = refbasis[offset + 1, 1]
    end
end

function get_basis(AT::Type{ON_CELLS}, ::Type{H1BR{3}}, EG::Type{<:Tetrahedron3D})
    refbasis_P1 = get_basis(AT, H1P1{3}, EG)
    offset = get_ndofs(AT, H1P1{3}, EG)
    return function closure(refbasis, xref)
        refbasis_P1(refbasis, xref)
        # add face bubbles to P1 basis
        refbasis[offset + 1, 1] = 60 * xref[1] * refbasis[1, 1] * xref[2]
        refbasis[offset + 2, 1] = 60 * refbasis[1, 1] * xref[1] * xref[3]
        refbasis[offset + 3, 1] = 60 * xref[1] * xref[2] * xref[3]
        refbasis[offset + 4, 1] = 60 * refbasis[1, 1] * xref[2] * xref[3]
        for j in 1:4, k in 2:3
            refbasis[offset + j, k] = refbasis[offset + j, 1]
        end
        return
    end
end

function get_basis(AT::Type{ON_CELLS}, ::Type{H1BR{3}}, EG::Type{<:Hexahedron3D})
    refbasis_P1 = get_basis(AT, H1Q1{3}, EG)
    offset = get_ndofs(AT, H1Q1{3}, EG)
    return function closure(refbasis, xref)
        refbasis_P1(refbasis, xref)
        # add face bubbles to Q1 basis
        refbasis[offset + 1, 1] = 36 * (1 - xref[1]) * (1 - xref[2]) * xref[1] * xref[2] * (1 - xref[3]) # bottom
        refbasis[offset + 2, 1] = 36 * (1 - xref[1]) * xref[1] * (1 - xref[3]) * xref[3] * (1 - xref[2]) # front
        refbasis[offset + 3, 1] = 36 * (1 - xref[1]) * (1 - xref[2]) * (1 - xref[3]) * xref[2] * xref[3] # left
        refbasis[offset + 4, 1] = 36 * (1 - xref[1]) * xref[1] * (1 - xref[3]) * xref[3] * xref[2]       # back
        refbasis[offset + 5, 1] = 36 * xref[1] * (1 - xref[2]) * (1 - xref[3]) * xref[2] * xref[3]       # right
        refbasis[offset + 6, 1] = 36 * (1 - xref[1]) * (1 - xref[2]) * xref[1] * xref[2] * xref[3]       # top
        for j in 1:6, k in 2:3
            refbasis[offset + j, k] = refbasis[offset + j, 1]
        end
        return
    end
end


function get_coefficients(::Type{ON_CELLS}, FE::FESpace{Tv, Ti, H1BR{3}, APT}, ::Type{<:Tetrahedron3D}, xgrid) where {Tv, Ti, APT}
    xFaceNormals::Array{Tv, 2} = xgrid[FaceNormals]
    xCellFaces::Adjacency{Ti} = xgrid[CellFaces]
    return function closure(coefficients::Array{<:Real, 2}, cell)
        # multiplication with normal vectors
        fill!(coefficients, 1.0)
        coefficients[1, 13] = xFaceNormals[1, xCellFaces[1, cell]]
        coefficients[2, 13] = xFaceNormals[2, xCellFaces[1, cell]]
        coefficients[3, 13] = xFaceNormals[3, xCellFaces[1, cell]]
        coefficients[1, 14] = xFaceNormals[1, xCellFaces[2, cell]]
        coefficients[2, 14] = xFaceNormals[2, xCellFaces[2, cell]]
        coefficients[3, 14] = xFaceNormals[3, xCellFaces[2, cell]]
        coefficients[1, 15] = xFaceNormals[1, xCellFaces[3, cell]]
        coefficients[2, 15] = xFaceNormals[2, xCellFaces[3, cell]]
        coefficients[3, 15] = xFaceNormals[3, xCellFaces[3, cell]]
        coefficients[1, 16] = xFaceNormals[1, xCellFaces[4, cell]]
        coefficients[2, 16] = xFaceNormals[2, xCellFaces[4, cell]]
        coefficients[3, 16] = xFaceNormals[3, xCellFaces[4, cell]]
        return nothing
    end
end


function get_coefficients(::Type{<:ON_FACES}, FE::FESpace{Tv, Ti, H1BR{3}, APT}, ::Type{<:Triangle2D}, xgrid) where {Tv, Ti, APT}
    xFaceNormals::Array{Tv, 2} = xgrid[FaceNormals]
    return function closure(coefficients::Array{<:Real, 2}, face)
        # multiplication of face bubble with normal vector of face
        fill!(coefficients, 1.0)
        coefficients[1, 10] = xFaceNormals[1, face]
        coefficients[2, 10] = xFaceNormals[2, face]
        return coefficients[3, 10] = xFaceNormals[3, face]
    end
end

function get_coefficients(::Type{ON_CELLS}, FE::FESpace{Tv, Ti, H1BR{3}, APT}, ::Type{<:Hexahedron3D}, xgrid) where {Tv, Ti, APT}
    xFaceNormals::Array{Tv, 2} = xgrid[FaceNormals]
    xCellFaces::Adjacency{Ti} = xgrid[CellFaces]
    return function closure(coefficients::Array{<:Real, 2}, cell)
        # multiplication with normal vectors
        fill!(coefficients, 1.0)
        coefficients[1, 25] = xFaceNormals[1, xCellFaces[1, cell]]
        coefficients[2, 25] = xFaceNormals[2, xCellFaces[1, cell]]
        coefficients[3, 25] = xFaceNormals[3, xCellFaces[1, cell]]
        coefficients[1, 26] = xFaceNormals[1, xCellFaces[2, cell]]
        coefficients[2, 26] = xFaceNormals[2, xCellFaces[2, cell]]
        coefficients[3, 26] = xFaceNormals[3, xCellFaces[2, cell]]
        coefficients[1, 27] = xFaceNormals[1, xCellFaces[3, cell]]
        coefficients[2, 27] = xFaceNormals[2, xCellFaces[3, cell]]
        coefficients[3, 27] = xFaceNormals[3, xCellFaces[3, cell]]
        coefficients[1, 28] = xFaceNormals[1, xCellFaces[4, cell]]
        coefficients[2, 28] = xFaceNormals[2, xCellFaces[4, cell]]
        coefficients[3, 28] = xFaceNormals[3, xCellFaces[4, cell]]
        coefficients[1, 29] = xFaceNormals[1, xCellFaces[5, cell]]
        coefficients[2, 29] = xFaceNormals[2, xCellFaces[5, cell]]
        coefficients[3, 29] = xFaceNormals[3, xCellFaces[5, cell]]
        coefficients[1, 30] = xFaceNormals[1, xCellFaces[6, cell]]
        coefficients[2, 30] = xFaceNormals[2, xCellFaces[6, cell]]
        coefficients[3, 30] = xFaceNormals[3, xCellFaces[6, cell]]
        return nothing
    end
end


function get_coefficients(::Type{<:ON_FACES}, FE::FESpace{Tv, Ti, H1BR{3}, APT}, ::Type{<:Quadrilateral2D}, xgrid) where {Tv, Ti, APT}
    xFaceNormals::Array{Tv, 2} = xgrid[FaceNormals]
    return function closure(coefficients::Array{<:Real, 2}, face)
        # multiplication of face bubble with normal vector of face
        fill!(coefficients, 1.0)
        coefficients[1, 13] = xFaceNormals[1, face]
        coefficients[2, 13] = xFaceNormals[2, face]
        coefficients[3, 13] = xFaceNormals[3, face]
        return nothing
    end
end


###########################
# RT0/BDM1 Reconstruction #
###########################


function get_reconstruction_coefficients!(xgrid::ExtendableGrid{Tv, Ti}, ::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FE::Type{<:H1BR{2}}, FER::Type{<:HDIVRT0{2}}, ::Type{<:Edge1D}) where {Tv, Ti}
    xFaceVolumes::Array{<:Real, 1} = xgrid[FaceVolumes]
    xFaceNormals::Array{<:Real, 2} = xgrid[FaceNormals]
    return function closure(coefficients::Array{<:Real, 2}, face::Int)
        coefficients[1, 1] = 1 // 2 * xFaceVolumes[face] * xFaceNormals[1, face]
        coefficients[2, 1] = 1 // 2 * xFaceVolumes[face] * xFaceNormals[1, face]
        coefficients[3, 1] = 1 // 2 * xFaceVolumes[face] * xFaceNormals[2, face]
        coefficients[4, 1] = 1 // 2 * xFaceVolumes[face] * xFaceNormals[2, face]
        coefficients[5, 1] = xFaceVolumes[face]
        return nothing
    end
end


function get_reconstruction_coefficients!(xgrid::ExtendableGrid{Tv, Ti}, ::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FE::Type{<:H1BR{2}}, FER::Type{<:HDIVBDM1{2}}, ::Type{<:Edge1D}) where {Tv, Ti}
    xFaceVolumes::Array{Tv, 1} = xgrid[FaceVolumes]
    xFaceNormals::Array{Tv, 2} = xgrid[FaceNormals]
    return function closure(coefficients::Array{<:Real, 2}, face)

        coefficients[1, 1] = 1 // 2 * xFaceVolumes[face] * xFaceNormals[1, face]
        coefficients[1, 2] = 1 // 12 * xFaceVolumes[face] * xFaceNormals[1, face]

        coefficients[2, 1] = 1 // 2 * xFaceVolumes[face] * xFaceNormals[1, face]
        coefficients[2, 2] = -1 // 12 * xFaceVolumes[face] * xFaceNormals[1, face]

        coefficients[3, 1] = 1 // 2 * xFaceVolumes[face] * xFaceNormals[2, face]
        coefficients[3, 2] = 1 // 12 * xFaceVolumes[face] * xFaceNormals[2, face]

        coefficients[4, 1] = 1 // 2 * xFaceVolumes[face] * xFaceNormals[2, face]
        coefficients[4, 2] = -1 // 12 * xFaceVolumes[face] * xFaceNormals[2, face]

        coefficients[5, 1] = xFaceVolumes[face]
        return nothing
    end
end


function get_reconstruction_coefficients!(xgrid::ExtendableGrid{Tv, Ti}, ::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FE::Type{<:H1BR{3}}, FER::Type{<:HDIVBDM1{3}}, EG::Type{<:Triangle2D}) where {Tv, Ti}
    xFaceVolumes::Array{Tv, 1} = xgrid[FaceVolumes]
    xFaceNormals::Array{Tv, 2} = xgrid[FaceNormals]
    nfacenodes::Int = num_nodes(EG)
    return function closure(coefficients::Array{<:Real, 2}, face::Int)
        for j in 1:nfacenodes, k in 1:3
            coefficients[(j - 1) * 3 + j, 1] = 1 // nfacenodes * xFaceVolumes[face] * xFaceNormals[k, face]
            coefficients[(j - 1) * 3 + j, 2] = -1 // 36 * xFaceVolumes[face] * xFaceNormals[k, face]
            coefficients[(j - 1) * 3 + j, 3] = -1 // 36 * xFaceVolumes[face] * xFaceNormals[k, face]
        end
        coefficients[nfacenodes * 3 + 1, 1] = xFaceVolumes[face]
        return nothing
    end
end
