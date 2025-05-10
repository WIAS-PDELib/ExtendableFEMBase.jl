"""
````
abstract type H1P3{ncomponents,edim} <: AbstractH1FiniteElement where {ncomponents<:Int,edim<:Int}
````

Continuous piecewise third-order polynomials.

allowed ElementGeometries:
- Edge1D
- Triangle2D
- Tetrahedron3D
"""
abstract type H1P3{ncomponents, edim} <: AbstractH1FiniteElement where {ncomponents <: Int, edim <: Int} end
H1P3(ncomponents::Int, edim = ncomponents) = H1P3{ncomponents, edim}

function Base.show(io::Core.IO, ::Type{<:H1P3{ncomponents, edim}}) where {ncomponents, edim}
    return print(io, "H1P3{$ncomponents,$edim}")
end

get_ncomponents(FEType::Type{<:H1P3}) = FEType.parameters[1]
get_edim(FEType::Type{<:H1P3}) = FEType.parameters[2]

get_ndofs(::Type{<:AssemblyType}, FEType::Type{<:H1P3}, EG::Type{<:AbstractElementGeometry0D}) = FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_EDGES}, Type{<:ON_BEDGES}}, FEType::Type{<:H1P3}, EG::Type{<:Union{AbstractElementGeometry1D}}) = FEType.parameters[1] * Int(factorial(FEType.parameters[2] + 1) / (6 * factorial(FEType.parameters[2] - 2)))
get_ndofs(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:H1P3}, EG::Type{<:Union{AbstractElementGeometry1D, Triangle2D}}) = FEType.parameters[1] * Int(factorial(FEType.parameters[2] + 2) / (6 * factorial(FEType.parameters[2] - 1)))
get_ndofs(::Type{<:ON_CELLS}, FEType::Type{<:H1P3}, EG::Type{<:Union{AbstractElementGeometry1D, Triangle2D, Tetrahedron3D}}) = FEType.parameters[1] * Int(factorial(FEType.parameters[2] + 3) / (6 * factorial(FEType.parameters[2])))

get_polynomialorder(::Type{<:H1P3}, ::Type{<:Edge1D}) = 3;
get_polynomialorder(::Type{<:H1P3}, ::Type{<:Triangle2D}) = 3;
get_polynomialorder(::Type{<:H1P3}, ::Type{<:Tetrahedron3D}) = 3;

get_dofmap_pattern(FEType::Type{<:H1P3}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry1D}) = "N1I2"
get_dofmap_pattern(FEType::Type{<:H1P3}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry2D}) = "N1F2I1"
get_dofmap_pattern(FEType::Type{<:H1P3}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry3D}) = "N1E2F1"
get_dofmap_pattern(FEType::Type{<:H1P3}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry0D}) = "N1"
get_dofmap_pattern(FEType::Type{<:H1P3}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry1D}) = "N1I2"
get_dofmap_pattern(FEType::Type{<:H1P3}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry2D}) = "N1E2I1"
get_dofmap_pattern(FEType::Type{<:H1P3}, ::Union{Type{EdgeDofs}, Type{BEdgeDofs}}, EG::Type{<:AbstractElementGeometry1D}) = "N1I2"

isdefined(FEType::Type{<:H1P3}, ::Type{<:AbstractElementGeometry1D}) = true
isdefined(FEType::Type{<:H1P3}, ::Type{<:Triangle2D}) = true
isdefined(FEType::Type{<:H1P3}, ::Type{<:Tetrahedron3D}) = true

interior_dofs_offset(::Type{<:AssemblyType}, ::Type{H1P3{ncomponents, edim}}, ::Type{Edge1D}) where {ncomponents, edim} = 2
interior_dofs_offset(::Type{<:AssemblyType}, ::Type{H1P3{ncomponents, edim}}, ::Type{Triangle2D}) where {ncomponents, edim} = 9
interior_dofs_offset(::Type{<:AssemblyType}, ::Type{H1P3{ncomponents, edim}}, ::Type{Tetrahedron3D}) where {ncomponents, edim} = 20

get_ref_cellmoments(::Type{<:H1P3}, ::Type{<:Triangle2D}) = [1 // 30, 1 // 30, 1 // 30, 3 // 40, 3 // 40, 3 // 40, 3 // 40, 3 // 40, 3 // 40, 1 // 1] # integrals of 1D basis functions over reference cell (divided by volume)

init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{AT_NODES}) where {Tv, Ti, FEType <: H1P3, APT} = NodalInterpolator(FES)
init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_EDGES}) where {Tv, Ti, FEType <: H1P3, APT} = MomentInterpolator(FES, ON_EDGES; order = 4 - get_edim(FEType))
init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}) where {Tv, Ti, FEType <: H1P3, APT} = MomentInterpolator(FES, ON_FACES; order = 3 - get_edim(FEType))
init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}) where {Tv, Ti, FEType <: H1P3, APT} = MomentInterpolator(FES, ON_CELLS; order = 2 - get_edim(FEType))

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{AT_NODES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1P3, APT}
    return get_interpolator(FE, AT_NODES).evaluate!(Target, exact_function!, items; kwargs...)
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_EDGES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1P3, APT}
    edim = get_edim(FEType)
    return if edim == 3
        # delegate edge nodes to node interpolation
        subitems = slice(FE.dofgrid[EdgeNodes], items)
        interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)

        # perform edge mean interpolation
        get_interpolator(FE, ON_EDGES).evaluate!(Target, exact_function!, items; kwargs...)
    end
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1P3, APT}
    edim = get_edim(FEType)
    return if edim == 2
        # delegate face nodes to node interpolation
        subitems = slice(FE.dofgrid[FaceNodes], items)
        interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)

        # perform face mean interpolation
        get_interpolator(FE, ON_FACES).evaluate!(Target, exact_function!, items; kwargs...)
    elseif edim == 3
        # delegate face edges to edge interpolation
        subitems = slice(FE.dofgrid[FaceEdges], items)
        interpolate!(Target, FE, ON_EDGES, exact_function!; items = subitems, kwargs...)

        # preserve face integral
        get_interpolator(FE, ON_FACES).evaluate!(Target, exact_function!, items; kwargs...)
    elseif edim == 1
        # delegate face nodes to node interpolation
        subitems = slice(FE.dofgrid[FaceNodes], items)
        interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)
    end
end


function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1P3, APT}
    edim = get_edim(FEType)
    return if edim == 2
        # delegate cell faces to face interpolation
        subitems = slice(FE.dofgrid[CellFaces], items)
        interpolate!(Target, FE, ON_FACES, exact_function!; items = subitems, kwargs...)

        # fix cell bubble value by preserving integral mean
        get_interpolator(FE, ON_CELLS).evaluate!(Target, exact_function!, items; kwargs...)
    elseif edim == 3
        # delegate cell faces to face interpolation
        subitems = slice(FE.dofgrid[CellFaces], items)
        interpolate!(Target, FE, ON_FACES, exact_function!; items = subitems, kwargs...)
    elseif edim == 1
        # delegate cell nodes to node interpolation
        subitems = slice(FE.dofgrid[CellNodes], items)
        interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)

        # preserve cell integral
        get_interpolator(FE, ON_CELLS).evaluate!(Target, exact_function!, items; kwargs...)
    end
end


function get_basis(::Type{<:AssemblyType}, ::Type{H1P3{ncomponents, edim}}, ::Type{<:Vertex0D}) where {ncomponents, edim}
    return function closure(refbasis, xref)
        for k in 1:ncomponents
            refbasis[k, k] = 1
        end
        return
    end
end

function get_basis(::Type{<:AssemblyType}, ::Type{H1P3{ncomponents, edim}}, ::Type{<:Edge1D}) where {ncomponents, edim}
    return function closure(refbasis, xref)
        refbasis[end] = 1 - xref[1]
        for k in 1:ncomponents
            refbasis[4 * k - 3, k] = 9 // 2 * refbasis[end] * (refbasis[end] - 1 // 3) * (refbasis[end] - 2 // 3)    # node 1 (scaled such that 1 at x = 0)
            refbasis[4 * k - 2, k] = 9 // 2 * xref[1] * (xref[1] - 1 // 3) * (xref[1] - 2 // 3)                      # node 2 (scaled such that 1 at x = 1)
            refbasis[4 * k - 1, k] = -27 / 2 * xref[1] * refbasis[end] * (xref[1] - 2 // 3)                                # face 1 (scaled such that 1 at x = 1//3)
            refbasis[4 * k, k] = 27 // 2 * xref[1] * refbasis[end] * (xref[1] - 1 // 3)                                  # face 2 (scaled such that 1 at x = 2//3)
        end
        return
    end
end

function get_basis(::Type{<:AssemblyType}, ::Type{H1P3{ncomponents, edim}}, ::Type{<:Triangle2D}) where {ncomponents, edim}
    function closure(refbasis, xref)
        refbasis[end] = 1 - xref[1] - xref[2]
        for k in 1:ncomponents
            refbasis[10 * k - 9, k] = 9 // 2 * refbasis[end] * (refbasis[end] - 1 // 3) * (refbasis[end] - 2 // 3)  # node 1
            refbasis[10 * k - 8, k] = 9 // 2 * xref[1] * (xref[1] - 1 // 3) * (xref[1] - 2 // 3)                    # node 2
            refbasis[10 * k - 7, k] = 9 // 2 * xref[2] * (xref[2] - 1 // 3) * (xref[2] - 2 // 3)                    # node 3
            refbasis[10 * k - 6, k] = 27 // 2 * xref[1] * refbasis[end] * (refbasis[end] - 1 // 3)                     # face 1.1
            refbasis[10 * k - 5, k] = 27 // 2 * xref[1] * refbasis[end] * (xref[1] - 1 // 3)                           # face 1.2
            refbasis[10 * k - 4, k] = 27 // 2 * xref[2] * xref[1] * (xref[1] - 1 // 3)                                 # face 2.1
            refbasis[10 * k - 3, k] = 27 // 2 * xref[2] * xref[1] * (xref[2] - 1 // 3)                                 # face 2.2
            refbasis[10 * k - 2, k] = 27 // 2 * refbasis[end] * xref[2] * (xref[2] - 1 // 3)                           # face 3.1
            refbasis[10 * k - 1, k] = 27 // 2 * refbasis[end] * xref[2] * (refbasis[end] - 1 // 3)                     # face 3.2
            refbasis[10 * k, k] = 27 * xref[1] * xref[2] * refbasis[end]                                             # cell (scaled such that 1 at x,y = 1//3)
        end
        return
    end
    return closure
end


function get_basis(::Type{<:AssemblyType}, ::Type{H1P3{ncomponents, edim}}, ::Type{<:Tetrahedron3D}) where {ncomponents, edim}
    function closure(refbasis, xref)
        refbasis[end] = 1 // 1 - xref[1] - xref[2] - xref[3]
        for k in 1:ncomponents
            ## assuming
            ## _local_cellfacenodes_Tetrahedron3D = [1 3 2; 1 2 4; 2 3 4; 1 4 3]'
            ## _local_celledgenodes_Tetrahedron3D = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]'
            refbasis[20 * k - 19, k] = 9 // 2 * refbasis[end] * (refbasis[end] - 1 // 3) * (refbasis[end] - 2 // 3)  # node 1
            refbasis[20 * k - 18, k] = 9 // 2 * xref[1] * (xref[1] - 1 // 3) * (xref[1] - 2 // 3)                    # node 2
            refbasis[20 * k - 17, k] = 9 // 2 * xref[2] * (xref[2] - 1 // 3) * (xref[2] - 2 // 3)                    # node 3
            refbasis[20 * k - 16, k] = 9 // 2 * xref[3] * (xref[3] - 1 // 3) * (xref[3] - 2 // 3)                    # node 4

            refbasis[20 * k - 15, k] = 27 // 2 * xref[1] * refbasis[end] * (refbasis[end] - 1 // 3)                     # edge 1.1
            refbasis[20 * k - 14, k] = 27 // 2 * xref[1] * refbasis[end] * (xref[1] - 1 // 3)                           # edge 1.2
            refbasis[20 * k - 13, k] = 27 // 2 * xref[2] * refbasis[end] * (refbasis[end] - 1 // 3)                     # edge 2.1
            refbasis[20 * k - 12, k] = 27 // 2 * xref[2] * refbasis[end] * (xref[2] - 1 // 3)                           # edge 2.2
            refbasis[20 * k - 11, k] = 27 // 2 * xref[3] * refbasis[end] * (refbasis[end] - 1 // 3)                     # edge 3.1
            refbasis[20 * k - 10, k] = 27 // 2 * xref[3] * refbasis[end] * (xref[3] - 1 // 3)                           # edge 3.2
            refbasis[20 * k - 9, k] = 27 // 2 * xref[1] * xref[2] * (xref[1] - 1 // 3)                                  # edge 4.1
            refbasis[20 * k - 8, k] = 27 // 2 * xref[1] * xref[2] * (xref[2] - 1 // 3)                                  # edge 4.2
            refbasis[20 * k - 7, k] = 27 // 2 * xref[1] * xref[3] * (xref[1] - 1 // 3)                                  # edge 5.1
            refbasis[20 * k - 6, k] = 27 // 2 * xref[1] * xref[3] * (xref[3] - 1 // 3)                                  # edge 5.2
            refbasis[20 * k - 5, k] = 27 // 2 * xref[2] * xref[3] * (xref[2] - 1 // 3)                                  # edge 6.1
            refbasis[20 * k - 4, k] = 27 // 2 * xref[2] * xref[3] * (xref[3] - 1 // 3)                                  # edge 6.2

            refbasis[20 * k - 3, k] = 27 // 1 * xref[1] * xref[2] * refbasis[end]                                            # face 1
            refbasis[20 * k - 2, k] = 27 // 1 * xref[1] * xref[3] * refbasis[end]                                            # face 2
            refbasis[20 * k - 1, k] = 27 // 1 * xref[1] * xref[2] * xref[3]                                                  # face 3
            refbasis[20 * k, k] = 27 // 1 * xref[2] * xref[3] * refbasis[end]                                              # face 4

        end
        return
    end
    return closure
end

# we need to change the ordering of the face dofs on faces that have a negative orientation sign
function get_basissubset(::Type{ON_CELLS}, FE::FESpace{Tv, Ti, H1P3{ncomponents, edim}, APT}, EG::Type{<:Triangle2D}, xgrid) where {ncomponents, edim, Tv, Ti, APT}
    xCellFaceSigns = xgrid[CellFaceSigns]
    nfaces::Int = num_faces(EG)
    return function closure(subset_ids::Array{Int, 1}, cell)
        for j in 1:nfaces
            if xCellFaceSigns[j, cell] != 1
                for c in 1:ncomponents
                    subset_ids[(c - 1) * 10 + 3 + 2 * j - 1] = (c - 1) * 10 + 3 + 2 * j
                    subset_ids[(c - 1) * 10 + 3 + 2 * j] = (c - 1) * 10 + 3 + 2 * j - 1
                end
            else
                for c in 1:ncomponents
                    subset_ids[(c - 1) * 10 + 3 + 2 * j - 1] = (c - 1) * 10 + 3 + 2 * j - 1
                    subset_ids[(c - 1) * 10 + 3 + 2 * j] = (c - 1) * 10 + 3 + 2 * j
                end
            end
        end
        return nothing
    end
end

# we need to change the ordering of the face dofs on faces that have a negative orientation sign
function get_basissubset(::Type{ON_CELLS}, FE::FESpace{Tv, Ti, H1P3{ncomponents, edim}, APT}, EG::Type{<:Tetrahedron3D}, xgrid) where {ncomponents, edim, Tv, Ti, APT}
    xCellEdgeSigns = xgrid[CellEdgeSigns]
    nedges::Int = num_edges(EG)
    return function closure(subset_ids::Array{Int, 1}, cell)
        for j in 1:nedges
            if xCellEdgeSigns[j, cell] != 1
                for c in 1:ncomponents
                    subset_ids[(c - 1) * 20 + 4 + 2 * j - 1] = (c - 1) * 20 + 4 + 2 * j
                    subset_ids[(c - 1) * 20 + 4 + 2 * j] = (c - 1) * 20 + 4 + 2 * j - 1
                end
            else
                for c in 1:ncomponents
                    subset_ids[(c - 1) * 20 + 4 + 2 * j - 1] = (c - 1) * 20 + 4 + 2 * j - 1
                    subset_ids[(c - 1) * 20 + 4 + 2 * j] = (c - 1) * 20 + 4 + 2 * j
                end
            end
        end
        return nothing
    end
end


# we need to change the ordering of the face dofs on faces that have a negative orientation sign
function get_basissubset(::Type{ON_FACES}, FE::FESpace{Tv, Ti, H1P3{ncomponents, edim}, APT}, EG::Type{<:Triangle2D}) where {ncomponents, edim, Tv, Ti, APT}
    xFaceEdgeSigns = FE.dofgrid[FaceEdgeSigns]
    nedges::Int = num_edges(EG)
    return function closure(subset_ids::Array{Int, 1}, face)
        for j in 1:nedges
            if xFaceEdgeSigns[j, face] != 1
                for c in 1:ncomponents
                    subset_ids[(c - 1) * 10 + 3 + 2 * j - 1] = (c - 1) * 10 + 3 + 2 * j
                    subset_ids[(c - 1) * 10 + 3 + 2 * j] = (c - 1) * 10 + 3 + 2 * j - 1
                end
            else
                for c in 1:ncomponents
                    subset_ids[(c - 1) * 10 + 3 + 2 * j - 1] = (c - 1) * 10 + 3 + 2 * j - 1
                    subset_ids[(c - 1) * 10 + 3 + 2 * j] = (c - 1) * 10 + 3 + 2 * j
                end
            end
        end
        return nothing
    end
end

function get_basissubset(::Type{ON_BFACES}, FE::FESpace{Tv, Ti, H1P3{ncomponents, edim}, APT}, EG::Type{<:Triangle2D}) where {ncomponents, edim, Tv, Ti, APT}
    xFaceEdgeSigns = FE.dofgrid[FaceEdgeSigns]
    xBFaceFaces = FE.dofgrid[BFaceFaces]
    nedges::Int = num_edges(EG)
    return function closure(subset_ids::Array{Int, 1}, bface)
        for j in 1:nedges
            if xFaceEdgeSigns[j, xBFaceFaces[bface]] != 1
                for c in 1:ncomponents
                    subset_ids[(c - 1) * 10 + 3 + 2 * j - 1] = (c - 1) * 10 + 3 + 2 * j
                    subset_ids[(c - 1) * 10 + 3 + 2 * j] = (c - 1) * 10 + 3 + 2 * j - 1
                end
            else
                for c in 1:ncomponents
                    subset_ids[(c - 1) * 10 + 3 + 2 * j - 1] = (c - 1) * 10 + 3 + 2 * j - 1
                    subset_ids[(c - 1) * 10 + 3 + 2 * j] = (c - 1) * 10 + 3 + 2 * j
                end
            end
        end
        return nothing
    end
end
