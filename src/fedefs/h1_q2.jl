"""
````
abstract type H1Q2{ncomponents,edim} <: AbstractH1FiniteElement where {ncomponents<:Int,edim<:Int}
````

Continuous piecewise second-order polynomials on simplices and quads. Can be used with mixed geometries (in 2D).

allowed ElementGeometries:
- Edge1D (P2 space)
- Triangle2D (P2 space)
- Quadrilateral2D (Q2 space with cell bubble)
- Tetrahedron3D (P2 space)
"""
abstract type H1Q2{ncomponents, edim} <: AbstractH1FiniteElement where {ncomponents <: Int, edim <: Int} end
H1Q2(ncomponents::Int, edim = ncomponents) = H1Q2{ncomponents, edim}

function Base.show(io::Core.IO, ::Type{<:H1Q2{ncomponents, edim}}) where {ncomponents, edim}
    return print(io, "H1Q2{$ncomponents,$edim}")
end

get_ncomponents(FEType::Type{<:H1Q2}) = FEType.parameters[1]
get_edim(FEType::Type{<:H1Q2}) = FEType.parameters[2]

get_ndofs(::Type{<:AssemblyType}, FEType::Type{<:H1Q2}, EG::Type{<:AbstractElementGeometry0D}) = FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_EDGES}, Type{<:ON_BEDGES}}, FEType::Type{<:H1Q2}, EG::Type{<:AbstractElementGeometry1D}) = 3 * FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:H1Q2}, EG::Type{<:Union{AbstractElementGeometry1D, Triangle2D, Tetrahedron3D}}) = Int((FEType.parameters[2]) * (FEType.parameters[2] + 1) / 2 * FEType.parameters[1])
get_ndofs(::Type{<:ON_CELLS}, FEType::Type{<:H1Q2}, EG::Type{<:Union{AbstractElementGeometry1D, Triangle2D, Tetrahedron3D}}) = Int((FEType.parameters[2] + 1) * (FEType.parameters[2] + 2) / 2 * FEType.parameters[1])
get_ndofs(::Type{<:ON_CELLS}, FEType::Type{<:H1Q2}, EG::Type{<:Quadrilateral2D}) = 9 * FEType.parameters[1]

get_polynomialorder(::Type{<:H1Q2}, ::Type{<:Edge1D}) = 2;
get_polynomialorder(::Type{<:H1Q2}, ::Type{<:Triangle2D}) = 2;
get_polynomialorder(::Type{<:H1Q2}, ::Type{<:Quadrilateral2D}) = 4;
get_polynomialorder(::Type{<:H1Q2}, ::Type{<:Tetrahedron3D}) = 2;

get_dofmap_pattern(FEType::Type{<:H1Q2}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry1D}) = "N1I1"
get_dofmap_pattern(FEType::Type{<:H1Q2}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry2D}) = "N1F1"
get_dofmap_pattern(FEType::Type{<:H1Q2}, ::Type{CellDofs}, EG::Type{<:Quadrilateral2D}) = "N1F1I1"
get_dofmap_pattern(FEType::Type{<:H1Q2}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry3D}) = "N1E1"
get_dofmap_pattern(FEType::Type{<:H1Q2}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry0D}) = "N1"
get_dofmap_pattern(FEType::Type{<:H1Q2}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry1D}) = "N1I1"
get_dofmap_pattern(FEType::Type{<:H1Q2}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry2D}) = "N1E1"
get_dofmap_pattern(FEType::Type{<:H1Q2}, ::Type{EdgeDofs}, EG::Type{<:AbstractElementGeometry1D}) = "N1I1"

isdefined(FEType::Type{<:H1Q2}, ::Type{<:AbstractElementGeometry1D}) = true
isdefined(FEType::Type{<:H1Q2}, ::Type{<:Triangle2D}) = true
isdefined(FEType::Type{<:H1Q2}, ::Type{<:Quadrilateral2D}) = true
isdefined(FEType::Type{<:H1Q2}, ::Type{<:Tetrahedron3D}) = true

interior_dofs_offset(::Type{<:AssemblyType}, ::Type{H1Q2{ncomponents, edim}}, ::Type{<:Edge1D}) where {ncomponents, edim} = 2
interior_dofs_offset(::Type{<:AssemblyType}, ::Type{H1Q2{ncomponents, edim}}, ::Type{<:Quadrilateral2D}) where {ncomponents, edim} = 8

init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{AT_NODES}) where {Tv, Ti, FEType <: H1Q2, APT} = NodalInterpolator(FES)
init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}) where {Tv, Ti, FEType <: H1Q2, APT} = MomentInterpolator(FES, ON_FACES)
init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_EDGES}) where {Tv, Ti, FEType <: H1Q2, APT} = MomentInterpolator(FES, ON_EDGES)
init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}) where {Tv, Ti, FEType <: H1Q2, APT} = MomentInterpolator(FES, ON_CELLS)

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{AT_NODES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1Q2, APT}
    return get_interpolator(FE, AT_NODES).evaluate!(Target, exact_function!, items; kwargs...)
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_EDGES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1Q2, APT}
    edim = get_edim(FEType)
    return if edim == 3
        # delegate edge nodes to node interpolation
        subitems = slice(FE.dofgrid[EdgeNodes], items)
        interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)

        # perform edge mean interpolation
        get_interpolator(FE, ON_EDGES).evaluate!(Target, exact_function!, items)
    end
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1Q2, APT}
    edim = get_edim(FEType)
    return if edim == 2
        # delegate face nodes to node interpolation
        subitems = slice(FE.dofgrid[FaceNodes], items)
        interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)

        # perform face mean interpolation
        get_interpolator(FE, ON_FACES).evaluate!(Target, exact_function!, items)
    elseif edim == 3
        # delegate face edges to edge interpolation
        subitems = slice(FE.dofgrid[FaceEdges], items)
        interpolate!(Target, FE, ON_EDGES, exact_function!; items = subitems, kwargs...)
    elseif edim == 1
        # delegate face nodes to node interpolation
        subitems = slice(FE.dofgrid[FaceNodes], items)
        interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)
    end
end


function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1Q2, APT}
    edim = get_edim(FEType)
    return if edim == 2
        # delegate cell faces to face interpolation
        subitems = slice(FE.dofgrid[CellFaces], items)
        interpolate!(Target, FE, ON_FACES, exact_function!; items = subitems, kwargs...)

        # preserve cell integral on quads
        if items == []
            items = 1:num_cells(FE.dofgrid)
        end
        items_quads = findall(i -> (i <: Quadrilateral2D), view(FE.dofgrid[CellGeometries], items))
        if length(items_quads) > 0
            get_interpolator(FE, ON_CELLS).evaluate!(Target, exact_function!, items)
        end

    elseif edim == 3
        # delegate cell edges to edge interpolation
        subitems = slice(FE.dofgrid[CellEdges], items)
        interpolate!(Target, FE, ON_EDGES, exact_function!; items = subitems, kwargs...)

        # todo: preserve face and cell integrals on quads
    elseif edim == 1
        # delegate cell nodes to node interpolation
        subitems = slice(FE.dofgrid[CellNodes], items)
        interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)

        # preserve cell integral
        get_interpolator(FE, ON_CELLS).evaluate!(Target, exact_function!, items)
    end
end


function get_basis(::Type{<:AssemblyType}, FEType::Type{H1Q2{ncomponents, edim}}, ::Type{<:Vertex0D}) where {ncomponents, edim}
    return function closure(refbasis, xref)
        for k in 1:ncomponents
            refbasis[k, k] = 1
        end
        return
    end
end

function get_basis(::Type{<:AssemblyType}, FEType::Type{H1Q2{ncomponents, edim}}, ::Type{<:Edge1D}) where {ncomponents, edim}
    return function closure(refbasis, xref)
        refbasis[end] = 1 - xref[1]
        for k in 1:ncomponents
            refbasis[3 * k - 2, k] = 2 * refbasis[end] * (refbasis[end] - 1 // 2)      # node 1
            refbasis[3 * k - 1, k] = 2 * xref[1] * (xref[1] - 1 // 2)                  # node 2
            refbasis[3 * k, k] = 4 * refbasis[end] * xref[1]                       # face 1
        end
        return
    end
end

function get_basis(::Type{<:AssemblyType}, FEType::Type{H1Q2{ncomponents, edim}}, ::Type{<:Triangle2D}) where {ncomponents, edim}
    return function closure(refbasis, xref)
        refbasis[end] = 1 - xref[1] - xref[2] # store last barycentric coordinate
        for k in 1:ncomponents
            refbasis[6 * k - 5, k] = 2 * refbasis[end] * (refbasis[end] - 1 // 2)      # node 1
            refbasis[6 * k - 4, k] = 2 * xref[1] * (xref[1] - 1 // 2)                  # node 2
            refbasis[6 * k - 3, k] = 2 * xref[2] * (xref[2] - 1 // 2)                  # node 3
            refbasis[6 * k - 2, k] = 4 * refbasis[end] * xref[1]                     # face 1
            refbasis[6 * k - 1, k] = 4 * xref[1] * xref[2]                           # face 2
            refbasis[6 * k, k] = 4 * xref[2] * refbasis[end]                       # face 3
        end
        return
    end
end


function get_basis(::Type{<:AssemblyType}, FEType::Type{H1Q2{ncomponents, edim}}, ::Type{<:Tetrahedron3D}) where {ncomponents, edim}
    return function closure(refbasis, xref)
        refbasis[end] = 1 - xref[1] - xref[2] - xref[3]
        for k in 1:ncomponents
            refbasis[10 * k - 9, k] = 2 * refbasis[end] * (refbasis[end] - 1 // 2)     # node 1
            refbasis[10 * k - 8, k] = 2 * xref[1] * (xref[1] - 1 // 2)                 # node 2
            refbasis[10 * k - 7, k] = 2 * xref[2] * (xref[2] - 1 // 2)                 # node 3
            refbasis[10 * k - 6, k] = 2 * xref[3] * (xref[3] - 1 // 2)                 # node 4
            refbasis[10 * k - 5, k] = 4 * refbasis[end] * xref[1]                    # edge 1
            refbasis[10 * k - 4, k] = 4 * refbasis[end] * xref[2]                    # edge 2
            refbasis[10 * k - 3, k] = 4 * refbasis[end] * xref[3]                    # edge 3
            refbasis[10 * k - 2, k] = 4 * xref[1] * xref[2]                          # edge 4
            refbasis[10 * k - 1, k] = 4 * xref[1] * xref[3]                          # edge 5
            refbasis[10 * k, k] = 4 * xref[2] * xref[3]                          # edge 6
        end
        return
    end
end


function get_basis(::Type{<:AssemblyType}, FEType::Type{H1Q2{ncomponents, edim}}, ::Type{<:Quadrilateral2D}) where {ncomponents, edim}
    return function closure(refbasis, xref)
        ## old Q2 functions without cell bubble
        #refbasis[1,1] = -2*(1 - xref[1])*(1 - xref[2])*(xref[1]+xref[2]-1//2);
        #refbasis[2,1] = -2*xref[1]*(1 - xref[2])*(xref[2]-xref[1]+1//2);
        #refbasis[3,1] = 2*xref[1]*xref[2]*(xref[1]+xref[2]-3//2);
        #refbasis[4,1] = -2*xref[2]*(1 - xref[1])*(xref[1]-xref[2]+1//2);
        #refbasis[5,1] = 4*xref[1]*(1 - xref[1])*(1 - xref[2])
        #refbasis[6,1] = 4*xref[2]*xref[1]*(1 - xref[2])
        #refbasis[7,1] = 4*xref[1]*xref[2]*(1 - xref[1])
        #refbasis[8,1] = 4*xref[2]*(1 - xref[1])*(1 - xref[2])

        # nodal basis functions
        refbasis[1, 1] = 4 * (xref[1] - 1 // 2) * (xref[1] - 1) * (xref[2] - 1 // 2) * (xref[2] - 1)
        refbasis[2, 1] = 4 * (xref[1] - 1 // 2) * (xref[1]) * (xref[2] - 1 // 2) * (xref[2] - 1)
        refbasis[3, 1] = 4 * (xref[1] - 1 // 2) * (xref[1]) * (xref[2] - 1 // 2) * xref[2]
        refbasis[4, 1] = 4 * (xref[1] - 1 // 2) * (xref[1] - 1) * (xref[2] - 1 // 2) * xref[2]
        # face basis functions
        refbasis[5, 1] = -8 * (xref[1]) * (xref[1] - 1) * (xref[2] - 1 // 2) * (xref[2] - 1)
        refbasis[6, 1] = -8 * (xref[1] - 1 // 2) * xref[1] * (xref[2]) * (xref[2] - 1)
        refbasis[7, 1] = -8 * (xref[1]) * (xref[1] - 1) * (xref[2] - 1 // 2) * xref[2]
        refbasis[8, 1] = -8 * (xref[1] - 1 // 2) * (xref[1] - 1) * (xref[2]) * (xref[2] - 1)
        # cell bubble
        refbasis[9, 1] = 16 * (1 - xref[1]) * (1 - xref[2]) * xref[1] * xref[2]
        for k in 2:ncomponents, j in 1:9
            refbasis[9 * (k - 1) + j, k] = refbasis[j, 1]
        end
        return
    end
end
