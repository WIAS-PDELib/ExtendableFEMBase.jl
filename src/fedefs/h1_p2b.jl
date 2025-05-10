"""
````
abstract type H1P2B{ncomponents,edim} <: AbstractH1FiniteElement where {ncomponents<:Int,edim<:Int}
````

Continuous piecewise second-order polynomials.

allowed ElementGeometries:
- Triangle2D
"""
abstract type H1P2B{ncomponents, edim} <: AbstractH1FiniteElement where {ncomponents <: Int, edim <: Int} end
H1P2B(ncomponents::Int, edim = ncomponents) = H1P2B{ncomponents, edim}

function Base.show(io::Core.IO, ::Type{<:H1P2B{ncomponents, edim}}) where {ncomponents, edim}
    return print(io, "H1P2B{$ncomponents,$edim}")
end

get_ncomponents(FEType::Type{<:H1P2B}) = FEType.parameters[1]
get_edim(FEType::Type{<:H1P2B}) = FEType.parameters[2]

get_ndofs(::Type{ON_CELLS}, FEType::Type{<:H1P2B}, EG::Type{<:Triangle2D}) = 7 * FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:H1P2B}, EG::Type{<:AbstractElementGeometry1D}) = 3 * FEType.parameters[1]

get_polynomialorder(::Type{<:H1P2B}, ::Type{<:Edge1D}) = 2;
get_polynomialorder(::Type{<:H1P2B}, ::Type{<:Triangle2D}) = 3;
get_polynomialorder(::Type{<:H1P2B}, ::Type{<:Tetrahedron3D}) = 4;

get_dofmap_pattern(FEType::Type{<:H1P2B}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry2D}) = "N1F1I1"
get_dofmap_pattern(FEType::Type{<:H1P2B}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry1D}) = "N1I1"

isdefined(FEType::Type{<:H1P2B}, ::Type{<:Triangle2D}) = true

interior_dofs_offset(::Type{<:ON_FACES}, ::Type{H1P2B{ncomponents, edim}}, ::Type{Edge1D}) where {ncomponents, edim} = 2
interior_dofs_offset(::Type{<:ON_CELLS}, ::Type{H1P2B{ncomponents, edim}}, ::Type{Triangle2D}) where {ncomponents, edim} = 6

get_ref_cellmoments(::Type{<:H1P2B}, ::Type{<:Triangle2D}) = [0 // 1, 0 // 1, 0 // 1, 1 // 3, 1 // 3, 1 // 3, 1 // 1] # integrals of 1D basis functions over reference cell (divided by volume)

init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{AT_NODES}) where {Tv, Ti, FEType <: H1P2B, APT} = NodalInterpolator(FES)
init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}) where {Tv, Ti, FEType <: H1P2B, APT} = MomentInterpolator(FES, ON_FACES)
init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}) where {Tv, Ti, FEType <: H1P2B, APT} = MomentInterpolator(FES, ON_CELLS)

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{AT_NODES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1P2B, APT}
    return get_interpolator(FE, AT_NODES).evaluate!(Target, exact_function!, items; kwargs...)
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_EDGES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1P2B, APT}
    edim = get_edim(FEType)
    return if edim == 3
        # delegate edge nodes to node interpolation
        subitems = slice(FE.dofgrid[EdgeNodes], items)
        interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)

        # perform edge mean interpolation
        get_interpolator(FE, ON_EDGES).evaluate!(Target, exact_function!, items)
    end
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1P2B, APT}
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

        # perform face mean interpolation
        # todo
    elseif edim == 1
        # delegate face nodes to node interpolation
        subitems = slice(FE.dofgrid[FaceNodes], items)
        interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)
    end
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1P2B, APT}
    edim = get_edim(FEType)
    if edim == 2
        # delegate cell faces to face interpolation
        subitems = slice(FE.dofgrid[CellFaces], items)
        interpolate!(Target, FE, ON_FACES, exact_function!; items = subitems, kwargs...)
    elseif edim == 3
        # delegate cell edges to edge interpolation
        subitems = slice(FE.dofgrid[CellEdges], items)
        interpolate!(Target, FE, ON_EDGES, exact_function!; items = subitems, kwargs...)
    elseif edim == 1
        # delegate cell nodes to node interpolation
        subitems = slice(FE.dofgrid[CellNodes], items)
        interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)
    end

    # fix cell bubble value by preserving integral mean
    get_interpolator(FE, ON_CELLS).evaluate!(Target, exact_function!, items; kwargs...)
end

function get_basis(AT::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, ::Type{H1P2B{ncomponents, edim}}, EG::Type{<:AbstractElementGeometry}) where {ncomponents, edim}
    # on faces same as P2
    return get_basis(AT, H1P2{ncomponents, edim}, EG)
end

function get_basis(AT::Type{ON_CELLS}, ::Type{H1P2B{ncomponents, edim}}, EG::Type{<:Triangle2D}) where {ncomponents, edim}
    refbasis_P2 = get_basis(AT, H1P2{1, edim}, EG)
    offset = get_ndofs(AT, H1P2{1, edim}, EG) + 1
    return function closure(refbasis, xref)
        refbasis_P2(refbasis, xref)
        # add cell bubbles to P2 basis
        refbasis[offset, 1] = 60 * (1 - xref[1] - xref[2]) * xref[1] * xref[2]
        for k in 1:(ncomponents - 1), j in 1:offset
            refbasis[k * offset + j, k + 1] = refbasis[j, 1]
        end
        return
    end
end
