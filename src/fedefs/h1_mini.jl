"""
````
abstract type H1MINI{ncomponents,edim} <: AbstractH1FiniteElement where {ncomponents<:Int,edim<:Int}
````

Mini finite element.

allowed element geometries:
- Triangle2D (linear polynomials + cubic cell bubble)
- Quadrilateral2D (Q1 space + quartic cell bubble)
- Tetrahedron3D (linear polynomials + cubic cell bubble)
"""
abstract type H1MINI{ncomponents, edim} <: AbstractH1FiniteElement where {ncomponents <: Int, edim <: Int} end
H1MINI(ncomponents::Int, edim = ncomponents) = H1MINI{ncomponents, edim}

function Base.show(io::Core.IO, ::Type{<:H1MINI{ncomponents, edim}}) where {ncomponents, edim}
    return print(io, "H1MINI{$ncomponents,$edim}")
end

get_ncomponents(FEType::Type{<:H1MINI}) = FEType.parameters[1]
get_edim(FEType::Type{<:H1MINI}) = FEType.parameters[2]
get_ndofs(::Union{Type{<:ON_FACES}, Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:H1MINI}, EG::Type{<:AbstractElementGeometry0D}) = FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:H1MINI}, EG::Type{<:AbstractElementGeometry}) = num_nodes(EG) * FEType.parameters[1]
get_ndofs(::Type{<:ON_CELLS}, FEType::Type{<:H1MINI}, EG::Type{<:AbstractElementGeometry0D}) = FEType.parameters[1]
get_ndofs(::Type{<:ON_CELLS}, FEType::Type{<:H1MINI}, EG::Type{<:AbstractElementGeometry}) = (1 + num_nodes(EG)) * FEType.parameters[1]

get_polynomialorder(FEType::Type{<:H1MINI}, ::Type{<:Edge1D}) = FEType.parameters[2] == 1 ? 2 : 1
get_polynomialorder(FEType::Type{<:H1MINI}, ::Type{<:Triangle2D}) = FEType.parameters[2] == 2 ? 3 : 1;
get_polynomialorder(FEType::Type{<:H1MINI}, ::Type{<:Quadrilateral2D}) = FEType.parameters[2] == 2 ? 4 : 2;
get_polynomialorder(FEType::Type{<:H1MINI}, ::Type{<:Tetrahedron3D}) = 4;

get_dofmap_pattern(FEType::Type{<:H1MINI}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry}) = "N1I1"
get_dofmap_pattern(FEType::Type{<:H1MINI}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry}) = "N1"

isdefined(FEType::Type{<:H1MINI}, ::Type{<:Triangle2D}) = true
isdefined(FEType::Type{<:H1MINI}, ::Type{<:Quadrilateral2D}) = true
isdefined(FEType::Type{<:H1MINI}, ::Type{<:Tetrahedron3D}) = true

get_ref_cellmoments(::Type{<:H1MINI}, ::Type{<:Triangle2D}) = [1 // 3, 1 // 3, 1 // 3, 1 // 1] # integrals of 1D basis functions over reference cell (divided by volume)
get_ref_cellmoments(::Type{<:H1MINI}, ::Type{<:Tetrahedron3D}) = [1 // 4, 1 // 4, 1 // 4, 1 // 4, 1 // 1] # integrals of 1D basis functions over reference cell (divided by volume)
get_ref_cellmoments(::Type{<:H1MINI}, ::Type{<:Quadrilateral2D}) = [1 // 4, 1 // 4, 1 // 4, 1 // 4, 1 // 1] # integrals of 1D basis functions over reference cell (divided by volume)

interior_dofs_offset(::Type{ON_CELLS}, ::Type{H1MINI{ncomponents, edim}}, EG::Type{<:AbstractElementGeometry}) where {ncomponents, edim} = num_nodes(EG)

init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{AT_NODES}) where {Tv, Ti, FEType <: H1MINI, APT} = NodalInterpolator(FES)
init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}) where {Tv, Ti, FEType <: H1MINI, APT} = MomentInterpolator(FES, ON_CELLS)

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{AT_NODES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1MINI, APT}
    return get_interpolator(FE, AT_NODES).evaluate!(Target, exact_function!, items; kwargs...)
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_EDGES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1MINI, APT}
    # delegate edge nodes to node interpolation
    subitems = slice(FE.dofgrid[EdgeNodes], items)
    return interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1MINI, APT}
    # delegate face nodes to node interpolation
    subitems = slice(FE.dofgrid[FaceNodes], items)
    return interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1MINI, APT}
    # delegate cell nodes to node interpolation
    subitems = slice(FE.dofgrid[CellNodes], items)
    interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)

    # fix cell bubble value by preserving integral mean
    return get_interpolator(FE, ON_CELLS).evaluate!(Target, exact_function!, items; kwargs...)
end

function nodevalues!(Target::AbstractArray{<:Real, 2}, Source::AbstractArray{<:Real, 1}, FE::FESpace{<:H1MINI})
    nnodes = num_sources(FE.dofgrid[Coordinates])
    ncells = num_sources(FE.dofgrid[CellNodes])
    FEType = eltype(FE)
    ncomponents = get_ncomponents(FEType)
    offset4component = 0:(nnodes + ncells):(ncomponents * (nnodes + ncells))
    for node in 1:nnodes
        for c in 1:ncomponents
            Target[c, node] = Source[offset4component[c] + node]
        end
    end
    return
end

function get_basis(AT::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, ::Type{H1MINI{ncomponents, edim}}, EG::Type{<:AbstractElementGeometry}) where {ncomponents, edim}
    # on faces same as P1
    return get_basis(AT, H1P1{ncomponents}, EG)
end

function get_basis(AT::Type{<:ON_CELLS}, ::Type{H1MINI{ncomponents, edim}}, EG::Type{<:Triangle2D}) where {ncomponents, edim}
    refbasis_P1 = get_basis(AT, H1P1{1}, EG)
    offset = get_ndofs(AT, H1P1{1}, EG) + 1
    return function closure(refbasis, xref)
        refbasis_P1(refbasis, xref)
        # add cell bubbles to P1 basis (scaled to have unit integral)
        refbasis[offset, 1] = 60 * (1 - xref[1] - xref[2]) * xref[1] * xref[2]
        for k in 1:(ncomponents - 1), j in 1:offset
            refbasis[k * offset + j, k + 1] = refbasis[j, 1]
        end
        return
    end
end

function get_basis(AT::Type{<:ON_CELLS}, ::Type{H1MINI{ncomponents, edim}}, EG::Type{<:Quadrilateral2D}) where {ncomponents, edim}
    refbasis_P1 = get_basis(AT, H1Q1{1}, EG)
    offset = get_ndofs(AT, H1Q1{1}, EG) + 1
    return function closure(refbasis, xref)
        refbasis_P1(refbasis, xref)
        # add cell bubbles to P1 basis (scaled to have unit integral)
        refbasis[offset, 1] = 36 * (1 - xref[1]) * (1 - xref[2]) * xref[1] * xref[2]
        for k in 1:(ncomponents - 1), j in 1:offset
            refbasis[k * offset + j, k + 1] = refbasis[j, 1]
        end
        return
    end
end

function get_basis(AT::Type{<:ON_CELLS}, ::Type{H1MINI{ncomponents, edim}}, EG::Type{<:Tetrahedron3D}) where {ncomponents, edim}
    refbasis_P1 = get_basis(AT, H1P1{1}, EG)
    offset = get_ndofs(AT, H1P1{1}, EG) + 1
    return function closure(refbasis, xref)
        refbasis_P1(refbasis, xref)
        # add cell bubbles to P1 basis (scaled to have unit integral)
        refbasis[offset, 1] = 840 * (1 - xref[1] - xref[2] - xref[3]) * xref[1] * xref[2] * xref[3]
        for k in 1:(ncomponents - 1), j in 1:offset
            refbasis[k * offset + j, k + 1] = refbasis[j, 1]
        end
        return
    end
end
