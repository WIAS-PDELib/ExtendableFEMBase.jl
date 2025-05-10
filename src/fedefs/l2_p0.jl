"""
````
abstract type L2P0{ncomponents} <: AbstractH1FiniteElement where {ncomponents<:Int}
````

Piecewise constant polynomials on cells.

allowed ElementGeometries:
- any
"""
abstract type L2P0{ncomponents} <: AbstractH1FiniteElement where {ncomponents <: Int} end
L2P0(ncomponents::Int) = L2P0{ncomponents}

function Base.show(io::Core.IO, ::Type{<:L2P0{ncomponents}}) where {ncomponents}
    return print(io, "L2P0{$ncomponents}")
end

get_ncomponents(FEType::Type{<:L2P0}) = FEType.parameters[1]
get_ndofs(::Type{<:AssemblyType}, FEType::Type{<:L2P0}, EG::Type{<:AbstractElementGeometry}) = FEType.parameters[1]

get_polynomialorder(::Type{<:L2P0}, ::Type{<:AbstractElementGeometry}) = 0

get_dofmap_pattern(FEType::Type{<:L2P0}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry}) = "I1"
get_dofmap_pattern(FEType::Type{<:L2P0}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry}) = "C1"

isdefined(FEType::Type{<:L2P0}, ::Type{<:AbstractElementGeometry}) = true

interior_dofs_offset(::Type{ON_CELLS}, ::Type{<:L2P0}, EG::Type{<:AbstractElementGeometry}) = 0

init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}) where {Tv, Ti, FEType <: L2P0, APT} = MomentInterpolator(FES, ON_CELLS)

function ExtendableGrids.interpolate!(Target::AbstractArray{T, 1}, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}, exact_function!; items = [], kwargs...) where {T, Tv, Ti, FEType <: L2P0, APT}
    get_interpolator(FE, ON_CELLS).evaluate!(Target, exact_function!, items; kwargs...)
    return
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: L2P0, APT}
    # delegate to node cell interpolation
    subitems = slice(FE.dofgrid[FaceCells], items)
    return interpolate!(Target, FE, ON_CELLS, exact_function!; items = subitems, kwargs...)
end

function nodevalues!(Target::AbstractArray{<:Real, 2}, Source::AbstractArray{<:Real, 1}, FE::FESpace{<:L2P0})
    xCoords = FE.dofgrid[Coordinates]
    xCellNodes = FE.dofgrid[CellNodes]
    xNodeCells = atranspose(xCellNodes)
    FEType = eltype(FE)
    ncomponents = get_ncomponents(FEType)
    value = 0.0
    nneighbours = 0
    for node in 1:num_sources(xCoords)
        for c in 1:ncomponents
            value = 0.0
            nneighbours = num_targets(xNodeCells, node)
            for n in 1:nneighbours
                value += Source[(xNodeCells[n, node] - 1) * ncomponents + c]
            end
            value /= nneighbours
            Target[c, node] = value
        end
    end
    return
end

function get_basis(::Type{<:AssemblyType}, FEType::Type{L2P0{ncomponents}}, ::Type{<:AbstractElementGeometry}) where {ncomponents}
    return function closure(refbasis, xref)
        for k in 1:ncomponents
            refbasis[k, k] = 1.0
        end
        return
    end
end
