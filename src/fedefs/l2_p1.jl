"""
````
abstract type L2P1{ncomponents} <: AbstractH1FiniteElement where {ncomponents<:Int}
````

Discontinuous piecewise first-order linear polynomials (same as H1P1 but enforces broken = true).

allowed ElementGeometries:
- Edge1D
- Triangle2D
- Tetrahedron3D
"""
abstract type L2P1{ncomponents} <: AbstractH1FiniteElement where {ncomponents <: Int} end
L2P1(ncomponents::Int) = L2P1{ncomponents}


function Base.show(io::Core.IO, ::Type{<:L2P1{ncomponents}}) where {ncomponents}
    return print(io, "L2P1{$ncomponents}")
end

get_ncomponents(FEType::Type{<:L2P1}) = FEType.parameters[1] # is this okay?
get_ndofs(::Type{<:AssemblyType}, FEType::Type{<:L2P1}, EG::Type{<:AbstractElementGeometry}) = (dim_element(EG) + 1) * FEType.parameters[1]

get_polynomialorder(::Type{<:L2P1}, ::Type{<:Edge1D}) = 1;
get_polynomialorder(::Type{<:L2P1}, ::Type{<:Triangle2D}) = 1;
get_polynomialorder(::Type{<:L2P1}, ::Type{<:Tetrahedron3D}) = 1;
get_polynomialorder(::Type{<:L2P1}, ::Type{<:Quadrilateral2D}) = 1;
get_polynomialorder(::Type{<:L2P1}, ::Type{<:Hexahedron3D}) = 1;

get_dofmap_pattern(FEType::Type{<:L2P1}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry1D}) = "I2"
get_dofmap_pattern(FEType::Type{<:L2P1}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry2D}) = "I3"
get_dofmap_pattern(FEType::Type{<:L2P1}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry3D}) = "I4"

isdefined(FEType::Type{<:L2P1}, ::Type{<:AbstractElementGeometry1D}) = true
isdefined(FEType::Type{<:L2P1}, ::Type{<:Triangle2D}) = true
isdefined(FEType::Type{<:L2P1}, ::Type{<:Tetrahedron3D}) = true

init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{AT_NODES}) where {Tv, Ti, FEType <: L2P1, APT} = NodalInterpolator(FES)
function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{AT_NODES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: L2P1, APT}
    return get_interpolator(FE, AT_NODES).evaluate!(Target, exact_function!, items; kwargs...)
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_EDGES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: L2P1, APT}
    # delegate edge nodes to node interpolation
    subitems = slice(FE.dofgrid[EdgeNodes], items)
    return interpolaste!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: L2P1, APT}
    # delegate face nodes to node interpolation
    subitems = slice(FE.dofgrid[FaceNodes], items)
    return interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: L2P1, APT}
    return if FE.broken == true
        # broken interpolation
        get_interpolator(FE, AT_NODES).evaluate!(Target, exact_function!, items; kwargs...)
    else
        # delegate cell nodes to node interpolation
        subitems = slice(FE.dofgrid[CellNodes], items)
        interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)
    end
end

function get_basis(::Type{<:AssemblyType}, FEType::Type{L2P1{ncomponents}}, ET::Type{<:Union{AbstractElementGeometry}}) where {ncomponents}
    edim::Int = dim_element(ET)
    return function closure(refbasis, xref)
        for k in 1:ncomponents
            refbasis[(edim + 1) * k - edim, k] = 1
            for j in 1:edim
                refbasis[(edim + 1) * k - edim, k] -= xref[j]
                refbasis[(edim + 1) * k - edim + j, k] = xref[j]
            end
        end
        return
    end
end
