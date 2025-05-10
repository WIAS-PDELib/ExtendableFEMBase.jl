"""
````
abstract type H1P1{ncomponents} <: AbstractH1FiniteElement where {ncomponents<:Int}
````

Continuous piecewise first-order linear polynomials.

allowed ElementGeometries:
- Edge1D
- Triangle2D
- Tetrahedron3D
"""
abstract type H1P1{ncomponents} <: AbstractH1FiniteElement where {ncomponents <: Int} end
H1P1(ncomponents::Int) = H1P1{ncomponents}

function Base.show(io::Core.IO, ::Type{<:H1P1{ncomponents}}) where {ncomponents}
    return print(io, "H1P1{$ncomponents}")
end

get_ncomponents(FEType::Type{<:H1P1}) = FEType.parameters[1] # is this okay?
get_ndofs(::Type{<:AssemblyType}, FEType::Type{<:H1P1}, EG::Type{<:AbstractElementGeometry}) = num_nodes(EG) * FEType.parameters[1]

get_polynomialorder(::Type{<:H1P1}, ::Type{<:Edge1D}) = 1;
get_polynomialorder(::Type{<:H1P1}, ::Type{<:Triangle2D}) = 1;
get_polynomialorder(::Type{<:H1P1}, ::Type{<:Tetrahedron3D}) = 1;

get_dofmap_pattern(FEType::Type{<:H1P1}, ::Type{CellDofs}, EG::Type{<:AbstractElementGeometry}) = "N1"
get_dofmap_pattern(FEType::Type{<:H1P1}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry}) = "N1"
get_dofmap_pattern(FEType::Type{<:H1P1}, ::Union{Type{EdgeDofs}, Type{BEdgeDofs}}, EG::Type{<:AbstractElementGeometry}) = "N1"

isdefined(FEType::Type{<:H1P1}, ::Type{<:AbstractElementGeometry1D}) = true
isdefined(FEType::Type{<:H1P1}, ::Type{<:Triangle2D}) = true
isdefined(FEType::Type{<:H1P1}, ::Type{<:Tetrahedron3D}) = true

init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{AT_NODES}) where {Tv, Ti, FEType <: H1P1, APT} = NodalInterpolator(FES)

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{AT_NODES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1P1, APT}
    return get_interpolator(FE, AT_NODES).evaluate!(Target, exact_function!, items; kwargs...)
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_EDGES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1P1, APT}
    # delegate edge nodes to node interpolation
    subitems = slice(FE.dofgrid[EdgeNodes], items)
    return interpolaste!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1P1, APT}
    # delegate face nodes to node interpolation
    subitems = slice(FE.dofgrid[FaceNodes], items)
    return interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)
end

function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}, exact_function!; items = [], kwargs...) where {Tv, Ti, FEType <: H1P1, APT}
    return if FE.broken == true
        # broken interpolation
        point_evaluation_broken!(Target, FE, ON_CELLS, exact_function!; items = items, time = time)
    else
        # delegate cell nodes to node interpolation
        subitems = slice(FE.dofgrid[CellNodes], items)
        interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)
    end
end

function get_basis(::Type{<:AssemblyType}, FEType::Type{H1P1{ncomponents}}, ET::Type{<:Union{AbstractElementGeometry}}) where {ncomponents}
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
