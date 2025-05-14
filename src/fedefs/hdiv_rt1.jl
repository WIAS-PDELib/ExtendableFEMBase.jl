"""
````
abstract type HDIVRT1{edim} <: AbstractHdivFiniteElement where {edim<:Int}
````

Hdiv-conforming vector-valued (ncomponents = edim) Raviart-Thomas space of order 1.

allowed ElementGeometries:
- Triangle2D
- Tetrahedron3D
"""
abstract type HDIVRT1{edim} <: AbstractHdivFiniteElement where {edim <: Int} end
HDIVRT1(edim::Int) = HDIVRT1{edim}

function Base.show(io::Core.IO, FEType::Type{<:HDIVRT1{edim}}) where {edim}
    return print(io, "HDIVRT1{$edim}")
end

get_ncomponents(FEType::Type{<:HDIVRT1}) = FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:HDIVRT1}, EG::Type{<:AbstractElementGeometry1D}) = 2
get_ndofs(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:HDIVRT1}, EG::Type{<:Triangle2D}) = 3
get_ndofs(::Type{ON_CELLS}, FEType::Type{<:HDIVRT1}, EG::Type{<:Triangle2D}) = 2 * num_faces(EG) + 2
get_ndofs(::Type{ON_CELLS}, FEType::Type{<:HDIVRT1}, EG::Type{<:Tetrahedron3D}) = 3 * num_faces(EG) + 3
get_ndofs_all(::Type{ON_CELLS}, FEType::Type{<:HDIVRT1}, EG::Type{<:Tetrahedron3D}) = 4 * num_faces(EG) + 3 # in 3D only 3 of 4 face dofs are used depending on orientation

get_polynomialorder(::Type{<:HDIVRT1{2}}, ::Type{<:AbstractElementGeometry1D}) = 1;
get_polynomialorder(::Type{<:HDIVRT1{3}}, ::Type{<:AbstractElementGeometry2D}) = 1;
get_polynomialorder(::Type{<:HDIVRT1{2}}, ::Type{<:AbstractElementGeometry2D}) = 2;
get_polynomialorder(::Type{<:HDIVRT1{3}}, ::Type{<:AbstractElementGeometry3D}) = 2;

get_dofmap_pattern(FEType::Type{<:HDIVRT1{2}}, ::Type{CellDofs}, EG::Type{<:Triangle2D}) = "f2i2"
get_dofmap_pattern(FEType::Type{<:HDIVRT1{2}}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry1D}) = "i2"

get_dofmap_pattern(FEType::Type{<:HDIVRT1{3}}, ::Type{CellDofs}, EG::Type{<:Tetrahedron3D}) = "f3i3"
get_dofmap_pattern(FEType::Type{<:HDIVRT1{3}}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:Triangle2D}) = "i3"

isdefined(FEType::Type{<:HDIVRT1}, ::Type{<:Triangle2D}) = true
isdefined(FEType::Type{<:HDIVRT1}, ::Type{<:Tetrahedron3D}) = true

interior_dofs_offset(::Type{<:ON_CELLS}, ::Type{<:HDIVRT1{2}}, ::Type{<:Triangle2D}) = 6
interior_dofs_offset(::Type{<:ON_CELLS}, ::Type{<:HDIVRT1{3}}, ::Type{<:Tetrahedron3D}) = 12

function RT1_normalflux_eval!(dim)
    function closure(result, f, qpinfo)
        result[1] = dot(f, qpinfo.normal)
        result[2] = result[1] * (qpinfo.xref[1] - 1 // dim)
        if dim == 3
            result[3] = result[1] * (qpinfo.xref[2] - 1 // dim)
        end
    end
end
init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}) where {Tv, Ti, FEType <: HDIVRT1, APT} = FunctionalInterpolator(RT1_normalflux_eval!(get_ncomponents(FEType)), FES, ON_FACES; bonus_quadorder = 1)
init_interpolator!(FES::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}) where {Tv, Ti, FEType <: HDIVRT1, APT} = MomentInterpolator(FES, ON_CELLS)

function ExtendableGrids.interpolate!(Target::AbstractArray{T, 1}, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_FACES}, exact_function!; items = [], kwargs...) where {T, Tv, Ti, FEType <: HDIVRT1, APT}
    get_interpolator(FE, ON_FACES).evaluate!(Target, exact_function!, items; kwargs...)
end

function ExtendableGrids.interpolate!(Target::AbstractArray{T, 1}, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{ON_CELLS}, exact_function!; items = [], kwargs...) where {T, Tv, Ti, FEType <: HDIVRT1, APT}
    # delegate cell faces to face interpolation
    subitems = slice(FE.dofgrid[CellFaces], items)
    interpolate!(Target, FE, ON_FACES, exact_function!; items = subitems, kwargs...)

    # set values of interior RT1 functions such that P0 moments are preserved
    get_interpolator(FE, ON_CELLS).evaluate!(Target, exact_function!, items; kwargs...)

    return nothing
end


# only normalfluxes on faces
function get_basis(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, ::Type{<:HDIVRT1{2}}, ::Type{<:AbstractElementGeometry1D})
    return function closure(refbasis, xref)
        refbasis[1, 1] = 1                # normal-flux of RT0 function on single face
        return refbasis[2, 1] = 12 * (xref[1] - 1 // 2) # linear normal-flux of RT1 function
    end
end

# only normalfluxes on faces
function get_basis(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, ::Type{<:HDIVRT1{3}}, ::Type{<:Triangle2D})
    return function closure(refbasis, xref)
        refbasis[1, 1] = 1                # normal-flux of RT0 function on single face
        refbasis[2, 1] = 12 * (2 * xref[1] + xref[2] - 1) # 1st linear normal-flux RT1 function (normal flux weighted with (phi_1 - 1/3))
        return refbasis[3, 1] = 12 * (2 * xref[2] + xref[1] - 1) # 2nd linear normal-flux RT1 function (normal flux weighted with (phi_2 - 1/3))
    end
end


function get_basis(::Type{ON_CELLS}, ::Type{HDIVRT1{2}}, ::Type{<:Triangle2D})
    return function closure(refbasis, xref)
        # RT0 basis
        refbasis[1, 1] = xref[1]
        refbasis[1, 2] = xref[2] - 1
        refbasis[3, 1] = xref[1]
        refbasis[3, 2] = xref[2]
        refbasis[5, 1] = xref[1] - 1
        refbasis[5, 2] = xref[2]

        for k in 1:2
            # additional RT1 face basis functions
            refbasis[2, k] = -12 * (1 // 2 - xref[1] - xref[2]) * refbasis[1, k]
            refbasis[4, k] = -(12 * (xref[1] - 1 // 2)) * refbasis[3, k]
            refbasis[6, k] = -(12 * (xref[2] - 1 // 2)) * refbasis[5, k]
            # interior functions
            refbasis[7, k] = 12 * xref[2] * refbasis[1, k]
            refbasis[8, k] = 12 * xref[1] * refbasis[5, k]
        end
        return
    end
end

function get_basis(::Type{ON_CELLS}, ::Type{HDIVRT1{3}}, ::Type{<:Tetrahedron3D})
    return function closure(refbasis, xref)
        refbasis[end] = 1 - xref[1] - xref[2] - xref[3]
        # RT0 basis
        refbasis[1, 1] = 2 * xref[1]
        refbasis[1, 2] = 2 * xref[2]
        refbasis[1, 3] = 2 * (xref[3] - 1)
        refbasis[5, 1] = 2 * xref[1]
        refbasis[5, 2] = 2 * (xref[2] - 1)
        refbasis[5, 3] = 2 * xref[3]
        refbasis[9, 1] = 2 * xref[1]
        refbasis[9, 2] = 2 * xref[2]
        refbasis[9, 3] = 2 * xref[3]
        refbasis[13, 1] = 2 * (xref[1] - 1)
        refbasis[13, 2] = 2 * xref[2]
        refbasis[13, 3] = 2 * xref[3]

        for k in 1:3
            # additional RT1 face basis functions (2 per face)          # Test with (phi_1-1/3,phi_3-1/3,phi_2-1/3)
            refbasis[2, k] = -12 * (2 * refbasis[end] + xref[2] - 1) * refbasis[1, k]     # [1,0,-1]
            refbasis[3, k] = -12 * (2 * xref[2] + xref[1] - 1) * refbasis[1, k]  # [-1,1,0]
            refbasis[4, k] = 12 * (2 * xref[2] + refbasis[end] - 1) * refbasis[1, k]      # [0,-1,1]

            refbasis[6, k] = -12 * (2 * refbasis[end] + xref[1] - 1) * refbasis[5, k]
            refbasis[7, k] = -12 * (2 * xref[1] + xref[3] - 1) * refbasis[5, k]
            refbasis[8, k] = 12 * (2 * xref[1] + refbasis[end] - 1) * refbasis[5, k]

            refbasis[10, k] = -12 * (2 * xref[1] + xref[2] - 1) * refbasis[9, k]
            refbasis[11, k] = -12 * (2 * xref[2] + xref[3] - 1) * refbasis[9, k]
            refbasis[12, k] = 12 * (2 * xref[2] + xref[1] - 1) * refbasis[9, k]

            refbasis[14, k] = -12 * (2 * refbasis[end] + xref[3] - 1) * refbasis[13, k]
            refbasis[15, k] = -12 * (2 * xref[3] + xref[2] - 1) * refbasis[13, k]
            refbasis[16, k] = 12 * (2 * xref[3] + refbasis[end] - 1) * refbasis[13, k]

            # interior functions
            refbasis[17, k] = 12 * xref[3] * refbasis[1, k]
            refbasis[18, k] = 12 * xref[2] * refbasis[5, k]
            refbasis[19, k] = 12 * xref[1] * refbasis[13, k]
        end
        return
    end
end

function get_coefficients(::Type{ON_CELLS}, FE::FESpace{Tv, Ti, <:HDIVRT1{2}, APT}, EG::Type{<:Triangle2D}, xgrid) where {Tv, Ti, APT}
    xCellFaceSigns = xgrid[CellFaceSigns]
    nfaces = num_faces(EG)
    return function closure(coefficients, cell)
        fill!(coefficients, 1.0)
        # multiplication with normal vector signs (only RT0)
        for j in 1:nfaces, k in 1:size(coefficients, 1)
            coefficients[k, 2 * j - 1] = xCellFaceSigns[j, cell]
        end
        return nothing
    end
end

function get_coefficients(::Type{ON_CELLS}, FE::FESpace{Tv, Ti, <:HDIVRT1{3}, APT}, EG::Type{<:Tetrahedron3D}, xgrid) where {Tv, Ti, APT}
    xCellFaceSigns = xgrid[CellFaceSigns]
    nfaces = num_faces(EG)
    return function closure(coefficients, cell)
        fill!(coefficients, 1.0)
        # multiplication with normal vector signs (only RT0)
        for j in 1:nfaces, k in 1:size(coefficients, 1)
            coefficients[k, 3 * j - 2] = xCellFaceSigns[j, cell] # RT0
            coefficients[k, 3 * j - 1] = -1
            coefficients[k, 3 * j] = 1
        end
        # @show coefficients
        return nothing
    end
end


# subset selector ensures that for every cell face
# the RT0 and those two BDM1 face functions are chosen
# such that they reflect the two moments with respect to the second and third node
# of the global face enumeration
function get_basissubset(::Type{ON_CELLS}, FE::FESpace{Tv, Ti, <:HDIVRT1{3}, APT}, EG::Type{<:Tetrahedron3D}, xgrid) where {Tv, Ti, APT}
    xCellFaceOrientations = xgrid[CellFaceOrientations]
    nfaces::Int = num_faces(EG)
    orientation = xCellFaceOrientations[1, 1]
    shift4orientation1::Array{Int, 1} = [1, 0, 1, 2]
    shift4orientation2::Array{Int, 1} = [2, 2, 0, 1]
    return function closure(subset_ids, cell)
        for j in 1:nfaces
            subset_ids[3 * j - 2] = 4 * j - 3 # always take the RT0 function
            orientation = xCellFaceOrientations[j, cell]
            subset_ids[3 * j - 1] = 4 * j - shift4orientation1[orientation]
            subset_ids[3 * j] = 4 * j - shift4orientation2[orientation]
        end
        for j in 1:3
            subset_ids[12 + j] = 16 + j # interior functions
        end
        # @show subset_ids
        return nothing
    end
end
