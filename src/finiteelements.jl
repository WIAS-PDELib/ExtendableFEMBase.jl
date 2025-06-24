#######################################################################################################
#######################################################################################################
### FFFFF II NN    N II TTTTTT EEEEEE     EEEEEE LL     EEEEEE M     M EEEEEE NN    N TTTTTT SSSSSS ###
### FF    II N N   N II   TT   EE         EE     LL     EE     MM   MM EE     N N   N   TT   SS     ###
### FFFF  II N  N  N II   TT   EEEEE      EEEEE  LL     EEEEE  M M M M EEEEE  N  N  N   TT    SSSS  ###
### FF    II N   N N II   TT   EE         EE     LL     EE     M  M  M EE     N   N N   TT       SS ###
### FF    II N    NN II   TT   EEEEEE     EEEEEE LLLLLL EEEEEE M     M EEEEEE N    NN   TT   SSSSSS ###
#######################################################################################################
#######################################################################################################

"""
$(TYPEDEF)

root type for finite element types
"""
abstract type AbstractFiniteElement end


#############################
# Finite Element Supertypes #
#############################
#
# they are used to steer the kind of local2global transformation
# below subtypes are defined that define basis functions on reference geometries
# and some other information like polyonomial degrees etc.

"""
$(TYPEDEF)

root type for Hdiv-conforming finite element types
"""
abstract type AbstractHdivFiniteElement <: AbstractFiniteElement end

"""
$(TYPEDEF)

root type for H1-conforming finite element types
"""
abstract type AbstractH1FiniteElement <: AbstractFiniteElement end

"""
$(TYPEDEF)

root type for H1-conforming finite element types with additional coefficients
"""
abstract type AbstractH1FiniteElementWithCoefficients <: AbstractH1FiniteElement end

"""
$(TYPEDEF)

root type for Hcurl-conforming finite element types with additional coefficients
"""
abstract type AbstractHcurlFiniteElement <: AbstractFiniteElement end

abstract type AbstractInterpolationOperator end

"""
````
struct FESpace{Tv, Ti, FEType<:AbstractFiniteElement,AT<:AssemblyType}
	name::String                          
	broken::Bool                         
	ndofs::Int                            
	coffset::Int                          
	xgrid::ExtendableGrid[Tv,Ti}           
	dofgrid::ExtendableGrid{Tv,Ti}	      
	dofmaps::Dict{Type{<:AbstractGridComponent},Any} 
    interpolators::Dict{Type{<:AssemblyType}, Any} 
end
````
A finite element space representing the global collection of degrees of freedom (DOFs) for a given finite element type and assembly type on a computational grid.

`FESpace` encapsulates the mapping between mesh entities (cells, faces, edges, etc.) and DOFs, as well as metadata and auxiliary structures needed for assembly, interpolation, and evaluation of finite element functions.

# Type Parameters
- `Tv`: Value type for grid coordinates (e.g., `Float64`).
- `Ti`: Integer type for grid indices (e.g., `Int64`).
- `FEType`: The finite element type (e.g., `H1P1`).
- `AT`: The assembly type (e.g., `ON_CELLS`).

# Fields
- `name::String`: Name of the finite element space.
- `broken::Bool`: Whether the space is "broken" (discontinuous across elements).
- `ndofs::Int64`: Total number of global degrees of freedom.
- `coffset::Int`: Offset for component DOFs (for vector-valued or mixed spaces).
- `xgrid::ExtendableGrid{Tv, Ti}`: Reference to the computational grid.
- `dofgrid::ExtendableGrid{Tv, Ti}`: Grid used for DOF numbering (may be a subgrid of `xgrid`).
- `dofmaps::Dict{Type{<:AbstractGridComponent}, Any}`: Dictionary of DOF maps for different grid components (cells, faces, edges, etc.).
- `interpolators::Dict{Type{<:AssemblyType}, Any}`: Dictionary of interpolation operators for different assembly types.

"""
struct FESpace{Tv, Ti, FEType <: AbstractFiniteElement, AT <: AssemblyType}
    name::String                          # full name of finite element space (used in messages)
    broken::Bool                          # if true, broken dofmaps are generated
    ndofs::Int64                          # total number of dofs
    coffset::Int                          # offset for component dofs
    xgrid::ExtendableGrid{Tv, Ti}         # link to (master/parent) grid
    dofgrid::ExtendableGrid{Tv, Ti}          # link to (sub) grid used for dof numbering (expected to be equal to or child grid of xgrid)
    dofmaps::Dict{Type{<:AbstractGridComponent}, Any} # backpack with dofmaps
    interpolators::Dict{Type{<:AssemblyType}, Any} # backpack with interpolators
end

function Base.copy(FES::FESpace{Tv, Ti, FEType, AT}) where {Tv, Ti, FEType, AT}
    return FESpace{Tv, Ti, FEType, AT}(deepcopy(FES.name), FES.broken, FES.ndofs, FES.coffset, FES.xgrid, FES.dofgrid, FES.dofmaps, FES.interpolators)
end

"""
$(TYPEDSIGNATURES)

returns the name of the finite element space.
"""
name(FES::FESpace) = FES.name

"""
$(TYPEDSIGNATURES)

returns the computational grid of the finite element space.
"""
xgrid(FES::FESpace) = FES.xgrid

"""
$(TYPEDSIGNATURES)

returns the dofgrid of the finite element space.
"""
dofgrid(FES::FESpace) = FES.dofgrid

"""
$(TYPEDSIGNATURES)

returns the total number of degrees of freedom of the finite element space.
"""
ndofs(FES::FESpace) = FES.ndofs

"""
$(TYPEDSIGNATURES)

returns the offset between the degrees of freedom of each component
(i.e. the number of scalar degrees of freedom that influence a component,
vector-valued degrees of freedom are stored at the end).
"""
coffset(FES::FESpace) = FES.coffset

"""
$(TYPEDSIGNATURES)

returns true if the finite element space is broken, false if not
"""
broken(FES::FESpace) = FES.broken

"""
$(TYPEDSIGNATURES)

returns the support of the finite element space
"""
get_AT(::FESpace{Tv, Ti, FEType, AT}) where {Tv, Ti, FEType, AT} = AT

"""
$(TYPEDSIGNATURES)

returns the finite element type of the finite element space
"""
get_FEType(::FESpace{Tv, Ti, FEType, AT}) where {Tv, Ti, FEType, AT} = FEType

include("dofmaps.jl")

"""
$(TYPEDSIGNATURES)
Set new dofmap
"""
Base.setindex!(FES::FESpace, v, DM::Type{<:DofMap}) = FES.dofmaps[DM] = v

"""
````
FESpace{FEType, AT}(xgrid::ExtendableGrid{Tv, Ti}; name = "", regions = nothing, broken::Bool = false)
````

Constructs a finite element space (`FESpace`) of the given finite element type (`FEType`) and assembly type (`AT`) on the provided computational grid.

- The `FESpace` represents the global collection of degrees of freedom (DOFs) for the specified finite element type and assembly type, and manages the mapping between mesh entities (cells, faces, edges, etc.) and DOFs.
- The `broken` switch allows the creation of a "broken" (discontinuous) finite element space, where basis functions are not required to be continuous across element boundaries.
- The `regions` argument can be used to restrict the space to a subset of the grid.
- If no `AT` is provided, the space is generated on cells (`ON_CELLS`).

# Arguments
- `xgrid::ExtendableGrid{Tv, Ti}`: The computational grid on which the finite element space is defined.

# Keyword Arguments
- `name::String: Name for the finite element space (for identification and debugging) (default: "", triggers automatic generation from FEType and broken).
- `regions`: Optional subset of the grid to restrict the space (default: all regions).
- `broken`: Whether to create a broken (discontinuous) space (default: false).

"""
function FESpace{FEType, AT}(
        xgrid::ExtendableGrid{Tv, Ti};
        name = "",
        regions = nothing,
        broken::Bool = false
    ) where {Tv, Ti, FEType <: AbstractFiniteElement, AT <: AssemblyType}

    # piecewise constants are always broken
    if FEType <: L2P0 || FEType <: L2P1
        broken = true
    end

    if AT == ON_FACES
        if isnothing(regions)
            regions = unique(xgrid[FaceRegions])
        end
        dofgrid = subgrid(xgrid, regions; support = ON_FACES, project = false)
    elseif AT == ON_BFACES
        if isnothing(regions)
            regions = unique(xgrid[BFaceRegions])
        end
        dofgrid = subgrid(xgrid, regions; support = ON_BFACES, project = false)
    elseif AT == ON_EDGES
        @assert false "not possible currently"
    end

    if isnothing(regions)
        regions = unique(xgrid[CellRegions])
    end

    if AT !== ON_BFACES && AT !== ON_FACES
        if regions != unique(xgrid[CellRegions])
            dofgrid = subgrid(xgrid, regions)
        else
            dofgrid = xgrid
        end
    end

    # first generate some empty FESpace
    if name == ""
        name = broken ? "$FEType (broken)" : "$FEType"
    end
    ndofs, coffset = count_ndofs(dofgrid, FEType, broken)
    FES = FESpace{Tv, Ti, FEType, AT}(name, broken, ndofs, coffset, xgrid, dofgrid, Dict{Type{<:AbstractGridComponent}, Any}(), Dict{Type{<:AbstractInterpolationOperator}, Any}())

    @debug "Generated FESpace $name ($AT, ndofs=$ndofs)"

    return FES
end

function FESpace{FEType}(
        xgrid::ExtendableGrid{Tv, Ti};
        kwargs...
    ) where {Tv, Ti, FEType <: AbstractFiniteElement}
    return FESpace{FEType, ON_CELLS}(xgrid; kwargs...)
end

"""
$(TYPEDSIGNATURES)

Custom `eltype` function for `FESpace` returns the finite element type parameter of the finite element space.
"""
Base.eltype(::FESpace{Tv, Ti, FEType, APT}) where {Tv, Ti, FEType <: AbstractFiniteElement, APT} = FEType


"""
$(TYPEDSIGNATURES)

returns the assembly type parameter of the finite element space, i.e. on which entities of the grid the finite element is defined.
"""
assemblytype(::FESpace{Tv, Ti, FEType, APT}) where {Tv, Ti, FEType <: AbstractFiniteElement, APT} = APT


"""
$(TYPEDSIGNATURES)

Custom `show` function for `FESpace` that prints some information and all available dofmaps.
"""
function Base.show(io::IO, FES::FESpace{Tv, Ti, FEType, APT}) where {Tv, Ti, FEType <: AbstractFiniteElement, APT}
    println(io, "\nFESpace information")
    println(io, "===================")
    println(io, "     name = $(FES.name)")
    println(io, "   FEType = $FEType ($(FES.broken ? "$APT, broken" : "$APT"))")
    println(io, "  FEClass = $(supertype(FEType))")
    println(io, "    ndofs = $(FES.ndofs) $(FES.coffset !== FES.ndofs ? "(coffset = $(FES.coffset))" : "")")
    println(io, "    xgrid = $(FES.xgrid)")
    println(io, "  dofgrid = $(FES.dofgrid !== FES.xgrid ? FES.dofgrid : "xgrid")")
    if !isempty(FES.dofmaps)
        println(io, "\nDofMaps:")
        for (k, v) in FES.dofmaps
            println(io, "  > $(k): $(typeof(v))")
        end
    end
    if !isempty(FES.interpolators)
        println(io, "\nInterpolators:")
        for (k, v) in FES.interpolators
            println(io, "  > $(k): $(typeof(v))")
        end
    end
    return
end

## used if no coefficient handler or subset handler is needed (to have a single Function type for all)
const NothingFunction = (x, y) -> nothing

#get_polynomialorder(::Type{<:AbstractFiniteElement}, ::Type{Vertex0D}) = 1
"""
$(TYPEDSIGNATURES)

returns the number of components of the FESpace (= number of components of its FEType)
"""
get_ncomponents(FES::FESpace) = get_ncomponents(get_FEType(FES))


"""
$(TYPEDSIGNATURES)

returns the coefficients for local evaluations of finite element functions
( see e.g. h1v_br.jl for a use-case)
"""
function get_coefficients(::Type{<:AssemblyType}, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{<:AbstractElementGeometry}, xgrid = FE.dofgrid) where {Tv, Ti, FEType <: AbstractFiniteElement, APT}
    return NothingFunction
end

"""
$(TYPEDSIGNATURES)

returns a closure function of the form
	
	closure(subset_ids::Array{Int, 1}, cell)

which returns the ids of the local basis functions needed on the cell.
Usually, subset_ids = 1:ndofs (meaning that all basis functions on the reference cells are used),

See e.g. the 3D implementation of BDM1 or H1P3
where different basis functions are chosen
depending on the face orientations (which in 3D is not just a sign)

"""
function get_basissubset(::Type{<:AssemblyType}, FE::FESpace{Tv, Ti, FEType, APT}, ::Type{<:AbstractElementGeometry}, xgrid = FE.dofgrid) where {Tv, Ti, FEType <: AbstractFiniteElement, APT}
    return NothingFunction
end

"""
$(TYPEDSIGNATURES)

returns the number of degrees of freedom for the given AssemblyType, FEType and geometry
"""
get_ndofs(AT::Type{<:AssemblyType}, FEType::Type{<:AbstractFiniteElement}, EG::Type{<:AbstractElementGeometry}) = 0 # element is undefined for this AT or EG

"""
$(TYPEDSIGNATURES)

returns the closure function of form

	closure(refbasis, xref)

that computes the evaluations of the basis functions on the reference geometry
"""
get_basis(AT::Type{<:AssemblyType}, FEType::Type{<:AbstractFiniteElement}, EG::Type{<:AbstractElementGeometry}) = (refbasis, xref) -> nothing


"""
$(TYPEDSIGNATURES)

returns the total number of degrees of freedom for the given AssemblyType, FEType and geometry
"""
get_ndofs_all(AT::Type{<:AssemblyType}, FEType::Type{<:AbstractFiniteElement}, EG::Type{<:AbstractElementGeometry}) = get_ndofs(AT, FEType, EG)

"""
$(TYPEDSIGNATURES)

tells if FEType is defined on this ElementGeometry
"""
isdefined(FEType::Type{<:AbstractFiniteElement}, EG::Type{<:AbstractElementGeometry}) = false
isdefined(FEType::Type{<:AbstractFiniteElement}, EG::Type{<:AbstractElementGeometry}, ::Bool) = isdefined(FEType, EG)


"""
$(TYPEDSIGNATURES)
To be called by getindex. This triggers lazy creation of 
non-existing dofmaps
"""
Base.get!(FES::FESpace, DM::Type{<:DofMap}) = get!(() -> init_dofmap!(FES, DM), FES.dofmaps, DM)

"""
````
Base.getindex(FES::FESpace,DM::Type{<:DofMap})
````
Generic method for obtaining dofmap.
This method is mutating in the sense that non-existing dofmaps
are created on demand.
Due to the fact that components are stored as Any the return
value triggers type instability.
"""
Base.getindex(FES::FESpace, DM::Type{<:DofMap}) = get!(FES, DM)

"""
$(TYPEDSIGNATURES)
counts the total number of degrees of freedom for the FEType
for the whole grid
"""
function count_ndofs(xgrid, FEType, broken::Bool)
    EG = xgrid[UniqueCellGeometries]
    xItemGeometries = xgrid[CellGeometries]
    if length(EG) == 1
        ncells4EG = [(EG[1], num_cells(xgrid))]
    else
        ncells4EG = [(i, count(==(i), xItemGeometries)) for i in EG] # allocations !!!
    end
    ncomponents::Int = get_ncomponents(FEType)
    totaldofs::Int = 0
    offset4component::Int = 0

    ## todo : below it is assumed that number of dofs on nodes/edges/faces is the same for all EG
    ##        (which usually makes sense for unbroken elements, but who knows...)
    for j in 1:length(EG)
        pattern = get_dofmap_pattern(FEType, CellDofs, EG[j])
        parsed_dofmap = ParsedDofMap(pattern, ncomponents, EG[j])
        if broken == true
            # if broken count all dofs here
            totaldofs += ncells4EG[j][2] * get_ndofs(parsed_dofmap)
        else
            # if not broken only count interior dofs on cell here
            totaldofs += ncells4EG[j][2] * get_ndofs(parsed_dofmap, DofTypeInterior)
            offset4component += ncells4EG[j][2] * get_ndofs4c(parsed_dofmap, DofTypeInterior)
        end

        if j == length(EG) && !broken
            # add continuous dofs here
            totaldofs += size(xgrid[Coordinates], 2) * get_ndofs(parsed_dofmap, DofTypeNode)
            if get_ndofs(parsed_dofmap, DofTypeFace) > 0
                totaldofs += num_sources(xgrid[FaceNodes]) * get_ndofs(parsed_dofmap, DofTypeFace)
            end
            if get_ndofs(parsed_dofmap, DofTypeEdge) > 0
                totaldofs += num_sources(xgrid[EdgeNodes]) * get_ndofs(parsed_dofmap, DofTypeEdge)
            end

            # compute also offset4component

            offset4component += num_nodes(xgrid) * get_ndofs4c(parsed_dofmap, DofTypeNode)
            if get_ndofs4c(parsed_dofmap, DofTypeFace) > 0
                nfaces = num_sources(xgrid[FaceNodes])
                offset4component += nfaces * get_ndofs4c(parsed_dofmap, DofTypeFace)
            end
            if get_ndofs4c(parsed_dofmap, DofTypeEdge) > 0
                nedges = num_sources(xgrid[EdgeNodes])
                offset4component += nedges * get_ndofs4c(parsed_dofmap, DofTypeEdge)
            end
            offset4component += num_cells(xgrid) * get_ndofs4c(parsed_dofmap, DofTypePCell)
            #nitems::Int = length(xItemGeometries)
            #offset4component += nitems * get_ndofs4c(parsed_dofmap, DofTypeInterior)
        end
    end

    return totaldofs, offset4component
end


include("fevector.jl");
include("fematrix.jl");
include("interpolations.jl")
include("interpolators.jl")


"""
$(TYPEDSIGNATURES)
returns the element dimension of the finite element
"""
get_edim(FEType::Type{<:AbstractFiniteElement}) = 0 # not defined
#get_polynomialorder(::Type{<:AbstractFiniteElement}, ::Type{<:Vertex0D}) = 0

"""
$(TYPEDSIGNATURES)
returns the offset for the interior degrees of freedom (all vertex, face and edge
dofs are assumed to come first and before that number)
"""
interior_dofs_offset(AT::Type{<:AssemblyType}, FE::Type{<:AbstractFiniteElement}, EG::Type{<:AbstractElementGeometry}) = -1 # should be specified by element if needed

"""
$(TYPEDSIGNATURES)
returns the expected polynomial order of the basis functions
(used to determine default quadrature orders)
"""
get_polynomialorder(AT::Type{<:AssemblyType}, FE::Type{<:AbstractFiniteElement}, EG::Type{<:AbstractElementGeometry}) = -1 # should be specified by element if needed


###########################
# Finite Element Subtypes #
###########################
#
# subtypes of the finite element supertypes above
# each defined in its own file

# Hdiv-conforming elements (only vector-valued)
# lowest order
include("fedefs/hdiv_rt0.jl");
include("fedefs/hdiv_bdm1.jl");
include("fedefs/hdiv_rt1.jl");
include("fedefs/hdiv_bdm2.jl");
include("fedefs/hdiv_rtk.jl");
include("fedefs/hdiv_rtk_enrich.jl");

# H1 conforming elements (also Crouzeix-Raviart)
# no order (just for certain purposes)
include("fedefs/h1_bubble.jl");
# lowest order
include("fedefs/l2_p0.jl");
include("fedefs/l2_p1.jl");
include("fedefs/h1_p1.jl");
include("fedefs/h1_q1.jl");
include("fedefs/h1_mini.jl");
include("fedefs/h1nc_cr.jl");
include("fedefs/h1v_br.jl"); # Bernardi--Raugel (only vector-valued, with coefficients)
include("fedefs/h1v_p1teb.jl"); # P1 + tangential edge bubbles (only vector-valued, with coefficients)
# second order
include("fedefs/h1_p2.jl");
include("fedefs/h1_q2.jl");
include("fedefs/h1_p2b.jl");
# third order
include("fedefs/h1_p3.jl");
# arbitrary order
include("fedefs/h1_pk.jl");

# Hcurl-conforming elements
include("fedefs/hcurl_n0.jl");
include("fedefs/hcurl_n1.jl");


function get_coefficients(::Type{ON_BFACES}, FE::FESpace{Tv, Ti, FEType, APT}, EG::Type{<:AbstractElementGeometry}, xgrid) where {Tv, Ti, FEType <: AbstractFiniteElement, APT}
    get_coeffs_on_face = get_coefficients(ON_FACES, FE, EG, xgrid)
    xBFaceFaces = FE.xgrid[BFaceFaces]
    return function closure(coefficients, bface)
        get_coeffs_on_face(coefficients, xBFaceFaces[bface])
        return nothing
    end
end


function get_reconstruction_matrix(T::Type{<:Real}, FE::FESpace, FER::FESpace)
    xgrid = FE.xgrid
    xCellGeometries = xgrid[CellGeometries]
    EG = xgrid[UniqueCellGeometries]

    FEType = eltype(FE)
    FETypeReconst = eltype(FER)

    ncells = num_sources(xgrid[CellNodes])
    rhandlers = [get_reconstruction_coefficient(ON_CELLS, FE, FER, EG[1])]
    chandlers = [get_coefficients(ON_CELLS, FER, EG[1], xgrid)]
    shandlers = [get_basissubset(ON_CELLS, FER, EG[1])]
    for j in 2:length(EG)
        append!(rhandlers, [get_reconstruction_coefficients(ON_CELLS, FE, FER, EG[j])])
        append!(chandlers, [get_coefficients(ON_CELLS, FER, EG[j], xgrid)])
        append!(shandlers, [get_basissubset(ON_CELLS, FER, EG[j])])
    end

    ndofs_FE = zeros(Int, length(EG))
    ndofs_FER = zeros(Int, length(EG))
    for j in 1:length(EG)
        ndofs_FE[j] = get_ndofs(ON_CELLS, FEType, EG[j])
        ndofs_FER[j] = get_ndofs(ON_CELLS, FETypeReconst, EG[j])
    end

    xCellDofs = FE[CellDofs]
    xCellDofsR = FER[CellDofs]

    ## generate matrix
    A = ExtendableSparseMatrix{T, Int64}(FER.ndofs, FE.ndofs)

    iEG = 1
    cellEG = EG[1]
    ncomponents = get_ncomponents(FEType)
    coefficients = zeros(T, ncomponents, maximum(ndofs_FER))
    basissubset = zeros(Int, maximum(ndofs_FER))
    rcoefficients = zeros(T, maximum(ndofs_FE), maximum(ndofs_FER))
    dof::Int = 0
    dofR::Int = 0
    for cell in 1:ncells
        if length(EG) > 1
            cellEG = xCellGeometries[cell]
            for j in 1:length(EG)
                if cellEG == EG[j]
                    iEG = j
                    break
                end
            end
        end
        # get Hdiv coefficients and subset
        chandlers[iEG](coefficients, cell)
        shandlers[iEG](basissubset, cell)

        # get reconstruction coefficients
        rhandlers[iEG](rcoefficients, cell)

        for dof_i in 1:ndofs_FE[iEG]
            dof = xCellDofs[dof_i, cell]
            for dof_j in 1:ndofs_FER[iEG]
                if rcoefficients[dof_i, dof_j] != 0
                    dofR = xCellDofsR[dof_j, cell]
                    A[dofR, dof] = rcoefficients[dof_i, dof_j]
                end
            end
        end
    end
    flush!(A)
    return A
end

function get_local_coffsets(FEType, AT, EG)
    ## get dofmap pattern
    pattern = get_dofmap_pattern(FEType, Dofmap4AssemblyType(AT), EG)
    ncomponents = get_ncomponents(FEType)
    parsed_dofmap = ParsedDofMap(pattern, ncomponents, EG)
    local_coffset = 0
    local_coffset += num_nodes(EG) * get_ndofs4c(parsed_dofmap, DofTypeNode)
    if get_ndofs4c(parsed_dofmap, DofTypeFace) > 0
        local_coffset += num_faces(EG) * get_ndofs4c(parsed_dofmap, DofTypeFace)
    end
    if get_ndofs4c(parsed_dofmap, DofTypeEdge) > 0
        local_coffset += num_edges(EG) * get_ndofs4c(parsed_dofmap, DofTypeEdge)
    end
    local_coffset += get_ndofs4c(parsed_dofmap, DofTypeInterior)
    if local_coffset == 0
        return Int[]
    else
        return 0:local_coffset:(ncomponents * local_coffset)
    end
end
