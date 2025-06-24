"""
$(TYPEDEF)

A `DofMap` encodes the association of global DoF indices to mesh entities of a given type (e.g., cells, faces, edges) for a particular finite element space. This mapping is essential for assembling system matrices and vectors, applying boundary conditions, and extracting solution values.
DofMaps are typically not constructed directly by users. Instead, they are generated/managed internally by the finite element space and accessed as needed for assembly or evaluation tasks.

# Implementation

- All concrete `DofMap` subtypes (e.g., `CellDofs`, `FaceDofs`, `EdgeDofs`, etc.) specify the mesh entity type to which DoFs are attached.
- DofMaps are stored as `ExtendableGrids.AbstractGridAdjacency` (usually VariableTargetAdjacency, SerialVariableTargetAdjacency or AbstractGridIntegerArray2D) objects within the finite element space and are generated automatically as needed.
- The appropriate DofMap for a given assembly type can be accessed via `FESpace[DofMapSubtype]`.

# Example

```julia
cell_dofs = FESpace[CellDofs]   # Get the cell-based DoF map
face_dofs = FESpace[FaceDofs]   # Get the face-based DoF map
```
"""
abstract type DofMap <: AbstractGridAdjacency end


"""
	$(TYPEDEF)

Key type describing the dofs for each cell of the dofgrid
"""
abstract type CellDofs <: DofMap end

"""
	$(TYPEDEF)

Key type describing the dofs for each face of the dofgrid
"""
abstract type FaceDofs <: DofMap end

"""
	$(TYPEDEF)

Key type describing the dofs for each edge of the dofgrid
"""
abstract type EdgeDofs <: DofMap end

"""
	$(TYPEDEF)

Key type describing the dofs for each boundary face of the dofgrid
"""
abstract type BFaceDofs <: DofMap end

"""
	$(TYPEDEF)

Key type describing the dofs for each boundary edge of the dofgrid
"""
abstract type BEdgeDofs <: DofMap end


"""
	$(TYPEDEF)

Key type describing the dofs for each cell of the parentgrid
"""
abstract type CellDofsParent <: DofMap end

"""
	$(TYPEDEF)

Key type describing the dofs for each face of the parentgrid
"""
abstract type FaceDofsParent <: DofMap end

"""
	$(TYPEDEF)

Key type describing the dofs for each edge of the parentgrid
"""
abstract type EdgeDofsParent <: DofMap end

"""
	$(TYPEDEF)

Key type describing the dofs for each boundary face of the parentgrid
"""
abstract type BFaceDofsParent <: DofMap end

"""
	$(TYPEDEF)

Key type describing the dofs for each boundary edge of the parentgrid
"""
abstract type BEdgeDofsParent <: DofMap end

"""
	$(TYPEDEF)

Union type for all dofmap adjacency types.
"""
const DofMapTypes{Ti} = Union{VariableTargetAdjacency{Ti}, SerialVariableTargetAdjacency{Ti}, Array{Ti, 2}}

"""
$(TYPEDSIGNATURES)
Unique geometry grid component for dofmaps
"""
UCG4DofMap(::Type{CellDofs}) = UniqueCellGeometries
UCG4DofMap(::Type{FaceDofs}) = UniqueFaceGeometries
UCG4DofMap(::Type{EdgeDofs}) = UniqueEdgeGeometries
UCG4DofMap(::Type{BFaceDofs}) = UniqueBFaceGeometries
UCG4DofMap(::Type{BEdgeDofs}) = UniqueBEdgeGeometries


"""
$(TYPEDSIGNATURES)
ItemNodes grid component for dofmaps
"""
SuperItemNodes4DofMap(::Type{CellDofs}) = CellNodes
SuperItemNodes4DofMap(::Type{FaceDofs}) = FaceNodes
SuperItemNodes4DofMap(::Type{EdgeDofs}) = EdgeNodes
SuperItemNodes4DofMap(::Type{BFaceDofs}) = FaceNodes
SuperItemNodes4DofMap(::Type{BEdgeDofs}) = EdgeNodes

"""
$(TYPEDSIGNATURES)
ItemGeomtries grid component for dofmaps
"""
ItemGeometries4DofMap(::Type{CellDofs}) = CellGeometries
ItemGeometries4DofMap(::Type{FaceDofs}) = FaceGeometries
ItemGeometries4DofMap(::Type{EdgeDofs}) = EdgeGeometries
ItemGeometries4DofMap(::Type{BFaceDofs}) = BFaceGeometries
ItemGeometries4DofMap(::Type{BEdgeDofs}) = BEdgeGeometries

"""
$(TYPEDSIGNATURES)
ItemEdges grid component for dofmaps
"""
ItemEdges4DofMap(::Type{CellDofs}) = CellEdges
ItemEdges4DofMap(::Type{FaceDofs}) = FaceEdges
ItemEdges4DofMap(::Type{BFaceDofs}) = FaceEdges

"""
$(TYPEDSIGNATURES)
SubItemEdges grid component for dofmaps
"""
Sub2Sup4DofMap(::Type{<:DofMap}) = nothing
Sub2Sup4DofMap(::Type{BFaceDofs}) = BFaceFaces
Sub2Sup4DofMap(::Type{BEdgeDofs}) = BEdgeEdges


"""
$(TYPEDSIGNATURES)
Default Dofmap for AssemblyType
"""
Dofmap4AssemblyType(::Type{ON_CELLS}) = CellDofs
Dofmap4AssemblyType(::Type{<:ON_FACES}) = FaceDofs
Dofmap4AssemblyType(::Type{ON_BFACES}) = BFaceDofs
Dofmap4AssemblyType(::Type{<:ON_EDGES}) = EdgeDofs
Dofmap4AssemblyType(::Type{ON_BEDGES}) = BEdgeDofs

"""
$(TYPEDSIGNATURES)
Parent Dofmap for Dofmap
"""
ParentDofmap4Dofmap(::Type{CellDofs}) = CellDofsParent
ParentDofmap4Dofmap(::Type{FaceDofs}) = FaceDofsParent
ParentDofmap4Dofmap(::Type{EdgeDofs}) = EdgeDofsParent
ParentDofmap4Dofmap(::Type{BFaceDofs}) = BFaceDofsParent
ParentDofmap4Dofmap(::Type{BEdgeDofs}) = BEdgeDofsParent

"""
$(TYPEDSIGNATURES)
Effective AssemblyType (on the subgrid)
for two AssemblyTypes where the first one
is related to where the finite element functions live and the second one
to where something should be assembled both with respect to the common parent grid
(e.g. face-based finite elements live on a subgrid of all faces, where the faces
are the cells in this subgrid, and they cannot be evaluated over the cells of the parentgrid,
but on the faces of the parengrid, which are the cells in the subgrid)
"""
EffAT4AssemblyType(::Type{ON_CELLS}, ::Type{ON_CELLS}) = ON_CELLS
EffAT4AssemblyType(::Type{ON_CELLS}, ::Type{<:ON_FACES}) = ON_FACES
EffAT4AssemblyType(::Type{ON_CELLS}, ::Type{ON_BFACES}) = ON_BFACES
EffAT4AssemblyType(::Type{ON_CELLS}, ::Type{<:ON_EDGES}) = ON_EDGES
EffAT4AssemblyType(::Type{ON_CELLS}, ::Type{ON_BEDGES}) = ON_BEDGES

EffAT4AssemblyType(::Type{ON_FACES}, ::Type{ON_CELLS}) = nothing
EffAT4AssemblyType(::Type{ON_FACES}, ::Type{<:ON_FACES}) = ON_CELLS
EffAT4AssemblyType(::Type{ON_FACES}, ::Type{<:ON_EDGES}) = ON_FACES
EffAT4AssemblyType(::Type{ON_FACES}, ::Type{<:ON_BEDGES}) = ON_BFACES

EffAT4AssemblyType(::Type{ON_BFACES}, ::Type{ON_CELLS}) = nothing
EffAT4AssemblyType(::Type{ON_BFACES}, ::Type{ON_BFACES}) = ON_CELLS
EffAT4AssemblyType(::Type{ON_BFACES}, ::Type{<:ON_FACES}) = nothing
EffAT4AssemblyType(::Type{ON_BFACES}, ::Type{<:ON_EDGES}) = nothing
EffAT4AssemblyType(::Type{ON_BFACES}, ::Type{<:ON_BEDGES}) = nothing

EffAT4AssemblyType(::Type{ON_EDGES}, ::Type{ON_CELLS}) = nothing
EffAT4AssemblyType(::Type{ON_EDGES}, ::Type{<:ON_FACES}) = nothing
EffAT4AssemblyType(::Type{ON_EDGES}, ::Type{<:ON_EDGES}) = ON_CELLS


"""
	$(TYPEDEF)

Abstract type for all dof types
"""
abstract type DofType end

"""
	$(TYPEDEF)

Dof connected to a single vertex
"""
abstract type DofTypeNode <: DofType end        # parsed from 'N' or 'n' (nodal continuous dof)

"""
	$(TYPEDEF)

Dof connected to a face
"""
abstract type DofTypeFace <: DofType end        # parsed from 'F' or 'f' (face continuous dof)

"""
	$(TYPEDEF)

Dof connected to an edge
"""
abstract type DofTypeEdge <: DofType end        # parsed from 'E' or 'e' (edge continuous dof)

"""
	$(TYPEDEF)

Dof connected to the interior of an item
"""
abstract type DofTypeInterior <: DofType end    # parsed from 'I' or 'i' (interior dof)

"""
	$(TYPEDEF)

Dof connected to a parent cell
"""
abstract type DofTypePCell <: DofType end       # parsed from 'C' or 'c' (parent cell dof, only needed by P0 element for BFaceDofs)

"""
$(TYPEDSIGNATURES)
parses a Char into a DofType
"""
function DofType(c::Char)
    if lowercase(c) == 'n'
        return DofTypeNode
    elseif lowercase(c) == 'f'
        return DofTypeFace
    elseif lowercase(c) == 'e'
        return DofTypeEdge
    elseif lowercase(c) == 'i'
        return DofTypeInterior
    elseif lowercase(c) == 'c'
        return DofTypePCell
    else
        @error "No DofType available to parse from $(lowercase(c))"
    end
end

const dofmap_type_chars = ['E', 'N', 'I', 'C', 'F', 'e', 'n', 'i', 'c', 'f']
const dofmap_number_chars = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']


"""
	$(TYPEDEF)

Pattern segment of a dofmap sequence
"""
struct DofMapPatternSegment
    type::Type{<:DofType}
    each_component::Bool    # single dof for all components or one for each component ?
    ndofs::Int              # how many dofs
end


"""
$(TYPEDSIGNATURES)
splits a pattern string into pairs of single chars and Ints
and generated a sequence of DofMapPatternSegments
"""
function parse_pattern(pattern::String)
    pairs = Array{DofMapPatternSegment, 1}([])
    for j in 1:length(pattern)
        if pattern[j] in dofmap_type_chars
            k = j + 1
            while (k < length(pattern))
                if (pattern[k + 1] in dofmap_number_chars)
                    k += 1
                else
                    break
                end
            end
            # @show pattern[j] SubString(pattern,j+1,k)
            push!(pairs, DofMapPatternSegment(DofType(pattern[j]), isuppercase(pattern[j]), tryparse(Int, SubString(pattern, j + 1, k))))
        end
    end
    return pairs
end


"""
	$(TYPEDEF)

Parsed dofmap, basically a sequence of DofMapPatternSegment
"""
struct ParsedDofMap
    segments::Array{DofMapPatternSegment, 1}
    ndofs_node4c::Int
    ndofs_face4c::Int
    ndofs_edge4c::Int
    ndofs_interior4c::Int
    ndofs_pcell4c::Int
    ndofs_node::Int
    ndofs_face::Int
    ndofs_edge::Int
    ndofs_interior::Int
    ndofs_pcell::Int
    ndofs_total::Int
end

"""
$(TYPEDSIGNATURES)
total number of dofs of the DofType for a single component
"""
get_ndofs4c(P::ParsedDofMap, ::Type{DofTypeNode}) = P.ndofs_node4c
get_ndofs4c(P::ParsedDofMap, ::Type{DofTypeEdge}) = P.ndofs_edge4c
get_ndofs4c(P::ParsedDofMap, ::Type{DofTypeFace}) = P.ndofs_face4c
get_ndofs4c(P::ParsedDofMap, ::Type{DofTypeInterior}) = P.ndofs_interior4c
get_ndofs4c(P::ParsedDofMap, ::Type{DofTypePCell}) = P.ndofs_pcell4c

"""
$(TYPEDSIGNATURES)
total number of dofs of the DofType
"""
get_ndofs(P::ParsedDofMap, ::Type{DofTypeNode}) = P.ndofs_node
get_ndofs(P::ParsedDofMap, ::Type{DofTypeEdge}) = P.ndofs_edge
get_ndofs(P::ParsedDofMap, ::Type{DofTypeFace}) = P.ndofs_face
get_ndofs(P::ParsedDofMap, ::Type{DofTypeInterior}) = P.ndofs_interior
get_ndofs(P::ParsedDofMap, ::Type{DofTypePCell}) = P.ndofs_pcell
get_ndofs(P::ParsedDofMap) = P.ndofs_total

"""
$(TYPEDSIGNATURES)
parses a dofmap string (defined in each finite element definition file)
and generated a ParsedDofMap for a certain number of components
and the given geometry
"""
function ParsedDofMap(pattern::String, ncomponents, EG::Type{<:AbstractElementGeometry})
    segments = parse_pattern(pattern)
    ndofs_node::Int = 0
    ndofs_face::Int = 0
    ndofs_edge::Int = 0
    ndofs_interior::Int = 0
    ndofs_pcell::Int = 0
    ndofs_node4c::Int = 0
    ndofs_face4c::Int = 0
    ndofs_edge4c::Int = 0
    ndofs_interior4c::Int = 0
    ndofs_pcell4c::Int = 0
    for j in 1:length(segments)
        if segments[j].type <: DofTypeNode
            ndofs_node += segments[j].ndofs * (segments[j].each_component ? ncomponents : 1)
            if segments[j].each_component
                ndofs_node4c += segments[j].ndofs
            end
        elseif segments[j].type <: DofTypeFace
            ndofs_face += segments[j].ndofs * (segments[j].each_component ? ncomponents : 1)
            if segments[j].each_component
                ndofs_face4c += segments[j].ndofs
            end
        elseif segments[j].type <: DofTypeEdge
            ndofs_edge += segments[j].ndofs * (segments[j].each_component ? ncomponents : 1)
            if segments[j].each_component
                ndofs_edge4c += segments[j].ndofs
            end
        elseif segments[j].type <: DofTypeInterior
            ndofs_interior += segments[j].ndofs * (segments[j].each_component ? ncomponents : 1)
            if segments[j].each_component
                ndofs_interior4c += segments[j].ndofs
            end
        elseif segments[j].type <: DofTypePCell
            ndofs_pcell += segments[j].ndofs * (segments[j].each_component ? ncomponents : 1)
            if segments[j].each_component
                ndofs_pcell4c += segments[j].ndofs
            end
        end
    end
    ndofs_total = num_nodes(EG) * ndofs_node + num_faces(EG) * ndofs_face + num_edges(EG) * ndofs_edge + ndofs_interior
    return ParsedDofMap(segments, ndofs_node4c, ndofs_face4c, ndofs_edge4c, ndofs_interior4c, ndofs_pcell4c, ndofs_node, ndofs_face, ndofs_edge, ndofs_interior, ndofs_pcell, ndofs_total)
end


"""
$(TYPEDSIGNATURES)
yields the coressponding dofmap of the FESpace for a given AssemblyType (assumed with respect to the
(parent)grid of the FESpace)
"""
function Dofmap4AssemblyType(FES::FESpace, AT::Type{<:AssemblyType})
    return FES[Dofmap4AssemblyType(EffAT4AssemblyType(assemblytype(FES), AT))]
end


"""
$(TYPEDSIGNATURES)
generates the DofMap for the FESpace based on the pattern of the FEType 
"""
function init_dofmap_from_pattern!(FES::FESpace{Tv, Ti, FEType, APT}, DM::Type{<:DofMap}) where {Tv, Ti, FEType <: AbstractFiniteElement, APT}
    ## Beware: Automatic broken DofMap generation currently only reliable for CellDofs

    @debug "Generating $DM for $(FES.name)"
    ## prepare dofmap patterns
    xgrid = FES.dofgrid
    EG = xgrid[UCG4DofMap(DM)]
    ncomponents::Int = get_ncomponents(FEType)
    need_faces::Bool = false
    need_edges::Bool = false
    maxdofs4item::Int = 0
    dofmap4EG::Array{ParsedDofMap, 1} = Array{ParsedDofMap, 1}(undef, length(EG))
    has_interiordofs = zeros(Bool, length(EG))
    for j in 1:length(EG)
        pattern = get_dofmap_pattern(FEType, DM, EG[j])
        dofmap4EG[j] = ParsedDofMap(pattern, ncomponents, EG[j])
        maxdofs4item = max(maxdofs4item, get_ndofs(dofmap4EG[j]))
        if get_ndofs(dofmap4EG[j], DofTypeFace) > 0
            need_faces = true
        end
        if get_ndofs(dofmap4EG[j], DofTypeEdge) > 0
            need_edges = true
        end
        if get_ndofs(dofmap4EG[j], DofTypeInterior) > 0
            has_interiordofs[j] = true
        end
    end

    ## prepare data
    nnodes::Int = size(xgrid[Coordinates], 2)
    xItemNodes::Adjacency{Ti} = xgrid[SuperItemNodes4DofMap(DM)]
    xItemGeometries::GridEGTypes = xgrid[ItemGeometries4DofMap(DM)]
    if need_faces
        xItemFaces::Adjacency{Ti} = xgrid[CellFaces]
        nfaces = num_sources(xgrid[FaceNodes])
    end
    if need_edges
        xItemEdges::Adjacency{Ti} = xgrid[ItemEdges4DofMap(DM)]
        nedges = num_sources(xgrid[EdgeNodes])
    end
    nitems::Int = num_sources(xItemNodes)
    offset4component = FES.coffset

    ## generate dofmap from patterns
    sub2sup = nothing
    nsubitems::Int = 0
    if Sub2Sup4DofMap(DM) !== nothing
        Sub2SupItems = xgrid[Sub2Sup4DofMap(DM)]
        sub2sup = (x) -> Sub2SupItems[x]
        nsubitems = length(Sub2SupItems)
    else
        sub2sup = (x) -> x
        nsubitems = nitems
    end

    nnodes4EG::Array{Int, 1} = zeros(Int, length(EG))
    nfaces4EG::Array{Int, 1} = zeros(Int, length(EG))
    nedges4EG::Array{Int, 1} = zeros(Int, length(EG))
    for j in 1:length(EG)
        nnodes4EG[j] = num_nodes(EG[j])
        nfaces4EG[j] = num_faces(EG[j])
        nedges4EG[j] = num_edges(EG[j])
    end

    itemEG = EG[1]
    iEG::Int = 1

    ## get total dofnumber for pre--allocation
    dofmap_totallength = 0
    colstarts = zeros(Ti, nsubitems + 1)
    for subitem in 1:nsubitems
        itemEG = xItemGeometries[subitem]
        if length(EG) > 1
            iEG = findfirst(isequal(itemEG), EG)
        end
        colstarts[subitem] = dofmap_totallength + 1
        dofmap_totallength += get_ndofs(dofmap4EG[iEG])
    end
    colstarts[end] = dofmap_totallength + 1

    if FES.broken
        xItemDofs = SerialVariableTargetAdjacency(colstarts)
    else

        colentries = zeros(Ti, dofmap_totallength)
        xItemDofs = VariableTargetAdjacency{Ti}(colentries, colstarts)

        cpattern::Array{DofMapPatternSegment, 1} = dofmap4EG[1].segments
        l::Int = 0
        offset::Int = 0
        pos::Int = 0
        q::Int = 0
        item_with_interiordofs::Int = 0
        for subitem in 1:nsubitems
            itemEG = xItemGeometries[subitem]
            item = sub2sup(subitem)
            if length(EG) > 1
                iEG = findfirst(isequal(itemEG), EG)
            end
            cpattern = dofmap4EG[iEG].segments
            l = length(cpattern)
            if has_interiordofs[iEG]
                item_with_interiordofs += 1
            end
            for c in 1:ncomponents
                offset = (c - 1) * offset4component
                for k in 1:l
                    q = cpattern[k].ndofs
                    if cpattern[k].type <: DofTypeNode && cpattern[k].each_component
                        for n in 1:nnodes4EG[iEG]
                            for m in 1:q
                                pos += 1
                                xItemDofs.colentries[pos] = xItemNodes[n, item] + offset + (m - 1) * nnodes
                            end
                        end
                        offset += nnodes * q
                    elseif cpattern[k].type <: DofTypeFace && cpattern[k].each_component
                        for n in 1:nfaces4EG[iEG]
                            for m in 1:q
                                pos += 1
                                xItemDofs.colentries[pos] = xItemFaces[n, item] + offset + (m - 1) * nfaces
                            end
                        end
                        offset += nfaces * q
                    elseif cpattern[k].type <: DofTypeEdge && cpattern[k].each_component
                        for n in 1:nedges4EG[iEG]
                            for m in 1:q
                                pos += 1
                                xItemDofs.colentries[pos] = xItemEdges[n, item] + offset + (m - 1) * nedges
                            end
                        end
                        offset += nedges * q
                    elseif cpattern[k].type <: DofTypeInterior && cpattern[k].each_component
                        for m in 1:q
                            pos += 1
                            xItemDofs.colentries[pos] = sub2sup(item_with_interiordofs) + offset
                            offset += nitems
                        end
                    end
                end
            end
            offset = ncomponents * offset4component
            for k in 1:l
                q = cpattern[k].ndofs
                if cpattern[k].type <: DofTypeFace && !cpattern[k].each_component
                    for n in 1:nfaces4EG[iEG]
                        for m in 1:q
                            pos += 1
                            xItemDofs.colentries[pos] = xItemFaces[n, item] + offset + (m - 1) * nfaces
                        end
                    end
                    offset += nfaces * q
                elseif cpattern[k].type <: DofTypeEdge && !cpattern[k].each_component
                    for n in 1:nedges4EG[iEG]
                        for m in 1:q
                            pos += 1
                            xItemDofs.colentries[pos] = xItemEdges[n, item] + offset + (m - 1) * nedges
                        end
                    end
                    offset += nedges * q
                elseif cpattern[k].type <: DofTypeInterior && !cpattern[k].each_component
                    for m in 1:q
                        pos += 1
                        xItemDofs.colentries[pos] = sub2sup(item_with_interiordofs) + offset
                        offset += nitems
                    end
                end
            end
        end
    end
    FES[DM] = xItemDofs


    if FES.dofgrid !== FES.xgrid
        ## assume parent relation between xgrid and dofgrid
        @assert FES.dofgrid[ParentGrid] == FES.xgrid "xgrid is not the parent grid of dofgrid !!!"
        @assert FES.dofgrid[ParentGridRelation] <: SubGrid "dofgrid needs to be a subgrid of xgrid"
        SubGridAssemblyType = FES.dofgrid[ParentGridRelation].parameters[1]
        ## lift subgrid dofmap to parent xgrid to allow assembly on common parent grid
        ## by constructing a variable target adjacency with empty dofs for parent items in xgrid
        ## that are not in the dofgrid --> this allows to use the dofmap in assembly over full xgrid
        xItemDofsLifted = VariableTargetAdjacency(Ti)
        if DM == CellDofs
            parentitems = FES.dofgrid[CellParents]
        elseif DM == FaceDofs
            parentitems = FES.dofgrid[FaceParents]
        elseif DM == EdgeDofs
            parentitems = FES.dofgrid[EdgeParents]
        elseif DM == BFaceDofs
            parentitems = FES.dofgrid[BFaceParents]
        elseif DM == BEdgeDofs
            parentitems = FES.dofgrid[BEdgeParents]
        end
        DMP = ParentDofmap4Dofmap(DM)
        subitem = 1
        for i in 1:num_sources(xItemDofs)
            t = parentitems[i]
            for j in subitem:(t - 1)
                append!(xItemDofsLifted, [])
                subitem += 1
            end
            append!(xItemDofsLifted, view(xItemDofs, :, i))
            subitem += 1
        end

        FES[DMP] = xItemDofsLifted
    end
    return FES[DM]
end


"""
$(TYPEDSIGNATURES)
generates the BFaceDofs or FaceDofs DofMap for a broken FESpace based on the pattern of the FEType 
"""
function init_broken_dofmap!(FES::FESpace{Tv, Ti, FEType, APT}, DM::Union{Type{BFaceDofs}, Type{FaceDofs}}) where {Tv, Ti, FEType <: AbstractFiniteElement, APT}

    ## prepare dofmap patterns
    xgrid = FES.dofgrid
    cellEG = xgrid[UniqueCellGeometries]
    EG = xgrid[UCG4DofMap(DM)]
    ncomponents::Int = get_ncomponents(FEType)
    need_nodes = false
    need_faces = false
    need_edges = false
    maxdofs4item::Int = 0
    dofmap4cellEG::Array{ParsedDofMap, 1} = Array{ParsedDofMap, 1}(undef, length(cellEG))
    dofmap4EG::Array{ParsedDofMap, 1} = Array{ParsedDofMap, 1}(undef, length(EG))
    for j in 1:length(cellEG)
        pattern = get_dofmap_pattern(FEType, CellDofs, cellEG[j])
        dofmap4cellEG[j] = ParsedDofMap(pattern, ncomponents, cellEG[j])
    end
    for j in 1:length(EG)
        pattern = get_dofmap_pattern(FEType, DM, EG[j])
        dofmap4EG[j] = ParsedDofMap(pattern, ncomponents, EG[j])
        maxdofs4item = max(maxdofs4item, get_ndofs(dofmap4EG[j]))
        if get_ndofs(dofmap4EG[j], DofTypeNode) > 0
            need_nodes = true
        end
        if get_ndofs(dofmap4EG[j], DofTypeInterior) > 0
            need_faces = true
        end
        if get_ndofs(dofmap4EG[j], DofTypeEdge) > 0
            need_edges = true
        end
    end

    xFaceCells = xgrid[FaceCells]
    xCellNodes = xgrid[CellNodes]
    xCellDofs = FES[CellDofs]
    xFaceNodes = xgrid[SuperItemNodes4DofMap(DM)]
    xCellGeometries = xgrid[CellGeometries]
    xFaceGeometries = xgrid[ItemGeometries4DofMap(DM)]
    xFaceDofs = VariableTargetAdjacency(Ti)
    if DM == BFaceDofs
        xRealFace = xgrid[BFaceFaces]
        nfaces = length(xRealFace)
    elseif DM == FaceDofs
        nfaces = num_sources(xFaceNodes)
        xRealFace = 1:nfaces
    end

    if need_faces
        xCellFaces = xgrid[CellFaces]
    end
    if need_edges
        xFaceEdges = xgrid[FaceEdges]
        xCellEdges = xgrid[CellEdges]
        local_edges = zeros(Int, max_num_targets_per_source(xFaceEdges))
    end

    cpattern::Array{DofMapPatternSegment, 1} = dofmap4EG[1].segments
    ccellpattern::Array{DofMapPatternSegment, 1} = dofmap4cellEG[1].segments
    iEG::Int = 1
    local_nodes = zeros(Int, max_num_targets_per_source(xFaceNodes))
    local_face::Int = 0
    local_dofs = zeros(Int, 2 * max_num_targets_per_source(xCellDofs))
    nldofs::Int = 0
    celldofoffset::Int = 0
    node::Int = 0
    face::Int = 0
    cell::Int = 0
    nfacenodes::Int = 0
    nfaceedges::Int = 0
    ncellnodes::Int = 0
    ncelledges::Int = 0
    ncellfaces::Int = 0
    pos::Int = 0
    cEG = xCellGeometries[1]
    for f in 1:nfaces
        face = xRealFace[f]
        cEG = xFaceGeometries[f]
        iEG = findfirst(isequal(cEG), EG)
        cpattern = dofmap4EG[iEG].segments
        nldofs = 0
        for icell in 1:2
            cell = xFaceCells[icell, face]
            if cell > 0
                cEG = xCellGeometries[cell]
                iEG = findfirst(isequal(cEG), cellEG)
                ccellpattern = dofmap4cellEG[iEG].segments

                ## get local nodes/faces/edges dofs for global face f
                if need_nodes
                    ncellnodes = num_targets(xCellNodes, cell)
                    nfacenodes = num_targets(xFaceNodes, face)
                    for k in 1:nfacenodes
                        node = xFaceNodes[k, face]
                        pos = 1
                        while xCellNodes[pos, cell] != node
                            pos += 1
                        end
                        local_nodes[k] = pos
                    end
                end
                if need_faces
                    ncellfaces = num_targets(xCellFaces, cell)
                    pos = 1
                    while xCellFaces[pos, cell] != face
                        pos += 1
                    end
                    local_face = pos
                end
                if need_edges
                    ncelledges = num_targets(xCellEdges, cell)
                    nfaceedges = num_targets(xFaceEdges, face)
                    for k in 1:nfaceedges
                        edge = xFaceEdges[k, face]
                        pos = 1
                        while xCellEdges[pos, cell] != edge
                            pos += 1
                        end
                        local_edges[k] = pos
                    end
                end

                ## get face dofs of cell and write them into local_dofs
                ## assuming that CellDofs and FaceDofs patterns are consistent (e.g. "N1F1" and "N1I1")

                ## get each-component dofs on cell
                celldofoffset = 0
                for c in 1:ncomponents
                    for s in 1:length(cpattern)
                        q = cpattern[s].ndofs
                        if cpattern[s].type <: DofTypePCell && cpattern[s].each_component
                            for dof in 1:q
                                nldofs += 1
                                local_dofs[nldofs] = xCellDofs[celldofoffset + dof, cell]
                            end
                            celldofoffset += q
                        elseif cpattern[s].type <: DofTypeNode && cpattern[s].each_component
                            for k in 1:nfacenodes
                                for dof in 1:q
                                    nldofs += 1
                                    local_dofs[nldofs] = xCellDofs[celldofoffset + (local_nodes[k] - 1) * q + dof, cell]
                                end
                            end
                            celldofoffset += q * ncellnodes
                        elseif cpattern[s].type <: DofTypeInterior && cpattern[s].each_component
                            for dof in 1:q
                                nldofs += 1
                                local_dofs[nldofs] = xCellDofs[celldofoffset + (local_face - 1) * q + dof, cell]
                            end
                            celldofoffset += q * ncellfaces
                        elseif cpattern[s].type <: DofTypeEdge && cpattern[s].each_component
                            for k in 1:nfaceedges
                                for dof in 1:q
                                    nldofs += 1
                                    local_dofs[nldofs] = xCellDofs[celldofoffset + (local_edges[k] - 1) * q + dof, cell]
                                end
                            end
                            celldofoffset += q * ncelledges
                        end
                    end
                    # increase celldofoffset for interior component dofs in cell pattern
                    for s in 1:length(ccellpattern)
                        if ccellpattern[s].type <: DofTypeInterior && ccellpattern[s].each_component
                            celldofoffset += ccellpattern[s].ndofs
                        end
                    end
                end

                ## add single component dofs on cell
                for s in 1:length(cpattern)
                    q = cpattern[s].ndofs
                    if cpattern[s].type <: DofTypeInterior && !cpattern[s].each_component
                        for dof in 1:q
                            nldofs += 1
                            local_dofs[nldofs] = xCellDofs[celldofoffset + (local_face - 1) * q + dof, cell]
                        end
                        celldofoffset += q
                    elseif cpattern[s].type <: DofTypePCell && !cpattern[s].each_component
                        for dof in 1:q
                            nldofs += 1
                            local_dofs[nldofs] = xCellDofs[celldofoffset + dof, cell]
                        end
                        celldofoffset += q
                    elseif cpattern[s].type <: DofTypeNode && !cpattern[s].each_component
                        for k in 1:nfacenodes
                            for dof in 1:q
                                nldofs += 1
                                local_dofs[nldofs] = xCellDofs[celldofoffset + (local_nodes[k] - 1) * q + dof, cell]
                            end
                        end
                        celldofoffset += q * ncellnodes
                    elseif cpattern[s].type <: DofTypeInterior && !cpattern[s].each_component
                        for dof in 1:q
                            nldofs += 1
                            local_dofs[nldofs] = xCellDofs[celldofoffset + (local_face - 1) * q + dof, cell]
                        end
                        celldofoffset += q
                    elseif cpattern[s].type <: DofTypeEdge && !cpattern[s].each_component
                        for k in 1:nfaceedges
                            for dof in 1:q
                                nldofs += 1
                                local_dofs[nldofs] = xCellDofs[celldofoffset + (local_edges[k] - 1) * q + dof, cell]
                            end
                        end
                        celldofoffset += q * ncelledges
                    end
                end
            end
        end
        append!(xFaceDofs, local_dofs[1:nldofs])
    end
    return FES[DM] = xFaceDofs
end


"""
$(TYPEDSIGNATURES)
generates the requested DofMap for the FESpace
"""
function init_dofmap!(FES::FESpace, DM::Type{<:DofMap})
    @debug "Generating dofmap $DM for FESpace $(FES.name)"

    return if (FES.broken == true) && DM == CellDofsParent
        init_dofmap_from_pattern!(FES, CellDofs)
    elseif (FES.broken == true) && (DM != CellDofs)
        ## dofmap needs to include all (e.g. face) dofs from neighbouring cells
        init_broken_dofmap!(FES, DM)
    else
        init_dofmap_from_pattern!(FES, DM)
    end
end


"""
$(TYPEDSIGNATURES)
returns an array with the number of all dofs associated to a boundary dofmap
(default: BFaceDofs) and certain boundary regions (default: all regions)
"""
function boundarydofs(FES; dofmap = BFaceDofs, regions = :all)
    bitemdofs = FES[dofmap]
    if regions === :all
        if typeof(bitemdofs) <: VariableTargetAdjacency
            return bitemdofs.colentries
        elseif typeof(bitemdofs) <: SerialVariableTargetAdjacency
            return 1:(bitemdofs.colstart[end] - 1)
        else
            @error "dofmap has unexpected type"
        end
    else
        if dofmap == BFaceDofs
            bitemregions = FES[BFaceRegions]
        elseif dofmap == BFaceEdges
            bitemregions = FES[BEdgeRegions]
            @error "do not know where to look for regions"
        end
        bdofs = []
        nbitems::Int = num_sources(bitemdofs)
        for bface in 1:nbitems
            for j in 1:num_targets(bitemdofs, 1)
                if bitemregions[bface] in regions
                    push!(bdofs, bitemdofs[j, bface])
                end
            end
        end
        return unique(bdofs)
    end
end
