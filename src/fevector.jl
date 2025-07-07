################
### FEVector ###
################
#
# used to store coefficients for FESpaces and can have several blocks of different FESpaces
# acts like an AbstractArray{T,1}

"""
$(TYPEDEF)

A block within an `FEVector` representing a contiguous segment of coefficients associated with a specific finite element space (`FESpace`).

Each `FEVectorBlock` provides array-like access to the degrees of freedom (DOFs) for its associated `FESpace`, mapping local indices to a shared global coefficient array. This enables efficient block-wise operations, assembly, and extraction of sub-vectors corresponding to different FE spaces.

# Type Parameters
- `T`: Value type of the vector entries (e.g., `Float64`).
- `Tv`: Value type for the associated `FESpace`.
- `Ti`: Integer type for the associated `FESpace`.
- `TVector`: Type of the entries vector.
- `FEType`: Type of the finite element.
- `APT`: Assembly type for the finite element.

# Fields
- `name::String`: Name or label for this block (for identification or debugging).
- `FES::FESpace{Tv, Ti, FEType, APT}`: The finite element space associated with this block.
- `offset::Int`: Global offset (start index in the global vector).
- `last_index::Int`: Global end index (inclusive).
- `entries::TVector`: Reference to the global coefficient array (shared with the parent `FEVector`).

# Usage
`FEVectorBlock` is typically created internally by `FEVector` constructors and provides efficient access to the coefficients for a particular FE space. Supports standard array operations (`getindex`, `setindex!`, `size`, `length`, etc.) and can be used for block-wise assembly, extraction, and manipulation.
"""
struct FEVectorBlock{T, Tv, Ti, TVector <: AbstractArray{T, 1}, FEType, APT} <: AbstractArray{T, 1}
    name::String
    FES::FESpace{Tv, Ti, FEType, APT}
    offset::Int
    last_index::Int
    entries::TVector # shares with parent object
end

function Base.copy(FEB::FEVectorBlock{T, Tv, Ti, TVector, FEType, APT}, entries) where {T, Tv, Ti, TVector, FEType, APT}
    return FEVectorBlock{T, Tv, Ti, TVector, FEType, APT}(deepcopy(FEB.name), copy(FEB.FES), FEB.offset, FEB.last_index, entries)
end

"""
$(TYPEDSIGNATURES)

Returns the number of components for the finite element in that block.
"""
get_ncomponents(FB::FEVectorBlock) = get_ncomponents(get_FEType(FB.FES))

"""
$(TYPEDSIGNATURES)

Custom `show` function for `FEVectorBlock` that prints some information and the view of that block.
"""
function Base.show(io::IO, ::MIME"text/plain", FEB::FEVectorBlock)
    @printf(io, "block %s [%d:%d] = ", FEB.name, FEB.offset + 1, FEB.last_index)
    return show(io, view(FEB))
end

"""
$(TYPEDEF)

A block-structured vector for storing coefficients associated with one or more finite element spaces (`FESpace`).

An `FEVector` consists of a global coefficient array subdivided into multiple `FEVectorBlock`s, each corresponding to a specific `FESpace`. This structure enables efficient block-wise access, assembly, and manipulation of solution vectors in finite element computations, especially for multi-field or mixed problems.

# Type Parameters
- `T`: Value type of the vector entries (e.g., `Float64`).
- `Tv`: Value type for the associated `FESpace`.
- `Ti`: Integer type for the associated `FESpace`.
- `TVector`: Type of the `entries` vector

# Fields
- `FEVectorBlocks::Array{FEVectorBlock{T, Tv, Ti}, 1}`: Array of blocks, each representing a segment of the global vector for a specific `FESpace`.
- `entries::TVector`: The global coefficient array, shared by all blocks.
- `tags::Vector{Any}`: Optional tags for identifying or accessing blocks (e.g., by name or symbol).

"""
struct FEVector{T, Tv, Ti, TVector <: AbstractArray{T, 1}} #<: AbstractVector{T}
    FEVectorBlocks::Array{FEVectorBlock{T, Tv, Ti, TVector}, 1}
    entries::TVector
    tags::Vector{Any}
end

function Base.copy(FEV::FEVector{T, Tv, Ti, TVector}) where {T, Tv, Ti, TVector}
    entries = deepcopy(FEV.entries)
    return FEVector{T, Tv, Ti, TVector}([copy(B, entries) for B in FEV.FEVectorBlocks], entries, [t for t in FEV.tags])
end

# overload stuff for AbstractArray{T,1} behaviour
Base.getindex(FEF::FEVector{T, Tv, Ti}, tag) where {T, Tv, Ti} = begin
    idx = findfirst(==(tag), FEF.tags)
    if idx === nothing
        error("Tag '$(tag)' not found in FEVector. Available tags: $(FEF.tags)")
    end
    FEF.FEVectorBlocks[idx]
end
Base.getindex(FEF::FEVector, i::Int) = begin
    if i < 1 || i > length(FEF.FEVectorBlocks)
        error("Index $(i) out of bounds for FEVector with $(length(FEF.FEVectorBlocks)) blocks.")
    end
    FEF.FEVectorBlocks[i]
end
Base.getindex(FEB::FEVectorBlock, i::Int) = FEB.entries[FEB.offset + i]
Base.getindex(FEB::FEVectorBlock, i::AbstractArray) = FEB.entries[FEB.offset .+ i]
Base.getindex(FEB::FEVectorBlock, ::Colon) = FEB.entries[(FEB.offset + 1):FEB.last_index]
Base.setindex!(FEB::FEVectorBlock, v, i::Int) = (FEB.entries[FEB.offset + i] = v)
Base.setindex!(FEB::FEVectorBlock, v, ::Colon) = (FEB.entries[(FEB.offset + 1):FEB.last_index] = v)
Base.setindex!(FEB::FEVectorBlock, v, i::AbstractArray) = (FEB.entries[FEB.offset .+ i] = v)
Base.eltype(::FEVector{T}) where {T} = T
Base.size(FEF::FEVector) = size(FEF.FEVectorBlocks)
Base.size(FEB::FEVectorBlock) = FEB.last_index - FEB.offset
Base.first(FEB::FEVectorBlock) = FEB.offset + 1
Base.last(FEB::FEVectorBlock) = FEB.last_index
Base.iterate(FEV::FEVector) = iterate(FEV.FEVectorBlocks)
Base.iterate(FEV::FEVector, state) = iterate(FEV.FEVectorBlocks, state)


"""
$(TYPEDSIGNATURES)

Returns a view of the part of the full FEVector that coressponds to the block. 
"""
Base.view(FEB::FEVectorBlock) = view(FEB.entries, (FEB.offset + 1):FEB.last_index)


"""
$(TYPEDSIGNATURES)

Returns a view of a slice of the FEVectorBlock, specified by local indices `inds` (which can be an integer, a range, or an array of indices).
The indices are relative to the block (i.e., `1` corresponds to the first entry of the block).

# Arguments
- `FEB::FEVectorBlock`: The FEVectorBlock to view.
- `inds`: Indices relative to the block (e.g., `1:30`, `[2,4,6]`).

# Returns
- A view into the underlying entries array for the specified slice.
"""
Base.view(FEB::FEVectorBlock, inds::Union{Integer, AbstractArray{<:Integer}, UnitRange{<:Integer}}) = view(FEB.entries, FEB.offset .+ inds)


function LinearAlgebra.norm(FEV::FEVector, p::Real = 2)
    return norm(FEV.entries, p)
end

function LinearAlgebra.norm(FEV::FEVectorBlock, p::Real = 2)
    return norm(view(FEV), p)
end

"""
$(TYPEDSIGNATURES)

Returns a vector with the individual norms of all blocks.
"""
function norms(FEV::FEVector{T}, p::Real = 2) where {T}
    norms = zeros(T, length(FEV))
    for j in 1:length(FEV)
        norms[j] = norm(view(FEV[j]), p)
    end
    return norms
end


"""
$(TYPEDSIGNATURES)

Returns the vector of FEspaces for the blocks of the given FEVector.
"""
function FESpaces(FEV::FEVector{T, Tv, Ti}) where {T, Tv, Ti}
    FEs = Array{FESpace{Tv, Ti}, 1}([])
    for j in 1:length(FEV.FEVectorBlocks)
        push!(FEs, FEV.FEVectorBlocks[j].FES)
    end
    return FEs
end


"""
$(TYPEDSIGNATURES)

Custom `length` function for `FEVector` that gives the number of defined FEMatrixBlocks in it
"""
Base.length(FEF::FEVector) = length(FEF.FEVectorBlocks)

"""
$(TYPEDSIGNATURES)

Custom `length` function for `FEVectorBlock` that gives the coressponding number of degrees of freedoms of the associated FESpace
"""
Base.length(FEB::FEVectorBlock) = FEB.last_index - FEB.offset

"""
````
FEVector{T}(FES; name = nothing, tags = nothing, kwargs...) where T <: Real
````

Constructs an `FEVector` for storing coefficients associated with one or more finite element spaces (`FESpace`).

- If `FES` is a single `FESpace`, the resulting `FEVector` contains one block.
- If `FES` is a vector of `FESpace` objects, the resulting `FEVector` is block-structured, with one block per space.

Optionally, you can assign a name (as a `String` for all blocks, or a vector of `String` for each block) and/or tags (as an array of any type) to the blocks for identification and access.

# Arguments
- `FES::FESpace` or `FES::Vector{<:FESpace}`: The finite element space(s) for which to create the vector.

# Keyword Arguments
- `entries`: Optional array of coefficients. If not provided, a zero vector of appropriate length is created.
- `name`: Name for the vector or for each block (default: `nothing` causes auto naming by index or tag).
- `tags`: Array of tags for the blocks (default: `[]`, i.e. block access only by index).

# Returns
- An `FEVector` object with one or more `FEVectorBlock`s, each corresponding to a given `FESpace`.

"""
function FEVector(FES::FESpace{Tv, Ti, FEType, APT}; kwargs...) where {Tv, Ti, FEType, APT}
    return FEVector{Float64}([FES]; kwargs...)
end
function FEVector{T}(FES::FESpace{Tv, Ti, FEType, APT}; kwargs...) where {T, Tv, Ti, FEType, APT}
    return FEVector{T}([FES]; kwargs...)
end
function FEVector(FES::Array{<:FESpace{Tv, Ti}, 1}; kwargs...) where {Tv, Ti}
    return FEVector{Float64}(FES; kwargs...)
end

# main constructor
function FEVector{T}(FES::Array{<:FESpace{Tv, Ti}, 1}; entries = nothing, name = nothing, tags = [], kwargs...) where {T, Tv, Ti}
    if name === nothing
        if length(tags) >= length(FES)
            names = ["$(tags[j])" for j in 1:length(FES)]
        else
            names = ["#$j" for j in 1:length(FES)]
        end
    elseif typeof(name) == String
        names = Array{String, 1}(undef, length(FES))
        for j in 1:length(FES)
            names[j] = name * (length(FES) == 1 ? "" : " [#$j]")
        end
    else
        names = name
    end
    @assert length(names) == length(FES)
    @debug "Creating FEVector mit blocks $((p -> p.name).(FES))"
    ndofs = 0
    for j in 1:length(FES)
        ndofs += FES[j].ndofs
    end
    if entries === nothing
        entries = zeros(T, ndofs)
    else
        @assert length(entries) == ndofs "length of given entries does not match number of dofs in given FESpace(s)"
    end
    Blocks = Array{FEVectorBlock{T, Tv, Ti}, 1}(undef, length(FES))
    offset = 0
    for j in 1:length(FES)
        Blocks[j] = FEVectorBlock{T, Tv, Ti, typeof(entries), eltype(FES[j]), assemblytype(FES[j])}(names[j], FES[j], offset, offset + FES[j].ndofs, entries)
        offset += FES[j].ndofs
    end
    return FEVector{T, Tv, Ti, typeof(entries)}(Blocks, entries, tags)
end


"""
$(TYPEDSIGNATURES)

Custom `show` function for `FEVector` that prints some information on its blocks.
"""
function Base.show(io::IO, FEF::FEVector)
    if length(FEF) == 0
        println(io, "FEVector is empty.")
        return
    end
    println(io, "\nFEVector information")
    println(io, "====================")
    @printf(io, "   block | starts |  ends  | length |     min      /     max      | FEType           | (name/tag)\n")
    for j in 1:length(FEF)
        ext = extrema(view(FEF[j]))
        name = FEF[j].FES.name
        tag = length(FEF.tags) >= j ? FEF.tags[j] : FEF[j].name
        @printf(
            io, " [%5d] | %6d | %6d | %6d | %12.4e / %12.4e | %-16s | %s\n",
            j, FEF[j].offset + 1, FEF[j].last_index, FEF[j].FES.ndofs, ext[1], ext[2], name, tag
        )
    end
    total_dofs = sum(FEF[j].FES.ndofs for j in 1:length(FEF))
    println(io, "\n total size = $total_dofs")
    return nothing
end


"""
$(TYPEDSIGNATURES)

Overloaded `append` function for `FEVector` that adds a FEVectorBlock at the end.
"""
function Base.append!(FEF::FEVector{T}, FES::FESpace{Tv, Ti, FEType, APT}; name = "", tag = nothing) where {T, Tv, Ti, FEType, APT}
    append!(FEF.entries, zeros(T, FES.ndofs))
    newBlock = FEVectorBlock{T, Tv, Ti, FEType, APT}(name, FES, FEF.FEVectorBlocks[end].last_index, FEF.FEVectorBlocks[end].last_index + FES.ndofs, FEF.entries)
    push!(FEF.FEVectorBlocks, newBlock)
    if tag !== nothing
        push!(FEF.tags, tag)
    end
    return length(FEF)
end

"""
$(TYPEDSIGNATURES)

Overloaded `fill` function for `FEVectorBlock` (only fills the block, not the complete FEVector).
"""
function Base.fill!(b::FEVectorBlock, value)
    fill!(view(b), value)
    return nothing
end


"""
$(TYPEDSIGNATURES)

Adds FEVectorBlock b to FEVectorBlock a.
"""
function addblock!(a::FEVectorBlock, b::FEVectorBlock; factor = 1)
    addblock!(a, b.entries; factor = factor, offset = b.offset)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Adds Array b to FEVectorBlock a.
"""
function addblock!(a::FEVectorBlock, b::AbstractVector; factor = 1, offset = 0)
    aoffset::Int = a.offset
    for j in 1:length(a)
        a.entries[aoffset + j] += b[j + offset] * factor
    end
    return nothing
end


"""
$(TYPEDSIGNATURES)

Scalar product between two FEVEctorBlocks.
"""
function LinearAlgebra.dot(a::FEVectorBlock{T}, b::FEVectorBlock{T}) where {T}
    return dot(view(a), view(b))
end

"""
$(TYPEDSIGNATURES)

Convert an `FEVector` to a standard Julia array containing all coefficients.
"""
Base.Array(FEV::FEVector) = copy(FEV.entries)

"""
$(TYPEDSIGNATURES)

Convert an `FEVectorBlock` to a standard Julia array containing the coefficients for that block.
"""
Base.Array(FEB::FEVectorBlock) = copy(FEB.entries[(FEB.offset + 1):FEB.last_index])
``````
