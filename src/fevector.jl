################
### FEVector ###
################
#
# used to store coefficients for FESpaces and can have several blocks of different FESpaces
# acts like an AbstractArray{T,1}

"""
$(TYPEDEF)

block of an FEVector that carries coefficients for an associated FESpace and can be assigned as an AbstractArray (getindex, setindex, size, length)
"""
struct FEVectorBlock{T, Tv, Ti, FEType, APT} <: AbstractArray{T, 1}
    name::String
    FES::FESpace{Tv, Ti, FEType, APT}
    offset::Int
    last_index::Int
    entries::Array{T, 1} # shares with parent object
end

function Base.copy(FEB::FEVectorBlock{T, Tv, Ti, FEType, APT}, entries) where {T, Tv, Ti, FEType, APT}
    return FEVectorBlock{T, Tv, Ti, FEType, APT}(deepcopy(FEB.name), copy(FEB.FES), FEB.offset, FEB.last_index, entries)
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
    @printf(io, "block %s [%d:%d] = ", FEB.name, FEB.offset+1, FEB.last_index)
    show(io, view(FEB))
end

"""
$(TYPEDEF)

a plain array but with an additional layer of several FEVectorBlock subdivisions each carrying coefficients for their associated FESpace.
The j-th block can be accessed by getindex(::FEVector, j) or by getindex(::FEVector, tag) if tags are associated.
The full vector can be accessed via FEVector.entries 
"""
struct FEVector{T, Tv, Ti} #<: AbstractVector{T}
    FEVectorBlocks::Array{FEVectorBlock{T, Tv, Ti}, 1}
    entries::Array{T, 1}
    tags::Vector{Any}
end

function Base.copy(FEV::FEVector{T, Tv, Ti}) where {T, Tv, Ti}
    entries = deepcopy(FEV.entries)
    return FEVector{T, Tv, Ti}([copy(B, entries) for B in FEV.FEVectorBlocks], entries, [t for t in FEV.tags])
end

# overload stuff for AbstractArray{T,1} behaviour
Base.IndexStyle(::Type{<:FEVector}) = LinearIndices
Base.IndexStyle(::Type{<:FEVectorBlock}) = LinearIndices
Base.getindex(FEF::FEVector{T, Tv, Ti}, tag) where {T, Tv, Ti} = FEF.FEVectorBlocks[findfirst(==(tag), FEF.tags)]
Base.getindex(FEF::FEVector, i::Int) = FEF.FEVectorBlocks[i]
Base.getindex(FEB::FEVectorBlock, i::Int) = FEB.entries[FEB.offset + i]
Base.getindex(FEB::FEVectorBlock, i::AbstractArray) = FEB.entries[FEB.offset .+ i]
Base.getindex(FEB::FEVectorBlock, ::Colon) = FEB.entries[(FEB.offset + 1):FEB.last_index]
Base.setindex!(FEB::FEVectorBlock, v, i::Int) = (FEB.entries[FEB.offset + i] = v)
Base.setindex!(FEB::FEVectorBlock, v, ::Colon) = (FEB.entries[(FEB.offset + 1):FEB.last_index] = v)
Base.setindex!(FEB::FEVectorBlock, v, i::AbstractArray) = (FEB.entries[FEB.offset .+ i] = v)
Base.eltype(::FEVector{T}) where {T} = T
Base.size(FEF::FEVector) = size(FEF.FEVectorBlocks)
Base.size(FEB::FEVectorBlock) = FEB.last_index - FEB.offset
Base.first(FEB::FEVectorBlock) = FEB.offset+1
Base.last(FEB::FEVectorBlock) = FEB.last_index


"""
$(TYPEDEF)

Returns a view of the part of the full FEVector that coressponds to the block. 
"""
Base.view(FEB::FEVectorBlock) = view(FEB.entries, (FEB.offset + 1):FEB.last_index)

function LinearAlgebra.norm(FEV::FEVector, p::Real = 2)
    return norm(FEV.entries, p)
end

function LinearAlgebra.norm(FEV::FEVectorBlock, p::Real = 2)
    return norm(view(FEV), p)
end

"""
$(TYPEDEF)

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

Creates FEVector that has one block if FES is a single FESpace, and a blockwise FEVector if FES is a vector of FESpaces.
Optionally a name for the vector (as a String) or each of the blocks (as a vector of Strings), or tags (as an Array{Any})
for the blocks can be specified.
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
        Blocks[j] = FEVectorBlock{T, Tv, Ti, eltype(FES[j]), assemblytype(FES[j])}(names[j], FES[j], offset, offset + FES[j].ndofs, entries)
        offset += FES[j].ndofs
    end
    return FEVector{T, Tv, Ti}(Blocks, entries, tags)
end


"""
$(TYPEDSIGNATURES)

Custom `show` function for `FEVector` that prints some information on its blocks.
"""
function Base.show(io::IO, FEF::FEVector)
    println(io, "\nFEVector information")
    println(io, "====================")
    print(io, "   block  |  starts  |   ends   |  length  |     min  /  max    \t| FEType \t\t (name/tag)")
    for j in 1:length(FEF)
        @printf(io, "\n [%5d]  |", j)
        @printf(io, "  %6d  |", FEF[j].offset+1)
        @printf(io, "  %6d  |", FEF[j].last_index)
        @printf(io, "  %6d  |", FEF[j].FES.ndofs)
        ext = extrema(view(FEF[j]))
        @printf(io, " %.2e/%.2e  \t|", ext[1], ext[2])
        if length(FEF.tags) >= j
            @printf(io, " %s  \t (%s)", FEF[j].FES.name, FEF.tags[j])
        else
            @printf(io, " %s  \t (%s)", FEF[j].FES.name, FEF[j].name)
        end
    end
    return
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
