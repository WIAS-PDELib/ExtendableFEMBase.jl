################
### FEMatrix ###
################
#
# used to store (sparse) matrix representations of PDEOperators for FESpaces
# and can have several matrix blocks (FEMatrixBlock) of different FESpaces
# each matrix block acts like an AbstractArray{T,2}

"""
$(TYPEDEF)

block of an FEMatrix that carries coefficients for an associated pair of FESpaces and can be assigned as an two-dimensional AbstractArray (getindex, setindex, size)
"""
struct FEMatrixBlock{TvM, TiM, TvG, TiG, FETypeX, FETypeY, APTX, APTY} <: AbstractArray{TvM, 2}
    name::String
    FES::FESpace{TvG, TiG, FETypeX, APTX}
    FESY::FESpace{TvG, TiG, FETypeY, APTY}
    offset::Int64
    offsetY::Int64
    last_index::Int64
    last_indexY::Int64
    entries::AbstractSparseMatrix{TvM, TiM} # shares with parent object
end

function Base.copy(FEMB::FEMatrixBlock{TvM, TiM, TvG, TiG, FETypeX, FETypeY, APTX, APTY}, entries) where {TvM, TiM, TvG, TiG, FETypeX, FETypeY, APTX, APTY}
    return FEMatrixBlock{TvM, TiM, TvG, TiG, FETypeX, FETypeY, APTX, APTY}(deepcopy(FEMB.name), copy(FEMB.FES), copy(FEMB.FESY), FEMB.offset, FEMB.offsetY, FEMB.last_index, FEMB.last_indexY, entries)
end

"""
$(TYPEDSIGNATURES)

Custom `show` function for `FEMatrixBlock` that prints its coordinates and the name.
"""
function Base.show(io::IO, FEB::FEMatrixBlock)
    @printf(io, "[%d:%d,%d:%d]: %s", FEB.offset+1, FEB.last_index, FEB.offsetY+1, FEB.last_indexY, FEB.name)
    return nothing
end


"""
$(TYPEDEF)

an AbstractMatrix (e.g. an ExtendableSparseMatrix) with an additional layer of several FEMatrixBlock subdivisions each carrying coefficients for their associated pair of FESpaces
"""
struct FEMatrix{TvM, TiM, TvG, TiG, nbrow, nbcol, nbtotal} <: AbstractSparseMatrix{TvM, TiM}
    FEMatrixBlocks::Array{FEMatrixBlock{TvM, TiM, TvG, TiG}, 1}
    entries::AbstractSparseMatrix{TvM, TiM}
    tags::Matrix{Any}
end

function Base.copy(FEV::FEMatrix{TvM, TiM, TvG, TiG, nbrow, nbcol, nbtotal}) where {TvM, TiM, TvG, TiG, nbrow, nbcol, nbtotal}
    entries = deepcopy(FEV.entries)
    return FEVector{TvM, TiM, TvG, TiG, nbrow, nbcol, nbtotal}([copy(B, entries) for B in FEV.FEMatrixBlocks], entries)
end


"""
$(TYPEDEF)

adds value v to matrix at position (i,j) if it is nonzero
"""
@inline function _addnz(ESM::AbstractExtendableSparseMatrixCSC, i, j, v::Tv, fac, part = 1) where {Tv}
    return if v != zero(Tv)
        rawupdateindex!(ESM, +, v * fac, i, j, part)
    end
end
@inline function _addnz(FEB::FEMatrixBlock, i, j, v::Tv, fac, part = 1) where {Tv}
    return if v != zero(Tv)
        rawupdateindex!(FEB.entries, +, v * fac, FEB.offset + i, FEB.offsetY + j, part)
    end
end

"""
$(TYPEDEF)

pre-allocates the expected pattern for the default dofmaps for the AssemblyType
"""
function apply_nonzero_pattern!(B::FEMatrixBlock, AT::Type{<:AssemblyType})
    dofmapX = Dofmap4AssemblyType(B.FESX, AT)
    dofmapY = Dofmap4AssemblyType(B.FESY, AT)
    @assert num_sources(dofmapX) == num_sources(dofmapY)
    for item in 1:num_sources(dofmapX)
        for j in 1:num_targets(dofmapX, item), k in 1:num_targets(dofmapY, item)
            rawupdateindex!(B.entries, +, 0, B.offset + dofmapX[j, item], B.offsetY + dofmapY[k, item])
        end
    end
    return
end

Base.getindex(FEF::FEMatrix, i) = FEF.FEMatrixBlocks[i]
Base.getindex(FEF::FEMatrix{TvM, TiM, TvG, TiG, nbrow, nbcol, nbtotal}, tagX, tagY) where {TvM, TiM, TvG, TiG, nbrow, nbcol, nbtotal} = (index = findfirst(==((tagX, tagY)), FEF.tags); return FEF.FEMatrixBlocks[(index[1] - 1) * nbcol + index[2]])
Base.getindex(FEF::FEMatrix{TvM, TiM, TvG, TiG, nbrow, nbcol, nbtotal}, i::Int, j::Int) where {TvM, TiM, TvG, TiG, nbrow, nbcol, nbtotal} = FEF.FEMatrixBlocks[(i - 1) * nbcol + j]
Base.getindex(FEB::FEMatrixBlock, i::Int, j::Int) = FEB.entries[FEB.offset + i, FEB.offsetY + j]
Base.getindex(FEB::FEMatrixBlock, i::Any, j::Any) = FEB.entries[FEB.offset .+ i, FEB.offsetY .+ j]
Base.setindex!(FEB::FEMatrixBlock, v, i::Int, j::Int) = setindex!(FEB.entries, v, FEB.offset + i, FEB.offsetY + j)
Base.first(FEB::FEMatrixBlock) = (FEB.offset+1, FEB.offsetY+1)
Base.last(FEB::FEMatrixBlock) = (FEB.last_index, FEB.last_indexY)



"""
$(TYPEDSIGNATURES)

Gives the number of FEMatrixBlocks in each column.
"""
nbrows(::FEMatrix{TvM, TiM, TvG, TiG, nbrow, nbcol, nbtotal}) where {TvM, TiM, TvG, TiG, nbrow, nbcol, nbtotal} = nbrow


"""
$(TYPEDSIGNATURES)

Gives the number of FEMatrixBlocks in each row.
"""
nbcols(::FEMatrix{TvM, TiM, TvG, TiG, nbrow, nbcol, nbtotal}) where {TvM, TiM, TvG, TiG, nbrow, nbcol, nbtotal} = nbcol

"""
$(TYPEDSIGNATURES)

Custom `size` function for `FEMatrix` that gives a tuple with the number of rows and columns of the FEBlock overlay
"""
Base.size(::FEMatrix{TvM, TiM, TvG, TiG, nbrow, nbcol, nbtotal}) where {TvM, TiM, TvG, TiG, nbrow, nbcol, nbtotal} = (nbrow, nbcol)


"""
$(TYPEDSIGNATURES)

Custom `length` function for `FEMatrix` that gives the total number of defined FEMatrixBlocks in it
"""
Base.length(::FEMatrix{TvM, TiM, TvG, TiG, nbrow, nbcol, nbtotal}) where {TvM, TiM, TvG, TiG, nbrow, nbcol, nbtotal} = nbtotal

"""
$(TYPEDSIGNATURES)

Custom `size` function for `FEMatrixBlock` that gives a tuple with the size of the block (that coressponds to the number of degrees of freedoms in X and Y)
"""
Base.size(FEB::FEMatrixBlock) = (FEB.last_index - FEB.offset, FEB.last_indexY - FEB.offsetY)

"""
$(TYPEDSIGNATURES)

Custom `show` function for `FEMatrix` that prints some information on its blocks.
"""
function Base.show(io::IO, FEM::FEMatrix{TvM, TiM, TvG, TiG, nbrow, nbcol, nbtotal}) where {TvM, TiM, TvG, TiG, nbrow, nbcol, nbtotal}
    println(io, "\nFEMatrix information")
    println(io, "====================")
    println(io, "  block   |      starts      |       ends       |       size       |    name ")
    for j in 1:length(FEM)
        n = mod(j - 1, nbrow) + 1
        m = Int(ceil(j / nbrow))
        @printf(io, " [%2d,%2d]  |", m, n)
        @printf(io, "  (%6d,%6d) |", FEM[j].offset+1, FEM[j].offsetY+1)
        @printf(io, "  (%6d,%6d) |", FEM[j].last_index, FEM[j].last_indexY)
        @printf(io, "  (%6d,%6d) |", FEM[j].FES.ndofs, FEM[j].FESY.ndofs)
        @printf(io, " %s \n", FEM[j].name)
    end
    return println(io, "\n    nnzvals = $(length(FEM.entries.cscmatrix.nzval))")
end

"""
````
FEMatrix{TvM,TiM}(name::String, FES::FESpace{TvG,TiG,FETypeX,APTX}) where {TvG,TiG,FETypeX,APTX}
````

Creates FEMatrix with one square block (FES,FES).
"""
function FEMatrix(FES::FESpace; name = :automatic)
    return FEMatrix{Float64, Int64}(FES; name = name)
end
function FEMatrix{TvM}(FES::FESpace; name = :automatic) where {TvM}
    return FEMatrix{TvM, Int64}(FES; name = name)
end
function FEMatrix{TvM, TiM}(FES::FESpace{TvG, TiG, FETypeX, APTX}; name = :automatic) where {TvM, TiM, TvG, TiG, FETypeX, APTX}
    return FEMatrix{TvM, TiM}([FES], [FES]; name = name)
end

"""
````
FEMatrix{TvM,TiM}(FESX, FESY; name = :automatic)
````

Creates FEMatrix with one rectangular block (FESX,FESY) if FESX and FESY are single FESpaces, or
a rectangular block matrix with blocks corresponding to the entries of the FESpace vectors FESX and FESY.
Optionally a name for the matrix can be given.
"""
function FEMatrix(FESX::FESpace, FESY::FESpace; kwargs...)
    return FEMatrix{Float64, Int64}(FESX, FESY; kwargs...)
end
function FEMatrix(FESX::Vector{<:FESpace}, FESY::Vector{<:FESpace}; kwargs...)
    return FEMatrix{Float64, Int64}(FESX, FESY; kwargs...)
end
function FEMatrix{TvM}(FESX::FESpace, FESY::FESpace; kwargs...) where {TvM}
    return FEMatrix{TvM, Int64}(FESX, FESY; kwargs...)
end
function FEMatrix{TvM, TiM}(FESX::FESpace, FESY::FESpace; kwargs...) where {TvM, TiM}
    return FEMatrix{TvM, TiM}([FESX], [FESY]; kwargs...)
end
function FEMatrix(FES::Array{<:FESpace{TvG, TiG}, 1}; kwargs...) where {TvG, TiG}
    return FEMatrix{Float64, Int64}(FES; kwargs...)
end
function FEMatrix{TvM}(FES::Array{<:FESpace{TvG, TiG}, 1}; kwargs...) where {TvM, TvG, TiG}
    return FEMatrix{TvM, Int64}(FES, FES; kwargs...)
end
function FEMatrix{TvM, TiM}(FES::Array{<:FESpace{TvG, TiG}, 1}; kwargs...) where {TvM, TiM, TvG, TiG}
    return FEMatrix{TvM, TiM}(FES, FES; kwargs...)
end

"""
````
FEMatrix{TvM,TiM}(FESX, FESY; name = :automatic)
````

Creates an FEMatrix with blocks coressponding to the ndofs of FESX (rows) and FESY (columns).
"""
function FEMatrix{TvM, TiM}(FESX::Array{<:FESpace{TvG, TiG}, 1}, FESY::Array{<:FESpace{TvG, TiG}, 1}; entries = nothing, name = :automatic, tags = nothing, tagsX = tags, tagsY = tagsX, npartitions = 1, kwargs...) where {TvM, TiM, TvG, TiG}
    ndofsX, ndofsY = 0, 0
    for j in 1:length(FESX)
        ndofsX += FESX[j].ndofs
    end
    for j in 1:length(FESY)
        ndofsY += FESY[j].ndofs
    end
    if entries === nothing
        if npartitions == 1
            entries = ExtendableSparseMatrixCSC{TvM, TiM}(ndofsX, ndofsY)
        elseif npartitions > 1
            entries = MTExtendableSparseMatrixCSC{TvM, TiM}(ndofsX, ndofsY, npartitions)
        end
    else
        @assert size(entries) == (ndofsX, ndofsY) "size of given entries not matching number of dofs in given FE space(s)"
    end

    if name === :automatic || name === nothing
        name = ""
    end

    if tagsX !== nothing
        @assert length(tagsX) == length(FESX)
    end
    if tagsY !== nothing
        @assert length(tagsY) == length(FESY)
    end

    Blocks = Array{FEMatrixBlock{TvM, TiM, TvG, TiG}, 1}(undef, length(FESX) * length(FESY))
    offset = 0
    offsetY = 0
    for j in 1:length(FESX)
        offsetY = 0
        for k in 1:length(FESY)
            if (tagsX !== nothing) && (tagsY !== nothing)
                blockname = name * " $(tagsX[j])[$(FESX[j].name)] x $(tagsY[k])[$(FESY[k].name)]"
            else
                blockname = name * " $(FESX[j].name) x $(FESY[k].name)"
            end
            Blocks[(j - 1) * length(FESY) + k] =
                FEMatrixBlock{TvM, TiM, TvG, TiG, eltype(FESX[j]), eltype(FESY[k]), assemblytype(FESX[j]), assemblytype(FESY[k])}(blockname, FESX[j], FESY[k], offset, offsetY, offset + FESX[j].ndofs, offsetY + FESY[k].ndofs, entries)
            offsetY += FESY[k].ndofs
        end
        offset += FESX[j].ndofs
    end

    if (tagsX !== nothing) && (tagsY !== nothing)
        tagmatrix = [(j, k) for j in tagsX, k in tagsY]
    else
        tagmatrix = zeros(Int, 0, 0)
    end
    return FEMatrix{TvM, TiM, TvG, TiG, length(FESX), length(FESY), length(FESX) * length(FESY)}(Blocks, entries, tagmatrix)
end

"""
$(TYPEDSIGNATURES)

Custom `fill` function for `FEMatrixBlock` (only fills the already present nzval in the block, not the complete FEMatrix).
"""
function Base.fill!(B::FEMatrixBlock{Tv, Ti}, value) where {Tv, Ti}
    cscmat::SparseMatrixCSC{Tv, Ti} = B.entries.cscmatrix
    rows::Array{Int, 1} = rowvals(cscmat)
    valsB::Array{Tv, 1} = cscmat.nzval
    for col in (B.offsetY + 1):B.last_indexY
        for r in nzrange(cscmat, col)
            if rows[r] > B.offset && rows[r] <= B.last_index
                valsB[r] = value
            end
        end
    end
    return nothing
end


"""
$(TYPEDSIGNATURES)

Adds FEMatrix/ExtendableSparseMatrix/CSCMatrix B to FEMatrix A.
"""
function add!(A::FEMatrix{Tv, Ti}, B::FEMatrix{Tv, Ti}; kwargs...) where {Tv, Ti}
    return add!(A.entries, B.entries; kwargs...)
end
function add!(AM::AbstractExtendableSparseMatrixCSC{Tv, Ti}, BM::AbstractExtendableSparseMatrixCSC{Tv, Ti}; kwargs...) where {Tv, Ti}
    return add!(AM, BM.cscmatrix; kwargs...)
end
function add!(AM::AbstractExtendableSparseMatrixCSC{Tv, Ti}, cscmat::SparseMatrixCSC{Tv, Ti}; factor = 1, rowoffset = 0, coloffset = 0, transpose::Bool = false) where {Tv, Ti}
    rows::Array{Ti, 1} = rowvals(cscmat)
    valsB::Array{Tv, 1} = cscmat.nzval
    ncols::Int = size(cscmat, 2)
    arow::Int = 0
    if transpose
        for col in 1:ncols
            for r in nzrange(cscmat, col)
                arow = rows[r]
                _addnz(AM, col + coloffset, arow + rowoffset, valsB[r], factor)
            end
        end
    else
        for col in 1:ncols
            for r in nzrange(cscmat, col)
                arow = rows[r]
                _addnz(AM, arow + rowoffset, col + coloffset, valsB[r], factor)
            end
        end
    end
    flush!(AM)
    return nothing
end


"""
$(TYPEDSIGNATURES)

Adds FEMatrixBlock B to FEMatrixBlock A.
"""
function addblock!(A::FEMatrixBlock{Tv, Ti}, B::FEMatrixBlock{Tv, Ti}; factor = 1, transpose::Bool = false) where {Tv, Ti}
    AM::AbstractExtendableSparseMatrixCSC{Tv, Ti} = A.entries
    BM::AbstractExtendableSparseMatrixCSC{Tv, Ti} = B.entries
    cscmat::SparseMatrixCSC{Tv, Ti} = BM.cscmatrix
    rows::Array{Ti, 1} = rowvals(cscmat)
    valsB::Array{Tv, 1} = cscmat.nzval
    arow::Int = 0
    acol::Int = 0
    if transpose
        for col in (B.offsetY + 1):B.last_indexY
            arow = col - B.offsetY + A.offset
            for r in nzrange(cscmat, col)
                if rows[r] > B.offset && rows[r] <= B.last_index
                    acol = rows[r] - B.offset + A.offsetY
                    ## add B[rows[r], col] to A[col, rows[r]]
                    _addnz(AM, arow, acol, valsB[r], factor)
                end
            end
        end
    else
        for col in (B.offsetY + 1):B.last_indexY
            acol = col - B.offsetY + A.offsetY
            for r in nzrange(cscmat, col)
                if rows[r] > B.offset && rows[r] <= B.last_index
                    arow = rows[r] - B.offset + A.offset
                    ## add B[rows[r], col] to A[rows[r], col]
                    _addnz(AM, arow, acol, valsB[r], factor)
                end
            end
        end
    end
    flush!(AM)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Adds ExtendableSparseMatrix B to FEMatrixBlock A.
"""
function addblock!(A::FEMatrixBlock{Tv}, B::AbstractExtendableSparseMatrixCSC{Tv, Ti}; factor = 1, transpose::Bool = false) where {Tv, Ti <: Integer}
    return addblock!(A, B.cscmatrix; factor = factor, transpose = transpose)
end


"""
$(TYPEDSIGNATURES)

Adds SparseMatrixCSC B to FEMatrixBlock A.
"""
function addblock!(A::FEMatrixBlock{Tv}, cscmat::SparseArrays.SparseMatrixCSC{Tv, Ti}; factor = 1, transpose::Bool = false) where {Tv, Ti <: Integer}
    AM::AbstractExtendableSparseMatrixCSC{Tv, Int64} = A.entries
    rows::Array{Int, 1} = rowvals(cscmat)
    valsB::Array{Tv, 1} = cscmat.nzval
    arow::Int = 0
    acol::Int = 0
    if transpose
        for col in 1:size(cscmat, 2)
            arow = col + A.offset
            for r in nzrange(cscmat, col)
                acol = rows[r] + A.offsetY
                _addnz(AM, arow, acol, valsB[r], factor)
                #A[rows[r],col] += B.cscmatrix.nzval[r] * factor
            end
        end
    else
        for col in 1:size(cscmat, 2)
            acol = col + A.offsetY
            for r in nzrange(cscmat, col)
                arow = rows[r] + A.offset
                _addnz(AM, arow, acol, valsB[r], factor)
                #A[rows[r],col] += B.cscmatrix.nzval[r] * factor
            end
        end
    end
    flush!(AM)
    return nothing
end


"""
$(TYPEDSIGNATURES)

sets penalty to the diagonal entries of fixed_dofs in A
"""
function apply_penalties!(A::AbstractExtendableSparseMatrixCSC, fixed_dofs, penalty)
    for dof in fixed_dofs
        A[dof, dof] = penalty
    end
    flush!(A)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Adds matrix-matrix product B times C to FEMatrixBlock A.
"""
function addblock_matmul!(A::FEMatrixBlock{Tv}, cscmatB::SparseMatrixCSC{Tv, Ti}, cscmatC::SparseMatrixCSC{Tv, Ti}; factor = 1, transposed::Bool = false) where {Tv, Ti}
    AM::AbstractExtendableSparseMatrixCSC{Tv, Int64} = A.entries
    rowsB::Array{Ti, 1} = rowvals(cscmatB)
    rowsC::Array{Ti, 1} = rowvals(cscmatC)
    valsB::Array{Tv, 1} = cscmatB.nzval
    valsC::Array{Tv, 1} = cscmatC.nzval
    arow::Int = 0
    acol::Int = 0
    if transposed # add (B*C)'= C'*B' to A
        for i in 1:size(cscmatC, 2)
            arow = i + A.offset
            for crow in nzrange(cscmatC, i)
                for j in nzrange(cscmatB, rowsC[crow])
                    acol = rowsB[j] + A.offsetY
                    # add b[j,crow]*c[crow,i] to a[i,j]
                    _addnz(AM, arow, acol, valsB[j] * valsC[crow], factor)
                end
            end
        end
    else # add B*C to A
        for j in 1:size(cscmatC, 2)
            acol = j + A.offsetY
            for crow in nzrange(cscmatC, j)
                for i in nzrange(cscmatB, rowsC[crow])
                    arow = rowsB[i] + A.offset
                    # add b[i,crow]*c[crow,j] to a[i,j]
                    _addnz(AM, arow, acol, valsB[i] * valsC[crow], factor)
                end
            end
        end
    end
    flush!(AM)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Adds matrix-vector product B times b (or B' times b if transposed = true) to FEVectorBlock a.
"""
function addblock_matmul!(a::FEVectorBlock{Tv}, B::FEMatrixBlock{Tv, Ti}, b::FEVectorBlock{Tv}; factor = 1, transposed::Bool = false) where {Tv, Ti}
    cscmat::SparseMatrixCSC{Tv, Ti} = B.entries.cscmatrix
    rows::Array{Ti, 1} = rowvals(cscmat)
    valsB::Array{Tv, 1} = cscmat.nzval
    row::Int = 0
    if transposed
        brow::Int = 0
        acol::Int = 0
        for col in (B.offsetY + 1):B.last_indexY
            acol = col - B.offsetY + a.offset
            for r in nzrange(cscmat, col)
                row = rows[r]
                if row > B.offset && row <= B.last_index
                    brow = row - B.offset + b.offset
                    a.entries[acol] += valsB[r] * b.entries[brow] * factor
                end
            end
        end
    else
        bcol::Int = 0
        arow::Int = 0
        for col in (B.offsetY + 1):B.last_indexY
            bcol = col - B.offsetY + b.offset
            for r in nzrange(cscmat, col)
                row = rows[r]
                if row > B.offset && row <= B.last_index
                    arow = row - B.offset + a.offset
                    a.entries[arow] += valsB[r] * b.entries[bcol] * factor
                end
            end
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Adds matrix-vector product B times b to FEVectorBlock a.
"""
function addblock_matmul!(a::AbstractVector{Tv}, B::FEMatrixBlock{Tv, Ti}, b::AbstractVector{Tv}; factor = 1, transposed::Bool = false) where {Tv, Ti}
    cscmat::SparseMatrixCSC{Tv, Ti} = B.entries.cscmatrix
    rows::Array{Ti, 1} = rowvals(cscmat)
    valsB::Array{Tv, 1} = cscmat.nzval
    bcol::Int = 0
    row::Int = 0
    arow::Int = 0
    if transposed
        for col in (B.offsetY + 1):B.last_indexY
            bcol = col - B.offsetY
            for r in nzrange(cscmat, col)
                row = rows[r]
                if row > B.offset && row <= B.last_index
                    arow = row - B.offset
                    a[bcol] += valsB[r] * b[arow] * factor
                end
            end
        end
    else
        for col in (B.offsetY + 1):B.last_indexY
            bcol = col - B.offsetY
            for r in nzrange(cscmat, col)
                row = rows[r]
                if row > B.offset && row <= B.last_index
                    arow = row - B.offset
                    a[arow] += valsB[r] * b[bcol] * factor
                end
            end
        end
    end
    return nothing
end


"""
$(TYPEDSIGNATURES)

Adds matrix-vector product B times b to FEVectorBlock a.
"""
function addblock_matmul!(a::FEVectorBlock{Tv}, B::AbstractExtendableSparseMatrixCSC{Tv, Ti}, b::FEVectorBlock{Tv}; factor = 1) where {Tv, Ti <: Integer}
    cscmat::SparseMatrixCSC{Tv, Ti} = B.cscmatrix
    rows::Array{Ti, 1} = rowvals(cscmat)
    valsB::Array{Tv, 1} = cscmat.nzval
    bcol::Int = 0
    arow::Int = 0
    for col in 1:size(B, 2)
        bcol = col + b.offset
        for r in nzrange(cscmat, col)
            arow = rows[r] + a.offset
            a.entries[arow] += valsB[r] * b.entries[bcol] * factor
        end
    end
    return nothing
end

function addblock_matmul!(a::FEVectorBlock{Tv}, cscmat::SparseMatrixCSC{Tv, Ti}, b::FEVectorBlock{Tv}; factor = 1) where {Tv, Ti <: Integer}
    rows::Array{Ti, 1} = rowvals(cscmat)
    valsB::Array{Tv, 1} = cscmat.nzval
    bcol::Int = 0
    arow::Int = 0
    for col in 1:size(cscmat, 2)
        bcol = col + b.offset
        for r in nzrange(cscmat, col)
            arow = rows[r] + a.offset
            a.entries[arow] += valsB[r] * b.entries[bcol] * factor
        end
    end
    return nothing
end


"""
$(TYPEDSIGNATURES)

Computes vector'-matrix-vector product a'*B*b.
"""
function lrmatmul(a::AbstractVector{Tv}, B::AbstractExtendableSparseMatrixCSC{Tv, Ti}, b::AbstractVector{Tv}; factor = 1) where {Tv, Ti <: Integer}
    cscmat::SparseMatrixCSC{Tv, Ti} = B.cscmatrix
    valsB::Array{Tv, 1} = cscmat.nzval
    rows::Array{Ti, 1} = rowvals(cscmat)
    result = 0.0
    for col in 1:size(B, 2)
        for r in nzrange(cscmat, col)
            result += valsB[r] * b[col] * factor * a[rows[r]]
        end
    end
    return result
end


"""
$(TYPEDSIGNATURES)

Computes vector'-matrix-vector product (a1-a2)'*B*(b1-b2).
"""
function ldrdmatmul(a1::AbstractVector{Tv}, a2::AbstractVector{Tv}, B::AbstractExtendableSparseMatrixCSC{Tv, Ti}, b1::AbstractVector{Tv}, b2::AbstractVector{Tv}; factor = 1) where {Tv, Ti <: Integer}
    cscmat::SparseMatrixCSC{Tv, Ti} = B.cscmatrix
    valsB::Array{Tv, 1} = cscmat.nzval
    rows::Array{Ti, 1} = rowvals(cscmat)
    result = 0.0
    for col in 1:size(B, 2)
        for r in nzrange(cscmat, col)
            result += valsB[r] * (b1[col] - b2[col]) * factor * (a1[rows[r]] - a2[rows[r]])
        end
    end
    return result
end


"""
$(TYPEDSIGNATURES)

Generates an ExtendableSparseMatrix from the submatrix
for the given row and col numbers
"""
function submatrix(A::AbstractExtendableSparseMatrixCSC{Tv, Ti}, srows, scols; factor = 1) where {Tv, Ti}
    cscmat::SparseMatrixCSC{Tv, Ti} = A.cscmatrix
    rows::Array{Ti, 1} = rowvals(cscmat)
    valsA = cscmat.nzval
    nrowsA = size(A,1)
    ncolsA = size(A,2)
    S = ExtendableSparseMatrix{Tv, Ti}(length(srows), length(scols))
    @assert maximum(srows) <= size(A, 1) "rows exceeds rowcount of A"
    @assert maximum(scols) <= size(A, 2) "cols exceeds colcount of A"
    ncols = length(scols)
    nrows = length(srows)
    newrows = zeros(Int, nrowsA)
    newcols = zeros(Int, ncolsA)
    newrows[srows] = 1:nrows
    newcols[scols] = 1:ncols
    minrow, maxrow = minimum(srows), maximum(srows)

    for scol in scols
        for r in nzrange(cscmat, scol)
            if newrows[rows[r]] > 0
                _addnz(S, newrows[rows[r]], newcols[scol], valsA[r], factor)
            else 
                break
            end
        end
    end
    flush!(S)
    return S
end

"""
$(TYPEDSIGNATURES)

Returns the FEMatrixBlock as an ExtendableSparseMatrix
"""
function submatrix(A::FEMatrixBlock{Tv, Ti}) where {Tv, Ti}
    return submatrix(A.entries, A.offset+1:A.last_index, A.offsetY+1:A.last_indexY)
end
