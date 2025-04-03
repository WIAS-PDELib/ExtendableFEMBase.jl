"""
````
abstract type H1PK{ncomponents,edim,order} <: AbstractH1FiniteElement where {ncomponents<:Int,edim<:Int,order<:Int}
````

Continuous piecewise polynomials of arbitrary order >= 1 with ncomponents components in edim space dimensions.

allowed ElementGeometries:
- Edge1D
- Triangle2D
"""
abstract type H1Pk{ncomponents, edim, order} <: AbstractH1FiniteElement where {ncomponents <: Int, edim <: Int, order <: Int} end
H1Pk(ncomponents::Int, edim = ncomponents, order = 2) = H1Pk{ncomponents, edim, order}
H1Pk(; ncomponents::Int, edim = ncomponents, order = 2) = H1Pk{ncomponents, edim, order}


function Base.show(io::Core.IO, ::Type{<:H1Pk{ncomponents, edim, order}}) where {ncomponents, edim, order}
    return print(io, "H1Pk{$ncomponents,$edim,$order}")
end

get_ncomponents(FEType::Type{<:H1Pk}) = FEType.parameters[1]
get_edim(FEType::Type{<:H1Pk}) = FEType.parameters[2]

## number of dofs on different geometries
get_ndofs(::Type{<:AssemblyType}, FEType::Type{H1Pk{n, e, order}}, EG::Type{<:AbstractElementGeometry0D}) where {n, e, order} = n
get_ndofs(::Type{<:AssemblyType}, FEType::Type{H1Pk{n, e, order}}, EG::Type{<:AbstractElementGeometry1D}) where {n, e, order} = (1 + order) * n
get_ndofs(::Type{<:AssemblyType}, FEType::Type{H1Pk{n, e, order}}, EG::Type{<:Triangle2D}) where {n, e, order} = Int(n * (2 + order) * (1 + order) / 2)

## polynomial order on different geometries
get_polynomialorder(::Type{H1Pk{n, e, order}}, ::Type{<:AbstractElementGeometry}) where {n, e, order} = order

## dofmap patterns on different geometries
get_dofmap_pattern(::Type{H1Pk{n, e, order}}, ::Type{<:CellDofs}, EG::Type{<:Triangle2D}) where {n, e, order} = (order == 1) ? "N1" : ((order == 2) ? "N1F$(order - 1)" : "N1F$(order - 1)I$(Int((order - 2) * (order - 1) / 2))")
get_dofmap_pattern(::Type{H1Pk{n, e, order}}, ::Type{<:CellDofs}, EG::Type{<:AbstractElementGeometry1D}) where {n, e, order} = (order == 1) ? "N1" : "N1I$(order - 1)"
get_dofmap_pattern(::Type{H1Pk{n, e, order}}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry0D}) where {n, e, order} = "N1"
get_dofmap_pattern(::Type{H1Pk{n, e, order}}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry1D}) where {n, e, order} = (order == 1) ? "N1" : "N1I$(order - 1)"

## on which geometries it is defined
isdefined(FEType::Type{<:H1Pk}, ::Type{<:AbstractElementGeometry1D}) = true
isdefined(FEType::Type{<:H1Pk}, ::Type{<:Triangle2D}) = true

## offset for interior dofs
interior_dofs_offset(::Type{<:AssemblyType}, ::Type{<:H1Pk}, ::Type{Edge1D}) = 2
interior_dofs_offset(::Type{<:AssemblyType}, ::Type{H1Pk{n, e, o}}, ::Type{Triangle2D}) where {n, e, o} = 3 * o

## standard interpolation at nodes = node evaluation
function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, H1Pk{ncomponents, edim, order}, APT}, ::Type{AT_NODES}, exact_function!; items = [], kwargs...) where {ncomponents, edim, order, Tv, Ti, APT}
    coffset = size(FE.dofgrid[Coordinates], 2)
    if edim == 1
        coffset += (order - 1) * num_sources(FE.dofgrid[CellNodes])
    elseif edim == 2
        if order > 1
            coffset += (order - 1) * num_sources(FE.dofgrid[FaceNodes])
            if order > 2
                coffset += Int(((order - 2) * (order - 1) / 2)) * num_sources(FE.dofgrid[CellNodes])
            end
        end
    elseif edim == 3
        # coffset += 2*num_sources(FE.dofgrid[EdgeNodes]) + num_sources(FE.dofgrid[FaceNodes]) + num_sources(FE.dofgrid[CellNodes])
    end

    return point_evaluation!(Target, FE, AT_NODES, exact_function!; items = items, component_offset = coffset, kwargs...)
end

## standard interpolation on edges
function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, H1Pk{ncomponents, edim, order}, APT}, ::Type{ON_EDGES}, exact_function!; items = [], kwargs...) where {ncomponents, edim, order, Tv, Ti, APT}
    # edim = get_edim(FEType)
    # if edim == 3
    #     # delegate edge nodes to node interpolation
    #     subitems = slice(FE.dofgrid[EdgeNodes], items)
    #     interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, time = time)

    #     # perform edge mean interpolation
    #     ensure_edge_moments!(Target, FE, ON_EDGES, exact_function!; order = 1, items = items, time = time)
    # end
end

## standard interpolation on faces = preserve face moments up to order-2
function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, H1Pk{ncomponents, edim, order}, APT}, ::Type{ON_FACES}, exact_function!; items = [], kwargs...) where {ncomponents, edim, order, Tv, Ti, APT}
    # edim = get_edim(FEType)
    return if edim == 2
        # delegate face nodes to node interpolation
        subitems = slice(FE.dofgrid[FaceNodes], items)
        interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)

        # preserve moments on face
        if order > 1
            ensure_moments!(Target, FE, ON_FACES, exact_function!; order = order - 2, items = items, kwargs...)
        end
        # elseif edim == 3
        #     # delegate face edges to edge interpolation
        #     subitems = slice(FE.dofgrid[FaceEdges], items)
        #     interpolate!(Target, FE, ON_EDGES, exact_function!; items = subitems, time = time)
    elseif edim == 1
        # delegate face nodes to node interpolation
        subitems = slice(FE.dofgrid[FaceNodes], items)
        interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)
    end
end

## standard interpolation on cells = preserve cell moments up to order-3
function ExtendableGrids.interpolate!(Target, FE::FESpace{Tv, Ti, H1Pk{ncomponents, edim, order}, APT}, ::Type{ON_CELLS}, exact_function!; items = [], kwargs...) where {ncomponents, edim, order, Tv, Ti, APT}
    return if edim == 2
        # delegate cell faces to face interpolation
        subitems = slice(FE.dofgrid[CellFaces], items)
        interpolate!(Target, FE, ON_FACES, exact_function!; items = subitems, kwargs...)

        # fix cell bubble value by preserving integral mean
        if order > 2
            ensure_moments!(Target, FE, ON_CELLS, exact_function!; order = order - 3, items = items, kwargs...)
        end
        # elseif edim == 3
        #     # todo
        #     # delegate cell edges to edge interpolation
        #     subitems = slice(FE.dofgrid[CellEdges], items)
        #     interpolate!(Target, FE, ON_EDGES, exact_function!; items = subitems, time = time)

        #     # fix face means

        #     # fix cell bubble value by preserving integral mean
        #     ensure_cell_moments!(Target, FE, exact_function!; facedofs = 1, edgedofs = 2, items = items, time = time)
    elseif edim == 1
        # delegate cell nodes to node interpolation
        subitems = slice(FE.dofgrid[CellNodes], items)
        interpolate!(Target, FE, AT_NODES, exact_function!; items = subitems, kwargs...)

        # preserve cell integral
        if order > 1
            ensure_moments!(Target, FE, ON_CELLS, exact_function!; order = order - 2, items = items, kwargs...)
        end
    end
end


## basis at Vertex0D
function get_basis(::Type{<:AssemblyType}, ::Type{H1Pk{ncomponents, edim, order}}, ::Type{<:Vertex0D}) where {ncomponents, edim, order}
    return function closure(refbasis, xref)
        for k in 1:ncomponents
            refbasis[k, k] = 1
        end
        return
    end
end

## basis on Edge1D
function get_basis(::Type{<:AssemblyType}, ::Type{H1Pk{ncomponents, edim, order}}, ::Type{<:Edge1D}) where {ncomponents, edim, order}
    coeffs::Array{Rational{Int}, 1} = (0 // 1):(1 // order):(1 // 1)
    # node functions first, then interior functions
    ordering::Array{Int, 1} = [1, order + 1]
    for j in 2:order
        push!(ordering, j)
    end
    # precalculate scaling factors
    factors::Array{Rational{Int}, 1} = ones(Rational{Int}, order + 1)
    for j in 1:length(ordering), k in 1:(order + 1)
        if k != ordering[j]
            factors[j] *= coeffs[ordering[j]] - coeffs[k]
        end
    end
    function closure(refbasis, xref)
        fill!(refbasis, 0)
        for j in 1:length(ordering)
            ## build basis function that is 1 at x = coeffs[order[j]] and 0 at other positions
            refbasis[j, 1] = 1
            for k in 1:(order + 1)
                if k != ordering[j]
                    refbasis[j, 1] *= (xref[1] - coeffs[k])
                end
            end
            refbasis[j, 1] /= factors[j]
        end
        # copy to other components
        for j in 1:(order + 1), k in 2:ncomponents
            refbasis[(order + 1) * (k - 1) + j, k] = refbasis[j, 1]
        end
        return
    end

    return closure
end

## basis on Triangle2D
function get_basis(::Type{<:AssemblyType}, FEType::Type{H1Pk{ncomponents, edim, order}}, EG::Type{<:Triangle2D}) where {ncomponents, edim, order}
    coeffs::Array{Rational{Int}, 1} = (0 // 1):(1 // order):(1 // 1)
    ndofs = get_ndofs(ON_CELLS, H1Pk{1, edim, order}, EG) # dofs of one component
    # precalculate scaling factors for node and face dofs
    factor_node::Rational{Int} = prod(coeffs[2:end])
    factors_face = ones(Rational{Int}, order - 1)
    for j in 2:order, k in 1:(order + 1)
        if k < j
            factors_face[j - 1] *= coeffs[j] - coeffs[k]
        elseif k > j
            factors_face[j - 1] *= coeffs[k] - coeffs[j]
        end
    end
    if order > 3 # use recursion to fill the interior dofs (+multiplication with cell bubble)
        interior_basis = get_basis(ON_CELLS, H1Pk{1, edim, order - 3}, Triangle2D)
        # todo: scaling factors for interior dofs (but may be omitted)
    end
    function closure(refbasis, xref)
        fill!(refbasis, 0)
        # store first nodal basis function (overwritten later by last basis function)
        refbasis[end] = 1 // 1 - xref[1] - xref[2]
        # nodal functions
        for j in 1:3
            refbasis[j, 1] = 1 // 1
            if j == 1
                for k in 1:order
                    refbasis[j, 1] *= (refbasis[end] - coeffs[k])
                end
            else
                for k in 1:order
                    refbasis[j, 1] *= (xref[j - 1] - coeffs[k])
                end
            end
            refbasis[j, 1] /= factor_node
        end
        # edge basis functions
        if order > 1
            for k in 1:(order - 1)
                # on each face find basis function that is 1 at s = k//order

                # first face (nodes [1,2])
                refbasis[3 + k, 1] = refbasis[end] * xref[1] / factors_face[k]
                if order > 2
                    for m in 1:(order - 1)
                        if m > k
                            refbasis[3 + k, 1] *= (refbasis[end] - (1 - coeffs[m + 1]))
                        elseif m < k
                            refbasis[3 + k, 1] *= (xref[1] - coeffs[m + 1])
                        end
                    end
                end

                # second face (nodes [2,3])
                refbasis[3 + (order - 1) + k, 1] = xref[1] * xref[2] / factors_face[k]
                if order > 2
                    for m in 1:(order - 1)
                        if m > k
                            refbasis[3 + (order - 1) + k, 1] *= (xref[1] - (1 - coeffs[m + 1]))
                        elseif m < k
                            refbasis[3 + (order - 1) + k, 1] *= (xref[2] - coeffs[m + 1])
                        end
                    end
                end

                # third face (nodes [3,1])
                refbasis[3 + 2 * (order - 1) + k, 1] = xref[2] * refbasis[end] / factors_face[k]
                if order > 2
                    for m in 1:(order - 1)
                        if m > k
                            refbasis[3 + 2 * (order - 1) + k, 1] *= (xref[2] - (1 - coeffs[m + 1]))
                        elseif m < k
                            refbasis[3 + 2 * (order - 1) + k, 1] *= (refbasis[end] - coeffs[m + 1])
                        end
                    end
                end
            end
        end
        # interior basis functions
        if order == 3
            refbasis[3 * order + 1, 1] = refbasis[end] * xref[1] * xref[2] * 27 // 1
        elseif order == 4
            refbasis[3 * order + 1, 1] = refbasis[end] * xref[1] * xref[2] * (refbasis[end] - 1 // 4) * 128 // 1
            refbasis[3 * order + 2, 1] = refbasis[end] * xref[1] * xref[2] * (xref[1] - 1 // 4) * 128 // 1
            refbasis[3 * order + 3, 1] = refbasis[end] * xref[1] * xref[2] * (xref[2] - 1 // 4) * 128 // 1
        elseif order == 5
            refbasis[3 * order + 1, 1] = refbasis[end] * xref[1] * xref[2] * (refbasis[end] - 1 // 5) * (refbasis[end] - 2 // 5) * 3125 // 6
            refbasis[3 * order + 2, 1] = refbasis[end] * xref[1] * xref[2] * (xref[1] - 1 // 5) * (xref[1] - 2 // 5) * 3125 // 6
            refbasis[3 * order + 3, 1] = refbasis[end] * xref[1] * xref[2] * (xref[2] - 1 // 5) * (xref[2] - 2 // 5) * 3125 // 6
            refbasis[3 * order + 4, 1] = refbasis[end] * xref[1] * xref[2] * (refbasis[end] - 1 // 5) * (xref[1] - 1 // 5) * 3125 // 4
            refbasis[3 * order + 5, 1] = refbasis[end] * xref[1] * xref[2] * (xref[1] - 1 // 5) * (xref[2] - 1 // 5) * 3125 // 4
            refbasis[3 * order + 6, 1] = refbasis[end] * xref[1] * xref[2] * (xref[2] - 1 // 5) * (refbasis[end] - 1 // 5) * 3125 // 4
        elseif order > 5
            interior_basis(view(refbasis, (3 * order + 1):(ncomponents * ndofs), :), xref)
            for k in (3 * order + 1):ndofs
                refbasis[k, 1] *= (1 - xref[1] - xref[2]) * xref[1] * xref[2] * 27 // 1
            end
        end

        # copy to other components
        for j in 1:ndofs, k in 2:ncomponents
            refbasis[(k - 1) * ndofs + j, k] = refbasis[j, 1]
        end
        return
    end

    return closure
end


## if order > 2 on Triangl2D we need to change the ordering of the face dofs on faces that have a negative orientation sign
function get_basissubset(::Type{ON_CELLS}, FE::FESpace{Tv, Ti, H1Pk{ncomponents, edim, order}, APT}, EG::Type{<:Triangle2D}, xgrid) where {ncomponents, edim, order, Tv, Ti, APT}
    if order < 3
        return NothingFunction # no reordering needed
    end
    xCellFaceSigns = xgrid[CellFaceSigns]
    nfaces::Int = num_faces(EG)
    ndofs_for_f::Int = order - 1
    ndofs_for_c = get_ndofs(ON_CELLS, H1Pk{1, edim, order}, EG)
    return function closure(subset_ids::Array{Int, 1}, cell)
        if order > 2
            for j in 1:nfaces
                if xCellFaceSigns[j, cell] != 1
                    for c in 1:ncomponents, k in 1:ndofs_for_f
                        subset_ids[(c - 1) * ndofs_for_c + 3 + (j - 1) * ndofs_for_f + k] = (c - 1) * ndofs_for_c + 3 + (j - 1) * ndofs_for_f + (1 + ndofs_for_f - k)
                    end
                else
                    for c in 1:ncomponents, k in 1:ndofs_for_f
                        subset_ids[(c - 1) * ndofs_for_c + 3 + (j - 1) * ndofs_for_f + k] = (c - 1) * ndofs_for_c + 3 + (j - 1) * ndofs_for_f + k
                    end
                end
            end
        end
        return nothing
    end
end
