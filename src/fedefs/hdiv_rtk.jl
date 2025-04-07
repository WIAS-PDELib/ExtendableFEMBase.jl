"""
````
abstract type HDIVRTk{edim, order} <: AbstractHdivFiniteElement where {edim<:Int}
````

Hdiv-conforming vector-valued (ncomponents = edim) Raviart-Thomas space of arbitrary order.

allowed ElementGeometries:
- Triangle2D
"""
abstract type HDIVRTk{edim, order} <: AbstractHdivFiniteElement where {edim <: Int, order <: Int} end
HDIVRTk(edim::Int, order::Int) = HDIVRTk{edim, order}
HDIVRTk(; edim::Int, order::Int) = HDIVRTk{edim, order}

function Base.show(io::Core.IO, FEType::Type{<:HDIVRTk{edim, order}}) where {edim, order}
    return print(io, "HDIVRTk{$edim, $order}")
end

get_ncomponents(FEType::Type{<:HDIVRTk}) = FEType.parameters[1]
get_ndofs(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, FEType::Type{<:HDIVRTk{edim, order}}, EG::Type{<:AbstractElementGeometry1D}) where {edim, order} = order + 1
get_ndofs(::Type{ON_CELLS}, FEType::Type{<:HDIVRTk{edim, order}}, EG::Type{<:Triangle2D}) where {edim, order} = (order + 1) * num_faces(EG) + order * (order + 1)

get_polynomialorder(::Type{<:HDIVRTk{2, order}}, ::Type{<:AbstractElementGeometry1D}) where {order} = order;
get_polynomialorder(::Type{<:HDIVRTk{2, order}}, ::Type{<:AbstractElementGeometry2D}) where {order} = order + 1;

get_dofmap_pattern(FEType::Type{<:HDIVRTk{2, order}}, ::Type{CellDofs}, EG::Type{<:Triangle2D}) where {order} = "f$(order + 1)" * (order > 0 ? "i$(Int((order + 1) * (order)))" : "")
get_dofmap_pattern(FEType::Type{<:HDIVRTk{2, order}}, ::Union{Type{FaceDofs}, Type{BFaceDofs}}, EG::Type{<:AbstractElementGeometry1D}) where {order} = "i$(order + 1)"

isdefined(FEType::Type{<:HDIVRTk}, ::Type{<:Triangle2D}) = true

interior_dofs_offset(::Type{<:ON_CELLS}, ::Type{<:HDIVRTk{2, order}}, ::Type{<:Triangle2D}) where {order} = 3 * (order + 1)


function ExtendableGrids.interpolate!(Target::AbstractArray{T, 1}, FE::FESpace{Tv, Ti, HDIVRTk{edim, order}, APT}, ::Type{ON_FACES}, data; items = [], bonus_quadorder = 0, kwargs...) where {T, Tv, Ti, edim, order, APT}
    ncomponents = edim
    xFaceNormals = FE.dofgrid[FaceNormals]
    nfaces = num_sources(xFaceNormals)
    if items == []
        items = 1:nfaces
    end

    # integrate normal flux of exact_function over edges
    data_eval = zeros(T, ncomponents)
    moments_weights = basis.(ShiftedLegendre, (0:order))
    offset = [k * nfaces for k in 0:order]
    function normalflux_eval(result, qpinfo)
        data(data_eval, qpinfo)
        result[1] = dot(data_eval, view(xFaceNormals, :, qpinfo.item))
        for j in 1:order
            result[j + 1] = result[1] * moments_weights[j + 1](qpinfo.xref[1])
        end
        return
    end
    return integrate!(Target, FE.dofgrid, ON_FACES, normalflux_eval; quadorder = 2 * order + bonus_quadorder, items = items, offset = offset, kwargs...)
end

function ExtendableGrids.interpolate!(Target::AbstractArray{T, 1}, FE::FESpace{Tv, Ti, HDIVRTk{edim, order}, APT}, ::Type{ON_CELLS}, data; items = [], time = 0, bonus_quadorder = 0, kwargs...) where {T, Tv, Ti, edim, order, APT}
    # delegate cell faces to face interpolation
    subitems = slice(FE.dofgrid[CellFaces], items)
    interpolate!(Target, FE, ON_FACES, data; items = subitems, bonus_quadorder = bonus_quadorder, kwargs...)

    if order == 0
        return nothing
    end

    # set values of interior functions as piecewise best-approximation
    FEType = HDIVRTk{edim, order}
    ncomponents = get_ncomponents(FEType)
    EG = (ncomponents == 2) ? Triangle2D : Tetrahedron3D
    ndofs = get_ndofs(ON_CELLS, FEType, EG)
    interior_offset::Int = interior_dofs_offset(ON_CELLS, HDIVRTk{edim, order}, Triangle2D)
    nidofs::Int = ndofs - interior_offset
    ncells = num_sources(FE.dofgrid[CellNodes])
    xCellVolumes::Array{Tv, 1} = FE.dofgrid[CellVolumes]
    xCellRegions = FE.dofgrid[CellRegions]
    xCellDofs::DofMapTypes{Ti} = FE[CellDofs]
    qf = QuadratureRule{T, EG}(max(2 * order, order + 1 + bonus_quadorder))
    FEB = FEEvaluator(FE, Identity, qf; T = T)
    QP = QPInfos(FE.dofgrid; kwargs...)

    # evaluation of gradient of P1 functions
    FE3 = order == 1 ? L2P0{ncomponents} : H1Pk{ncomponents, 2, order - 1}
    FES3 = FESpace{FE3, ON_CELLS}(FE.dofgrid)
    FEBPk = FEEvaluator(FES3, Identity, qf; T = T)

    if items == []
        items = 1:ncells
    end

    interiordofs = zeros(Int, nidofs)
    basisvals::Array{T, 3} = FEB.cvals
    basisvalsPk::Array{T, 3} = FEBPk.cvals
    IMM_face = zeros(T, nidofs, interior_offset)
    IMM = zeros(T, nidofs, nidofs)
    lb = zeros(T, nidofs)
    temp::T = 0
    data_eval = zeros(T, ncomponents)
    for cell in items
        # update basis
        update_basis!(FEB, cell)
        update_basis!(FEBPk, cell)
        fill!(IMM, 0)
        fill!(IMM_face, 0)
        fill!(lb, 0)

        QP.item = cell
        QP.cell = cell
        QP.region = xCellRegions[cell]

        # quadrature loop
        for i in 1:length(qf.w)
            # right-hand side : f times P1
            eval_trafo!(QP.x, FEB.L2G, FEB.xref[i])
            QP.xref = FEB.xref[i]
            data(data_eval, QP)
            data_eval .*= xCellVolumes[cell] * qf.w[i]
            for dof in 1:nidofs
                for k in 1:ncomponents
                    lb[dof] += data_eval[k] * basisvalsPk[k, dof, i]
                end

                # mass matrix of interior basis functions
                for dof2 in 1:nidofs
                    temp = 0
                    for k in 1:ncomponents
                        temp += basisvals[k, interior_offset + dof, i] * basisvalsPk[k, dof2, i]
                    end
                    IMM[dof2, dof] += temp * xCellVolumes[cell] * qf.w[i]
                end

                # mass matrix of face basis functions
                for dof2 in 1:interior_offset
                    temp = 0
                    for k in 1:ncomponents
                        temp += basisvalsPk[k, dof, i] * basisvals[k, dof2, i]
                    end
                    IMM_face[dof, dof2] += temp * xCellVolumes[cell] * qf.w[i]
                end
            end
        end

        # subtract face interpolation from right-hand side
        for dof in 1:nidofs, dof2 in 1:interior_offset
            lb[dof] -= Target[xCellDofs[dof2, cell]] * IMM_face[dof, dof2]
        end

        # solve local system
        for dof in 1:nidofs
            interiordofs[dof] = xCellDofs[interior_offset + dof, cell]
        end

        Target[interiordofs] = IMM \ lb
    end
    return
end

# only normalfluxes on faces
function get_basis(::Union{Type{<:ON_FACES}, Type{<:ON_BFACES}}, ::Type{<:HDIVRTk{2, order}}, ::Type{<:AbstractElementGeometry1D}) where {order}
    moments_weights = basis.(ShiftedLegendre, (0:order))
    return function closure(refbasis, xref)
        refbasis[1, 1] = 1
        for j in 1:order
            refbasis[j + 1, 1] = (2 * j + 1) * moments_weights[j + 1](xref[1])
        end
        return
    end
end

function get_basis(::Type{ON_CELLS}, ::Type{HDIVRTk{2, order}}, ::Type{<:Triangle2D}) where {order}
    if order == 1
        interior_basis = get_basis(ON_CELLS, L2P0{1}, Triangle2D)
        ninterior = 1
    elseif order > 0
        interior_basis = get_basis(ON_CELLS, H1Pk{1, 2, order - 1}, Triangle2D)
        ninterior = Int((order + 1) * (order) / 2)
    else
        ninterior = 0
    end
    interior_offset = 3 * (order + 1)
    moments_weights = basis.(ShiftedLegendre, (0:order))
    moment_factors = [convert(Polynomial, m).coeffs[end] for m in moments_weights]
    return function closure(refbasis, xref)
        # RT0 basis
        refbasis[1, 1] = xref[1]
        refbasis[1, 2] = xref[2] - 1
        refbasis[2 + order, 1] = xref[1]
        refbasis[2 + order, 2] = xref[2]
        refbasis[3 + 2 * order, 1] = xref[1] - 1
        refbasis[3 + 2 * order, 2] = xref[2]
        for k in 1:2
            for j in 1:order
                refbasis[1 + j, k] = (-1)^j * (2 * j + 1) * moments_weights[j + 1](1 - xref[1] - xref[2]) * refbasis[1, k]
                refbasis[2 + j + order, k] = (-1)^j * (2 * j + 1) * moments_weights[j + 1](xref[1]) * refbasis[2 + order, k]
                refbasis[3 + j + 2 * order, k] = (-1)^j * (2 * j + 1) * moments_weights[j + 1](xref[2]) * refbasis[3 + 2 * order, k]
            end

            # interior RT2 functions (RT1 interior functions times P1) = 6 dofs
            if order > 0
                interior_basis(view(refbasis, (interior_offset + ninterior + 1):(interior_offset + 2 * ninterior), 2), xref)
                for j in 1:ninterior
                    refbasis[interior_offset + j, k] = 12 * xref[2] * refbasis[1, k] * refbasis[interior_offset + ninterior + j, 2]
                    refbasis[interior_offset + ninterior + j, k] = 12 * xref[1] * refbasis[3 + 2 * order, k] * refbasis[interior_offset + ninterior + j, 2]
                end
            end
        end
        return
    end
end


function get_coefficients(::Type{ON_CELLS}, FE::FESpace{Tv, Ti, <:HDIVRTk{2, order}, APT}, EG::Type{<:Triangle2D}, xgrid) where {Tv, Ti, APT, order}
    xCellFaceSigns::Union{VariableTargetAdjacency{Int32}, Array{Int32, 2}} = xgrid[CellFaceSigns]
    nfaces::Int = num_faces(EG)
    dim::Int = dim_element(EG)
    return function closure(coefficients::Array{<:Real, 2}, cell::Int)
        fill!(coefficients, 1.0)
        # multiplication with normal vector signs (only RT0, RT2, RT4 etc.)
        for j in 1:nfaces, k in 1:dim
            coefficients[k, (order + 1) * j - order] = xCellFaceSigns[j, cell]
            for o in 2:order
                if iseven(o)
                    coefficients[k, (order + 1) * j - (order - o)] = xCellFaceSigns[j, cell]
                end
            end
        end
        return nothing
    end
end
