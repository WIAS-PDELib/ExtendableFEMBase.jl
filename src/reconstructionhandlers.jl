abstract type ReconstructionOperator{FETypeR, O} <: AbstractFunctionOperator end

"""
$(TYPEDEF)

reconstruction operator: evaluates a reconstructed version of the finite element function.

FETypeR specifies the reconstruction space (needs to be defined for the finite element that it is applied to).
O specifies the StandardFunctionOperator that shall be evaluated.
"""
abstract type Reconstruct{FETypeR, O} <: ReconstructionOperator{FETypeR, O} end


"""
$(TYPEDEF)

Weighted reconstruction operator: evaluates a reconstructed version of a finite element function, multiplied by a weight function.

# Parameters
- `FETypeR`: The reconstruction finite element space type (target space for reconstruction).
- `O`: The standard function operator to be evaluated (e.g., identity, gradient, etc.).
- `F`: The type of the weight function (should be callable, e.g., a function or functor).
"""
abstract type WeightedReconstruct{FETypeR, O, F} <: Reconstruct{FETypeR, O} end

weight_type(::Type{<:WeightedReconstruct{FETypeR, O, F}}) where {FETypeR, O, F} = F


################## SPECIAL INTERPOLATORS ####################

"""
````
function interpolator_matrix(::Type{<:HDIVRT0{ncomponents}}, V1::FESpace{Tv, Ti, H1Pk{ncomponents, edim, order}, ON_CELLS})
````

produces a matrix representation for the interpolation
of H1Pk (currently only up to order 4, and ncomponents = 2) 
into HDIVRT0
"""
function interpolator_matrix(::Type{<:HDIVRT0{ncomponents}}, V1::FESpace{Tv, Ti, H1Pk{ncomponents, edim, order}, ON_CELLS}) where {Tv, Ti, ncomponents, edim, order}

    @assert order in 1:4 "higher order H1PK->RT0 interpolator not yet implemented"
    @assert ncomponents == 2 && edim == ncomponents
    ## setup interpolation matrix F
    ## that maps a V1 function to a RT0 function
    ## with same negative divergence in ̂P_0
    ## we just need to evaluate normal fluxes of basis functions in V1

    xgrid = V1.xgrid
    facedofs_V1::VariableTargetAdjacency{Int32} = V1[FaceDofs]
    ndofs_V1 = max_num_targets_per_source(facedofs_V1)
    if edim == 2
        @assert ndofs_V1 == 2 * (2 + (order - 1))
        ## quadrature weights for Newton-Cotes formula associated to V1
        if order == 1
            factors = [1 // 2, 1 // 2]
        elseif order == 2
            factors = [1 // 6, 1 // 6, 2 // 3]
        elseif order == 3
            factors = [1 // 8, 1 // 8, 3 // 8, 3 // 8]
        elseif order == 4
            factors = [7 // 90, 7 // 90, 32 // 90, 12 // 90, 32 // 90]
        end
        append!(factors, factors)
        nc = ones(Int, ndofs_V1)
        nc[(3 + (order - 1)):end] .= 2
    end

    FE = ExtendableSparseMatrix{Float64, Int64}(V1.ndofs, size(xgrid[FaceNodes], 2))
    xFaceVolumes = xgrid[FaceVolumes]
    xFaceNormals = xgrid[FaceNormals]
    nfaces = num_sources(xFaceNormals)
    for face in 1:nfaces
        for dof1 in 1:ndofs_V1
            FE[facedofs_V1[dof1, face], face] = -factors[dof1] * xFaceNormals[nc[dof1], face] * xFaceVolumes[face]
        end
    end
    flush!(FE)
    return FE
end


"""
$(TYPEDEF)

abstract grid component type for storing reconstruction coefficients
"""
abstract type ReconstructionCoefficients{FE1, FE2, AT} <: AbstractGridFloatArray2D end

"""
$(TYPEDEF)

abstract grid component type for storing reconstruction dof numbers
"""
abstract type ReconstructionDofs{FE1, FE2, AT} <: AbstractGridIntegerArray2D end

"""
$(TYPEDEF)

abstract grid component type that can be used to funnel reconstruction weights into the operator
(default is 1)
"""
abstract type ReconstructionWeights{AT} <: AbstractGridFloatArray2D end

"""
$(TYPEDEF)

struct for storing information needed to evaluate a reconstruction operator
"""
struct ReconstructionHandler{Tv, Ti, RT, FE1, FE2, AT, EG, F}
    FES::FESpace{Tv, Ti, FE1, ON_CELLS}
    FER::FESpace{Tv, Ti, FE2, ON_CELLS}
    xCoordinates::Array{Tv, 2}
    xFaceNodes::Adjacency{Ti}
    xFaceVolumes::Array{Tv, 1}
    xFaceNormals::Array{Tv, 2}
    xCellFaceOrientations::Adjacency{Ti}
    xCellFaces::Adjacency{Ti}
    interior_offset::Int
    interior_ndofs::Int
    interior_coefficients::Matrix{Tv} # coefficients for interior basis functions are precomputed
    weight::F
end

function default_weight_function(x)
    return 1
end

"""
$(TYPEDSIGNATURES)

generates a reconstruction handler
returns the local coefficients need to evaluate a reconstruction operator
of one finite element space into another
"""
function ReconstructionHandler(FES::FESpace{Tv, Ti, FE1, APT}, FES_Reconst::FESpace{Tv, Ti, FE2, APT}, AT, EG, RT, weight = default_weight_function) where {Tv, Ti, FE1, FE2, APT}
    xgrid = FES.xgrid
    interior_offset = interior_dofs_offset(AT, FE2, EG)
    interior_ndofs = get_ndofs(AT, FE2, EG) - interior_offset
    if interior_offset != -1 && interior_ndofs > 0
        coeffs = xgrid[ReconstructionCoefficients{FE1, FE2, AT}]
    else
        interior_ndofs = 0
        coeffs = zeros(Tv, 0, 0)
    end

    xFaceVolumes = xgrid[FaceVolumes]
    xFaceNormals = xgrid[FaceNormals]
    xCoordinates = xgrid[Coordinates]
    xFaceNodes = xgrid[FaceNodes]
    xCellFaceOrientations = dim_element(EG) == 2 ? xgrid[CellFaceSigns] : xgrid[CellFaceOrientations]
    xCellFaces = xgrid[CellFaces]
    return ReconstructionHandler{Tv, Ti, RT, FE1, FE2, AT, EG, typeof(weight)}(FES, FES_Reconst, xCoordinates, xFaceNodes, xFaceVolumes, xFaceNormals, xCellFaceOrientations, xCellFaces, interior_offset, interior_ndofs, coeffs, weight)
end

"""
$(TYPEDSIGNATURES)

caller function to extract the coefficients of the reconstruction
on an item
"""
function get_rcoefficients!(coefficients, RH::ReconstructionHandler, item)
    boundary_coefficients!(coefficients, RH, item)
    for dof in 1:size(coefficients, 1), k in 1:RH.interior_ndofs
        coefficients[dof, RH.interior_offset + k] = RH.interior_coefficients[(dof - 1) * RH.interior_ndofs + k, item]
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

generates the interior coefficients for the RT1 reconstruction of P2B
"""
function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tv, Ti}, ::Type{ReconstructionCoefficients{FE1, FE2, AT}}) where {Tv, Ti, FE1 <: H1P2B{2, 2}, FE2 <: HDIVRT1{2}, AT <: ON_CELLS}
    @info "Computing interior reconstruction coefficients for $FE1 > $FE2 ($AT)"
    xCellFaces = xgrid[CellFaces]
    xCoordinates = xgrid[Coordinates]
    xCellNodes = xgrid[CellNodes]
    xFaceVolumes::Array{Tv, 1} = xgrid[FaceVolumes]
    xFaceNormals::Array{Tv, 2} = xgrid[FaceNormals]
    xCellVolumes::Array{Tv, 1} = xgrid[CellVolumes]
    xCellFaceSigns = xgrid[CellFaceSigns]
    EG = xgrid[UniqueCellGeometries]

    @assert EG == [Triangle2D]

    face_rule::Array{Int, 2} = local_cellfacenodes(EG[1])
    RT1_coeffs = _P1_INTO_BDM1_COEFFS
    nnf::Int = size(face_rule, 2)
    ndofs4component::Int = 2 * nnf + 1
    ndofs1::Int = get_ndofs(AT, FE1, EG[1])
    ncells::Int = num_sources(xCellFaces)
    coefficients::Array{Tv, 2} = zeros(Tv, ndofs1, 6)
    interior_coefficients::Array{Tv, 2} = zeros(Tv, 2 * ndofs1, ncells)

    C = zeros(Tv, 2, 3)  # vertices
    E = zeros(Tv, 2, 3)  # edge midpoints
    M = zeros(Tv, 2)    # midpoint of current cell
    A = zeros(Tv, 2, 8)  # integral means of RT1 functions (from analytic formulas)
    b = zeros(Tv, 2)    # right-hand side for integral mean
    dof::Int = 0
    det::Tv = 0
    face::Int = 0
    node::Int = 0
    for cell in 1:ncells

        for f in 1:nnf
            face = xCellFaces[f, cell]
            for n in 1:2
                node = face_rule[n, f]
                for k in 1:2
                    # RT0 reconstruction coefficients for node P2 functions on reference element
                    coefficients[ndofs4component * (k - 1) + node, 2 * (f - 1) + 1] = 1 // 6 * xFaceVolumes[face] * xFaceNormals[k, face]

                    # RT1 reconstruction coefficients for node P2 functions on reference element
                    coefficients[ndofs4component * (k - 1) + node, 2 * (f - 1) + 2] = RT1_coeffs[n] * xFaceVolumes[face] * xFaceNormals[k, face] * xCellFaceSigns[f, cell]
                end
            end
            for k in 1:2
                # RT0 reconstruction coefficients for face P2 functions (=face bubbles) on reference element
                coefficients[ndofs4component * (k - 1) + f + nnf, 2 * (f - 1) + 1] = 2 // 3 * xFaceVolumes[face] * xFaceNormals[k, face]
            end
        end

        # get coordinates of cells
        fill!(M, 0)
        fill!(E, 0)
        for n in 1:3, k in 1:2
            C[k, n] = xCoordinates[k, xCellNodes[n, cell]]
            M[k] += C[k, n] / 3
        end

        # get edge midpoints
        for f in 1:nnf
            for n in 1:2, k in 1:2
                E[k, f] += C[k, face_rule[n, f]] / 2
            end
        end

        # compute integral means of RT1 functions
        for k in 1:2
            A[k, 1] = (M[k] - C[k, 3]) / 2 * xCellFaceSigns[1, cell]
            A[k, 2] = C[k, 2] - E[k, 2]
            A[k, 3] = (M[k] - C[k, 1]) / 2 * xCellFaceSigns[2, cell]
            A[k, 4] = C[k, 3] - E[k, 3]
            A[k, 5] = (M[k] - C[k, 2]) / 2 * xCellFaceSigns[3, cell]
            A[k, 6] = C[k, 1] - E[k, 1]
        end
        # directly assign inverted A[1:2,7:8] for faster solve of local systems
        A[2, 8] = (E[1, 1] - C[1, 3]) # A[1,7]
        A[2, 7] = -(E[2, 1] - C[2, 3]) # A[2,7]
        A[1, 8] = -(E[1, 3] - C[1, 2]) # A[1,8]
        A[1, 7] = (E[2, 3] - C[2, 2]) # A[2,8]

        det = A[1, 7] * A[2, 8] - A[2, 7] * A[1, 8]
        A[1:2, 7:8] ./= det

        # correct integral means with interior RT1 functions
        for k in 1:2
            for n in 1:3
                # nodal P2 functions have integral mean zero
                dof = (ndofs4component * (k - 1) + n)
                fill!(b, 0)
                for c in 1:2, j in 1:6
                    b[c] -= coefficients[dof, j] * A[c, j]
                end
                for k in 1:2
                    interior_coefficients[(dof - 1) * 2 + k, cell] = A[k, 7] * b[1] + A[k, 8] * b[2]
                end

                # face P2 functions have integral mean 1//3
                dof = (ndofs4component * (k - 1) + n + nnf)
                fill!(b, 0)
                b[k] = xCellVolumes[cell] / 3
                for c in 1:2, j in 1:6
                    b[c] -= coefficients[dof, j] * A[c, j]
                end
                for k in 1:2
                    interior_coefficients[(dof - 1) * 2 + k, cell] = A[k, 7] * b[1] + A[k, 8] * b[2]
                end
            end

            # cell bubbles have integral mean 1
            dof = ndofs4component * k
            fill!(b, 0)
            b[k] = xCellVolumes[cell]
            for k in 1:2
                interior_coefficients[(dof - 1) * 2 + k, cell] = A[k, 7] * b[1] + A[k, 8] * b[2]
            end
        end
    end
    return interior_coefficients
end

"""
$(TYPEDSIGNATURES)

generates the interior coefficients for the BDM2 reconstruction of P2B
"""
function ExtendableGrids.instantiate(xgrid::ExtendableGrid{Tv, Ti}, ::Type{ReconstructionCoefficients{FE1, FE2, AT}}) where {Tv, Ti, FE1 <: H1P2B{2, 2}, FE2 <: HDIVBDM2{2}, AT <: ON_CELLS}
    @info "Computing interior reconstruction coefficients for $FE1 > $FE2 ($AT)"
    xCellFaces = xgrid[CellFaces]
    xCellVolumes::Array{Tv, 1} = xgrid[CellVolumes]
    xCellFaceSigns = xgrid[CellFaceSigns]
    xFaceVolumes::Array{Tv, 1} = xgrid[FaceVolumes]
    xFaceNormals::Array{Tv, 2} = xgrid[FaceNormals]
    EG = xgrid[UniqueCellGeometries]

    @assert EG == [Triangle2D]

    face_rule::Array{Int, 2} = local_cellfacenodes(EG[1])
    coeffs1 = _P1_INTO_BDM1_COEFFS
    nnf = size(face_rule, 2)
    ndofs4component = 2 * nnf + 1

    ndofs1::Int = get_ndofs(AT, FE1, EG[1])
    ncells::Int = num_sources(xCellFaces)
    interior_offset::Int = 9
    interior_ndofs::Int = 3
    coefficients::Array{Tv, 2} = zeros(Tv, ndofs1, interior_offset)
    interior_coefficients::Array{Tv, 2} = zeros(Tv, interior_ndofs * ndofs1, ncells)

    qf = QuadratureRule{Tv, EG[1]}(4)
    weights::Array{Tv, 1} = qf.w
    # evaluation of FE1 and FE2 basis
    FES1 = FESpace{FE1, ON_CELLS}(xgrid)
    FES2 = FESpace{FE2, ON_CELLS}(xgrid)
    FEB1 = FEEvaluator(FES1, Identity, qf; T = Tv)
    FEB2 = FEEvaluator(FES2, Identity, qf; T = Tv)
    # evaluation of gradient of P1 functions
    FE3 = H1P1{1}
    FES3 = FESpace{FE3, ON_CELLS}(xgrid)
    FEB3 = FEEvaluator(FES3, Gradient, qf; T = Tv)
    # evaluation of curl of bubble functions
    FE4 = H1BUBBLE{1}
    FES4 = FESpace{FE4, ON_CELLS}(xgrid)
    FEB4 = FEEvaluator(FES4, CurlScalar, qf; T = Tv)

    basisvals1::Array{Tv, 3} = FEB1.cvals
    basisvals2::Array{Tv, 3} = FEB2.cvals
    basisvals3::Array{Tv, 3} = FEB3.cvals
    basisvals4::Array{Tv, 3} = FEB4.cvals
    IMM_face = zeros(Tv, interior_ndofs, interior_offset)
    IMM = zeros(Tv, interior_ndofs, interior_ndofs)
    for k in 1:interior_ndofs
        IMM[k, k] = 1
    end
    lb = zeros(Tv, interior_ndofs)
    lx = zeros(Tv, interior_ndofs)
    temp::Tv = 0
    face = 0
    node = 0
    offset::Int = 0
    IMMfact = lu(IMM)
    for cell in 1:ncells

        # get reconstruction coefficients for boundary dofs
        for f in 1:nnf
            face = xCellFaces[f, cell]
            for n in 1:2
                node = face_rule[n, f]
                for k in 1:2
                    # RT0 reconstruction coefficients for node P2 functions on reference element
                    coefficients[ndofs4component * (k - 1) + node, 3 * (f - 1) + 1] = 1 // 6 * xFaceVolumes[face] * xFaceNormals[k, face]

                    # 1st BDM2 reconstruction coefficients for node P2 functions on reference element
                    coefficients[ndofs4component * (k - 1) + node, 3 * (f - 1) + 2] = coeffs1[n] * xFaceVolumes[face] * xFaceNormals[k, face] * xCellFaceSigns[f, cell]

                    # 2nd BDM2 reconstruction coefficients for node P2 functions on reference element
                    coefficients[ndofs4component * (k - 1) + node, 3 * (f - 1) + 3] = 1 // 90 * xFaceVolumes[face] * xFaceNormals[k, face]
                end
            end
            for k in 1:2
                # RT0 reconstruction coefficients for face P2 functions (=face bubbles) on reference element
                coefficients[ndofs4component * (k - 1) + f + nnf, 3 * (f - 1) + 1] = 2 // 3 * xFaceVolumes[face] * xFaceNormals[k, face]

                # 2nd BDM2 reconstruction coefficients for face P2 functions on reference element
                coefficients[ndofs4component * (k - 1) + f + nnf, 3 * (f - 1) + 3] = -1 // 45 * xFaceVolumes[face] * xFaceNormals[k, face]
            end
        end

        # update basis
        update_basis!(FEB1, cell)
        update_basis!(FEB2, cell)
        update_basis!(FEB3, cell)
        update_basis!(FEB4, cell)

        # compute local mass matrices
        fill!(IMM, 0)
        fill!(IMM_face, 0)
        for i in eachindex(weights)
            for dof in 1:interior_ndofs
                # interior FE2 basis functions times grad(P1) of first two P1 functions
                for dof2 in 1:(interior_ndofs - 1)
                    temp = 0
                    for k in 1:2
                        temp += basisvals2[k, interior_offset + dof, i] * basisvals3[k, dof2, i]
                    end
                    IMM[dof2, dof] += temp * xCellVolumes[cell] * weights[i]
                end
                # interior FE2 basis functions times curl(bubble)
                temp = 0
                for k in 1:2
                    temp += basisvals2[k, interior_offset + dof, i] * basisvals4[k, 1, i]
                end
                IMM[3, dof] += temp * xCellVolumes[cell] * weights[i]

                # mass matrix of face basis functions x grad(P1) and curl(bubble)
                if dof < 3
                    for dof2 in 1:interior_offset
                        temp = 0
                        for k in 1:2
                            temp += basisvals3[k, dof, i] * basisvals2[k, dof2, i]
                        end
                        IMM_face[dof, dof2] += temp * xCellVolumes[cell] * weights[i]
                    end
                    # mass matrix of face basis functions x interior basis functions
                elseif dof == 3
                    for dof2 in 1:interior_offset
                        temp = 0
                        for k in 1:2
                            temp += basisvals4[k, 1, i] * basisvals2[k, dof2, i]
                        end
                        IMM_face[dof, dof2] += temp * xCellVolumes[cell] * weights[i]
                    end
                end
            end
        end

        # solve local systems
        IMMfact = lu(IMM)
        for dof1 in 1:ndofs1
            # right-hand side
            fill!(lb, 0)
            for i in eachindex(weights)
                for idof in 1:(interior_ndofs - 1)
                    temp = 0
                    for k in 1:2
                        temp += basisvals1[k, dof1, i] * basisvals3[k, idof, i]
                    end
                    lb[idof] += temp * xCellVolumes[cell] * weights[i]
                end
                temp = 0
                for k in 1:2
                    temp += basisvals1[k, dof1, i] * basisvals4[k, 1, i]
                end
                lb[3] += temp * xCellVolumes[cell] * weights[i]
            end

            # subtract face interpolation from right-hand side
            for idof in 1:interior_ndofs, dof2 in 1:interior_offset
                lb[idof] -= coefficients[dof1, dof2] * IMM_face[idof, dof2]
            end

            # solve local system
            ldiv!(lx, IMMfact, lb)
            offset = interior_ndofs * (dof1 - 1)
            for idof in 1:interior_ndofs
                interior_coefficients[offset + idof, cell] = lx[idof]
            end
        end
    end
    return interior_coefficients
end


##### BOUNDARY COEFFICIENTS #####

"""
````
boundary_coefficients!(coefficients, RH::ReconstructionHandler, cell)
````

generates the coefficients for the facial dofs of the reconstruction operator on the cell
"""
function boundary_coefficients!(coefficients, RH::ReconstructionHandler{Tv, Ti, <:H1CR{ncomponents}, <:HDIVRT0{ncomponents}, <:ON_CELLS, EG}, cell) where {Tv, Ti, ncomponents, EG}
    xFaceVolumes = RH.xFaceVolumes
    xFaceNormals = RH.xFaceNormals
    xCellFaces = RH.xCellFaces
    face_rule = local_cellfacenodes(EG)
    nfaces = size(face_rule, 2)
    face = 0
    for f in 1:nfaces
        face = xCellFaces[f, cell]
        for k in 1:ncomponents
            coefficients[nfaces * (k - 1) + f, f] = xFaceVolumes[face] * xFaceNormals[k, face]
        end
    end
    return nothing
end

const _P1_INTO_BDM1_COEFFS = [-1 // 12, 1 // 12]

function boundary_coefficients!(coefficients, RH::ReconstructionHandler{Tv, Ti, RT, FE1, FE2, AT, EG}, cell) where {Tv, Ti, RT <: Reconstruct, FE1 <: H1BR{2}, FE2 <: HDIVBDM1{2}, AT <: ON_CELLS, EG <: Union{Triangle2D, Quadrilateral2D}}
    xFaceVolumes = RH.xFaceVolumes
    xFaceNormals = RH.xFaceNormals
    xCellFaceSigns = RH.xCellFaceOrientations
    xCellFaces = RH.xCellFaces
    xFaceNodes = RH.xFaceNodes
    xCoordinates = RH.xCoordinates
    face_rule = local_cellfacenodes(EG)
    nnodes = size(face_rule, 1)
    nfaces = size(face_rule, 2)
    node = 0
    face = 0
    BDM1_coeffs = _P1_INTO_BDM1_COEFFS
    weight = RH.weight
    xmid = zeros(Tv, 2)
    w = ones(Tv, 2)
    for f in 1:nfaces
        face = xCellFaces[f, cell]
        x1 = view(xCoordinates, :, xFaceNodes[1, face])
        x2 = view(xCoordinates, :, xFaceNodes[2, face])
        xmid .= 0.5 * (x1 .+ x2)
        w[1] = weight(x1)
        w[2] = weight(x2)
        wmid = weight(xmid)
        for n in 1:nnodes
            node = face_rule[n, f]
            for k in 1:2
                # RT0 reconstruction coefficients for P1 functions on reference element
                coefficients[nfaces * (k - 1) + node, 2 * (f - 1) + 1] = xFaceVolumes[face] * (1 // 6 * w[n] + 1 // 3 * wmid) * xFaceNormals[k, face]
                # BDM1 reconstruction coefficients for P1 functions on reference element
                coefficients[nfaces * (k - 1) + node, 2 * (f - 1) + 2] = xFaceVolumes[face] * xFaceNormals[k, face] * xCellFaceSigns[f, cell] * (BDM1_coeffs[n] * w[n] + 1 // 3 * wmid)
            end
        end
        # RT0 reconstruction coefficients for face bubbles on reference element
        coefficients[nfaces * 2 + f, 2 * (f - 1) + 1] = xFaceVolumes[face] * wmid
    end

    return nothing
end

function boundary_coefficients!(coefficients, RH::ReconstructionHandler{Tv, Ti, RT, FE1, FE2, AT, EG}, cell) where {Tv, Ti, RT <: Reconstruct, FE1 <: H1BR{2}, FE2 <: HDIVRT0{2}, AT <: ON_CELLS, EG <: Union{Triangle2D, Quadrilateral2D}}
    xFaceVolumes = RH.xFaceVolumes
    xFaceNormals = RH.xFaceNormals
    xCellFaces = RH.xCellFaces
    xCellFaceSigns = RH.xCellFaceOrientations
    xFaceNodes = RH.xFaceNodes
    xCoordinates = RH.xCoordinates
    face_rule = local_cellfacenodes(EG)
    nnodes = size(face_rule, 1)
    nfaces = size(face_rule, 2)
    node = 0
    face = 0
    weight = RH.weight
    xmid = zeros(Tv, 2)
    w = ones(Tv, 2)
    for f in 1:nfaces
        face = xCellFaces[f, cell]
        x1 = view(xCoordinates, :, xFaceNodes[1, face])
        x2 = view(xCoordinates, :, xFaceNodes[2, face])
        xmid .= 0.5 * (x1 .+ x2)
        if xCellFaceSigns[f, cell] == -1
            w[1] = weight(x2)
            w[2] = weight(x1)
        else
            w[1] = weight(x1)
            w[2] = weight(x2)
        end
        wmid = weight(xmid)
        # reconstruction coefficients for P1 functions on reference element
        for n in 1:nnodes
            node = face_rule[n, f]
            for k in 1:2
                coefficients[nfaces * (k - 1) + node, f] = xFaceVolumes[face] * (1 // 6 * w[n] + 1 // 3 * wmid) * xFaceNormals[k, face]
            end
        end
        # reconstruction coefficients for face bubbles on reference element
        coefficients[2 * nfaces + f, f] = xFaceVolumes[face] * wmid
    end
    return nothing
end

function boundary_coefficients!(coefficients, RH::ReconstructionHandler{Tv, Ti, RT, FE1, FE2, AT, EG}, cell) where {Tv, Ti, RT <: Reconstruct, FE1 <: H1P2B{2, 2}, FE2 <: HDIVRT1{2}, AT <: ON_CELLS, EG <: Triangle2D}
    xFaceVolumes = RH.xFaceVolumes
    xFaceNormals = RH.xFaceNormals
    xCellFaceSigns = RH.xCellFaceOrientations
    xCellFaces = RH.xCellFaces
    face_rule = local_cellfacenodes(EG)
    node = 0
    face = 0
    nnf = size(face_rule, 2)
    ndofs4component = 2 * nnf + 1
    RT1_coeffs = _P1_INTO_BDM1_COEFFS
    for f in 1:nnf
        face = xCellFaces[f, cell]
        for n in 1:2
            node = face_rule[n, f]
            for k in 1:2
                # RT0 reconstruction coefficients for node P2 functions on reference element
                coefficients[ndofs4component * (k - 1) + node, 2 * (f - 1) + 1] = 1 // 6 * xFaceVolumes[face] * xFaceNormals[k, face]

                # RT1 reconstruction coefficients for node P2 functions on reference element
                coefficients[ndofs4component * (k - 1) + node, 2 * (f - 1) + 2] = RT1_coeffs[n] * xFaceVolumes[face] * xFaceNormals[k, face] * xCellFaceSigns[f, cell]
            end
        end
        for k in 1:2
            # RT0 reconstruction coefficients for face P2 functions (=face bubbles) on reference element
            coefficients[ndofs4component * (k - 1) + f + nnf, 2 * (f - 1) + 1] = 2 // 3 * xFaceVolumes[face] * xFaceNormals[k, face]
        end
    end
    return nothing
end

function boundary_coefficients!(coefficients, RH::ReconstructionHandler{Tv, Ti, RT, FE1, FE2, AT, EG}, cell) where {Tv, Ti, RT <: Reconstruct, FE1 <: H1P2B{2, 2}, FE2 <: HDIVBDM2{2}, AT <: ON_CELLS, EG <: Triangle2D}
    xFaceVolumes = RH.xFaceVolumes
    xFaceNormals = RH.xFaceNormals
    xCellFaceSigns = RH.xCellFaceOrientations
    xCellFaces = RH.xCellFaces
    face_rule = local_cellfacenodes(EG)
    node = 0
    face = 0
    nnf = size(face_rule, 2)
    ndofs4component = 2 * nnf + 1
    coeffs1 = _P1_INTO_BDM1_COEFFS
    for f in 1:nnf
        face = xCellFaces[f, cell]
        for n in 1:2
            node = face_rule[n, f]
            for k in 1:2
                # RT0 reconstruction coefficients for node P2 functions on reference element
                coefficients[ndofs4component * (k - 1) + node, 3 * (f - 1) + 1] = 1 // 6 * xFaceVolumes[face] * xFaceNormals[k, face]

                # 1st BDM2 reconstruction coefficients for node P2 functions on reference element
                coefficients[ndofs4component * (k - 1) + node, 3 * (f - 1) + 2] = coeffs1[n] * xFaceVolumes[face] * xFaceNormals[k, face] * xCellFaceSigns[f, cell]

                # 2nd BDM2 reconstruction coefficients for node P2 functions on reference element
                coefficients[ndofs4component * (k - 1) + node, 3 * (f - 1) + 3] = 1 // 90 * xFaceVolumes[face] * xFaceNormals[k, face]
            end
        end
        for k in 1:2
            # RT0 reconstruction coefficients for face P2 functions (=face bubbles) on reference element
            coefficients[ndofs4component * (k - 1) + f + nnf, 3 * (f - 1) + 1] = 2 // 3 * xFaceVolumes[face] * xFaceNormals[k, face]

            # 2nd BDM2 reconstruction coefficients for face P2 functions on reference element
            coefficients[ndofs4component * (k - 1) + f + nnf, 3 * (f - 1) + 3] = -1 // 45 * xFaceVolumes[face] * xFaceNormals[k, face]
        end
    end
    return nothing
end


function boundary_coefficients!(coefficients, RH::ReconstructionHandler{Tv, Ti, RT, FE1, FE2, AT, EG}, cell) where {Tv, Ti, RT <: Reconstruct, FE1 <: H1BR{3}, FE2 <: HDIVRT0{3}, AT <: ON_CELLS, EG <: Tetrahedron3D}
    xFaceVolumes = RH.xFaceVolumes
    xFaceNormals = RH.xFaceNormals
    xCellFaces = RH.xCellFaces
    face_rule = local_cellfacenodes(EG)
    node = 0
    face = 0
    # fill!(coefficients,0.0)
    for f in 1:4
        face = xCellFaces[f, cell]
        # reconstruction coefficients for P1 functions on reference element
        for n in 1:3
            node = face_rule[n, f]
            for k in 1:3
                coefficients[4 * (k - 1) + node, f] = 1 // 3 * xFaceVolumes[face] * xFaceNormals[k, face]
            end
        end
        # reconstruction coefficients for face bubbles on reference element
        coefficients[12 + f, f] = xFaceVolumes[face]
    end
    return nothing
end

const _P1_INTO_BDM1_COEFFS_3D = [-1 // 36 -1 // 36 1 // 18; -1 // 36 1 // 18 -1 // 36; 1 // 18 -1 // 36 -1 // 36]

function boundary_coefficients!(coefficients, RH::ReconstructionHandler{Tv, Ti, RT, FE1, FE2, AT, EG}, cell) where {Tv, Ti, RT <: Reconstruct, FE1 <: H1BR{3}, FE2 <: HDIVBDM1{3}, AT <: ON_CELLS, EG <: Tetrahedron3D}
    xFaceVolumes = RH.xFaceVolumes
    xFaceNormals = RH.xFaceNormals
    xCellFaces = RH.xCellFaces
    xCellFaceOrientations = RH.xCellFaceOrientations
    face_rule = local_cellfacenodes(EG)
    node = 0
    face = 0
    face_rule = local_cellfacenodes(EG)
    BDM1_coeffs = _P1_INTO_BDM1_COEFFS_3D
    orientation = 0
    index1 = 0
    index2 = 0
    row4orientation1::Array{Int, 1} = [2, 2, 3, 1]
    row4orientation2::Array{Int, 1} = [1, 3, 1, 2]
    # fill!(coefficients,0.0)
    for f in 1:4
        face = xCellFaces[f, cell]
        index1 = 0
        for k in 1:3
            for n in 1:3
                node = face_rule[n, f]
                # RT0 reconstruction coefficients for P1 functions on reference element
                coefficients[index1 + node, index2 + 1] = 1 // 3 * xFaceNormals[k, face] * xFaceVolumes[face]
                orientation = xCellFaceOrientations[f, cell]
                # BDM1 reconstruction coefficients for P1 functions on reference element
                coefficients[index1 + node, index2 + 2] = BDM1_coeffs[n, row4orientation1[orientation]] * xFaceNormals[k, face] * xFaceVolumes[face]
                coefficients[index1 + node, index2 + 3] = BDM1_coeffs[n, row4orientation2[orientation]] * xFaceNormals[k, face] * xFaceVolumes[face]
            end
            index1 += 4
        end
        # RT0 reconstruction coefficients for face bubbles on reference element
        coefficients[index1 + f, index2 + 1] = xFaceVolumes[face]
        index2 += 3
    end
    return nothing
end

