StandardFunctionOperator(::Type{<:ReconstructionOperator{FETypeR, O}}) where {FETypeR, O} = O
ReconstructionSpace(::Type{<:ReconstructionOperator{FETypeR, O}}) where {FETypeR, O} = FETypeR

NeededDerivative4Operator(::Type{<:ReconstructionOperator{FETypeR, O}}) where {FETypeR, O} = NeededDerivative4Operator(O)
Length4Operator(::Type{<:ReconstructionOperator{FETypeR, O}}, xdim, nc) where {FETypeR, O} = Length4Operator(O, xdim, nc)
DefaultName4Operator(::Type{Reconstruct{FETypeR, O}}) where {FETypeR, O} = "R(" * DefaultName4Operator(O) * ")"
DefaultName4Operator(::Type{WeightedReconstruct{FETypeR, O}}) where {FETypeR, O} = "R(r" * DefaultName4Operator(O) * ")"

struct FEReconstEvaluator{T, TvG, TiG, FEType, FEType2, stdop, O <: ReconstructionOperator} <: FEEvaluator{T, TvG, TiG}
    citem::Base.RefValue{Int}                       # current item
    FE::FESpace{TvG, TiG, FEType}                     # link to full FE (e.g. for coefficients)
    FEB::SingleFEEvaluator{T, TvG, TiG, stdop, FEType2}                # FEBasisEvaluator for stdop in reconstruction space
    cvals::Array{T, 3}                               # current operator vals on item (reconstruction)
    coefficients::Array{TvG, 2}                      # additional coefficients for reconstruction
    reconst_handler::ReconstructionHandler{T, TiG, O}  # handler for reconstruction coefficients
end

# constructor for reconstruction operators
function FEEvaluator(
        FE::FESpace{TvG, TiG, FEType, FEAPT},
        operator::Type{<:ReconstructionOperator{FETypeReconst, stdop}},
        qrule::QuadratureRule{TvR, EG},
        xgrid = FE.xgrid;
        L2G = nothing,
        T = Float64,
        AT = ON_CELLS
    ) where {TvG, TiG, TvR, FEType <: AbstractFiniteElement, stdop <: StandardFunctionOperator, FETypeReconst <: AbstractFiniteElement, EG <: AbstractElementGeometry, FEAPT <: AssemblyType}

    @debug "Creating FEBasisEvaluator for reconstruction of $stdop operator of $FEType into $FETypeReconst on $EG"

    ## generate FESpace for reconstruction
    FE2 = FESpace{FETypeReconst}(xgrid)

    ## collect dimension information
    ncomponents::Int = get_ncomponents(FEType)
    ncomponents2::Int = get_ncomponents(FETypeReconst)
    if AT <: Union{ON_BFACES, <:ON_FACES}
        if FETypeReconst <: AbstractHdivFiniteElement
            ncomponents2 = 1
        end
    end
    ndofs4item = get_ndofs(AT, FEType, EG)
    ndofs4item2 = get_ndofs(AT, FETypeReconst, EG)
    coefficients2 = zeros(T, ndofs4item, ndofs4item2)
    edim = dim_element(EG)
    resultdim = Int(Length4Operator(operator, edim, ncomponents))
    cvals = zeros(T, resultdim, ndofs4item, length(qrule.xref))

    ## FEEvaluator for reconstruction space
    if L2G === nothing
        L2G = L2GTransformer(EG, xgrid, AT)
    end
    FEB = FEEvaluator(FE2, stdop, qrule, xgrid; L2G = L2G, T = T, AT = AT)

    ## reconstruction coefficient handler
    reconst_handler = ReconstructionHandler(FE, FE2, AT, EG, operator)

    return FEReconstEvaluator{T, TvG, TiG, FEType, FETypeReconst, stdop, operator}(FEB.citem, FE, FEB, cvals, coefficients2, reconst_handler)
end

function update_basis!(FEBE::FEReconstEvaluator)
    ## evaluate standard operator in reconstruction basis (->cvals_reconst)
    update_basis!(FEBE.FEB)
    cvals_reconst = FEBE.FEB.cvals

    # get local reconstruction coefficients
    rcoeffs = FEBE.coefficients
    get_rcoefficients!(rcoeffs, FEBE.reconst_handler, FEBE.citem[])

    # accumulate
    cvals = FEBE.cvals
    fill!(cvals, 0)
    for dof_i in 1:size(rcoeffs, 1), dof_j in 1:size(rcoeffs, 2)
        if rcoeffs[dof_i, dof_j] != 0
            for i in 1:size(cvals, 3), k in 1:size(cvals_reconst, 1)
                cvals[k, dof_i, i] += rcoeffs[dof_i, dof_j] * cvals_reconst[k, dof_j, i]
            end
        end
    end
    return nothing
end

function Base.show(io::IO, FEB::FEReconstEvaluator)
    println(io, "FEReconstEvaluator (reconstruction operator)")
    println(io, "-------------------------------------------")
    println(io, "Reconstruction FE space: ", typeof(FEB.FEB.FE))
    println(io, "Reconstruction operator: ", typeof(FEB.FEB).parameters[4])
    println(io, "Underlying SingleFEEvaluator:")
    show(io, FEB.FEB)
    return
end
