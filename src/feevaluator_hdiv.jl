# IDENTITY HDIV
function update_basis!(FEBE::SingleFEEvaluator{<:Real, <:Real, <:Integer, <:Identity, <:AbstractHdivFiniteElement})
    L2GM = _update_piola!(FEBE)
    subset = _update_subset!(FEBE)
    coefficients = _update_coefficients!(FEBE)
    det = FEBE.L2G.det # 1 alloc
    cvals = FEBE.cvals
    refbasisvals = FEBE.refbasisvals
    fill!(cvals, 0)
    for i in 1:size(cvals, 3), dof_i in 1:size(cvals, 2)
        for k in 1:size(L2GM, 1)
            for l in 1:size(L2GM, 2)
                cvals[k, dof_i, i] += L2GM[k, l] * refbasisvals[i][subset[dof_i], l]
            end
            cvals[k, dof_i, i] *= coefficients[k, dof_i] / det
        end
    end
    return nothing
end


# IDENTITYCOMPONENT HDIV
function update_basis!(FEBE::SingleFEEvaluator{<:Real, <:Real, <:Integer, <:IdentityComponent{c}, <:AbstractHdivFiniteElement}) where {c}
    L2GM = _update_piola!(FEBE)
    subset = _update_subset!(FEBE)
    coefficients = _update_coefficients!(FEBE)
    det = FEBE.L2G.det # 1 alloc
    cvals = FEBE.cvals
    refbasisvals = FEBE.refbasisvals
    fill!(cvals, 0)
    for i in 1:size(cvals, 3), dof_i in 1:size(cvals, 2)
        for l in 1:size(L2GM, 2)
            cvals[1, dof_i, i] += L2GM[c, l] * refbasisvals[i][subset[dof_i], l]
        end
        cvals[1, dof_i, i] *= coefficients[c, dof_i] / det
    end
    return nothing
end


# NORMALFLUX HDIV
function update_basis!(FEBE::SingleFEEvaluator{<:Real, <:Real, <:Integer, <:NormalFlux, <:AbstractHdivFiniteElement})
    xItemVolumes = FEBE.L2G.ItemVolumes
    cvals = FEBE.cvals
    refbasisvals = FEBE.refbasisvals
    for i in 1:size(cvals, 3), dof_i in 1:size(cvals, 2), k in 1:size(cvals, 1)
        cvals[k, dof_i, i] = refbasisvals[i][dof_i, k] / xItemVolumes[FEBE.citem[]]
    end
    return nothing
end


# function update_basis!(FEBE::SingleFEEvaluator{<:Real,<:Real,<:Integer,<:Jump{NormalFlux},<:AbstractHdivFiniteElement})
#     xItemVolumes = FEBE.L2G.ItemVolumes
#     cvals = FEBE.cvals
#     refbasisvals = FEBE.refbasisvals
#     for i = 1 : size(cvals,3), dof_i = 1 : size(refbasisvals[1],1), k = 1 : size(cvals,1)
#         cvals[k,dof_i,i] = refbasisvals[i][dof_i,k] / xItemVolumes[FEBE.citem[]]
#         cvals[k,size(refbasisvals[1],1)+dof_i,i] = -cvals[k,dof_i,i]
#     end
#     return nothing
# end


# DIVERGENCE HDIV
function update_basis!(FEBE::SingleFEEvaluator{<:Real, <:Real, <:Integer, <:Divergence, <:AbstractHdivFiniteElement})
    # update transformation
    _update_piola!(FEBE)
    subset = _update_subset!(FEBE)
    coefficients = _update_coefficients!(FEBE)
    det = FEBE.L2G.det # 1 alloc
    cvals = FEBE.cvals
    offsets2 = FEBE.offsets2
    refbasisderivvals = FEBE.refbasisderivvals
    fill!(cvals, 0)
    for i in 1:size(cvals, 3), dof_i in 1:size(cvals, 2)
        for j in 1:size(refbasisderivvals, 2)
            cvals[1, dof_i, i] += refbasisderivvals[subset[dof_i] + offsets2[j], j, i]
        end
        cvals[1, dof_i, i] *= coefficients[1, dof_i] / det
    end
    return nothing
end


# GRADIENT HDIV
function update_basis!(FEBE::SingleFEEvaluator{<:Real, <:Real, <:Integer, <:Gradient, <:AbstractHdivFiniteElement})
    L2GAinv = _update_trafo!(FEBE)
    L2GM = _update_piola!(FEBE)
    subset = _update_subset!(FEBE)
    coefficients = _update_coefficients!(FEBE)
    cvals = FEBE.cvals
    offsets = FEBE.offsets
    offsets2 = FEBE.offsets2
    refbasisderivvals = FEBE.refbasisderivvals
    fill!(cvals, 0)
    det = FEBE.L2G.det # 1 alloc
    for i in 1:size(cvals, 3), dof_i in 1:size(cvals, 2)
        for c in 1:size(L2GM, 1), k in 1:size(L2GAinv, 1)
            # compute duc/dxk
            for j in 1:size(L2GM, 2), m in 1:size(L2GAinv, 2)
                cvals[k + offsets[c], dof_i, i] += L2GAinv[k, m] * L2GM[c, j] * refbasisderivvals[subset[dof_i] + offsets2[j], m, i]
            end
            cvals[k + offsets[c], dof_i, i] *= coefficients[c, dof_i] / det
        end
    end
    return nothing
end


# CURL2D HDIV
function update_basis!(FEBE::SingleFEEvaluator{<:Real, <:Real, <:Integer, <:Curl2D, <:AbstractHdivFiniteElement})
    L2GAinv = _update_trafo!(FEBE)
    L2GM = _update_piola!(FEBE)
    subset = _update_subset!(FEBE)
    coefficients = _update_coefficients!(FEBE)
    cvals = FEBE.cvals
    offsets2 = FEBE.offsets2
    refbasisderivvals = FEBE.refbasisderivvals
    fill!(cvals, 0)
    det = FEBE.L2G.det # 1 alloc
    for i in 1:size(cvals, 3), dof_i in 1:size(cvals, 2)
        for j in 1:size(L2GM, 2), m in 1:size(L2GAinv, 2)
            cvals[1, dof_i, i] -= L2GAinv[2, m] * L2GM[1, j] * refbasisderivvals[subset[dof_i] + offsets2[j], m, i] * coefficients[1, dof_i] / det
            cvals[1, dof_i, i] += L2GAinv[1, m] * L2GM[2, j] * refbasisderivvals[subset[dof_i] + offsets2[j], m, i] * coefficients[2, dof_i] / det
        end
    end
    return nothing
end


# CURL3D HDIV
function update_basis!(FEBE::SingleFEEvaluator{<:Real, <:Real, <:Integer, <:Curl3D, <:AbstractHdivFiniteElement})
    L2GAinv = _update_trafo!(FEBE)
    L2GM = _update_piola!(FEBE)
    subset = _update_subset!(FEBE)
    coefficients = _update_coefficients!(FEBE)
    cvals = FEBE.cvals
    offsets2 = FEBE.offsets2
    refbasisderivvals = FEBE.refbasisderivvals
    fill!(cvals, 0)
    det = FEBE.L2G.det # 1 alloc
    for i in 1:size(cvals, 3), dof_i in 1:size(cvals, 2)
        for j in 1:size(L2GM, 2), m in 1:size(L2GAinv, 2)
            cvals[1, dof_i, i] -= L2GAinv[3, m] * L2GM[2, j] * refbasisderivvals[subset[dof_i] + offsets2[j], m, i] * coefficients[2, dof_i] / det
            cvals[1, dof_i, i] += L2GAinv[2, m] * L2GM[3, j] * refbasisderivvals[subset[dof_i] + offsets2[j], m, i] * coefficients[3, dof_i] / det
            cvals[2, dof_i, i] -= L2GAinv[1, m] * L2GM[3, j] * refbasisderivvals[subset[dof_i] + offsets2[j], m, i] * coefficients[3, dof_i] / det
            cvals[2, dof_i, i] += L2GAinv[3, m] * L2GM[1, j] * refbasisderivvals[subset[dof_i] + offsets2[j], m, i] * coefficients[1, dof_i] / det
            cvals[3, dof_i, i] -= L2GAinv[2, m] * L2GM[1, j] * refbasisderivvals[subset[dof_i] + offsets2[j], m, i] * coefficients[1, dof_i] / det
            cvals[3, dof_i, i] += L2GAinv[1, m] * L2GM[2, j] * refbasisderivvals[subset[dof_i] + offsets2[j], m, i] * coefficients[2, dof_i] / det
        end
    end
    return nothing
end
