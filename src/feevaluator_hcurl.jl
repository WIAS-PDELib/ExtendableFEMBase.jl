# IDENTITY HCURL
function update_basis!(FEBE::SingleFEEvaluator{<:Real, <:Real, <:Integer, <:Identity, <:AbstractHcurlFiniteElement})
    L2GAinv = _update_trafo!(FEBE)
    subset = _update_subset!(FEBE)
    coefficients = _update_coefficients!(FEBE)
    refbasisvals = FEBE.refbasisvals
    cvals = FEBE.cvals
    for i in 1:size(cvals, 3), dof_i in 1:size(cvals, 2), k in 1:size(cvals, 1)
        cvals[k, dof_i, i] = dot(view(L2GAinv, k, :), view(refbasisvals[i], subset[dof_i], :)) * coefficients[k, dof_i]
    end
    return nothing
end


# TANGENTLFLUX HCURL
# (just divide by face volume)
function update_basis!(FEBE::SingleFEEvaluator{<:Real, <:Real, <:Integer, <:TangentFlux, <:AbstractHcurlFiniteElement})
    subset = _update_subset!(FEBE)
    refbasisvals = FEBE.refbasisvals
    xItemVolumes = FEBE.L2G.ItemVolumes
    item = FEBE.citem[]
    cvals = FEBE.cvals
    for i in 1:size(cvals, 3), dof_i in 1:size(cvals, 2), k in 1:size(cvals, 1)
        cvals[k, dof_i, i] = refbasisvals[i][subset[dof_i], k] / xItemVolumes[item]
    end
    return nothing
end


# CURL2D HCURL
# covariant Piola transformation preserves curl2D (up to a factor 1/det(A))
function update_basis!(FEBE::SingleFEEvaluator{<:Real, <:Real, <:Integer, <:Curl2D, <:AbstractHcurlFiniteElement})
    _update_trafo!(FEBE)
    subset = _update_subset!(FEBE)
    coefficients = _update_coefficients!(FEBE)
    refbasisderivvals = FEBE.refbasisderivvals
    offsets2 = FEBE.offsets2
    cvals = FEBE.cvals
    fill!(cvals, 0)
    det = FEBE.L2G.det
    for i in 1:size(cvals, 3), dof_i in 1:size(cvals, 2)
        cvals[1, dof_i, i] = refbasisderivvals[subset[dof_i] + offsets2[2], 1, i] # du2/dx
        cvals[1, dof_i, i] -= refbasisderivvals[subset[dof_i] + offsets2[1], 2, i] # -du1/dy
        cvals[1, dof_i, i] *= coefficients[1, dof_i] / det
    end
    return nothing
end


# CURL3D HCURL
# covariant Piola transformation preserves curl3D (up to a factor 1/det(A))
function update_basis!(FEBE::SingleFEEvaluator{<:Real, <:Real, <:Integer, <:Curl3D, <:AbstractHcurlFiniteElement})
    L2GM = _update_piola!(FEBE)
    subset = _update_subset!(FEBE)
    coefficients = _update_coefficients!(FEBE)
    refbasisderivvals = FEBE.refbasisderivvals
    offsets2 = FEBE.offsets2
    cvals = FEBE.cvals
    fill!(cvals, 0)
    det = FEBE.L2G.det
    for i in 1:size(cvals, 3), dof_i in 1:size(cvals, 2)
        for k in 1:3
            cvals[k, dof_i, i] += L2GM[k, 1] * refbasisderivvals[subset[dof_i] + offsets2[3], 2, i] # du3/dx2
            cvals[k, dof_i, i] -= L2GM[k, 1] * refbasisderivvals[subset[dof_i] + offsets2[2], 3, i] # - du2/dx3
            cvals[k, dof_i, i] += L2GM[k, 2] * refbasisderivvals[subset[dof_i] + offsets2[1], 3, i] # du3/dx1
            cvals[k, dof_i, i] -= L2GM[k, 2] * refbasisderivvals[subset[dof_i] + offsets2[3], 1, i] # - du1/dx3
            cvals[k, dof_i, i] += L2GM[k, 3] * refbasisderivvals[subset[dof_i] + offsets2[2], 1, i] # du2/dx1
            cvals[k, dof_i, i] -= L2GM[k, 3] * refbasisderivvals[subset[dof_i] + offsets2[1], 2, i] # - du1/dx2
            cvals[k, dof_i, i] *= coefficients[k, dof_i] / det
        end
    end
    return nothing
end
