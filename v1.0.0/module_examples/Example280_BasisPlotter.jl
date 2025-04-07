#=

# 280 : Basis-Plotter
([source code](@__SOURCE_URL__))

This example plots all the basis functions of a H1 finite element on Edge1D or Triangle2D
as unicode plots. This is the result with the default parameters (dim = 1, order = 3):

![](https://github.com/chmerdon/ExtendableFEMBase.jl/blob/master/docs/src/assets/example280.png?raw=true")

=#

module Example280_BasisPlotter

using ExtendableFEMBase
using ExtendableGrids

## everything is wrapped in a main function
function main(; dim = 1, order = 3)

	## generate two grids
	@assert dim in [1, 2] "dim must be 1 or 2"
	refgeom = dim == 1 ? Edge1D : Triangle2D
	xgrid = reference_domain(refgeom)

	## set finite element type and get some information
	FEType = H1Pk{1, dim, order}
	ndofs = get_ndofs(ON_CELLS, FEType, refgeom)
	FEType = H1Pk{ndofs, dim, order}

	## generate FEVector with ncomponents = ndofs
	## that will carry one basis function in each component
	FEFunc = FEVector(FESpace{FEType}(xgrid))
	coffsets = ExtendableFEMBase.get_local_coffsets(FEType, ON_CELLS, refgeom)
	for j ∈ 1:ndofs
		FEFunc[1][j+coffsets[j]] = 1
	end

    ## plot
	println(stdout, unicode_scalarplot(FEFunc[1]; title = "φ", ylim = (-0.5, 1), resolution = dim == 1 ? (40, 10) : (20, 15), nrows = order))
end
end
