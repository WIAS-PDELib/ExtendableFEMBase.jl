### A Pluto.jl notebook ###
# v0.19.35

using Markdown
using InteractiveUtils

# ╔═╡ 014c3380-b361-11ed-273e-37541f1ed74f
begin
	using ExtendableFEMBase
	using ForwardDiff
	using DiffResults
	using LinearAlgebra
	using ExtendableGrids
	using ExtendableSparse
	using GridVisualize
	using PlutoVista
end

# ╔═╡ 3e8fd4f5-b6be-4019-9c76-a80cc985b70e
md"""
# Tutorial notebook: Navier--Stokes problem

Consider the Navier-Stokes problem that seeks ``u`` and ``p`` such that
```math
\begin{aligned}
	- \mu \Delta u + (u \cdot \nabla) u + \nabla p &= f\\
			\mathrm{div}(u) & = 0.
\end{aligned}
```

The weak formulation seeks ``u \in V := H^1_0(\Omega)`` and ``p \in Q := L^2_0(\Omega)`` such that
```math
\begin{aligned}
	\mu (\nabla u, \nabla v) + ((u \cdot \nabla) u, v) - (p, \mathrm{div}(v)) & = (f, v)
	& \text{for all } v \in V\\
	(q, \mathrm{div}(u)) & = 0
	& \text{for all } q \in Q\\
\end{aligned}
```
This tutorial notebook compute a planar lattice flow with inhomogeneous Dirichlet boundary conditions (which requires some modification above). Newton's method with automatic differentation is used to handle the nonlinear convection term.
"""

# ╔═╡ b5119656-dbf0-4249-a8bc-eac47deb8a43
md"""
This is a plot of the computed velocity and pressure:
"""

# ╔═╡ 6160364e-6ab2-475d-9345-3436e6f2b3e3
begin
	## PDE data
	const μ = 1e-2
	function f!(fval, x, t) # right-hand side
		fval[1] = 8.0*π*π*μ * exp(-8.0*π*π*μ*t) * sin(2.0*π*x[1])*sin(2.0*π*x[2])
		fval[2] = 8.0*π*π*μ * exp(-8.0*π*π*μ*t) * cos(2.0*π*x[1])*cos(2.0*π*x[2])
		return nothing
	end

	# exact velocity (for boundary data and error calculation)
	function u!(uval, qpinfo) 
		x = qpinfo.x
		t = qpinfo.time
		uval[1] = exp(-8.0*π*π*μ*t) * sin(2.0*π*x[1])*sin(2.0*π*x[2])
		uval[2] = exp(-8.0*π*π*μ*t) * cos(2.0*π*x[1])*cos(2.0*π*x[2])
		return nothing
	end
	
	## discretization parameters
	const nref = 5
	const teval = 0
	const order = 2
	
	## prepare error calculation
	function p!(pval,x,t) # exact pressure (for error calculation)
		pval[1] = exp(-16*pi*pi*μ*t)*(cos(4*pi*x[1])-cos(4*pi*x[2]))/4
		return nothing
	end
end

# ╔═╡ 8e0a11e7-7998-403a-8484-7d0180ea40b2
begin
	## create grid
	X = LinRange(0,1,2^nref+1)
	Y = LinRange(0,1,2^nref+1)
	println("Creating grid...")
	@time xgrid = simplexgrid(X,Y)
	println("Preparing FaceNodes...")
	@time xgrid[FaceNodes]
	println("Preparing CellVolumes...")
    @time xgrid[CellVolumes]
	xgrid
end

# ╔═╡ 7e2a082f-5457-47ac-96d5-b688a70a5506
begin
	## create finite element space (Taylor--Hood)
    FETypes = [H1Pk{2,2,order}, H1Pk{1,2,order-1}]

	## prepare finite element space and dofmaps
	println("Creating FESpace...")
    @time FES = [FESpace{FETypes[1]}(xgrid; name = "velocity space"),                  	             FESpace{FETypes[2]}(xgrid; name = "pressure space")]
	FES
end

# ╔═╡ 6dcca265-b8fd-4043-a2a7-4ff9bf579dd5
function prepare_assembly!(A, b, FESu, FESp, Solution, f, μ = 1)

	A = A.entries
	b = b.entries
	Solution = Solution.entries
    xgrid = FESu.xgrid
    EG = xgrid[UniqueCellGeometries][1]
    FEType_u = eltype(FESu)
    FEType_p = eltype(FESp)
    L2G = L2GTransformer(EG, xgrid, ON_CELLS)
    cellvolumes = xgrid[CellVolumes]
    ncells::Int = num_cells(xgrid)

    ## dofmap
    CellDofs_u = FESu[ExtendableFEMBase.CellDofs]
    CellDofs_p = FESp[ExtendableFEMBase.CellDofs]
	offset_p = FESu.ndofs
	
    ## quadrature formula
	qf = QuadratureRule{Float64, EG}(3*get_polynomialorder(FEType_u, EG)-1)
	weights::Vector{Float64} = qf.w
	xref::Vector{Vector{Float64}} = qf.xref
	nweights::Int = length(weights)
	
	## FE basis evaluator
	FEBasis_∇u = FEEvaluator(FESu, Gradient, qf)
	∇uvals = FEBasis_∇u.cvals
	FEBasis_idu = FEEvaluator(FESu, Identity, qf)
	iduvals = FEBasis_idu.cvals
	FEBasis_idp = FEEvaluator(FESp, Identity, qf)
	idpvals = FEBasis_idp.cvals

	## prepare automatic differentation of convection operator
	function operator!(result, input)
		# result = (u ⋅ ∇)u
		result[1] = input[1]*input[3]+input[2]*input[4]
		result[2] = input[1]*input[5]+input[2]*input[6]
	end
    result = Vector{Float64}(undef,2)
    input = Vector{Float64}(undef,6)
    tempV = zeros(Float64, 2)
    Dresult = DiffResults.JacobianResult(result, input)
    cfg = ForwardDiff.JacobianConfig(operator!, result, input, ForwardDiff.Chunk{6}())
    jac = DiffResults.jacobian(Dresult)
    value = DiffResults.value(Dresult)
	
   
    ## ASSEMBLY LOOP
    function barrier(EG, L2G::L2GTransformer, linear::Bool, nonlinear::Bool)
		## barrier function to avoid allocations caused by L2G
		
    	ndofs4cell_u::Int = get_ndofs(ON_CELLS, FEType_u, EG)
    	ndofs4cell_p::Int = get_ndofs(ON_CELLS, FEType_p, EG)
    	Aloc = zeros(Float64, ndofs4cell_u, ndofs4cell_u)
    	Bloc = zeros(Float64, ndofs4cell_u, ndofs4cell_p)
    	dof_j::Int, dof_k::Int = 0, 0
		fval::Vector{Float64} = zeros(Float64,2)
        x::Vector{Float64} = zeros(Float64, 2)
		
        for cell = 1 : ncells
			## update FE basis evaluators
	        update_basis!(FEBasis_∇u, cell)
	        update_basis!(FEBasis_idu, cell)
	        update_basis!(FEBasis_idp, cell) 
	
			## assemble local stiffness matrix (symmetric)
			if (linear)
		        for j = 1 : ndofs4cell_u, k = 1 : ndofs4cell_u
					temp = 0
					for qp = 1 : nweights
						temp += weights[qp] * dot(view(∇uvals,:,j,qp), view(∇uvals,:,k,qp))
					end
					Aloc[k,j] = μ * temp
		        end

				## assemble div-pressure coupling
		        for j = 1 : ndofs4cell_u, k = 1 : ndofs4cell_p
					temp = 0
					for qp = 1 : nweights
						temp -= weights[qp] * (∇uvals[1,j,qp] + ∇uvals[4,j,qp]) * 
						idpvals[1,k,qp]
					end
					Bloc[j,k] = temp
		        end
		        Bloc .*= cellvolumes[cell]
				
				## assemble right-hand side
	            update_trafo!(L2G, cell)
	            for j = 1 : ndofs4cell_u
	                ## right-hand side
	                temp = 0
	                for qp = 1 : nweights
	                    ## get global x for quadrature point
	                    eval_trafo!(x, L2G, xref[qp])
	                    ## evaluate (f(x), v_j(x))
						f!(fval, x, teval)
	                    temp += weights[qp] * dot(view(iduvals,: , j, qp), fval)
	                end
					## write into global vector
	                dof_j = CellDofs_u[j, cell]
	                b[dof_j] += temp * cellvolumes[cell]
	            end
			end

			## assemble nonlinear term
			if (nonlinear)
				for qp = 1 : nweights
					fill!(input,0)
					for j = 1 : ndofs4cell_u
						dof_j = CellDofs_u[j, cell]
						for d = 1 : 2
							input[d] += Solution[dof_j] * iduvals[d,j,qp]
						end
						for d = 1 : 4
							input[2+d] += Solution[dof_j] * ∇uvals[d,j,qp]
						end
					end
					
                	## evaluate jacobian
					ForwardDiff.chunk_mode_jacobian!(Dresult, operator!, result, input, cfg)
					
	                # update matrix
	                for j = 1 : ndofs4cell_u
	                    # multiply ansatz function with local jacobian
						fill!(tempV,0)
						for d = 1 : 2
							tempV[1] += jac[1,d] * iduvals[d,j,qp]
							tempV[2] += jac[2,d] * iduvals[d,j,qp]
						end
						for d = 1 : 4
							tempV[1] += jac[1,2+d] * ∇uvals[d,j,qp]
							tempV[2] += jac[2,2+d] * ∇uvals[d,j,qp]
						end
	
	                    # multiply test function operator evaluation
	                    for k = 1 : ndofs4cell_u
	                        Aloc[k,j] += dot(tempV,view(iduvals,:,k,qp)) * weights[qp]
	                    end
	                end 
	
	                # update rhs
	                mul!(tempV, jac, input)
	                tempV .-= value
	                for j = 1 : ndofs4cell_u
                		dof_j = CellDofs_u[j, cell]
                		b[dof_j] += dot(tempV, view(iduvals,:,j,qp)) * weights[qp] * cellvolumes[cell]
	                end
				end
			end
			
			## add local matrices to global matrix
	        Aloc .*= cellvolumes[cell]
	        for j = 1 : ndofs4cell_u
	            dof_j = CellDofs_u[j, cell]
	            for k = 1 : ndofs4cell_u
	                dof_k = CellDofs_u[k, cell]
	                rawupdateindex!(A, +, Aloc[j,k], dof_j, dof_k)
	            end
				if (linear)
					for k = 1 : ndofs4cell_p
						dof_k = CellDofs_p[k, cell] + offset_p
		                rawupdateindex!(A, +, Bloc[j,k], dof_j, dof_k) 
		                rawupdateindex!(A, +, Bloc[j,k], dof_k, dof_j)
					end
				end
	        end
	        fill!(Aloc, 0)
	        fill!(Bloc, 0)
        end
    end

	function update_system!(linear::Bool, nonlinear::Bool)
	    barrier(EG, L2G, linear, nonlinear)
	    flush!(A)
	end
	update_system!
end

# ╔═╡ 7299586b-6859-41af-bcd0-1a0daa454c81
function solve_stokes_lowlevel(FES, μ, f!)
	
	println("Initializing system...")
	Solution = FEVector(FES)
	A = FEMatrix(FES)
	b = FEVector(FES)
	@time update_system! = prepare_assembly!(A, b, FES[1], FES[2], Solution, f!, μ)
	@time update_system!(true, false)
	Alin = deepcopy(A) # = keep linear part of system matrix
	blin = deepcopy(b) # = keep linear part of right-hand side
	
	println("Pepare boundary conditions...")
	@time begin
		u_init = FEVector(FES)
		interpolate!(u_init[1], u!; time = teval)
		
		fixed_dofs = [size(A.entries,1)] # fix one pressure dof = last dof
		BFaceDofs::Adjacency{Int32} = FES[1][ExtendableFEMBase.BFaceDofs]
		nbfaces::Int = num_sources(BFaceDofs)
		AM::ExtendableSparseMatrix{Float64,Int64} = A.entries
		dof_j::Int = 0
		for bface = 1 : nbfaces
			for j = 1 : num_targets(BFaceDofs,1)
				dof_j = BFaceDofs[j, bface]
				push!(fixed_dofs, dof_j)
			end
		end
	end

	
	for it = 1 : 20
	    ## solve
		println("\nITERATION $it\n=============")
		println("Solving linear system...")
	    @time copyto!(Solution.entries, A.entries \ b.entries)
		res = A.entries.cscmatrix * Solution.entries .- b.entries
		for dof in fixed_dofs
			res[dof] = 0
		end
		linres = norm(res)
		println("linear residual = $linres")
		
		fill!(A.entries.cscmatrix.nzval,0)
		fill!(b.entries,0)
		println("Updating linear system...")
		@time begin
			update_system!(false,true)
			A.entries.cscmatrix += Alin.entries.cscmatrix
			b.entries .+= blin.entries
		end
		
	    ## fix boundary dofs
		for dof in fixed_dofs
			AM[dof,dof] = 1e60
			b.entries[dof] = 1e60 * u_init.entries[dof]
		end
		ExtendableSparse.flush!(A.entries)

		## calculate nonlinear residual
		res = A.entries.cscmatrix * Solution.entries .- b.entries
		for dof in fixed_dofs
			res[dof] = 0
		end
		nlres = norm(res)
		println("nonlinear residual = $nlres")
		if nlres < max(1e-12, 20*linres)
			break
		end
	end

	return Solution, u_init
end

# ╔═╡ 0785624c-e95a-4e76-8071-2eea591091e0
begin
	## call low level solver
	sol, u_init = solve_stokes_lowlevel(FES, μ, f!)
	sol
end

# ╔═╡ 6f7a1407-dfb3-497b-8c57-8efac6592194
[tricontour(xgrid[Coordinates],xgrid[CellNodes],nodevalues(sol[1]; abs = true)[:]; levels = 5),
tricontour(xgrid[Coordinates],xgrid[CellNodes],nodevalues(sol[2])[:]; levels = 5)
]


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DiffResults = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
ExtendableFEMBase = "12fb9182-3d4c-4424-8fd1-727a0899810c"
ExtendableGrids = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
ExtendableSparse = "95c220a8-a1cf-11e9-0c77-dbfce5f500b3"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
GridVisualize = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoVista = "646e1f28-b900-46d7-9d87-d554eb38a413"

[compat]
DiffResults = "~1.1.0"
ExtendableFEMBase = "~0.0.16"
ExtendableGrids = "~1.2.2"
ExtendableSparse = "~1.2.1"
ForwardDiff = "~0.10.36"
GridVisualize = "~1.5.0"
PlutoVista = "~1.0.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.1"
manifest_format = "2.0"
project_hash = "4784243049be372dfc5617ed3f1f0585bc622af2"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "793501dcd3fa7ce8d375a2c878dca2296232686e"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cde29ddf7e5726c9fb511f340244ea3481267608"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.7.2"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bijections]]
git-tree-sha1 = "c9b163bd832e023571e86d0b90d9de92a9879088"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.6"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "2118cb2765f8197b08e5958cdd17c165427425ee"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.19.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "c0216e792f518b39b22212127d4a84dc31e4e386"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.5"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "cd67fc487743b2f0fd4380d4cbd3a24660d0eec8"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.3"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "886826d76ea9e72b35fcd000e535588f7b60f21d"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.10.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c53fc348ca4d40d7b371e71fd52251839080cbc9"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.4"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.ElasticArrays]]
deps = ["Adapt"]
git-tree-sha1 = "e1c40d78de68e9a2be565f0202693a158ec9ad85"
uuid = "fdbdab4c-e67f-52f5-8c3f-e7b388dad3d4"
version = "1.2.11"

[[deps.ExtendableFEMBase]]
deps = ["DiffResults", "DocStringExtensions", "ExtendableGrids", "ExtendableSparse", "ForwardDiff", "LinearAlgebra", "Printf", "SparseArrays", "Term", "UnicodePlots"]
git-tree-sha1 = "c355ee4adc8187513e654c8778610832f6a8ad55"
uuid = "12fb9182-3d4c-4424-8fd1-727a0899810c"
version = "0.0.16"

[[deps.ExtendableGrids]]
deps = ["AbstractTrees", "Bijections", "Dates", "DocStringExtensions", "ElasticArrays", "InteractiveUtils", "LinearAlgebra", "Printf", "Random", "Requires", "SparseArrays", "StaticArrays", "StatsBase", "Test", "UUIDs", "WriteVTK"]
git-tree-sha1 = "25c9892fe1f0cc5386b947713a30a46686f7add7"
uuid = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
version = "1.2.2"

    [deps.ExtendableGrids.extensions]
    ExtendableGridsGmshExt = "Gmsh"

    [deps.ExtendableGrids.weakdeps]
    Gmsh = "705231aa-382f-11e9-3f0c-b7cb4346fdeb"

[[deps.ExtendableSparse]]
deps = ["DocStringExtensions", "ILUZero", "LinearAlgebra", "Printf", "Requires", "SparseArrays", "Sparspak", "StaticArrays", "SuiteSparse", "Test"]
git-tree-sha1 = "c5bef30b6a553bc5130129dee9d6b9dfab7cfd46"
uuid = "95c220a8-a1cf-11e9-0c77-dbfce5f500b3"
version = "1.2.1"

    [deps.ExtendableSparse.extensions]
    ExtendableSparseAlgebraicMultigridExt = "AlgebraicMultigrid"
    ExtendableSparseIncompleteLUExt = "IncompleteLU"
    ExtendableSparsePardisoExt = "Pardiso"

    [deps.ExtendableSparse.weakdeps]
    AlgebraicMultigrid = "2169fc97-5a83-5252-b627-83903c6c433c"
    IncompleteLU = "40713840-3770-5561-ab4c-a76e7d0d7895"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"

[[deps.Extents]]
git-tree-sha1 = "2140cd04483da90b2da7f99b2add0750504fc39c"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.2"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "5b93957f6dcd33fc343044af3d48c215be2562f1"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.9.3"

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

    [deps.FillArrays.weakdeps]
    PDMats = "90014a1f-27ba-587c-ab20-58faa44d9150"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "2d6ca471a6c7b536127afccfa7564b5b39227fe0"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.5"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "d53480c0793b13341c40199190f92c611aa2e93c"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.2"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "424a5a6ce7c5d97cca7bcc4eac551b97294c54af"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.9"

[[deps.GridVisualize]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "ElasticArrays", "ExtendableGrids", "GeometryBasics", "GridVisualizeTools", "HypertextLiteral", "Interpolations", "IntervalSets", "LinearAlgebra", "Observables", "OrderedCollections", "Printf", "StaticArrays"]
git-tree-sha1 = "f88733a32e49542e3237d7e03ddc77d7c79a1825"
uuid = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
version = "1.5.0"

    [deps.GridVisualize.weakdeps]
    CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
    GLMakie = "e9467ef8-e4e7-5192-8a1a-b1aee30e663a"
    Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
    PlutoVista = "646e1f28-b900-46d7-9d87-d554eb38a413"
    PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee"

[[deps.GridVisualizeTools]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "StaticArraysCore"]
git-tree-sha1 = "e111f256aa000c4e4662d1119281b751aa66dc37"
uuid = "5573ae12-3b76-41d9-b48c-81d0b6e61cc5"
version = "1.1.0"

[[deps.Highlights]]
deps = ["DocStringExtensions", "InteractiveUtils", "REPL"]
git-tree-sha1 = "0341077e8a6b9fc1c2ea5edc1e93a956d2aec0c7"
uuid = "eafb193a-b7ab-5a9e-9068-77385905fa72"
version = "0.5.2"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.ILUZero]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "b007cfc7f9bee9a958992d2301e9c5b63f332a90"
uuid = "88f59080-6952-5380-9ea5-54057fb9a43f"
version = "0.2.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "274ad8005db8b4ef6cc46d1392927083405813c2"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.0"

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

    [deps.Interpolations.weakdeps]
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.IntervalSets]]
deps = ["Dates", "Random"]
git-tree-sha1 = "3d8866c029dd6b16e69e0d4a939c4dfcb98fac47"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.8"
weakdeps = ["Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsStatisticsExt = "Statistics"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "274c38bd733f9d29036d0a73658fff1dc1d3a065"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.9.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.LightXML]]
deps = ["Libdl", "XML2_jll"]
git-tree-sha1 = "3a994404d3f6709610701c7dabfc03fed87a81f8"
uuid = "9c8b4983-aa76-5018-a973-4c85ecc9e179"
version = "0.9.1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.MarchingCubes]]
deps = ["PrecompileTools", "StaticArrays"]
git-tree-sha1 = "27d162f37cc29de047b527dab11a826dd3a650ad"
uuid = "299715c1-40a9-479a-aaf9-4a633d36f717"
version = "0.1.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.MyterialColors]]
git-tree-sha1 = "01d8466fb449436348999d7c6ad740f8f853a579"
uuid = "1c23619d-4212-4747-83aa-717207fae70f"
version = "0.3.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "2ac17d29c523ce1cd38e27785a7d23024853a4bb"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.10"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[deps.PlutoVista]]
deps = ["AbstractPlutoDingetjes", "ColorSchemes", "Colors", "DocStringExtensions", "GridVisualizeTools", "HypertextLiteral", "UUIDs"]
git-tree-sha1 = "5be7548065d668761814809e2c7ee33310a3d82f"
uuid = "646e1f28-b900-46d7-9d87-d554eb38a413"
version = "1.0.1"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "5165dfb9fd131cf0c6957a3a7605dede376e7b63"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "342cf4b449c299d8d1ceaf00b7a49f4fbc7940e7"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.9"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "fba11dbe2562eecdfcac49a05246af09ee64d055"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.8.1"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "1d77abd07f617c4868c33d4f5b9e1dbb2643c9cf"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.2"

[[deps.StructArrays]]
deps = ["Adapt", "ConstructionBase", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "0a3db38e4cce3c54fe7a71f831cd7b6194a54213"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.16"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "cb76cf677714c095e535e3501ac7954732aeea2d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Term]]
deps = ["AbstractTrees", "CodeTracking", "Dates", "Highlights", "InteractiveUtils", "Logging", "Markdown", "MyterialColors", "OrderedCollections", "Parameters", "PrecompileTools", "ProgressLogging", "REPL", "Tables", "UUIDs", "Unicode", "UnicodeFun"]
git-tree-sha1 = "ffac67f6fbcbb32027d924b93ba91b7633af9220"
uuid = "22787eb5-b846-44ae-b979-8e399b8463ab"
version = "2.0.5"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.UnicodePlots]]
deps = ["ColorSchemes", "ColorTypes", "Contour", "Crayons", "Dates", "LinearAlgebra", "MarchingCubes", "NaNMath", "PrecompileTools", "Printf", "Requires", "SparseArrays", "StaticArrays", "StatsBase"]
git-tree-sha1 = "b96de03092fe4b18ac7e4786bee55578d4b75ae8"
uuid = "b8865327-cd53-5732-bb35-84acbb429228"
version = "3.6.0"

    [deps.UnicodePlots.extensions]
    FreeTypeExt = ["FileIO", "FreeType"]
    ImageInTerminalExt = "ImageInTerminal"
    IntervalSetsExt = "IntervalSets"
    TermExt = "Term"
    UnitfulExt = "Unitful"

    [deps.UnicodePlots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    FreeType = "b38be410-82b0-50bf-ab77-7b57e271db43"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    Term = "22787eb5-b846-44ae-b979-8e399b8463ab"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.VTKBase]]
git-tree-sha1 = "c2d0db3ef09f1942d08ea455a9e252594be5f3b6"
uuid = "4004b06d-e244-455f-a6ce-a5f9919cc534"
version = "1.0.1"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "5f24e158cf4cee437052371455fe361f526da062"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.6"

[[deps.WriteVTK]]
deps = ["Base64", "CodecZlib", "FillArrays", "LightXML", "TranscodingStreams", "VTKBase"]
git-tree-sha1 = "41f0dc2a8f6fd860c266b91fd5cdf4fead65ae69"
uuid = "64499a7a-5c06-52f2-abe2-ccb03c286192"
version = "1.18.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "801cbe47eae69adc50f36c3caec4758d2650741b"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.2+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─014c3380-b361-11ed-273e-37541f1ed74f
# ╟─3e8fd4f5-b6be-4019-9c76-a80cc985b70e
# ╟─b5119656-dbf0-4249-a8bc-eac47deb8a43
# ╟─6f7a1407-dfb3-497b-8c57-8efac6592194
# ╠═6160364e-6ab2-475d-9345-3436e6f2b3e3
# ╠═8e0a11e7-7998-403a-8484-7d0180ea40b2
# ╠═7e2a082f-5457-47ac-96d5-b688a70a5506
# ╠═0785624c-e95a-4e76-8071-2eea591091e0
# ╠═7299586b-6859-41af-bcd0-1a0daa454c81
# ╠═6dcca265-b8fd-4043-a2a7-4ff9bf579dd5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
