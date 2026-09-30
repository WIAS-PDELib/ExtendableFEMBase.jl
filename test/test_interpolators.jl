function run_interpolator_tests()

    ## function that computes errors at enough quadrature points for polynomial of degree order
    function compute_error(uh::FEVectorBlock, u::Function, order = get_polynomialorder(get_FEType(uh), uh.FES.xgrid[CellGeometries][1]))
        xgrid = uh.FES.xgrid
        FES = uh.FES
        EGs = xgrid[UniqueCellGeometries]
        ncomponents = get_ncomponents(uh)
        cells4eg = xgrid[ExtendableGrids.CellAssemblyGroups]
        celldofs = FES[CellDofs]
        error = zeros(Float64, ncomponents, num_cells(xgrid))
        uhval = zeros(Float64, ncomponents)
        uval = zeros(Float64, ncomponents)
        for (j, EG) in enumerate(EGs)
            cells = view(cells4eg, :, j)
            L2G = L2GTransformer(EG, xgrid, ON_CELLS)
            QP = QPInfos(xgrid)
            qf = VertexRule(EG, order)
            FEB = FEEvaluator(FES, Identity, qf)
            show(devnull, FEB)
            for cell::Int in cells
                update_trafo!(L2G, cell)
                update_basis!(FEB, cell)
                for (qp, weight) in enumerate(qf.w)
                    ## evaluate uh
                    fill!(uhval, 0)
                    eval_febe!(uhval, FEB, view(uh.entries, view(celldofs, :, cell)), qp)

                    ## evaluate u
                    fill!(uval, 0)
                    eval_trafo!(QP.x, L2G, qf.xref[qp])
                    u(uval, QP)

                    ## evaluate error
                    view(error, :, cell) .+= abs.(uval - uhval)
                end
            end
        end
        return error
    end

    function test_interpolation(xgrid, FEType, order, broken::Bool = false)

        u, ~ = exact_function(Val(size(xgrid[Coordinates], 1)), order)

        # choose FE and generate FESpace
        FES = FESpace{FEType}(xgrid; broken = broken)
        AT = ON_CELLS

        # interpolate
        Solution = FEVector(FES)
        interpolate!(Solution[1], u; bonus_quadorder = order)
        show(devnull, Solution)

        # compute error
        error = compute_error(Solution[1], u, order + 1)
        println("FEType = $FEType $(broken ? "broken" : "") $AT | ndofs = $(FES.ndofs) | order = $order | error = $(norm(error, Inf))")
        return @test norm(error) < tolerance
    end

    @testset "Interpolations" begin
        println("\n")
        println("============================")
        println("Testing Interpolations in 1D")
        println("============================")
        for_each_catalog([Edge1D], TestCatalog1D) do ~, xgrid, FEType, order, broken
            test_interpolation(xgrid, FEType, order, broken)
        end
        println("\n")
        println("============================")
        println("Testing Interpolations in 2D")
        println("============================")
        for_each_catalog([Triangle2D, Parallelogram2D], TestCatalog2D; grid = testgrid_refdomain) do ~, xgrid, FEType, order, broken
            test_interpolation(xgrid, FEType, order, broken)
        end
        println("\n")
        println("============================")
        println("Testing Interpolations in 3D")
        println("============================")
        for_each_catalog([Tetrahedron3D, Parallelepiped3D], TestCatalog3D; grid = testgrid_refdomain) do ~, xgrid, FEType, order, broken
            test_interpolation(xgrid, FEType, order, broken)
        end
    end
    return println("")
end
