function run_interpolator_tests()

    # list of FETypes that should be tested
    TestCatalog1D = [
        L2P0{1},
        H1P1{1},
        H1P2{1, 1},
        H1P3{1, 1},
        H1Pk{1, 1, 3},
        H1Pk{1, 1, 4},
        H1Pk{1, 1, 5},
    ]
    ExpectedOrders1D = [0, 1, 2, 3, 3, 4, 5]

    TestCatalog2D = [
        HCURLN0{2},
        HCURLN1{2},
        HDIVRT0{2},
        HDIVRTk{2, 0},
        HDIVBDM1{2},
        HDIVRT1{2},
        HDIVRTk{2, 1},
        HDIVBDM2{2},
        HDIVRTk{2, 2},
        HDIVRTk{2, 3},
        HDIVRTk{2, 4},
        L2P0{2},
        L2P1{2},
        H1P1{2},
        H1Q1{2},
        H1CR{2},
        H1MINI{2, 2},
        H1P1TEB{2},
        H1BR{2},
        H1P2{2, 2},
        H1P2B{2, 2},
        H1Q2{2, 2},
        H1P3{2, 2},
        H1Pk{2, 2, 3},
        H1Pk{2, 2, 4},
        H1Pk{2, 2, 5},
    ]
    ExpectedOrders2D = [0, 1, 0, 0, 1, 1, 1, 2, 2, 3, 4, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 5]

    TestCatalog3D = [
        HCURLN0{3},
        HDIVRT0{3},
        HDIVBDM1{3},
        HDIVRT1{3},
        L2P0{3},
        H1P1{3},
        H1Q1{3},
        H1CR{3},
        H1MINI{3, 3},
        H1P1TEB{3},
        H1BR{3},
        H1P2{3, 3},
        H1P3{3, 3},
    ]
    ExpectedOrders3D = [0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 2, 3]

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
        print("FEType = $FEType $(broken ? "broken" : "") $AT | ndofs = $(FES.ndofs) | order = $order")

        # interpolate
        Solution = FEVector(FES)
        interpolate!(Solution[1], u; bonus_quadorder = order)

        # compute error
        error = compute_error(Solution[1], u, order)
        println(" | error = $(norm(error, Inf))")
        return @test norm(error) < tolerance
    end

    @testset "Interpolations" begin
        println("\n")
        println("============================")
        println("Testing Interpolations in 1D")
        println("============================")
        xgrid = testgrid(Edge1D)
        for n in 1:length(TestCatalog1D)
            test_interpolation(xgrid, TestCatalog1D[n], ExpectedOrders1D[n])
            test_interpolation(xgrid, TestCatalog1D[n], ExpectedOrders1D[n], true)
        end
        println("\n")
        println("============================")
        println("Testing Interpolations in 2D")
        println("============================")
        for EG in [Triangle2D, Parallelogram2D]
            xgrid = uniform_refine(reference_domain(EG), 1)
            println("EG = $EG")
            for n in 1:length(TestCatalog2D), broken in (false, true)
                if ExtendableFEMBase.isdefined(TestCatalog2D[n], EG, broken)
                    test_interpolation(xgrid, TestCatalog2D[n], ExpectedOrders2D[n], broken)
                else
                    @warn "$(TestCatalog2D[n]) (broken = $broken) not defined on $EG (skipping test case)"
                end
            end
        end
        println("\n")
        println("============================")
        println("Testing Interpolations in 3D")
        println("============================")
        for EG in [Tetrahedron3D, Parallelepiped3D]
            xgrid = uniform_refine(reference_domain(EG), 1)
            for n in 1:length(TestCatalog3D), broken in (false, true)
                if ExtendableFEMBase.isdefined(TestCatalog3D[n], EG, broken)
                    test_interpolation(xgrid, TestCatalog3D[n], ExpectedOrders3D[n], broken)
                else
                    @warn "$(TestCatalog3D[n]) (broken = $broken) not defined on $EG (skipping test case)"
                end
            end
        end
    end
    return println("")
end
