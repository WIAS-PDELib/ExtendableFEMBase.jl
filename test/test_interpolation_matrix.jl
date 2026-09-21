function run_grid_interpolation_matrix_tests()
    # test interpolation of same space between refined grids
    function test_grid_matrix_computation(xgrid, FEType, order; broken::Bool = false, use_cellparents::Bool = false)
        u, ~ = exact_function(Val(dim_grid(xgrid)), order)

        source_FES = FESpace{FEType}(xgrid; broken)
        target_FES = FESpace{FEType}(uniform_refine(xgrid); broken)

        source_vector = FEVector(source_FES)
        interpolate!(source_vector[1], u; bonus_quadorder = order)

        target_vector = FEVector(target_FES)
        interpolate!(target_vector[1], u; bonus_quadorder = order)

        interpolation_matrix = compute_lazy_interpolation_jacobian(target_FES, source_FES; use_cellparents)
        matrix_interpolated_entries = interpolation_matrix * source_vector.entries

        return @test norm(target_vector.entries - matrix_interpolated_entries) < tolerance
    end

    @testset "Grid Interpolation Matrix Tests" begin
        println("\n")
        println("============================")
        println("Testing Grid Interpolation Matrices in 1D")
        println("============================")
        for_each_catalog([Edge1D], TestCatalog1D) do ~, xgrid, FEType, order, broken
            @info "Element: ($(FEType), $(order), broken = $(broken)) \n"
            test_grid_matrix_computation(xgrid, FEType, order; broken)
        end

        println("\n")
        println("============================")
        println("Testing Grid Interpolation Matrices in 2D")
        println("============================")
        for_each_catalog([Triangle2D, Parallelogram2D], TestCatalog2D; grid = testgrid_refdomain) do EG, xgrid, FEType, order, broken
            @info "Element: ($(EG), $(FEType), $(order), broken = $(broken)) \n"
            test_grid_matrix_computation(xgrid, FEType, order; broken)
        end

        println("\n")
        println("============================")
        println("Testing Grid Interpolation Matrices in 3D")
        println("============================")
        for_each_catalog([Tetrahedron3D, Parallelepiped3D], TestCatalog3D; grid = testgrid_refdomain) do EG, xgrid, FEType, order, broken
            @info "Element: ($(EG), $(FEType), $(order), broken = $(broken)) \n"
            test_grid_matrix_computation(xgrid, FEType, order; broken)
        end
    end

    return println("")
end

function run_space_interpolation_matrix_tests()
    # list of space pairs that should be interpolated into one another

    PairTestCatalog1D = [
        (L2P0{1}, H1P1{1}) => 0,
        (H1P1{1}, H1P2{1, 1}) => 1,
        (H1P2{1, 1}, H1P3{1, 1}) => 2,
        (H1Pk{1, 1, 3}, H1Pk{1, 1, 5}) => 3,
    ]

    PairTestCatalog2D = [
        (H1P1{2}, HDIVRT0{2}) => 0,
        (H1P2{2, 2}, HDIVRT0{2}) => 0,
        (H1P2{2, 2}, HDIVRT1{2}) => 1,
        (H1P2{2, 2}, HDIVRTk{2, 2}) => 2,
        (HDIVRT1{2}, HDIVBDM1{2}) => 1,
        (L2P0{2}, L2P1{2}) => 0,
        (H1P2B{2, 2}, H1BR{2}) => 1,
    ]

    PairTestCatalog3D = [
        (H1P1{3}, HDIVRT0{3}) => 0,
        (H1P2{3, 3}, HDIVRT0{3}) => 0,
        (H1P2{3, 3}, HDIVRT1{3}) => 1,
        (HDIVRT1{3}, HDIVBDM1{3}) => 1,
        (L2P0{3}, L2P1{3}) => 0,
        (H1P2B{3, 3}, H1BR{3}) => 1,
    ]

    # test interpolation for different elements on same grid
    function test_space_matrix_computation(xgrid, source_FEType, target_FEType, order; broken::Bool = false, use_cellparents::Bool = false)
        u, ~ = exact_function(Val(dim_grid(xgrid)), order)

        source_FES = FESpace{source_FEType}(xgrid)
        target_FES = FESpace{target_FEType}(xgrid; broken)

        source_vector = FEVector(source_FES)
        interpolate!(source_vector[1], u; bonus_quadorder = order)

        target_vector = FEVector(target_FES)
        interpolate!(target_vector[1], u; bonus_quadorder = order)

        interpolation_matrix = compute_lazy_interpolation_jacobian(target_FES, source_FES; use_cellparents)
        matrix_interpolated_entries = interpolation_matrix * source_vector.entries

        return @test norm(target_vector.entries - matrix_interpolated_entries) < tolerance
    end

    @testset "Space Interpolation Matrix Tests" begin
        println("\n")
        println("============================")
        println("Testing Space Interpolation Matrices in 1D")
        println("============================")
        xgrid = testgrid(Edge1D)
        for ((source_element, target_element), order) in PairTestCatalog1D
            @info "Element pair: ($(source_element), $(target_element)), order: $(order) \n"
            test_space_matrix_computation(xgrid, source_element, target_element, order; broken = false)
            test_space_matrix_computation(xgrid, source_element, target_element, order; broken = true)
        end

        println("\n")
        println("============================")
        println("Testing Space Interpolation Matrices in 2D")
        println("============================")
        for EG in [Triangle2D, Parallelogram2D]
            xgrid = testgrid_refdomain(EG)
            for ((source_element, target_element), order) in PairTestCatalog2D, broken in (false, true)
                @info "Element pair: ($(EG), $(source_element), $(target_element)), order: $(order), broken = $(broken) \n"
                if ExtendableFEMBase.isdefined(target_element, EG, broken) && ExtendableFEMBase.isdefined(source_element, EG, broken)
                    test_space_matrix_computation(xgrid, source_element, target_element, order; broken)

                    if (source_element, target_element) == (H1P1{2}, HDIVRT0{2}) && !broken
                        source_FES = FESpace{source_element}(xgrid; broken)
                        target_FES = FESpace{target_element}(xgrid; broken)
                        autodiff_matrix = compute_lazy_interpolation_jacobian(target_FES, source_FES; use_cellparents = false)
                        RT0_matrix = H1Pk_to_HDIVRT0_interpolator(target_FES, source_FES)

                        @test norm(autodiff_matrix' - RT0_matrix[1]) < tolerance
                    end
                else
                    @warn "($(target_element),$(order)) (broken = $(broken)) not defined on $(EG) (skipping test case)"
                end
            end
        end

        println("\n")
        println("============================")
        println("Testing Space Interpolation Matrices in 3D")
        println("============================")
        for EG in [Tetrahedron3D, Parallelepiped3D]
            xgrid = testgrid_refdomain(EG)
            for ((source_element, target_element), order) in PairTestCatalog3D, broken in (false, true)
                @info "Element pair: ($(EG), $(source_element), $(target_element)), order: $(order), broken = $(broken) \n"
                if ExtendableFEMBase.isdefined(target_element, EG, broken) && ExtendableFEMBase.isdefined(source_element, EG, broken)
                    test_space_matrix_computation(xgrid, source_element, target_element, order; broken)
                else
                    @warn "($(target_element),$(order)) (broken = $(broken)) not defined on $(EG) (skipping test case)"
                end
            end
        end
    end

    return println("")
end
