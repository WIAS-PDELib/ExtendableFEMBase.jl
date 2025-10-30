function run_grid_interpolation_matrix_tests()
    # list of FETypes that should be tested
    TestCatalog1D = [
        L2P0{1} => 0,
        H1P1{1} => 1,
        H1P2{1, 1} => 2,
        H1P3{1, 1} => 3,
        H1Pk{1, 1, 3} => 3,
        H1Pk{1, 1, 4} => 4,
        H1Pk{1, 1, 5} => 5,
    ]

    PairTestCatalog1D = [
        (L2P0{1}, H1P1{1}) => 0,
        (H1P1{1}, H1P2{1, 1}) => 1,
        (H1P2{1, 1}, H1P3{1, 1}) => 2,
        (H1Pk{1, 1, 3}, H1Pk{1, 1, 5}) => 3,
    ]

    TestCatalog2D = [
        HCURLN0{2} => 0,
        HCURLN1{2} => 1,
        HDIVRT0{2} => 0,
        HDIVRTk{2, 0} => 0,
        HDIVBDM1{2} => 1,
        HDIVRT1{2} => 1,
        HDIVRTk{2, 1} => 1,
        HDIVBDM2{2} => 2,
        HDIVRTk{2, 2} => 2,
        HDIVRTk{2, 3} => 3,
        HDIVRTk{2, 4} => 4,
        L2P0{2} => 0,
        L2P1{2} => 1,
        H1P1{2} => 1,
        H1Q1{2} => 1,
        H1CR{2} => 1,
        H1MINI{2, 2} => 1,
        H1P1TEB{2} => 1,
        H1BR{2} => 1,
        H1P2{2, 2} => 2,
        H1P2B{2, 2} => 2,
        H1Q2{2, 2} => 2,
        H1P3{2, 2} => 3,
        H1Pk{2, 2, 3} => 3,
        H1Pk{2, 2, 4} => 4,
        H1Pk{2, 2, 5} => 5,
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

    TestCatalog3D = [
        HCURLN0{3} => 0,
        HDIVRT0{3} => 0,
        HDIVBDM1{3} => 1,
        HDIVRT1{3} => 1,
        L2P0{3} => 0,
        H1P1{3} => 1,
        H1Q1{3} => 1,
        H1CR{3} => 1,
        H1MINI{3, 3} => 1,
        H1P1TEB{3} => 1,
        H1BR{3} => 1,
        H1P2{3, 3} => 2,
        H1P3{3, 3} => 3,
    ]

    PairTestCatalog3D = [
        (H1P1{3}, HDIVRT0{3}) => 0,
        (H1P2{3, 3}, HDIVRT0{3}) => 0,
        (H1P2{3, 3}, HDIVRT1{3}) => 1,
        (HDIVRT1{3}, HDIVBDM1{3}) => 1,
        (L2P0{3}, L2P1{3}) => 0,
        (H1P2B{3, 3}, H1BR{3}) => 1,
    ]

    # test interpolation of same space between refined grids
    function test_grid_matrix_computation(xgrid, FEType, order; broken::Bool = false, use_cellparents::Bool = false)
        u, ~ = exact_function(Val(dim_grid(xgrid)), order)

        source_FES = FESpace{FEType}(xgrid; broken)
        target_FES = FESpace{FEType}(uniform_refine(xgrid); broken)

        source_vector = FEVector(source_FES)
        interpolate!(source_vector[1], u; bonus_quadorder = order)

        target_vector = FEVector(target_FES)
        interpolate!(target_vector[1], u; bonus_quadorder = order)

        interpolation_matrix = compute_interpolation_jacobian(target_FES, source_FES; use_cellparents)
        matrix_interpolated_entries = interpolation_matrix * source_vector.entries

        return @test norm(target_vector.entries - matrix_interpolated_entries) < tolerance
    end

    @testset "Grid Interpolation Matrix Tests" begin
        println("\n")
        println("============================")
        println("Testing Grid Interpolation Matrices in 1D")
        println("============================")
        xgrid = testgrid(Edge1D)
        for (element, order) in TestCatalog1D
            @info "Element: ($(element), $(order)) \n"
            test_grid_matrix_computation(xgrid, element, order; broken = false)
            test_grid_matrix_computation(xgrid, element, order; broken = true)
        end

        println("\n")
        println("============================")
        println("Testing Grid Interpolation Matrices in 2D")
        println("============================")
        for EG in [Triangle2D, Parallelogram2D]
            xgrid = uniform_refine(reference_domain(EG), 1)
            for (element, order) in TestCatalog2D, broken in (false, true)
                @info "Element: ($(EG), $(element), $(order), broken = $(broken)) \n"
                if ExtendableFEMBase.isdefined(element, EG, broken)
                    test_grid_matrix_computation(xgrid, element, order; broken)
                else
                    @warn "($(element),$(order)) (broken = $(broken)) not defined on $(EG) (skipping test case)"
                end
            end
        end

        println("\n")
        println("============================")
        println("Testing Grid Interpolation Matrices in 3D")
        println("============================")
        for EG in [Tetrahedron3D, Parallelepiped3D]
            xgrid = uniform_refine(reference_domain(EG), 1)
            for (element, order) in TestCatalog3D, broken in (false, true)
                @info "Element: ($(EG), $(element), $(order), broken = $(broken)) \n"
                if ExtendableFEMBase.isdefined(element, EG, broken)
                    test_grid_matrix_computation(xgrid, element, order; broken)
                else
                    @warn "($(element),$(order)) (broken = $(broken)) not defined on $(EG) (skipping test case)"
                end
            end
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


    # V1: H1-conforming space, computes matrix into RT0 space
    function RT0_interpolator(V1::FESpace, VRT0::FESpace)
        FEType = get_FEType(V1)
        facedofs_V1::VariableTargetAdjacency{Int32} = V1[FaceDofs]
        dim = get_ncomponents(V1)
        EGF = dim == 2 ? Edge1D : Triangle2D
        ndofs = get_ndofs(ON_FACES, FEType, EGF)
        order = get_polynomialorder(FEType, EGF)
        face_basis = get_basis(ON_FACES, FEType, EGF)
        result = zeros(ndofs, dim)
        ref_integrate!(result, EGF, order, face_basis)

        F0 = FEMatrix(V1, VRT0)
        FE::ExtendableSparseMatrix{Float64, Int64} = F0.entries
        xgrid = V1.xgrid
        xFaceVolumes = xgrid[FaceVolumes]
        xFaceNormals = xgrid[FaceNormals]
        nfaces = num_sources(xFaceNormals)
        for face in 1:nfaces
            for dof1 in 1:ndofs
                FE[facedofs_V1[dof1, face], face] = dot(view(result, dof1, :), view(xFaceNormals, :, face)) * xFaceVolumes[face]
            end
        end
        flush!(FE)
        return F0
    end

    # test interpolation for different elements on same grid
    function test_space_matrix_computation(xgrid, source_FEType, target_FEType, order; broken::Bool = false, use_cellparents::Bool = false)
        u, ~ = exact_function(Val(dim_grid(xgrid)), order)

        source_FES = FESpace{source_FEType}(xgrid)
        target_FES = FESpace{target_FEType}(xgrid; broken)

        source_vector = FEVector(source_FES)
        interpolate!(source_vector[1], u; bonus_quadorder = order)

        target_vector = FEVector(target_FES)
        interpolate!(target_vector[1], u; bonus_quadorder = order)

        interpolation_matrix = compute_interpolation_jacobian(target_FES, source_FES; use_cellparents)
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
            xgrid = uniform_refine(reference_domain(EG), 1)
            for ((source_element, target_element), order) in PairTestCatalog2D, broken in (false, true)
                @info "Element pair: ($(EG), $(source_element), $(target_element)), order: $(order), broken = $(broken) \n"
                if ExtendableFEMBase.isdefined(target_element, EG, broken) && ExtendableFEMBase.isdefined(source_element, EG, broken)
                    test_space_matrix_computation(xgrid, source_element, target_element, order; broken)

                    if (source_element, target_element) == (H1P1{2}, HDIVRT0{2}) && !broken
                        source_FES = FESpace{source_element}(xgrid; broken)
                        target_FES = FESpace{target_element}(xgrid; broken)
                        autodiff_matrix = compute_interpolation_jacobian(target_FES, source_FES; use_cellparents = false)
                        RT0_matrix = RT0_interpolator(source_FES, target_FES)

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
            xgrid = uniform_refine(reference_domain(EG), 1)
            for ((source_element, target_element), order) in PairTestCatalog3D, broken in (false, true)
                @info "Element pair: ($(EG), $(source_element), $(target_element)), order: $(order), broken =$(broken) \n"
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
