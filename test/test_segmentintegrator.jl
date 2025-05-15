function test_segmentintegrator_nokernel()
    ## initial grid
    xgrid = grid_unitsquare(Triangle2D)

    ## Taylor--Hood FESpace
    FES = FESpace{H1P2{2, 2}}(xgrid)

    ## Hagen-Poiseuille flow
    function u(result, qpinfo)
        x = qpinfo.x
        result[1] = x[2] * (1.0 - x[2])
        return result[2] = 0.0
    end

    ## interpolate
    uh = FEVector(FES)
    interpolate!(uh[1], u)

    ## init segment integrator
    SI = SegmentIntegrator(Edge1D, [(1, Identity)])
    initialize!(SI, uh)

    ## integrate along line [1/4,1/4] to [3/4,1/4] in first triangle
    ## exact integral should be [3//32,0]
    result = zeros(Float64, 2)
    world = Array{Array{Float64, 1}, 1}([[1 // 4, 1 // 4], [3 // 4, 1 // 4]])
    bary = Array{Array{Float64, 1}, 1}([[1 // 4, 1 // 2], [3 // 4, 1 // 2]])
    integrate_segment!(result, SI, world, bary, 1)
    error1 = sqrt((result[1] - 3 // 32)^2 + result[2]^2)
    println("error1 without kernel = $error1 (result = $result)")

    ## integrate along line [1/2, 0] to [1/2, 1/2]
    ## exact integral should be [1//12, 0]
    world = Array{Array{Float64, 1}, 1}([[1 // 2, 0], [1 // 2, 1 // 2]])
    bary = Array{Array{Float64, 1}, 1}([[1 // 2, 0], [1 // 2, 1 // 1]])
    integrate_segment!(result, SI, world, bary, 1)
    error2 = sqrt((result[1] - 1 // 12)^2 + result[2]^2)
    println("error2 without kernel = $error2 (result = $result)")

    return max(error1, error2) < 1e-15
end

function test_segmentintegrator_withkernel()
    ## initial grid
    xgrid = grid_unitsquare(Triangle2D)

    ## Taylor--Hood FESpace
    FES = FESpace{H1P2{2, 2}}(xgrid)

    ## stagnation flow
    function u(result, qpinfo)
        x = qpinfo.x
        result[1] = x[1]
        return result[2] = -2 * x[2]
    end

    function multiply_r!(result, input, qpinfo)
        return result .= input * qpinfo.x[1]
    end

    ## interpolate
    uh = FEVector(FES)
    interpolate!(uh[1], u)

    ## init segment integrator
    SI = SegmentIntegrator(Edge1D, multiply_r!, [(1, Identity)]; bonus_quadorder = 1)
    initialize!(SI, uh)

    @show xgrid[Coordinates], xgrid[CellNodes]

    L2G = L2GTransformer(Triangle2D, xgrid, ON_CELLS)
    update_trafo!(L2G, 1)

    b1 = [0 // 1, 1 // 2]
    b2 = [1 // 2, 1 // 2]
    x1 = zeros(Float64, 2)
    x2 = zeros(Float64, 2)
    eval_trafo!(x1, L2G, b1)
    eval_trafo!(x2, L2G, b2)

    ## integrate along line [1/4,1/4] to [3/4,1/4] in first triangle
    ## exact integral should be [13//96, -1//8]
    result = zeros(Float64, 2)
    world = Array{Array{Float64, 1}, 1}([x1, x2])
    bary = Array{Array{Float64, 1}, 1}([b1, b2])
    integrate_segment!(result, SI, world, bary, 1)
    error1 = sqrt((result[1] - 13 // 96)^2 + (result[2] + 1 // 8)^2)
    println("error with kernel = $error1 (result = $result)")
    return error1 < 1.0e-15
end

function run_segmentintegrator_tests()

    return @testset "SegmentIntegrator" begin
        println("\n")
        println("=========================")
        println("Testing SegmentIntegrator")
        println("=========================")

        @test test_segmentintegrator_nokernel()
        @test test_segmentintegrator_withkernel()
    end
end
