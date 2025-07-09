function run_pointevaluator_tests()
    return @testset "PointEvaluator" begin
        println("\n")
        println("======================")
        println("Testing PointEvaluator")
        println("======================")
        test_pointevaluation2D()
        test_pointevaluation3D()
    end
end

function test_pointevaluation2D()
    ## check if a linear function is interpolated correctly
    ## and can be evaluated correctly within the cell
    xgrid = uniform_refine(grid_unitsquare(Triangle2D), 1)
    FES = FESpace{H1P1{1}}(xgrid)
    Iu = FEVector(FES)
    interpolate!(Iu[1], (result, qpinfo) -> (result[1] = qpinfo.x[1] + qpinfo.x[2]))
    PE = PointEvaluator([(1, Identity)])
    @show PE
    initialize!(PE, Iu)
    CF = CellFinder(xgrid)
    eval = zeros(Float64, 1)
    x = [0.15234, 0.2234] # point inside the cell
    xref = [0.0, 0.0]
    cell = gFindLocal!(xref, CF, x; icellstart = 1)
    evaluate_bary!(eval, PE, xref, cell)
    @test abs(eval[1] - sum(x)) < 1.0e-15

    evaluate!(eval, PE, x)
    return @test abs(eval[1] - sum(x)) < 1.0e-15
end

function test_pointevaluation3D()
    ## check if a linear function is interpolated correctly
    ## and can be evaluated correctly within the cell
    xgrid = uniform_refine(grid_unitcube(Tetrahedron3D), 1)
    FES = FESpace{H1P1{1}}(xgrid)
    Iu = FEVector(FES)
    interpolate!(Iu[1], (result, qpinfo) -> (result[1] = qpinfo.x[1] + qpinfo.x[2] + qpinfo.x[3]))
    PE = PointEvaluator([(1, Identity)])
    @show PE
    initialize!(PE, Iu)
    CF = CellFinder(xgrid)
    eval = zeros(Float64, 1)
    x = [0.15234, 0.2234, 0.5674] # point inside the cell
    xref = [0.0, 0.0, 0.0]
    cell = gFindLocal!(xref, CF, x; icellstart = 1)
    evaluate_bary!(eval, PE, xref, cell)
    @test abs(eval[1] - sum(x)) < 1.0e-15

    evaluate!(eval, PE, x)
    return @test abs(eval[1] - sum(x)) < 1.0e-15
end
