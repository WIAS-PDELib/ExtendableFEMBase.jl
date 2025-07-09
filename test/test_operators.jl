function run_operator_tests()
    return @testset "Operators" begin
        println("\n")
        println("============================")
        println("Testing Operator Evaluations")
        println("============================")
        error = test_derivatives2D()
        @test error < 1.0e-14
        error = test_derivatives3D()
        @test error < 1.0e-14
    end
end

function test_derivatives2D()
    ## define test function and expected operator evals
    function testf(result, qpinfo)
        x = qpinfo.x
        result[1] = x[1]^2
        return result[2] = 3 * x[2]^2 + x[1] * x[2]
    end

    ## expected values of operators in cell midpoint
    expected_L = [2, 6] # expected Laplacian
    expected_H = [2, 0, 0, 0, 0, 1, 1, 6] # expected Hessian
    expected_symH = [2, 0, 0, 0, 6, 1] # expected symmetric Hessian
    expected_symH2 = [2, 0, 0, 0, 6, sqrt(2)] # expected symmetric Hessian
    expected_curl2 = [1 / 3]

    ## define grid = a single non-refenrece triangle
    xgrid = grid_triangle([-1.0 0.0; 1.0 0.0; 0.0 1.0]')

    ## define P2-Courant finite element space
    FEType = H1P2{2, 2}
    FES = FESpace{FEType}(xgrid)
    @show FES

    ## get midpoint quadrature rule for constants
    qf = QuadratureRule{Float64, Triangle2D}(0)

    ## define FE basis Evaluator for Hessian
    FEBE_curl2 = FEEvaluator(FES, Curl2D, qf)
    FEBE_L = FEEvaluator(FES, Laplacian, qf)
    FEBE_H = FEEvaluator(FES, Hessian, qf)
    FEBE_symH = FEEvaluator(FES, SymmetricHessian{1}, qf)
    FEBE_symH2 = FEEvaluator(FES, SymmetricHessian{sqrt(2)}, qf)

    ## update on cell 1
    update_basis!(FEBE_L, 1)
    update_basis!(FEBE_H, 1)
    update_basis!(FEBE_symH, 1)
    update_basis!(FEBE_symH2, 1)
    update_basis!(FEBE_curl2, 1)

    ## interpolate quadratic testfunction
    Iu = FEVector(FES)
    interpolate!(Iu[1], testf)

    ## check if operator evals have the correct length
    @assert size(FEBE_L.cvals, 1) == length(expected_L)
    @assert size(FEBE_H.cvals, 1) == length(expected_H)
    @assert size(FEBE_symH.cvals, 1) == length(expected_symH)
    @assert size(FEBE_symH2.cvals, 1) == length(expected_symH2)
    @assert size(FEBE_curl2.cvals, 1) == length(expected_curl2)

    ## eval 2nd order derivatives at only quadrature point 1
    ## since function is quadratic this should be constant
    H = zeros(Float64, 8)
    symH = zeros(Float64, 6)
    symH2 = zeros(Float64, 6)
    L = zeros(Float64, 2)
    curl2 = zeros(Float64, 1)
    eval_febe!(L, FEBE_L, Iu.entries[FES[CellDofs][:, 1]], 1)
    eval_febe!(H, FEBE_H, Iu.entries[FES[CellDofs][:, 1]], 1)
    eval_febe!(symH, FEBE_symH, Iu.entries[FES[CellDofs][:, 1]], 1)
    eval_febe!(symH2, FEBE_symH2, Iu.entries[FES[CellDofs][:, 1]], 1)
    eval_febe!(curl2, FEBE_curl2, Iu.entries[FES[CellDofs][:, 1]], 1)

    ## compute errors to expected values
    error_L = sqrt(sum((L - expected_L) .^ 2))
    error_H = sqrt(sum((H - expected_H) .^ 2))
    error_symH = sqrt(sum((symH - expected_symH) .^ 2))
    error_symH2 = sqrt(sum((symH2 - expected_symH2) .^ 2))
    error_curl2 = sqrt(sum((curl2 - expected_curl2) .^ 2))
    println("EG = Triangle2D | operator = Curl2 | error = $error_curl2")
    println("EG = Triangle2D | operator = Laplacian | error = $error_L")
    println("EG = Triangle2D | operator = Hessian | error = $error_H")
    println("EG = Triangle2D | operator = SymmetricHessian{1} | error = $error_symH")
    println("EG = Triangle2D | operator = SymmetricHessian{√2} | error = $error_symH2")

    return maximum([error_curl2, error_L, error_H, error_symH, error_symH2])
end


function test_derivatives3D()
    ## define test function and expected operator evals
    function testf(result, qpinfo)
        x = qpinfo.x
        result[1] = x[1]^2 + x[3] * x[2]
        result[2] = 3 * x[3]^2 + x[1] * x[2]
        return result[3] = x[1] * x[2]
    end

    ## expected values of operators in cell midpoint
    expected_L = [2, 6, 0] # expected Laplacian
    expected_H = [2, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 6, 0, 1, 0, 1, 0, 0, 0, 0, 0] # expected Hessian
    expected_symH = [2, 0, 0, 1, 0, 0, 0, 0, 6, 0, 0, 1, 0, 0, 0, 0, 0, 1] # expected symmetric Hessian
    expected_symH2 = [2, 0, 0, sqrt(2), 0, 0, 0, 0, 6, 0, 0, sqrt(2), 0, 0, 0, 0, 0, sqrt(2)] # expected symmetric Hessian
    expected_curl3 = [1 / 2 - 6 / 4, 0, 1 / 4 - 1 / 4]

    ## define grid = a single non-refenrece triangle
    xgrid = reference_domain(Tetrahedron3D)
    xgrid[Coordinates][:, 2] = [2, 0, 0]

    ## define P2-Courant finite element space
    FEType = H1P2{3, 3}
    FES = FESpace{FEType}(xgrid)
    @show FES

    ## get midpoint quadrature rule for constants
    qf = QuadratureRule{Float64, Tetrahedron3D}(0)

    ## define FE basis Evaluator for Hessian
    FEBE_curl3 = FEEvaluator(FES, Curl3D, qf)
    FEBE_L = FEEvaluator(FES, Laplacian, qf)
    FEBE_H = FEEvaluator(FES, Hessian, qf)
    FEBE_symH = FEEvaluator(FES, SymmetricHessian{1}, qf)
    FEBE_symH2 = FEEvaluator(FES, SymmetricHessian{sqrt(2)}, qf)

    ## update on cell 1
    update_basis!(FEBE_L, 1)
    update_basis!(FEBE_H, 1)
    update_basis!(FEBE_symH, 1)
    update_basis!(FEBE_symH2, 1)
    update_basis!(FEBE_curl3, 1)

    ## interpolate quadratic testfunction
    Iu = FEVector(FES)
    interpolate!(Iu[1], testf)

    ## check if operator evals have the correct length
    @assert size(FEBE_L.cvals, 1) == length(expected_L)
    @assert size(FEBE_H.cvals, 1) == length(expected_H)
    @assert size(FEBE_symH.cvals, 1) == length(expected_symH)
    @assert size(FEBE_symH2.cvals, 1) == length(expected_symH2)
    @assert size(FEBE_curl3.cvals, 1) == length(expected_curl3)

    ## eval 2nd order derivatives at only quadrature point 1
    ## since function is quadratic this should be constant
    H = zeros(Float64, 27)
    symH = zeros(Float64, 18)
    symH2 = zeros(Float64, 18)
    L = zeros(Float64, 3)
    curl3 = zeros(Float64, 3)
    eval_febe!(L, FEBE_L, Iu.entries[FES[CellDofs][:, 1]], 1)
    eval_febe!(H, FEBE_H, Iu.entries[FES[CellDofs][:, 1]], 1)
    eval_febe!(symH, FEBE_symH, Iu.entries[FES[CellDofs][:, 1]], 1)
    eval_febe!(symH2, FEBE_symH2, Iu.entries[FES[CellDofs][:, 1]], 1)
    eval_febe!(curl3, FEBE_curl3, Iu.entries[FES[CellDofs][:, 1]], 1)

    ## compute errors to expected values
    error_L = sqrt(sum((L - expected_L) .^ 2))
    error_H = sqrt(sum((H - expected_H) .^ 2))
    error_symH = sqrt(sum((symH - expected_symH) .^ 2))
    error_symH2 = sqrt(sum((symH2 - expected_symH2) .^ 2))
    error_curl3 = sqrt(sum((curl3 - expected_curl3) .^ 2))
    println("EG = Tetrahedron3D | operator = Curl3 | error = $error_curl3")
    println("EG = Tetrahedron3D | operator = Laplacian | error = $error_L")
    println("EG = Tetrahedron3D | operator = Hessian | error = $error_H")
    println("EG = Tetrahedron3D | operator = SymmetricHessian{1} | error = $error_symH")
    println("EG = Tetrahedron3D | operator = SymmetricHessian{√2} | error = $error_symH2")

    return maximum([error_curl3, error_L, error_H, error_symH, error_symH2])
end
