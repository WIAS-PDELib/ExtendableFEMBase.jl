function run_operator_tests()
    return @testset "Operators" begin
        println("\n")
        println("============================")
        println("Testing Operator Evaluations")
        println("============================")
        error = test_derivatives2D(H1P2{2, 2}, 2)
        @test error < 1.0e-14
        error = test_derivatives2D(HDIVBDM2{2}, 2)
        @test error < 1.0e-14
        error = test_derivatives2D(HCURLN1{2}, 1)
        @test error < 1.0e-14
        error = test_derivatives3D(H1P2{3, 3}, 2)
        @test error < 1.0e-14
        error = test_derivatives3D(HDIVRT1{3}, 1)
        @test error < 1.0e-14
        error = test_derivatives3D(HCURLN0{3}, 1)
        @test error < 1.0e-14
        test_reconstructions()
    end
end

function test_reconstructions()
    ## divergence-free axisymmetric velocity field u(r,z) = (r,-2z) in cylindrical coordinates
    function u!(result, qpinfo)
        x = qpinfo.x
        result[1] = x[1]
        return result[2] = -2 * x[2]
    end

    ## vector field times r, should have zero divergence div(ru) = (d/dr, d/dz) \cdot (ru) = 0
    function ru!(result, qpinfo)
        x = qpinfo.x
        result[1] = x[1]^2
        return result[2] = -2 * x[1] * x[2]
    end

    xgrid = testgrid(Triangle2D)

    ## interpolate u into H1BR{2} (inf-sup stable Stokes element)
    FES = FESpace{H1BR{2}}(xgrid)
    uh = FEVector(FES)
    interpolate!(uh[1], u!; bonus_quadorder = 2)

    for FETypeR in [HDIVRT0{2}, HDIVBDM1{2}]
        ## interpolate ru into HDIVRTO{2} or HDIVBDM1{2}
        FES2 = FESpace{FETypeR}(xgrid)
        Πur = FEVector(FES2)
        interpolate!(Πur[1], ru!; bonus_quadorder = 2)

        ## test if interpolate of ru is divergence-free by interpolating into P0 function and checking its coefficients
        FES3 = FESpace{L2P0{1}}(xgrid)
        divΠur = FEVector(FES3)
        lazy_interpolate!(divΠur[1], Πur, [(1, Divergence)])
        @test sqrt(sum((divΠur.entries .^ 2))) < 1.0e-14

        ## test if r-weighted reconstruction of uh is divergence-free by interpolating into P0 function and checking its coefficients
        weight = (x) -> (x[1])
        lazy_interpolate!(divΠur[1], uh, [(1, WeightedReconstruct{FETypeR, Divergence, typeof(weight)})])
        @test sqrt(sum((divΠur.entries .^ 2))) < 1.0e-14

        ## test if weighted reconstruction of uh and interpolation of ru are identical
        FES4 = FESpace{L2P1{1}}(xgrid)
        diff = FEVector(FES4)
        lazy_interpolate!(diff[1], [uh[1], Πur[1]], [(1, WeightedReconstruct{FETypeR, Identity, typeof(weight)}), (2, Identity)]; postprocess = (result, args, qpinfo) -> (result[1] = (args[1] - args[3])^2 + (args[2] - args[4])^2))
        @test sqrt(sum((diff.entries))) < 1.0e-14
    end
    return
end

function test_derivatives2D(fetype, order)
    ## define test function and expected operator evals
    function testf(result, qpinfo)
        x = qpinfo.x
        if order == 2
            result[1] = x[1]^2
            result[2] = 3 * x[2]^2 + x[1] * x[2]
        elseif order == 1
            result[1] = x[1] + x[2] + 1
            result[2] = 3 * x[2] - x[1]
        end
        return nothing
    end

    ## define grid = a single non-reference triangle
    xgrid = grid_triangle([-1.0 0.0; 1.0 0.0; 0.0 1.0]') # midpoint = [0.0, 1 / 3]

    ## expected values of operators in cell midpoint
    if order == 2
        expected_id = [0, 1 / 3]
        expected_L = [2, 6] # expected Laplacian
        expected_H = [2, 0, 0, 0, 0, 1, 1, 6] # expected Hessian
        expected_symH = [2, 0, 0, 0, 6, 1] # expected symmetric Hessian
        expected_symH2 = [2, 0, 0, 0, 6, sqrt(2)] # expected symmetric Hessian
        expected_curl2 = [1 / 3]
        expected_grad = [0, 0, 1 / 3, 2]
        expected_div = [2]
    elseif order == 1
        expected_id = [4 / 3, 1]
        expected_curl2 = [-2]
    end

    ## define P2-Courant finite element space
    FES = FESpace{fetype}(xgrid)
    show(devnull, FES)

    ## interpolate quadratic testfunction
    Iu = FEVector(FES)
    interpolate!(Iu[1], testf)

    ## get midpoint quadrature rule for constants
    qf = QuadratureRule{Float64, Triangle2D}(0)

    FEBE_id = FEEvaluator(FES, Identity, qf)
    FEBE_curl2 = FEEvaluator(FES, Curl2D, qf)
    update_basis!(FEBE_id, 1)
    update_basis!(FEBE_curl2, 1)
    ## check if operator evals have the correct length
    @assert size(FEBE_id.cvals, 1) == length(expected_id)
    @assert size(FEBE_curl2.cvals, 1) == length(expected_curl2)
    # evaluate at quadrature points = cell midpoint
    id = zeros(Float64, 2)
    curl2 = zeros(Float64, 1)
    eval_febe!(id, FEBE_id, Iu.entries[FES[CellDofs][:, 1]], 1)
    eval_febe!(curl2, FEBE_curl2, Iu.entries[FES[CellDofs][:, 1]], 1)
    ## compute errors to expected values
    error_id = sqrt(sum((id - expected_id) .^ 2))
    error_curl2 = sqrt(sum((curl2 - expected_curl2) .^ 2))
    println("EG = Triangle2D | $fetype | operator = Identity | error = $error_id")
    println("EG = Triangle2D | $fetype | operator = Curl2 | error = $error_curl2")
    if fetype <: AbstractHcurlFiniteElement
        return maximum([error_id, error_curl2])
    else
        grad = zeros(Float64, 4)
        div = zeros(Float64, 1)
        FEBE_grad = FEEvaluator(FES, Gradient, qf)
        FEBE_div = FEEvaluator(FES, Divergence, qf)
        update_basis!(FEBE_grad, 1)
        update_basis!(FEBE_div, 1)
        @assert size(FEBE_grad.cvals, 1) == length(expected_grad)
        @assert size(FEBE_div.cvals, 1) == length(expected_div)
        eval_febe!(grad, FEBE_grad, Iu.entries[FES[CellDofs][:, 1]], 1)
        eval_febe!(div, FEBE_div, Iu.entries[FES[CellDofs][:, 1]], 1)
        error_grad = sqrt(sum((grad - expected_grad) .^ 2))
        error_div = sqrt(sum((div - expected_div) .^ 2))
        println("EG = Triangle2D | $fetype | operator = Gradient | error = $error_grad")
        println("EG = Triangle2D | $fetype | operator = Divergence | error = $error_div")

        if fetype <: AbstractH1FiniteElement
            # do the same for second order derivatives
            FEBE_L = FEEvaluator(FES, Laplacian, qf)
            FEBE_H = FEEvaluator(FES, Hessian, qf)
            FEBE_symH = FEEvaluator(FES, SymmetricHessian{1}, qf)
            FEBE_symH2 = FEEvaluator(FES, SymmetricHessian{sqrt(2)}, qf)
            update_basis!(FEBE_L, 1)
            update_basis!(FEBE_H, 1)
            update_basis!(FEBE_symH, 1)
            update_basis!(FEBE_symH2, 1)
            @assert size(FEBE_L.cvals, 1) == length(expected_L)
            @assert size(FEBE_H.cvals, 1) == length(expected_H)
            @assert size(FEBE_symH.cvals, 1) == length(expected_symH)
            @assert size(FEBE_symH2.cvals, 1) == length(expected_symH2)
            H = zeros(Float64, 8)
            symH = zeros(Float64, 6)
            symH2 = zeros(Float64, 6)
            L = zeros(Float64, 2)
            eval_febe!(L, FEBE_L, Iu.entries[FES[CellDofs][:, 1]], 1)
            eval_febe!(H, FEBE_H, Iu.entries[FES[CellDofs][:, 1]], 1)
            eval_febe!(symH, FEBE_symH, Iu.entries[FES[CellDofs][:, 1]], 1)
            eval_febe!(symH2, FEBE_symH2, Iu.entries[FES[CellDofs][:, 1]], 1)
            error_L = sqrt(sum((L - expected_L) .^ 2))
            error_H = sqrt(sum((H - expected_H) .^ 2))
            error_symH = sqrt(sum((symH - expected_symH) .^ 2))
            error_symH2 = sqrt(sum((symH2 - expected_symH2) .^ 2))
            println("EG = Triangle2D | $fetype | operator = Laplacian | error = $error_L")
            println("EG = Triangle2D | $fetype | operator = Hessian | error = $error_H")
            println("EG = Triangle2D | $fetype | operator = SymmetricHessian{1} | error = $error_symH")
            println("EG = Triangle2D | $fetype | operator = SymmetricHessian{√2} | error = $error_symH2")
            return maximum([error_id, error_curl2, error_L, error_H, error_symH, error_symH2, error_grad, error_div])
        else
            return maximum([error_id, error_curl2, error_grad])
        end
    end


end

function test_derivatives3D(fetype, order)
    ## define test function and expected operator evals
    function testf(result, qpinfo)
        x = qpinfo.x
        if order == 2
            result[1] = x[1]^2 + x[3] * x[2]
            result[2] = 3 * x[3]^2 + x[1] * x[2]
            result[3] = x[1] * x[2]
        elseif order == 1
            result[1] = x[2]
            result[2] = 3 * x[3]
            result[3] = x[2]
        end
        return nothing
    end

    ## define grid = a single non-refenrece triangle
    xgrid = reference_domain(Tetrahedron3D)
    xgrid[Coordinates][:, 2] = [2, 0, 0] # midpoint = [0.5, 0.25, 0.25]

    ## define P2-Courant finite element space
    FES = FESpace{fetype}(xgrid)
    show(devnull, FES)

    ## interpolate quadratic testfunction
    Iu = FEVector(FES)
    interpolate!(Iu[1], testf)

    ## expected values of operators in cell midpoint
    if order == 2
        expected_L = [2, 6, 0] # expected Laplacian
        expected_H = [2, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 6, 0, 1, 0, 1, 0, 0, 0, 0, 0] # expected Hessian
        expected_symH = [2, 0, 0, 1, 0, 0, 0, 0, 6, 0, 0, 1, 0, 0, 0, 0, 0, 1] # expected symmetric Hessian
        expected_symH2 = [2, 0, 0, sqrt(2), 0, 0, 0, 0, 6, 0, 0, sqrt(2), 0, 0, 0, 0, 0, sqrt(2)] # expected symmetric Hessian
        expected_curl3 = [1 / 2 - 6 / 4, 0, 1 / 4 - 1 / 4]
        expected_grad = [1, 0.25, 0.25, 0.25, 0.5, 1.5, 0.25, 0.5, 0] # expected Gradient
        expected_div = [1.5]
    elseif order == 1
        expected_curl3 = [1 - 3, 0, -1]
        expected_grad = [0, 1, 0, 0, 0, 3, 0, 1, 0] # expected Gradient
        expected_div = [0]
    end

    ## get midpoint quadrature rule for constants
    qf = QuadratureRule{Float64, Tetrahedron3D}(0)

    FEBE_curl3 = FEEvaluator(FES, Curl3D, qf)
    update_basis!(FEBE_curl3, 1)
    @assert size(FEBE_curl3.cvals, 1) == length(expected_curl3)

    FEBE_grad = FEEvaluator(FES, Gradient, qf)
    FEBE_div = FEEvaluator(FES, Divergence, qf)
    curl3 = zeros(Float64, 3)
    eval_febe!(curl3, FEBE_curl3, Iu.entries[FES[CellDofs][:, 1]], 1)
    error_curl3 = sqrt(sum((curl3 - expected_curl3) .^ 2))
    println("EG = Tetrahedron3D | $fetype | operator = Curl3 | error = $error_curl3")

    if fetype <: AbstractHcurlFiniteElement
        return maximum([error_curl3])
    else

        update_basis!(FEBE_grad, 1)
        update_basis!(FEBE_div, 1)

        @assert size(FEBE_grad.cvals, 1) == length(expected_grad)
        @assert size(FEBE_div.cvals, 1) == length(expected_div)
        grad = zeros(Float64, 9)
        div = zeros(Float64, 1)
        eval_febe!(grad, FEBE_grad, Iu.entries[FES[CellDofs][:, 1]], 1)
        eval_febe!(div, FEBE_div, Iu.entries[FES[CellDofs][:, 1]], 1)
        error_grad = sqrt(sum((grad - expected_grad) .^ 2))
        error_div = sqrt(sum((div - expected_div) .^ 2))
        println("EG = Tetrahedron3D | $fetype | operator = Gradient | error = $error_grad")
        println("EG = Tetrahedron3D | $fetype | operator = Divergence | error = $error_div")

        if fetype <: AbstractH1FiniteElement
            FEBE_L = FEEvaluator(FES, Laplacian, qf)
            FEBE_H = FEEvaluator(FES, Hessian, qf)
            FEBE_symH = FEEvaluator(FES, SymmetricHessian{1}, qf)
            FEBE_symH2 = FEEvaluator(FES, SymmetricHessian{sqrt(2)}, qf)

            ## update on cell 1
            update_basis!(FEBE_L, 1)
            update_basis!(FEBE_H, 1)
            update_basis!(FEBE_symH, 1)
            update_basis!(FEBE_symH2, 1)

            ## check if operator evals have the correct length
            @assert size(FEBE_L.cvals, 1) == length(expected_L)
            @assert size(FEBE_H.cvals, 1) == length(expected_H)
            @assert size(FEBE_symH.cvals, 1) == length(expected_symH)
            @assert size(FEBE_symH2.cvals, 1) == length(expected_symH2)

            ## eval 2nd order derivatives at only quadrature point 1
            ## since function is quadratic this should be constant
            H = zeros(Float64, 27)
            symH = zeros(Float64, 18)
            symH2 = zeros(Float64, 18)
            L = zeros(Float64, 3)
            eval_febe!(L, FEBE_L, Iu.entries[FES[CellDofs][:, 1]], 1)
            eval_febe!(H, FEBE_H, Iu.entries[FES[CellDofs][:, 1]], 1)
            eval_febe!(symH, FEBE_symH, Iu.entries[FES[CellDofs][:, 1]], 1)
            eval_febe!(symH2, FEBE_symH2, Iu.entries[FES[CellDofs][:, 1]], 1)

            ## compute errors to expected values
            error_L = sqrt(sum((L - expected_L) .^ 2))
            error_H = sqrt(sum((H - expected_H) .^ 2))
            error_symH = sqrt(sum((symH - expected_symH) .^ 2))
            error_symH2 = sqrt(sum((symH2 - expected_symH2) .^ 2))
            println("EG = Tetrahedron3D | $fetype | operator = Laplacian | error = $error_L")
            println("EG = Tetrahedron3D | $fetype | operator = Hessian | error = $error_H")
            println("EG = Tetrahedron3D | $fetype | operator = SymmetricHessian{1} | error = $error_symH")
            println("EG = Tetrahedron3D | $fetype | operator = SymmetricHessian{√2} | error = $error_symH2")
            return maximum([error_curl3, error_L, error_H, error_symH, error_symH2, error_grad, error_div])
        else
            return maximum([error_curl3, error_grad, error_div])
        end
    end

end
