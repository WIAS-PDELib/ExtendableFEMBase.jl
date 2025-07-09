function run_quadrature_tests()
    ############################
    # TESTSET QUADRATURE RULES #
    ############################
    maxorder1D = 20
    EG2D = [Triangle2D, Parallelogram2D]
    maxorder2D = [20, 20]
    EG3D = [Parallelepiped3D, Tetrahedron3D]
    maxorder3D = [12, 8]

    return @testset "QuadratureRules" begin
        println("\n")
        println("=============================")
        println("Testing QuadratureRules in 1D")
        println("=============================")
        xgrid = testgrid(Edge1D)
        for order in 1:maxorder1D
            integrand, exactvalue = exact_function(Val(1), order)
            qf = QuadratureRule{Float64, Edge1D}(order)
            show(devnull, qf)
            quadvalue = integrate(xgrid, ON_CELLS, integrand, length(exactvalue); force_quadrature_rule = qf)
            println("EG = Edge1D | order = $order ($(qf.name), $(length(qf.w)) points) | error = $(quadvalue - exactvalue)")
            @test isapprox(quadvalue, exactvalue)
        end
        println("\n")
        println("=============================")
        println("Testing QuadratureRules in 2D")
        println("=============================")
        for (j, EG) in enumerate(EG2D)
            xgrid = testgrid(EG)
            for order in 1:maxorder2D[j]
                integrand, exactvalue = exact_function(Val(2), order)
                qf = QuadratureRule{Float64, EG}(order)
                show(devnull, qf)
                quadvalue = integrate(xgrid, ON_CELLS, integrand, length(exactvalue); force_quadrature_rule = qf)
                println("EG = $EG | order = $order ($(qf.name), $(length(qf.w)) points) | error = $(quadvalue - exactvalue)")
                @test isapprox(quadvalue, exactvalue)
            end
        end
        println("\n")
        println("=============================")
        println("Testing QuadratureRules in 3D")
        println("=============================")
        for (j, EG) in enumerate(EG3D)
            xgrid = testgrid(EG)
            for order in 1:maxorder3D[j]
                integrand, exactvalue = exact_function(Val(3), order)
                qf = QuadratureRule{Float64, EG}(order)
                show(devnull, qf)
                quadvalue = integrate(xgrid, ON_CELLS, integrand, length(exactvalue); force_quadrature_rule = qf)
                println("EG = $EG | order = $order ($(qf.name), $(length(qf.w)) points) | error = $(quadvalue - exactvalue)")
                @test isapprox(quadvalue, exactvalue)
            end
        end
        println("")
    end

end
