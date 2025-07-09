function run_fematrix_tests()

    @testset "FEMatrixVector" begin
        println("\n")
        println("===========================")
        println("Testing FEMatrix & FEVector")
        println("===========================")
        xgrid = simplexgrid(0:0.1:1, 0:0.1:1)
        FES1 = FESpace{H1Pk{1, 1, 1}}(xgrid)
        show(devnull, FES1)
        FES2 = FESpace{H1Pk{1, 1, 2}}(xgrid)
        show(devnull, FES2)
        A = FEMatrix(FES1, FES2)
        @test size(A.entries) == (FES1.ndofs, FES2.ndofs)
        @test size(A[1, 1]) == (FES1.ndofs, FES2.ndofs)
        show(devnull, A)

        B = FEMatrix([FES1, FES2])
        @test length(B) == 4
        @test size(B.entries) == (FES1.ndofs + FES2.ndofs, FES1.ndofs + FES2.ndofs)
        @test size(B[1, 2]) == (FES1.ndofs, FES2.ndofs)
        show(devnull, B)

        C = FEMatrix([FES2, FES2], [FES1, FES1])
        @test length(C) == 4
        @test size(C.entries) == (2 * FES2.ndofs, 2 * FES1.ndofs)
        @test size(C[1, 2]) == (FES2.ndofs, FES1.ndofs)
        C.entries.cscmatrix = sprand(2 * FES2.ndofs, 2 * FES1.ndofs, 0.5)
        show(devnull, C)

        b = FEVector([FES1, FES2])
        b.entries .= rand(FES1.ndofs + FES2.ndofs)

        addblock!(A[1, 1], C[1, 1]; transpose = true)
        @test norm(A[1, 1]) ≈ norm(C[1, 1])

        fill!(b[2], 0)
        addblock_matmul!(b[2], C[1, 1], b[1])
        @test all(view(b[2]) .== view(C.entries, 1:FES2.ndofs, 1:FES1.ndofs) * view(b[1]))
    end
    return println("")
end
