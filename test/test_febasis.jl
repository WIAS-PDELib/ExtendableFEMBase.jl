function run_febasis_tests()

    @testset "FEBasis" begin
        println("\n")
        println("============================")
        println("Testing Interpolations in 1D")
        println("============================")
        for order in 1:6
            test_vertex_values(Edge1D, H1Pk{1, 1, order}, order)
        end
        for order in 1:5
            test_vertex_values(Triangle2D, H1Pk{1, 2, order}, order)
        end
        test_vertex_values(Parallelogram2D, H1Q1{1}, 1)
        test_vertex_values(Parallelogram2D, H1Q2{1, 2}, 2)
        test_vertex_values(Tetrahedron3D, H1P1{1}, 1)
        test_vertex_values(Tetrahedron3D, H1P2{1, 3}, 2)
        test_vertex_values(Tetrahedron3D, H1P3{1, 3}, 3)
        test_vertex_values(Parallelepiped3D, H1Q1{1}, 1)
    end
    return println("")
end

function test_vertex_values(EG, FEType, order)
    qf = VertexRule(EG, order; T = Rational{Int})
    show(devnull, qf)
    basis = get_basis(ON_CELLS, FEType, EG)
    ndofs = get_ndofs(ON_CELLS, FEType, EG)
    basis_vals = zeros(Rational{Int}, ndofs, 1)
    @info "Testing vertex values of $FEType on $EG"
    @assert ndofs == length(qf.xref)
    expected = zeros(Int, ndofs)
    error = zeros(Float64, ndofs)
    for (j, xref) in enumerate(qf.xref)
        basis(basis_vals, xref)
        expected[j] = 1
        error[j] = norm(basis_vals[:, 1] - expected, Inf)
        if error[j] > 0
            @info "dof $j does not have the expected values"
            @show basis_vals, expected
        end
        expected[j] = 0
    end
    J = findall(>(1.0e-15), error)
    return @test length(J) == 0
end
