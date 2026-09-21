using Test
using ExtendableGrids
using ExtendableFEMBase
using ExtendableSparse
using ExplicitImports
using ExampleJuggler
using ForwardDiff
using LinearAlgebra
using SparseArrays
using Aqua

## shared test utilities: catalogs, test grids, exact functions, loop helpers
include("test_utils.jl")

@testset "Aqua.jl" begin
    Aqua.test_all(
        ExtendableFEMBase;
        ambiguities = false,
    )
    Aqua.test_ambiguities(ExtendableFEMBase)
end

@testset "ExplicitImports" begin
    @test ExplicitImports.check_no_implicit_imports(ExtendableFEMBase) === nothing
    @test ExplicitImports.check_no_stale_explicit_imports(ExtendableFEMBase) === nothing
end

if isdefined(Docs, :undocumented_names) # >=1.11
    @testset "UndocumentedNames" begin
        @test isempty(Docs.undocumented_names(ExtendableFEMBase))
    end
end


include("test_quadrature.jl")
include("test_interpolators.jl")
include("test_interpolation_matrix.jl")
include("test_operators.jl")
include("test_febasis.jl")
include("test_segmentintegrator.jl")
include("test_pointevaluator.jl")
include("test_fematrix_and_vector.jl")

function run_examples()
    ExampleJuggler.verbose!(true)

    example_dir = joinpath(@__DIR__, "..", "examples")

    modules = [
        "Example200_LowLevelPoisson.jl",
        "Example210_LowLevelNavierStokes.jl",
        "Example220_LowLevelHeatEquation.jl",
        "Example290_InterpolationBetweenMeshes.jl",
    ]

    return @testset "module examples" begin
        @testmodules(example_dir, modules)
    end
end

function run_all_tests()
    return begin
        run_fematrix_tests()
        run_examples()
        run_febasis_tests()
        run_operator_tests()
        run_quadrature_tests()
        run_interpolator_tests()
        run_grid_interpolation_matrix_tests()
        run_space_interpolation_matrix_tests()
        run_segmentintegrator_tests()
        run_pointevaluator_tests()
    end
end

run_all_tests()
