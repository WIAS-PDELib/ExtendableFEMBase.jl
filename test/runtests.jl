using Test
using ExtendableGrids
using ExtendableFEMBase
using ExplicitImports
using ExampleJuggler
using SparseArrays
using Aqua

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
    ]

    return @testset "module examples" begin
        @testmodules(example_dir, modules)
    end
end

function testgrid(::Type{Edge1D})
    return uniform_refine(simplexgrid([0.0, 1 // 4, 2 // 3, 1.0]), 1)
end
function testgrid(EG::Type{<:AbstractElementGeometry2D})
    return uniform_refine(grid_unitsquare(EG), 1)
end
function testgrid(EG::Type{<:AbstractElementGeometry3D})
    return uniform_refine(grid_unitcube(EG), 1)
end
function testgrid(::Type{Triangle2D}, ::Type{Parallelogram2D})
    return uniform_refine(grid_unitsquare_mixedgeometries(), 1)
end

tolerance = 6.0e-12

function exact_function(::Val{1}, polyorder)
    function polynomial(result, qpinfo)
        x = qpinfo.x
        return result[1] = x[1]^polyorder + 1
    end
    function gradient(result, qpinfo)
        x = qpinfo.x
        return result[1] = polyorder * x[1]^(polyorder - 1)
    end
    function hessian(result, qpinfo)
        x = qpinfo.x
        return result[1] = polyorder * (polyorder - 1) * x[1]^(polyorder - 2)
    end
    exact_integral = 1 // (polyorder + 1) + 1
    return polynomial, exact_integral, gradient, hessian
end

function exact_function(::Val{2}, polyorder)
    function polynomial(result, qpinfo)
        x = qpinfo.x
        result[1] = x[1]^polyorder + 2 * x[2]^polyorder + 1
        return result[2] = 3 * x[1]^polyorder - x[2]^polyorder - 1
    end
    function gradient(result, qpinfo)
        x = qpinfo.x
        result[1] = polyorder * x[1]^(polyorder - 1)
        result[2] = 2 * polyorder * x[2]^(polyorder - 1)
        result[3] = 3 * polyorder * x[1]^(polyorder - 1)
        return result[4] = -polyorder * x[2]^(polyorder - 1)
    end
    function hessian(result, qpinfo)
        x = qpinfo.x
        result[1] = polyorder * (polyorder - 1) * x[1]^(polyorder - 2)
        result[2] = 0
        result[3] = 0
        result[4] = 2 * polyorder * (polyorder - 1) * x[2]^(polyorder - 2)
        result[5] = 3 * polyorder * (polyorder - 1) * x[1]^(polyorder - 2)
        result[6] = 0
        result[7] = 0
        return result[8] = -polyorder * (polyorder - 1) * x[2]^(polyorder - 2)
    end
    exact_integral = [3 // (polyorder + 1) + 1, 2 // (polyorder + 1) - 1]
    return polynomial, exact_integral, gradient, hessian
end

function exact_function(::Val{3}, polyorder)
    function polynomial(result, qpinfo)
        x = qpinfo.x
        result[1] = 2 * x[3]^polyorder - x[2]^polyorder - 1
        result[2] = x[1]^polyorder + 2 * x[2]^polyorder + 1
        return result[3] = 3 * x[1]^polyorder - x[2]^polyorder - 1
    end
    function gradient(result, qpinfo)
        x = qpinfo.x
        result[1] = 0
        result[2] = -polyorder * x[2]^(polyorder - 1)
        result[3] = 2 * polyorder * x[3]^(polyorder - 1)
        result[4] = polyorder * x[1]^(polyorder - 1)
        result[5] = 2 * polyorder * x[2]^(polyorder - 1)
        result[6] = 0
        result[7] = 3 * polyorder * x[2]^(polyorder - 1)
        result[8] = -polyorder * x[2]^(polyorder - 1)
        return result[9] = 0
    end
    function hessian(result, qpinfo)
        x = qpinfo.x
        fill!(result, 0)
        result[5] = -polyorder * (polyorder - 1) * x[2]^(polyorder - 2)
        result[9] = 2 * polyorder * (polyorder - 1) * x[3]^(polyorder - 2)
        result[10] = polyorder * (polyorder - 1) * x[1]^(polyorder - 2)
        result[14] = 2 * polyorder * (polyorder - 1) * x[2]^(polyorder - 2)
        result[19] = 3 * polyorder * (polyorder - 1) * x[1]^(polyorder - 2)
        return result[23] = -polyorder * (polyorder - 1) * x[2]^(polyorder - 2)
    end
    exact_integral = [1 // (polyorder + 1) - 1, 3 // (polyorder + 1) + 1, 2 // (polyorder + 1) - 1]
    return polynomial, exact_integral, gradient, hessian
end


function run_all_tests()
    return begin
        run_fematrix_tests()
        run_examples()
        run_febasis_tests()
        run_operator_tests()
        run_quadrature_tests()
        run_interpolator_tests()
        run_segmentintegrator_tests()
        run_pointevaluator_tests()
    end
end

run_all_tests()
