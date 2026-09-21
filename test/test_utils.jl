############################
# SHARED TEST UTILITIES    #
# catalogs, test grids and #
# helpers used by several  #
# test files               #
############################

## tolerance for polynomial-precision tests
tolerance = 6.0e-12

## list of FETypes that should be tested
## (shared between interpolation and interpolation matrix tests)
const TestCatalog1D = [
    L2P0{1} => 0,
    H1P1{1} => 1,
    H1P2{1, 1} => 2,
    H1P3{1, 1} => 3,
    H1Pk{1, 1, 3} => 3,
    H1Pk{1, 1, 4} => 4,
    H1Pk{1, 1, 5} => 5,
]

const TestCatalog2D = [
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

const TestCatalog3D = [
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

## test grids: nonuniform 1D simplex, unit square, unit cube
function testgrid(::Type{Edge1D})
    return uniform_refine(simplexgrid([0.0, 1 // 4, 2 // 3, 1.0]), 1)
end
function testgrid(EG::Type{<:AbstractElementGeometry2D})
    return uniform_refine(grid_unitsquare(EG), 1)
end
function testgrid(EG::Type{<:AbstractElementGeometry3D})
    return uniform_refine(grid_unitcube(EG), 1)
end

## grids used by the catalog tests: one refinement of the reference domain
function testgrid_refdomain(EG::Type{<:AbstractElementGeometry2D})
    return uniform_refine(reference_domain(EG), 1)
end
function testgrid_refdomain(EG::Type{<:AbstractElementGeometry3D})
    return uniform_refine(reference_domain(EG), 1)
end

## generic loop over element geometries x catalog entries x broken flag
## the callback f (do-block convention, like map(f, xs)) is invoked as
## f(EG, xgrid, FEType, order, broken) for every catalog entry defined on
## EG (skipped with a warning otherwise)
function for_each_catalog(f, EGs, catalog; grid = testgrid, broken_values = (false, true))
    for EG in EGs
        println("EG = $EG")
        xgrid = grid(EG)
        for (FEType, order) in catalog
            for broken in broken_values
                if ExtendableFEMBase.isdefined(FEType, EG, broken)
                    f(EG, xgrid, FEType, order, broken)
                else
                    @warn "$(FEType) (broken = $(broken)) not defined on $(EG) (skipping test case)"
                end
            end
        end
    end
    return nothing
end

## exact polynomial test functions used by quadrature and
## interpolation error tests; return (function, exact integral, gradient, hessian)
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

## closure around a PointEvaluator that caches one evaluator per number type
## (used by the automatic differentiation tests in test_pointevaluator.jl)
function point_evaluator_closure(coeffs = [1.0, 0.0, 0.0])
    input_types = Dict{DataType, Any}()
    result_types = Dict{DataType, Any}()

    grid_types = Dict{DataType, Any}()
    FES_types = Dict{DataType, Any}()
    RVec_types = Dict{DataType, Any}()
    PE_types = Dict{DataType, Any}()

    function closure(x::Vector{T}) where {T}
        if !haskey(PE_types, T)
            input_types[T] = zeros(T, 2)
            result_types[T] = zeros(T, 1)

            grid_types[T] = reference_domain(Triangle2D, T)
            FES_types[T] = FESpace{H1P1{1}}(grid_types[T])
            RVec_types[T] = FEVector(FES_types[T])
            RVec_types[T].entries .= coeffs
            PE_types[T] = PointEvaluator([(1, Identity)], [RVec_types[T][1]]; Tv = T, TCoeff = T)
        end

        input_types[T][1] = log(x[1])
        input_types[T][2] = log(x[2])

        ExtendableFEMBase.evaluate!(result_types[T], PE_types[T], input_types[T])

        return result_types[T]
    end

    return closure
end
