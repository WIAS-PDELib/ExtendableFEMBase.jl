"""
$(TYPEDEF)

A mutable struct that encapsulates information about the current quadrature point during finite element assembly or evaluation.

# Fields
- `item::Ti`: Index of the current item (with respect to the current assembly type).
- `cell::Ti`: Index of the current cell, if applicable.
- `region::Ti`: Index of the current region (with respect to the current assembly type).
- `volume::TvG`: Volume associated with the current item.
- `normal::Vector{TvG}`: Outward normal vector at the quadrature point, if applicable (e.g., for face assembly/integration).
- `time::Ttime`: Current time value, useful for time-dependent problems.
- `x::Vector{Tx}`: World (physical) coordinates of the quadrature point.
- `xref::Vector{Txref}`: Reference (local) coordinates of the quadrature point within the reference element.
- `grid::ExtendableGrid{TvG, TiG}`: Reference to the underlying grid or mesh structure.
- `params::PT`: Additional user- or problem-specific parameters.

# Description
`QPInfos` is used to pass contextual information about the current quadrature point to finite element kernels, integrators, and evaluators. It provides access to geometric, topological, and problem-specific data required for local assembly, evaluation, or postprocessing routines.

# Notes
- Not all fields are always meaningful for every assembly type (e.g., `normal` may be unused for volume integrals).

"""
mutable struct QPInfos{Ti, Tv, Ttime, Tx, Txref, TvG, TiG, PT}
    item::Ti
    cell::Ti
    region::Ti
    volume::TvG
    normal::Vector{TvG}
    time::Ttime
    x::Vector{Tx}
    xref::Vector{Txref}
    grid::ExtendableGrid{TvG, TiG}
    params::PT
end


"""
$(TYPEDSIGNATURES)

constructor for QPInfos
"""
function QPInfos(xgrid::ExtendableGrid{Tv, Ti}; time = 1.0, dim = size(xgrid[Coordinates], 1), T = Tv, x = ones(T, dim), params = [], kwargs...) where {Tv, Ti}
    return QPInfos{Ti, Tv, typeof(time), T, T, Tv, Ti, typeof(params)}(Ti(1), Ti(1), Ti(1), Tv(1.0), zeros(Tv, dim), time, x, ones(T, dim), xgrid, params)
end


"""
$(TYPEDSIGNATURES)

standard kernel that just copies the input to the result
"""
function standard_kernel(result, input, qpinfo)
    result .= input
    return nothing
end

"""
$(TYPEDSIGNATURES)

a kernel that acts as the constant function one
"""
function constant_one_kernel(result, qpinfo)
    result .= 1
    return nothing
end
