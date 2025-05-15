"""
	ExtendableFEMBase

$(read(joinpath(@__DIR__, "..", "README.md"), String))
"""
module ExtendableFEMBase

using DocStringExtensions: DocStringExtensions, TYPEDEF, TYPEDSIGNATURES
using ExtendableGrids: ExtendableGrids, AT_NODES, AbstractElementGeometry,
    AbstractElementGeometry0D, AbstractElementGeometry1D,
    AbstractElementGeometry2D, AbstractElementGeometry3D,
    AbstractGridAdjacency, AbstractGridComponent,
    AbstractGridFloatArray2D, AbstractGridIntegerArray2D,
    Adjacency, AssemblyType, BEdgeEdges, BEdgeGeometries,
    BEdgeRegions, BEdgeVolumes, BFaceEdges, BFaceFaces,
    BFaceGeometries, BFaceParents, BFaceRegions,
    BFaceVolumes, CellEdgeSigns, CellEdges,
    CellFaceOrientations, CellFaceSigns, CellFaces,
    CellFinder, CellGeometries, CellNodes, CellParents,
    CellRegions, CellVolumes, Coordinates, Edge1D,
    EdgeCells, EdgeGeometries, EdgeNodes, EdgeTangents,
    EdgeParents, EdgeVolumes, ElementGeometries, ExtendableGrid,
    FaceCells, FaceEdgeSigns, FaceEdges, FaceGeometries,
    FaceNodes, FaceNormals, FaceParents, FaceRegions,
    FaceVolumes, GridComponentGeometries4AssemblyType,
    NodeCells,
    GridComponentNodes4AssemblyType,
    GridComponentRegions4AssemblyType,
    GridComponentUniqueGeometries4AssemblyType,
    GridComponentVolumes4AssemblyType, GridEGTypes,
    GridRegionTypes, Hexahedron3D, L2GTransformer,
    ON_BEDGES, ON_BFACES, ON_CELLS, ON_EDGES, ON_FACES,
    ON_IFACES, Parallelepiped3D, Parallelogram2D,
    ParentGrid, ParentGridRelation, Quadrilateral2D,
    SerialVariableTargetAdjacency, SubGrid, Tetrahedron3D,
    Triangle2D, UniqueBEdgeGeometries,
    UniqueBFaceGeometries, UniqueCellGeometries,
    UniqueEdgeGeometries, UniqueFaceGeometries,
    VariableTargetAdjacency, Vertex0D, append!, atranspose,
    dim_element, eval_trafo!, gFindLocal!, interpolate!,
    local_celledgenodes, local_cellfacenodes, mapderiv!,
    max_num_targets_per_source, num_cells, num_edges,
    num_faces, num_nodes, num_sources, num_targets,
    reference_domain, subgrid, unique,
    update_trafo!
using ExtendableSparse: ExtendableSparse, ExtendableSparseMatrix, flush!,
    AbstractExtendableSparseMatrixCSC, ExtendableSparseMatrixCSC, MTExtendableSparseMatrixCSC,
    rawupdateindex!
using ForwardDiff: ForwardDiff, DiffResults
using LinearAlgebra: LinearAlgebra, convert, det, diagm, dot, eigen, ldiv!, lu,
    mul!, norm, transpose
using Polynomials: Polynomials, Polynomial, coeffs
using Printf: Printf, @printf
using SparseArrays: SparseArrays, AbstractSparseArray, AbstractSparseMatrix,
    SparseMatrixCSC, nzrange, rowvals
using SpecialPolynomials: SpecialPolynomials, ShiftedLegendre, basis

include("functionoperators.jl")
export AbstractFunctionOperator
export StandardFunctionOperator
export Identity, IdentityComponent
export NormalFlux, TangentFlux
export Gradient
export SymmetricGradient, TangentialGradient
export Divergence
export CurlScalar, Curl2D, Curl3D
export Laplacian, Hessian, SymmetricHessian
export Trace, Deviator
export NeededDerivative4Operator, Length4Operator, QuadratureOrderShift4Operator, DefaultName4Operator
export OperatorPair, OperatorTriple

include("qpinfos.jl")
#include("feview.jl")
export QPInfos

include("quadrature.jl")
export QuadratureRule
export VertexRule
export integrate!, integrate, ref_integrate!

include("finiteelements.jl") # also includes dofmaps.jl and feevaluator*.jl
export DofMap
export CellDofs, FaceDofs, EdgeDofs, BFaceDofs, BEdgeDofs
export CellDofsParent, FaceDofsParent, EdgeDofsParent, BFaceDofsParent, BEdgeDofsParent
export DofMapTypes
export Dofmap4AssemblyType, ItemGeometries4DofMap, EffAT4AssemblyType, ParentDofmap4Dofmap
export AbstractFiniteElement
export FESpace, FESpaces, get_AT, get_FEType
export boundarydofs, ndofs, broken

export AbstractH1FiniteElement
export H1BUBBLE, L2P0, H1P1, H1P2, H1P2B, H1MINI, H1CR, H1P3, H1Pk
export L2P1
export H1Q1, H1Q2

export AbstractH1FiniteElementWithCoefficients
export H1BR, H1P1TEB

export AbstractHdivFiniteElement
export HDIVRT0, HDIVBDM1, HDIVRT1, HDIVBDM2, HDIVRTk
export HDIVRTkENRICH

export AbstractHcurlFiniteElement
export HCURLN0, HCURLN1

export get_AT
export get_polynomialorder, get_ndofs, get_ndofs_all
export get_ncomponents, get_edim
export get_basis, get_coefficients, get_basissubset

export NodalInterpolator

export interpolate! # must be defined separately by each FEdefinition
export nodevalues, continuify
export nodevalues!, nodevalues_subset!
export nodevalues_view

export interpolator_matrix

export FEVectorBlock, FEVector
export dot, norm, norms
export FEMatrixBlock, FEMatrix, _addnz
export fill!, addblock!, addblock_matmul!, lrmatmul, mul!, add!, apply_penalties!, coffset
export submatrix

export displace_mesh, displace_mesh!

include("reconstructionhandlers.jl")
export ReconstructionHandler, get_rcoefficients!

include("feevaluator.jl")
export FEEvaluator, update_basis!, eval_febe!

include("reconstructionoperators.jl")
export Reconstruct

include("accumvector.jl")
export AccumulatingVector


#
# Print default dict for solver parameters into docstrings
#
function _myprint(dict::Dict{Symbol, Tuple{Any, String}})
    lines_out = IOBuffer()
    for (k, v) in dict
        if typeof(v[1]) <: String
            println(lines_out, "  - $(k): $(v[2]). Default: ''$(v[1])''\n")
        else
            println(lines_out, "  - $(k): $(v[2]). Default: $(v[1])\n")
        end
    end
    return String(take!(lines_out))
end
#
# Update solver params from dict
#
function _update_params!(parameters, kwargs)
    for (k, v) in kwargs
        parameters[Symbol(k)] = v
    end
    return nothing
end

include("segment_integrator.jl")
export SegmentIntegrator, initialize!, integrate_segment!

include("point_evaluator.jl")
export PointEvaluator, evaluate!, evaluate_bary!, eval_func, eval_func_bary

include("lazy_interpolate.jl")
export lazy_interpolate!


# ExtendableFEMBaseUnicodePlotsExt extension


"""
````
function unicode_gridplot(
	xgrid::ExtendableGrid;
	title = "gridplot",
	resolution = (40,20),
	color = (200,200,200),
	bface_color = (255,0,0),
	CanvasType = BrailleCanvas,
	plot_based = ON_CELLS,   # or ON_FACES/ON_EDGES
	kwargs...
````

(via extension that requires UnicodePlots)

Plots the grid on a UnicodePlots canvas (default: BrailleCanvas) by drawing all edges in the triangulation.
"""
function unicode_gridplot end

"""
````
function unicode_scalarplot(
	u::FEVectorBlock; 
	components = 1:get_ncomponents(u),
	abs = false,
	resolution = (30,30),
	colormap = :viridis,
	title = u.name,
	kwargs...)
````

(via extension that requires UnicodePlots)

Plots components of the finite element function in the FEVectorBlock u by
using a lazy_interpolate! onto a coarse uniform mesh and UnicodePlots.jl
lineplot or heatmap for 1D or 2D, respectively.

In 1D all components all plotted in the same lineplot, while
in 2D all components are plotted in a separate heatmap.

If abs = true, only the absolute value over the components is plotted.
"""
function unicode_scalarplot end
export unicode_gridplot, unicode_scalarplot


end # module ExtendableFEMBase.
