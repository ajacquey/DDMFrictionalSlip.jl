module DDMFrictionalSlip

using LinearAlgebra
using SparseArrays
using IterativeSolvers
using Preconditioners
using IncompleteLU
# using AlgebraicMultigrid
using TimerOutputs
using UnPack
using DelimitedFiles

# Mesh
include("mesh/point.jl")
include("mesh/elem.jl")
include("mesh/mesh.jl")
export Mesh1D

# Collocation
include("collocation/utils.jl")
include("collocation/collocation_points.jl")
include("collocation/elastic_kernels.jl")
include("collocation/collocation_matrix.jl")

include("time_stepper.jl")
export TimeSequence

include("variable.jl")

include("problem.jl")
export Problem, addVariable!, addAuxVariable!

include("solver.jl")
export IterativeSolver
# export solve!

include("assembly.jl")

include("output.jl")
export DomainOutput, MaximumOutput

include("executioner.jl")
export run!

end
