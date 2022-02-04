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

include("kernels/kernel.jl")
include("kernels/function_kernel.jl")
include("kernels/local_elastic.jl")
include("kernels/coupled_force.jl")
include("kernels/local_elasto_plastic.jl")
export FunctionKernel, LocalElasticKernel, LocalElastoPlasticKernel, CoupledForceKernel

include("aux_kernels/aux_kernel.jl")
include("aux_kernels/function_aux_kernel.jl")
include("aux_kernels/local_elastic.jl")
include("aux_kernels/local_elasto_plastic.jl")
export FunctionAuxKernel, LocalElasticAuxKernel, LocalElastoPlasticAuxKernel

include("problem.jl")
export Problem
export addVariable!, addKernel!
export addAuxVariable!, addAuxKernel!

include("solver.jl")
export IterativeSolver

include("assembly.jl")

include("output.jl")
export DomainOutput, MaximumOutput

include("executioner.jl")
export run!

end
