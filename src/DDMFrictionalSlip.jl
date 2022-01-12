module DDMFrictionalSlip

using LinearAlgebra
using DDMMesh
using DDMCollocationMatrix
# using SparseArrays
using IterativeSolvers
# using Preconditioners
# using AlgebraicMultigrid
using TimerOutputs
using UnPack

include("problem.jl")
export MechanicalProblem, HydroMechanicalProblem

include("solver.jl")
export IterativeSolver
# export solve!

include("assembly.jl")

include("time_stepper.jl")
export TimeSequence

include("executioner.jl")
export run!

end
