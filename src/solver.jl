abstract type Solver{T<:Real} end

" Must be overriden by Solver implementations to actually do the solve "
function solve!(solver::Solver)
    throw(MethodError(solve!, solver))
end

# mutable struct NLSolver{T <: Real } <: Solver{T}
#     "The Problem this solver will act on"
#     problem::Problem{T}

#     "The jacobian matrix"
#     mat::Matrix{T}

#     "The residual vector"
#     rhs::Vector{T}

#     "The solution vector"
#     solution::Vector{T}

#     "Whether or not this Solver has been initialized"
#     initialized::Bool

#     "Constructor"
#     NLSolver(pro::Problem) = new(
#         pro, 
#         Matrix(T, problem.n_dofs, problem.n_dofs), 
#         Vector{T}(undef, problem.n_dofs), 
#         Vector{T}(undef, problem.n_dofs), 
#         false
#     )
# end

# "Initializes the jacobian and residuals to the correct size"
# function initialize!(solver::NLSolver{T}) where T <: Real
#     @assert solver.problem.initialized

#     n_dofs = solver.problem.n_dofs

#     solver.mat = Matrix(Float64, problem.n_dofs, problem.n_dofs)
#     solver.rhs = zeros(Float64, n_dofs)
#     solver.solution = zeros(Float64, n_dofs)
#     solver.initialized = true
# end

# "Solve the system using the NLSolve package"
# function solve!(solver::NLSolver{T}; assemble=true) where T <: Real
#     if !solver.initialized
#         initialize!(solver)
#     end

#     if assemble
#         assembleResidualAndJacobian!(solver, solver.system)
#     end

#     # Use of the NLSolve package
#     sol = nlsolve(df, solver.solution)
#     solver.solution = sol.zero
# end

mutable struct IterativeSolver{T<:Real} <: Solver{T}
    "The Problem this solver will act on"
    problem::Problem{T}

    "The jacobian matrix"
    mat::AbstractMatrix{T}

    "The residual vector"
    rhs::Vector{T}

    "The solution vector"
    solution::Vector{T}

    "Whether or not this Solver has been initialized"
    initialized::Bool

    "Maximum number of nonlinear iterations"
    nl_max_iters::Int64

    "Nonlinear absolute tolerance"
    nl_abs_tol::T

    "Nonlinear relative tolerance"
    nl_rel_tol::T

    # "Constructor"
    # IterativeSolver(pro::Problem{T}; nl_max_iters::Int64=100, nl_abs_tol::T=1.0e-10, nl_rel_tol::T=1.0e-10) where T <: Real = 
    # new{T}(
    #     pro, 
    #     Matrix{T}(undef, pro.n_dofs, pro.n_dofs), 
    #     Vector{T}(undef, pro.n_dofs), 
    #     Vector{T}(undef, pro.n_dofs), 
    #     false,
    #     nl_max_iters,
    #     nl_abs_tol,
    #     nl_rel_tol,
    # )
end

"""
Constructor for MechanicalProblem
"""
function IterativeSolver(problem::MechanicalProblem{T}; nl_max_iters::Int64 = 100, nl_abs_tol::T = 1.0e-10, nl_rel_tol::T = 1.0e-10) where {T<:Real}
    return IterativeSolver{T}(
        problem,
        Matrix{T}(undef, problem.n_dofs, problem.n_dofs),
        Vector{T}(undef, problem.n_dofs),
        Vector{T}(undef, problem.n_dofs),
        false,
        nl_max_iters,
        nl_abs_tol,
        nl_rel_tol,
    )
end

"""
Custom `show` function for `IterativeSolver` that prints some information.
"""
function Base.show(io::IO, solv::IterativeSolver{T}) where {T<:Real}
    println("Solver information:")
    # println("  ", show(solv.pro))
    println("     res size: $(length(solv.rhs))")
    println("     Jac size: $(size(solv.mat))")
end

"Initializes the jacobian and residuals to the correct size"
function initialize!(solver::IterativeSolver{T}) where {T<:Real}
    # @assert solver.system.initialized

    n_dofs = solver.problem.n_dofs

    solver.mat = zeros(T, n_dofs, n_dofs)
    solver.rhs = zeros(T, n_dofs)
    solver.solution = zeros(T, n_dofs)
    solver.initialized = true
    return
end

function reinit!(solver::IterativeSolver{T}) where {T<:Real}
    fill!(solver.rhs, 0.0)
    fill!(solver.mat, 0.0)
    fill!(solver.solution, 0.0)
    return
end

"Solve the problem using the IterativeSolvers package"
function solve!(solver::IterativeSolver{T}, timer::TimerOutput; log::Bool = true) where {T<:Real}
    ##### Newton loop #####
    # Non-linear iterations
    nl_iter = 0
    # Initial residual
    assembleResidualAndJacobian!(solver, solver.problem, timer)
    r = norm(solver.rhs)
    r0 = r
    # Preconditioner 
    # @timeit timer "Preconditionning" precond = CholeskyPreconditioner(solver.mat, 2)
    if log
        println("  ", 0, " Nonlinear Iteration: |R| = ", r0)
    end
    # Main loop
    while (nl_iter <= solver.nl_max_iters)
        # Check convergence
        if (r <= solver.nl_abs_tol)
            if log
                println("Nonlinear Solve converged with absolute tolerance!")
                println()
            end
            return
        end
        if (r / r0 <= solver.nl_rel_tol)
            if log
                println("Nonlinear Solve converged with relative tolerance!")
                println()
            end
            return
        end
        # Linear Solve
        @timeit timer "Solve" dx, ch = cg(solver.mat, -solver.rhs; log = true, verbose = false)
        if log
            println("    -> Linear Solve converged after ", ch.iters, " iterations")
        end

        # Update solution
        solver.solution += dx
        # Update residuals and jacobian
        assembleResidualAndJacobian!(solver, solver.problem, timer)
        r = norm(solver.rhs)
        # # Preconditioner 
        # @timeit timer "Preconditionning" UpdatePreconditioner!(precond, -solver.mat)

        nl_iter += 1
        if log
            println("  ", nl_iter, " Nonlinear Iteration: |R| = ", norm(r))
        end
    end
    # Error if exceeded maximum number of iterations
    if (nl_iter > solver.nl_max_iters)
        throw(ErrorException("Exceeded the maximum number of nonlinear iterations!"))
    end
end