abstract type Solver{T<:Real} end

mutable struct IterativeSolver{T<:Real} <: Solver{T}
    "The Problem this solver will act on"
    problem::AbstractProblem{T}

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
end

"""
Constructor for Solver
"""
function IterativeSolver(problem::AbstractProblem{T}; nl_max_iters::Int64 = 100, nl_abs_tol::T = 1.0e-10, nl_rel_tol::T = 1.0e-10) where {T<:Real}
    return IterativeSolver{T}(
        problem,
        spzeros(T, problem.n_dofs, problem.n_dofs),
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
    # To call after initializing the problem
    n_dofs = solver.problem.n_dofs
    solver.mat = spzeros(T, n_dofs, n_dofs)
    reinitSparsityPattern!(solver.mat, solver.problem; first_run = true)
    solver.rhs = zeros(T, n_dofs)
    solver.solution = zeros(T, n_dofs)
    solver.initialized = true
    return
end

function reinit!(solver::IterativeSolver{T}) where {T<:Real}
    fill!(solver.rhs, 0.0)
    reinitSparsityPattern!(solver.mat, solver.problem)
    fill!(solver.solution, 0.0)
    return
end

function reinitSparsityPattern!(mat::AbstractMatrix{T}, problem::AbstractProblem{T}; first_run = false) where {T<:Real}
    # Loop over variables
    for var_i in problem.vars
        # Index range
        idx_i = ((var_i.id-1)*problem.n_cps+1):(var_i.id*problem.n_cps)
        for var_j in problem.vars
            # Index range
            idx_j = ((var_j.id-1)*problem.n_cps+1):(var_j.id*problem.n_cps)
            if first_run
                if var_i.id == var_j.id # Diagonal blocks
                    view(mat, idx_i, idx_j) .= sprand(T, problem.n_cps, problem.n_cps, 1.0)
                end
            end
            if var_i.id != var_j.id # Off-diagonal blocks
                view(mat, idx_i, idx_j) .= spdiagm(0 => zeros(T, problem.n_cps))
            end
        end
    end

    return
end

"Solve the problem using the IterativeSolvers package"
function solve!(solver::IterativeSolver{T}, timer::TimerOutput; log::Bool = true) where {T<:Real}
    ##### Newton loop #####
    # Non-linear iterations
    nl_iter = 0
    # Declare solution
    dx = zeros(T, length(solver.rhs))
    # Initial residual
    assembleResidualAndJacobian!(solver, solver.problem, timer)
    r = norm(solver.rhs)
    r0 = r
    # Pardiso
    # ps = MKLPardisoSolver()
    # Preconditioner 
    # @timeit timer "Preconditionning" precond = ilu(solver.mat, τ = 0.01)
    # @timeit timer "Preconditionning" precond = JacobiPreconditioner(solver.mat)
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
        @timeit timer "Solve" dx, ch = bicgstabl!(dx, solver.mat, -solver.rhs; log = true, verbose = false, abstol = 1.0e-10, reltol = 1.0e-10)
        # @timeit timer "Solve" dx = jacobi!(dx, solver.mat, -solver.rhs; maxiter = 200)
        # @timeit timer "Solve" dx = Pardiso.solve(ps, solver.mat, -solver.rhs)
        if log
            if ch.isconverged
                println("    -> Linear Solve converged after ", ch.iters, " iterations")
            else
                println("    -> Linear Solve did NOT converge after ", ch.iters, " iterations")
            end
        end

        # Update solution
        solver.solution .+= dx
        # Update problem
        @timeit timer "Update problem" update!(solver.problem, solver)
        # Update residuals and jacobian
        assembleResidualAndJacobian!(solver, solver.problem, timer)
        r = norm(solver.rhs)
        # Preconditioner 
        # @timeit timer "Preconditionning" precond = ilu(solver.mat, τ = 0.01)

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