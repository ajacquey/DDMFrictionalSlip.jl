# function assembleInitialJacobian(problem::AbstractProblem{T})::AbstractMatrix{T} where {T<:Real}
#     # Number of dof per variable
#     n_cps = problem.n_cps
#     # Number of dofs
#     n_dofs = problem.n_dofs
#     # Numver of vars
#     n_vars = length(problem.vars)

#     # Diagonal blocks - Collocation matrix 
#     E = sprand(T, n_cps, n_cps, 1.0)
#     collocationMatrix!(E, problem.mesh, problem.order; μ = problem.μ)
    
#     # Off-digonal blocks
#     I = spdiagm(0 => zeros(T, n_cps))

#     # List of blocks - CHECK IF WE CAN AVOID DECLARING THAT!!!
#     blocks = Vector{AbstractMatrix{T}}(undef, 0)

#     # Build matrix 
#     for i in 1:n_vars
#         for j in 1:n_vars
#             if i == j
#                 push!(blocks, E)
#             else
#                 push!(blocks, I)
#             end
#         end
#     end

#     return hvcat(Tuple(n_vars for _ in 1:n_vars), blocks...)
# end

function assembleInitialJacobian(problem::AbstractProblem{T})::AbstractMatrix{T} where {T<:Real}
    # Number of dof per variable
    n_cps = problem.n_cps
    # Number of dofs
    n_dofs = problem.n_dofs
    # Numver of vars
    n_vars = length(problem.vars)

    E = zeros(T, n_dofs, n_dofs)
    # Diagonal blocks - Collocation matrix 
    for var_i in problem.vars
        # Index range
        idx_i = ((var_i.id-1)*problem.n_cps+1):(var_i.id*problem.n_cps)
        for var_j in problem.vars
            # Index range
            idx_j = ((var_j.id-1)*problem.n_cps+1):(var_j.id*problem.n_cps)
            if (var_i == var_j)
                collocationMatrix!(view(E, idx_i, idx_j), problem.mesh, problem.order, μ = problem.μ)
            end
        end
    end
    return E
end

function assembleResidualAndJacobian!(solver::Solver{T}, problem::AbstractProblem{T}, timer::TimerOutput) where {T<:Real}
    @timeit timer "Assembly" begin
        # Jacobian = Collocation matrix
        @timeit timer "Jacobian" begin
            @timeit timer "Collocation" begin
                for var_i in problem.vars
                    # Index range
                    idx_i = ((var_i.id-1)*problem.n_cps+1):(var_i.id*problem.n_cps)
                    for var_j in problem.vars
                        # Index range
                        idx_j = ((var_j.id-1)*problem.n_cps+1):(var_j.id*problem.n_cps)
                        if var_i == var_j
                            # Collocation matrix
                            # Start building matrix element-wise
                            n_cp = problem.order + 1
                            Threads.@threads for elem in problem.mesh.elems
                                idx = (elem.id-1)*n_cp+1:elem.id*n_cp
                                localCollocationMatrix!(view(view(solver.mat, idx_i, idx_j), idx, idx), elem, elem, problem.order, problem.μ)
                                # collocationMatrix!(view(solver.mat, idx_i, idx_j), problem.mesh, problem.order; μ = problem.μ)
                            end
                        else
                            # # Off diagonal blocks
                            # view(solver.mat, idx_i, idx_j) .= spdiagm(0 => zeros(T, problem.n_cps))
                            # Off diagonal blocks
                            view(solver.mat, idx_i, idx_j)[diagind(view(solver.mat, idx_i, idx_j))] .= 0.0
                        end

                    end
                end
            end
        end

        # Residuals = collocation stress - imposed stress
        @timeit timer "Residuals" begin
            # Collocation stress
            @timeit timer "Collocation" mul!(solver.rhs, solver.mat, solver.solution)
            # Imposed stress
            @timeit timer "Stress residuals" kernelResiduals!(solver, problem)
        end

        @timeit timer "Jacobian" begin

            @timeit timer "Stress jacobian" kernelJacobian!(solver, problem)
        end
    end
    return
end

function kernelResiduals!(solver::Solver{T}, problem::AbstractProblem{T}) where {T<:Real}
    # Number of collocation points
    n_cp = problem.order + 1

    # Loop over elements
    Threads.@threads for elem in problem.mesh.elems
        for kernel in problem.kernels
            for i in 1:n_cp
                # Effective idx
                idx = (kernel.u.id - 1) * problem.n_cps + (elem.id - 1) * n_cp + i
                idx_var = (elem.id - 1) * n_cp + i

                solver.rhs[idx] -= computeCpResidual(kernel, problem.time, problem.x[idx_var], idx_var)
            end
        end
    end

    return
end

function kernelJacobian!(solver::Solver{T}, problem::AbstractProblem{T}) where {T<:Real}
    # Number of collocation points
    n_cp = problem.order + 1

    # Loop over elements
    # THIS IS NOT THREAD SAFE BECAUSE WE MODIFY nzval
    cond = Base.Threads.Condition()
    Threads.@threads for elem in problem.mesh.elems
        lock(cond)
        for kernel in problem.kernels
            for i in 1:n_cp
                # Effective idx
                idx_i = (kernel.u.id - 1) * problem.n_cps + (elem.id - 1) * n_cp + i
                idx_var = (elem.id - 1) * n_cp + i

                for j_var in problem.vars
                    # Effective idx
                    idx_j = (j_var.id - 1) * problem.n_cps + (elem.id - 1) * n_cp + i
                    
                    solver.mat[idx_i, idx_j] -= computeCpJacobian(kernel, j_var, problem.time, problem.x[idx_var], idx_var)
                end
            end
        end
        unlock(cond)
    end

    return
end