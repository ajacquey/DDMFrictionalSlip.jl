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
                            collocationMatrix!(view(solver.mat, idx_i, idx_j), problem.mesh, problem.order; μ = problem.μ)
                        else
                            # Off diagonal blocks
                            view(solver.mat, idx_i, idx_j) .= spdiagm(0 => zeros(T, problem.n_cps))
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