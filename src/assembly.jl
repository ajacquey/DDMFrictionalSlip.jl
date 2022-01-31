function assembleResidualAndJacobian!(solver::Solver{T}, problem::AbstractProblem{T}, timer::TimerOutput) where {T<:Real}
    @timeit timer "Assembly" begin
        # Jacobian = Collocation matrix
        @timeit timer "Jacobian" begin
            @timeit timer "Collocation" begin
                for var in problem.vars
                    # Index range
                    idx = ((var.id-1)*problem.n_cps+1):(var.id*problem.n_cps)
                    collocationMatrix!(view(solver.mat, idx, idx), problem.mesh, problem.order; μ = problem.μ)
                end
            end
        end

        # Residuals = collocation stress - imposed stress
        @timeit timer "Residuals" begin
            # Collocation stress
            @timeit timer "Collocation" mul!(solver.rhs, solver.mat, solver.solution)
            # Imposed stress
            # @timeit timer "Imposed stress" solver.rhs .-= problem.stress_function.(problem.x)
            @timeit timer "Stress residuals" stressResiduals!(solver, problem)
        end

        @timeit timer "Jacobian" begin

            @timeit timer "Stress jacobian" stressJacobian!(solver, problem)
        end
    end
    return
end

function stressResiduals!(solver::Solver{T}, problem::SteadyProblem{T}) where {T<:Real}
    # Number of collocation points
    n_cp = problem.order + 1

    # Loop over elements
    Threads.@threads for elem in problem.mesh.elems
        for i in 1:n_cp
            for var in problem.vars
                # Effective idx
                idx = (var.id - 1) * problem.n_cps + (elem.id - 1) * n_cp + i
                idx_var = (elem.id - 1) * n_cp + i

                for aux_var in problem.aux_vars
                    if aux_var.ass_var == var.sym
                        solver.rhs[idx] -= aux_var.func(problem.x[idx_var], solver.solution[idx])
                    end
                end
            end
        end
    end

    return
end

function stressResiduals!(solver::Solver{T}, problem::TransientProblem{T}) where {T<:Real}
    # Number of collocation points
    n_cp = problem.order + 1

    # Loop over elements
    Threads.@threads for elem in problem.mesh.elems
        for i in 1:n_cp
            for var in problem.vars
                # Effective idx
                idx = (var.id - 1) * problem.n_cps + (elem.id - 1) * n_cp + i
                idx_var = (elem.id - 1) * n_cp + i

                for aux_var in problem.aux_vars
                    if aux_var.ass_var == var.sym
                        solver.rhs[idx] -= (aux_var.func(aux_var.u_old[idx_var], problem.x[idx_var], problem.time, solver.solution[idx]) - aux_var.u_old[idx_var])
                    end
                end
            end
        end
    end

    return
end

function stressJacobian!(solver::Solver{T}, problem::SteadyProblem{T}) where {T<:Real}
    # Number of collocation points
    n_cp = problem.order + 1

    # Loop over elements
    Threads.@threads for elem in problem.mesh.elems
        for i in 1:n_cp
            for var in problem.vars
                # Effective idx
                idx = (var.id - 1) * problem.n_cps + (elem.id - 1) * n_cp + i
                idx_var = (elem.id - 1) * n_cp + i

                for aux_var in problem.aux_vars
                    if aux_var.ass_var == var.sym
                        solver.mat[idx, idx] -= aux_var.dfunc(problem.x[idx_var], solver.solution[idx])
                    end
                end
            end
        end
    end

    return
end

function stressJacobian!(solver::Solver{T}, problem::TransientProblem{T}) where {T<:Real}
    # Number of collocation points
    n_cp = problem.order + 1

    # Loop over elements
    Threads.@threads for elem in problem.mesh.elems
        for i in 1:n_cp
            for var in problem.vars
                # Effective idx
                idx = (var.id - 1) * problem.n_cps + (elem.id - 1) * n_cp + i
                idx_var = (elem.id - 1) * n_cp + i

                for aux_var in problem.aux_vars
                    if aux_var.ass_var == var.sym
                        solver.mat[idx, idx] -= aux_var.dfunc(aux_var.u_old[idx_var], problem.x[idx_var], problem.time, solver.solution[idx])
                    end
                end
            end
        end
    end

    return
end