function assembleResidualAndJacobian!(solver::Solver{T}, problem::MechanicalProblem{T}, timer::TimerOutput) where T <: Real
    @timeit timer "Assembly" begin
        # Jacobian = Collocation matrix
        @timeit timer "Jacobian" begin
            @timeit timer "Collocation" collocationMatrix!(solver.mat, problem.mesh, problem.order)
        end

        # Residuals = collocation stress - imposed stress
        @timeit timer "Residuals" begin
            # Collocation stress
            @timeit timer "Collocation" mul!(solver.rhs, solver.mat, solver.solution)
            # Imposed stress
            @timeit timer "Imposed stress" solver.rhs .-= problem.stress_function.(problem.x)
        end
    end
    return
end

function assembleResidualAndJacobian!(solver::Solver{T}, problem::HydroMechanicalProblem{T}, timer::TimerOutput) where T <: Real
    @timeit timer "Assembly" begin
        # Jacobian = Collocation matrix
        @timeit timer "Jacobian" begin
            @unpack μ = problem.par
            # Index range
            idx = 1:Int(problem.n_dofs / 2)
            @timeit timer "Collocation normal" collocationMatrix!(view(solver.mat, idx, idx), problem.mesh, problem.order; μ = μ)
            # Index range
            idx = (Int(problem.n_dofs / 2) + 1):problem.n_dofs
            @timeit timer "Collocation shear" collocationMatrix!(view(solver.mat, idx, idx), problem.mesh, problem.order; μ = μ)
        end

        # Residuals = collocation stress - frictional stress
        @timeit timer "Residuals" begin
            # Collocation stress
            @timeit timer "Collocation" mul!(solver.rhs, solver.mat, solver.solution)
        end
        
        # Frictional stress
        @timeit timer "Frictional stress" frictionalConstraints(solver, problem)
    end
    return
end

function frictionalConstraints(solver::Solver{T}, problem::HydroMechanicalProblem{T}) where T <: Real
    # Number of collocation points
    n_cp = problem.order + 1
    # Number of dofs per variable
    n_dofs_per_var = problem.mesh.n_elems * n_cp
    # Unpack parameters
    @unpack f, λ, μ, h = problem.par
    # Loop over elements
    for elem in problem.mesh.elems
        for i in 1:n_cp
            # Effective idx
            i_loc = (elem.id - 1) * n_cp + i
            i_glo = n_dofs_per_var + (elem.id - 1) * n_cp + i

            # Compute trial elastic stresses
            # Normal stress
            problem.stress[0][i_loc] = problem.stress_old[0][i_loc] + λ / h * solver.solution[i_loc]
            # Shear stress
            problem.stress[1][i_loc] = problem.stress_old[1][i_loc] + μ / h * solver.solution[i_glo]

            # # Yield stress
            # F = problem.stress[1][i_loc] - f * problem.stress[0][i_loc]

            # # Check plastic yield
            # if (F > 0.0) # plastic update
            #     # Correct stress
            #     problem.stress[1][i_loc] -= F

            #     # Jacobian
            #     # Normal stress
            #     solver.mat[i_loc, i_loc] += λ / h
            #     # Shear stress
            #     solver.mat[i_glo, i_loc] += f * λ / h
            # else # only elastic
                # Jacobian
                # Normal stress
                solver.mat[i_loc, i_loc] += λ / h
                # Shear stress
                solver.mat[i_glo, i_glo] += μ / h
            # end
            
            # Residuals
            # Normal stress
            solver.rhs[i_loc] -= (problem.p - problem.p_old) + (problem.stress[0][i_loc] -problem.stress_old[0][i_loc]) 
            # Shear stress
            solver.rhs[i_glo] -= (problem.stress[1][i_loc] -problem.stress_old[1][i_loc]) 
        end
    end

end