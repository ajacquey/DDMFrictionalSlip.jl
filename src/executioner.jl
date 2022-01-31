function run!(problem::SteadyProblem{T}, solver::Solver{T}; log::Bool = true) where {T<:Real}
    # Timer
    timer = TimerOutput()

    # Initialize problem
    @timeit timer "Initialize Problem" initialize!(problem)

    # Initialize solver
    @timeit timer "Initialize Solver" initialize!(solver)

    # Display some information about simulation

    # Apply ICs

    # Steady state problem
    solve!(solver, timer; log)

    # Update problem
    @timeit timer "Update problem" update!(problem, solver)

    # End of simulation information - TimerOutputs
    if log
        show(timer, title = "Performance graph")
        println()
    end

    return
end

function run!(problem::TransientProblem{T}, solver::Solver{T}, time_stepper::TimeStepper; log::Bool = true, outputs::Vector{AbstractOutput} = Vector{AbstractOutput}(undef, 0)) where {T<:Real}
    # Timer
    timer = TimerOutput()

    # Initialize problem
    @timeit timer "Initialize Problem" initialize!(problem)

    # Initialize solver
    @timeit timer "Initialize Solver" initialize!(solver)

    # Initialize outputs
    if ~isempty(outputs)
        @timeit timer "Initialize Outputs" initializeOutputs!(outputs, problem)
    end

    # Display some information about simulation

    # Apply ICs

    # Transient problem
    problem.time = time_stepper.start_time
    problem.time_old = time_stepper.start_time
    if log
        println("Time Step ", problem.time_step, ": time = ", problem.time, " dt = ", 0.0)
        println()
    end

    while problem.time < (time_stepper.end_time - time_stepper.tol)
        # Update iteration number
        problem.time_step = problem.time_step + 1

        # Save old state
        @timeit timer "Reinitialize problem" reinit!(problem, time_stepper)

        # Print time step information
        if log
            println("Time Step ", problem.time_step, ": time = ", problem.time, " dt = ", problem.dt)
        end

        # Actual solve
        solve!(solver, timer; log)

        # Update problem
        @timeit timer "Update problem" update!(problem, solver)

        # Output
        if ~isempty(outputs)
            @timeit timer "Outputs" outputResults!(outputs, problem)
        end

        # Reinit solver
        @timeit timer "Reinitialize Solver" reinit!(solver)
    end

    # End of simulation information - TimerOutputs
    if log
        show(timer, title = "Performance graph")
        println()
    end

    return
end

function update!(problem::SteadyProblem{T}, solver::Solver{T}) where {T<:Real}
    n_cp = problem.order + 1
    # Update displacements
    Threads.@threads for elem in problem.mesh.elems
        for i in 1:n_cp
            for var in problem.vars
                # Effective idx
                idx = (var.id - 1) * problem.n_cps + (elem.id - 1) * n_cp + i
                idx_var = (elem.id - 1) * n_cp + i

                # Update displacements
                var.u[idx_var] = solver.solution[idx]

                # Update stress
                for aux_var in problem.aux_vars
                    if aux_var.ass_var == var.sym
                        aux_var.u[idx_var] = aux_var.func(problem.x[idx_var], solver.solution[idx])
                    end
                end
            end
        end
    end
    return
end

function update!(problem::TransientProblem{T}, solver::Solver{T}) where {T<:Real}
    n_cp = problem.order + 1
    # Update displacements
    Threads.@threads for elem in problem.mesh.elems
        for var in problem.vars
            for i in 1:n_cp
                # Effective idx
                idx = (var.id - 1) * problem.n_cps + (elem.id - 1) * n_cp + i
                idx_var = (elem.id - 1) * n_cp + i

                # Update displacements
                var.u[idx_var] = var.u_old[idx_var] + solver.solution[idx]

                # Update stress
                for aux_var in problem.aux_vars
                    if aux_var.ass_var == var.sym
                        aux_var.u[idx_var] = aux_var.func(aux_var.u_old[idx_var], problem.x[idx_var], problem.time, solver.solution[idx])
                    end
                end
            end
        end
    end
    # Time
    problem.time_old = problem.time
    return
end