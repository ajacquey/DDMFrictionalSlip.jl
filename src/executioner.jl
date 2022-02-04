function run!(problem::SteadyProblem{T}, solver::Solver{T}; log::Bool = true) where {T<:Real}
    # Timer
    timer = TimerOutput()

    # Initialize problem
    @timeit timer "Initialize Problem" initialize!(problem)

    # Initialize solver
    @timeit timer "Initialize Solver" initialize!(solver)

    # Display some information about simulation

    # Apply ICs
    @timeit timer "Apply Initial Conditions" applyIC!(problem)

    # Steady state problem
    solve!(solver, timer; log)

    # Final update
    @timeit timer "Final update" final_update!(problem)

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
    @timeit timer "Apply Initial Conditions" applyIC!(problem)

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

        # Final update
        @timeit timer "Final update" final_update!(problem)

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

function update!(problem::AbstractProblem{T}, solver::Solver{T}) where {T<:Real}
    n_cp = problem.order + 1
    # Update problem
    Threads.@threads for elem in problem.mesh.elems
        for i in 1:n_cp
            # Effective idx
            idx = (elem.id - 1) * n_cp + i
            # Update main variables
            for var in problem.vars
                # Effective idx
                idx_sol = (var.id - 1) * problem.n_cps + (elem.id - 1) * n_cp + i

                # Update vars
                var.value[idx] = var.value_old[idx] + solver.solution[idx_sol]

            end
            # Update aux variables
            for aux_kernel in problem.aux_kernels

                # Update auxiliary vars
                if aux_kernel.execute_on == :non_linear
                    aux_kernel.u.value[idx] = computeCpValue(aux_kernel, problem.time, problem.x[idx], idx)
                end
            end
        end
    end
    return
end

function final_update!(problem::AbstractProblem{T}) where {T<:Real}
    n_cp = problem.order + 1
    # Update auxiliary variables on time step end
    Threads.@threads for elem in problem.mesh.elems
        for i in 1:n_cp
            # Effective idx
            idx = (elem.id - 1) * n_cp + i
            # Update aux variables
            for aux_kernel in problem.aux_kernels

                # Update auxiliary vars
                if aux_kernel.execute_on == :time_step_end
                    aux_kernel.u.value[idx] = computeCpValue(aux_kernel, problem.time, problem.x[idx], idx)
                end
            end
        end
    end
    return
end