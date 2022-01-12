function run!(problem::MechanicalProblem{T}, solver::Solver{T}) where T <: Real
    # Timer
    timer = TimerOutput()

    # Display some information about simulation

    # Initialize solver
    @timeit timer "Initialize Solver" initialize!(solver)

    # Steady state problem
    solve!(solver, timer)

    # End of simulation information - TimerOutputs
    show(timer, title="Performance graph")
    println()

    return problem.x, solver.solution
end

function run!(problem::HydroMechanicalProblem{T}, solver::Solver{T}, time_stepper::TimeStepper) where T <: Real
    # Timer
    timer = TimerOutput()

    # Display some information about simulation

    # Initialize solver
    @timeit timer "Initialize Solver" initialize!(solver)

    # Apply ICs

    # Transient problem
    it = 0 # time step
    while time_stepper.time < time_stepper.end_time
        # Save old state
        @timeit timer "Reinitialize problem" reinit!(problem)

        # Pressure update
        @timeit timer "Calculate fluid pressure" problem.pressure .= problem.pressure_function.(problem.x, time_stepper.time)

        # Actual solve
        solve!(solver, timer)

        # Update problem
        @timeit timer "Update problem" update!(problem, solver)
        
        # Update time stepper
        it += 1
        time_stepper.time = time_stepper.time_seq[it]

        # Reinit solver
        @timeit timer "Reinitialize Solver" reinit!(solver)
    end

    # End of simulation information - TimerOutputs
    show(timer, title="Performance graph")
    println()

    return
end

function update!(problem::HydroMechanicalProblem{T}, solver::Solver{T}) where T <: Real
    # Update displacements
    for elem in problem.mesh.elems
        for i in 1:n_cp
            # Effective idx
            i_loc = (elem.id - 1) * n_cp + i
            i_glo = n_dofs_per_var + (elem.id - 1) * n_cp + i

            # Update displacements
            problem.disp[1][i_loc] = problem.disp_old[1][i_loc] + solver.solution[i_loc]
            problem.disp[2][i_loc] = problem.disp_old[2][i_loc] + solver.solution[i_glo]
        end
    end
    return
end

