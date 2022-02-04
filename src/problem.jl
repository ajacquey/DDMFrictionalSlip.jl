abstract type AbstractProblem{T<:Real} end

mutable struct SteadyProblem{T<:Real} <: AbstractProblem{T}
    " The mesh of the problem"
    mesh::Mesh{T}

    " The order of the basis functions"
    order::Int64

    " The elastic modulus"
    μ::Float64

    " The total number of collocation points"
    n_cps::Int64

    " The total number of degrees of freedom"
    n_dofs::Int64

    " The collocation points coordinates"
    x::Vector{T}

    " The variables"
    vars::Vector{Variable{T}}

    " The auxiliary variables"
    aux_vars::Vector{AuxVariable{T}}

    " The kernels"
    kernels::Vector{AbstractKernel{T}}

    " The aux_kernels"
    aux_kernels::Vector{AbstractAuxKernel{T}}

    " The current time"
    time::T

    " A boolean to check if system is initialized"
    initialized::Bool
end

mutable struct TransientProblem{T<:Real} <: AbstractProblem{T}
    " The mesh of the problem"
    mesh::Mesh{T}

    " The order of the basis functions"
    order::Int64

    " The elastic modulus"
    μ::Float64

    " The total number of collocation points"
    n_cps::Int64

    " The total number of degrees of freedom"
    n_dofs::Int64

    " The collocation points coordinates"
    x::Vector{T}

    " The variables"
    vars::Vector{Variable{T}}

    " The auxiliary variables"
    aux_vars::Vector{AuxVariable{T}}

    " The kernels"
    kernels::Vector{AbstractKernel{T}}

    " The aux_kernels"
    aux_kernels::Vector{AbstractAuxKernel{T}}

    " The current time"
    time::T

    " The old time"
    time_old::T

    " The current time step"
    dt::T

    " The current time step"
    time_step::Int64

    " A boolean to check if system is initialized"
    initialized::Bool
end

function Problem(mesh::Mesh{T}; μ::T = 1.0, order::Int64 = 0, transient::Bool = false) where {T<:Real}
    # Collocation point coordinates
    x = collocationPointsCoordinates(mesh, order)

    # Total number of cps
    n_cps = mesh.n_elems * (order + 1)

    if (transient)
        return TransientProblem{T}(mesh, order, μ, n_cps, 0, x, Vector{Variable}(undef, 0), Vector{AuxVariable}(undef, 0), Vector{AbstractKernel}(undef, 0), Vector{AbstractAuxKernel}(undef, 0), 0.0, 0.0, 0.0, 0, false)
    else
        return SteadyProblem{T}(mesh, order, μ, n_cps, 0, x, Vector{Variable}(undef, 0), Vector{AuxVariable}(undef, 0), Vector{AbstractKernel}(undef, 0), Vector{AbstractAuxKernel}(undef, 0), 0.0, false)
    end
end

function reinit!(problem::TransientProblem{T}, time_stepper::TimeStepper{T}) where {T<:Real}
    # Time
    problem.time_old = problem.time
    problem.time = time_stepper.time_seq[problem.time_step]
    problem.dt = problem.time - problem.time_old

    n_cp = problem.order + 1

    # Move current values as old
    Threads.@threads for elem in problem.mesh.elems
        for i in 1:n_cp
            # Effective idx
            idx = (elem.id - 1) * n_cp + i
            # Variables
            for var in problem.vars
                var.value_old[idx] = var.value[idx]
            end
            # Auxiliary Variables
            for aux_var in problem.aux_vars
                aux_var.value_old[idx] = aux_var.value[idx]
            end
            # Update auxiliary variable on time step begin
            for aux_kernel in problem.aux_kernels
                if aux_kernel.execute_on == :time_step_begin
                    aux_kernel.u.value[idx] = computeCpValue(aux_kernel, problem.time, problem.x[idx], idx)
                end
            end
        end
    end

    return
end

function default_ic(x::Vector{T}) where {T<:Real}
    return zero(x)
end

function addVariable!(problem::AbstractProblem{T}, sym::Symbol; func_ic::Function = default_ic)::Variable{T} where {T<:Real}
    # Check if symbol is not already used
    for var in problem.vars
        if var.sym == sym
            throw(DomainError(sym, "The symbol $(sym) is already in use by this system!"))
        end
    end

    # Current number of var
    n_var = length(problem.vars)

    # ID for new variable
    id = n_var + 1

    # Add variable to problem
    if problem isa SteadyProblem{T}
        var = Variable{T}(id, sym, zeros(T, problem.n_cps), zeros(T, 0), func_ic)
        push!(problem.vars, var)
    elseif problem isa TransientProblem{T}
        var = Variable{T}(id, sym, zeros(T, problem.n_cps), zeros(T, 0), func_ic)
        push!(problem.vars, var)
    end

    return var
end

function addAuxVariable!(problem::AbstractProblem{T}, sym::Symbol; func_ic::Function = default_ic)::AuxVariable{T} where {T<:Real}
    # Check if symbol is not already used
    for aux_var in problem.aux_vars
        if aux_var.sym == sym
            throw(DomainError(sym, "The symbol $(sym) is already in use by this system!"))
        end
    end
    for var in problem.vars
        if var.sym == sym
            throw(DomainError(sym, "The symbol $(sym) is already in use by this system!"))
        end
    end

    # Current number of auxiliary vars
    n_aux_var = length(problem.aux_vars)

    # ID for new variable
    id = n_aux_var + 1

    # Add aux variable to problem
    if problem isa SteadyProblem{T}
        aux_var = AuxVariable{T}(id, sym, zeros(T, problem.n_cps), zeros(T, 0), func_ic)
        push!(problem.aux_vars, aux_var)
    elseif problem isa TransientProblem{T}
        aux_var = AuxVariable{T}(id, sym, zeros(T, problem.n_cps), zeros(T, problem.n_cps), func_ic)
        push!(problem.aux_vars, aux_var)
    end

    return aux_var
end

function addKernel!(problem::AbstractProblem{T}, kernel::AbstractKernel{T}) where {T<:Real}
    push!(problem.kernels, kernel)
    return
end

function addAuxKernel!(problem::AbstractProblem{T}, aux_kernel::AbstractAuxKernel{T}) where {T<:Real}
    push!(problem.aux_kernels, aux_kernel)
    return
end

function initialize!(problem::AbstractProblem{T}) where {T<:Real}
    # Check if problem has variables
    if length(problem.vars) == 0
        throw(ErrorException("No variable has been added to the problem!"))
    end

    # Update number of dofs
    problem.n_dofs = problem.n_cps * length(problem.vars)

    # Check that each variable has an associated kernel
    for var in problem.vars
        exist = false
        for kernel in problem.kernels
            if kernel.u == var
                exist = true
            end
        end
        if ~exist
            throw(ErrorException("Variable $(string(var.sym)) has no associated kernel!"))
        end
    end

    # Check that each aux var has a single associated aux kernel
    for aux_var in problem.aux_vars
        exist = false
        i = 0
        for aux_kernel in problem.aux_kernels
            if aux_kernel.u == aux_var
                exist = true
                i += 1
            end
        end
        if ~exist
            throw(ErrorException("Auxiliary variable $(string(aux_var.sym)) has no associated auxiliary kernel!"))
        end
        if (i > 1)
            throw(ErrorException("Auxiliary variable $(string(aux_var.sym)) has more than one auxiliary kernel!"))
        end
    end

    return
end

function applyIC!(problem::AbstractProblem{T}) where {T<:Real}
    for var in problem.vars
        var.value = var.func_ic(problem.x)
        var.value_old = copy(var.value)
    end
    for aux_var in problem.aux_vars
        aux_var.value = aux_var.func_ic(problem.x)
        aux_var.value_old = copy(aux_var.value)
    end
    return
end