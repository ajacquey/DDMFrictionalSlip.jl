# abstract type Problem{T<:Real} end

# abstract type MechanicalProblem{T<:Real} <: Problem{T} end

abstract type AbstractProblem{T<:Real} end

mutable struct SteadyProblem{T<:Real} <: AbstractProblem{T}
    " The mesh of the problem"
    mesh::Mesh{T}

    " The order of the basis functions"
    order::Int64

    " The elastic modulus"
    μ::Float64

    # " Functions describing the stress residuals and jacobian"
    # stress_residuals::Function
    # stress_jacobian::Function

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

    " A boolean to check if system is initialized"
    initialized::Bool

    # " The displacements"
    # disp::Vector{T}

    # " The stress"
    # stress::Vector{T}
end

mutable struct TransientProblem{T<:Real} <: AbstractProblem{T}
    " The mesh of the problem"
    mesh::Mesh{T}

    " The order of the basis functions"
    order::Int64

    " The elastic modulus"
    μ::Float64

    # " Functions describing the stress residuals and jacobian"
    # stress_residuals::Function
    # stress_jacobian::Function

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

    # " The displacements"
    # disp::Vector{T}
    # disp_old::Vector{T}

    # " The stress"
    # stress::Vector{T}
    # stress_old::Vector{T}
end

# mutable struct HydroMechanicalProblem{T<:Real} <: Problem{T}
#     " The mesh of the problem"
#     mesh::Mesh{T}

#     " The order of the basis functions"
#     order::Int64

#     " Functiom describing the fluid pressure"
#     pressure_function::Function

#     " The total number of degrees of freedom"
#     n_dofs::Int64

#     " The collocation points coordinates"
#     x::Vector{T}

#     " The displacements"
#     disp::Vector{Vector{T}}
#     disp_old::Vector{Vector{T}}

#     " The stresses"
#     stress::Vector{Vector{T}}
#     stress_old::Vector{Vector{T}}

#     " The fluid pressure"
#     p::Vector{T}
#     p_old::Vector{T}

#     " Physical parameters"
#     par::NamedTuple{}
# end

function Problem(mesh::Mesh{T}; μ::T = 1.0, order::Int64 = 0, transient::Bool = false) where {T<:Real}
    # Collocation point coordinates
    x = collocationPointsCoordinates(mesh, order)

    # Total number of cps
    n_cps = mesh.n_elems * (order + 1)

    if (transient)
        return TransientProblem{T}(mesh, order, μ, n_cps, 0, x, Vector{Variable}(undef, 0), Vector{AuxVariable}(undef, 0), 0.0, 0.0, 0.0, 0, false)
    else
        return SteadyProblem{T}(mesh, order, μ, n_cps, 0, x, Vector{Variable}(undef, 0), Vector{AuxVariable}(undef, 0), false)
    end
end

# function Problem(mesh::Mesh{T}, stress_res::Function, stress_jac::Function; μ::T = 1.0, order::Int64 = 0, transient::Bool = false) where {T<:Real}
#     # Collocation point coordinates
#     x = collocationPointsCoordinates(mesh, order)
#     # Total number of dofs
#     n_dofs = mesh.n_elems * (order + 1)

#     # Displacements and stress
#     disp = zeros(T, n_dofs)
#     stress = zeros(T, n_dofs)

#     if (transient)
#         disp_old = zeros(T, n_dofs)
#         stress_old = zeros(T, n_dofs)
#         return TransientProblem{T}(mesh, order, μ, stress_res, stress_jac, n_dofs, x, 0.0, 0.0, 0.0, 0, disp, disp_old, stress, stress_old)
#     else
#         return SteadyProblem{T}(mesh, order, μ, stress_res, stress_jac, n_dofs, x, disp, stress)
#     end
# end

# function HydroMechanicalProblem(mesh::Mesh{T}, pressure_function::Function, param::NamedTuple{}; order::Int64 = 0) where {T<:Real}
#     # Collocation point coordinates
#     x = collocationPointsCoordinates(mesh, order)
#     # Number of dofs per variable
#     n_dofs_per_var = mesh.n_elems * (order + 1)
#     # Total number of dofs
#     n_dofs = 2 * n_dofs_per_var
#     # Displacement Discontinuities
#     disp = Vector{Vector{T}}(undef, 2)
#     for i in 1:length(disp)
#         disp[i] = zeros(T, n_dofs_per_var)
#     end
#     # Stresses
#     stress = Vector{Vector{T}}(undef, 2)
#     # for i in 1:length(stress)
#     #     stress[i] = zeros(T, n_dofs_per_var)
#     # end
#     stress[1] = ones(T, n_dofs_per_var)
#     stress[2] = 12.0 * ones(T, n_dofs_per_var)
#     # Fluid pressure
#     p = zeros(T, n_dofs_per_var)

#     return HydroMechanicalProblem{T}(mesh, order, pressure_function, n_dofs, x, disp, copy(disp), stress, copy(stress), p, copy(p), param)
# end

function reinit!(problem::TransientProblem{T}, time_stepper::TimeStepper{T}) where {T<:Real}
    n_cp = problem.order + 1
    # Move current values as old
    Threads.@threads for elem in problem.mesh.elems
        for i in 1:n_cp
            for var in problem.vars
                # Effective idx
                idx = (elem.id - 1) * n_cp + i

                # Variables
                var.u_old[idx] = var.u[idx]

                # Aux variables
                for aux_var in problem.aux_vars
                    if aux_var.ass_var == var.sym
                        aux_var.u_old[idx] = aux_var.u[idx]
                    end
                end
            end
        end
    end

    # Time
    problem.time = time_stepper.time_seq[problem.time_step]
    problem.dt = problem.time - problem.time_old
    return
end

# function reinit!(problem::HydroMechanicalProblem{T}) where {T<:Real}
#     n_cp = problem.order + 1
#     # Move current values as old
#     for elem in problem.mesh.elems
#         for i in 1:n_cp
#             # Effective idx
#             idx = (elem.id - 1) * n_cp + i

#             # Displacements
#             problem.disp_old[1][idx] = problem.disp[1][idx]
#             problem.disp_old[2][idx] = problem.disp[2][idx]

#             # Stress
#             problem.stress_old[1][idx] = problem.stress[1][idx]
#             problem.stress_old[2][idx] = problem.stress[2][idx]

#             # Fluid pressure
#             problem.p_old[idx] = problem.p[idx]
#         end
#     end
#     return
# end

function addVariable!(problem::AbstractProblem{T}, sym::Symbol) where {T<:Real}
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
        push!(problem.vars, Variable{T}(id, sym, zeros(T, problem.n_cps), zeros(T, 0)))
    elseif problem isa TransientProblem{T}
        push!(problem.vars, Variable{T}(id, sym, zeros(T, problem.n_cps), zeros(T, problem.n_cps)))
    end
    
    return
end

function addAuxVariable!(problem::AbstractProblem{T}, sym::Symbol, ass_var::Symbol, func::Function, dfunc::Function) where {T<:Real}
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
        push!(problem.aux_vars, AuxVariable{T}(id, sym, ass_var, zeros(T, problem.n_cps), zeros(T, 0), func, dfunc))
    elseif problem isa TransientProblem{T}
        push!(problem.aux_vars, AuxVariable{T}(id, sym, ass_var, zeros(T, problem.n_cps), zeros(T, problem.n_cps), func, dfunc))
    end

    return
end

function initialize!(problem::AbstractProblem{T}) where {T<:Real}
    # Check if problem has variables
    if length(problem.vars) == 0
        throw(ErrorException("No variable has been added to the problem!"))
    end

    # Update number of dofs
    problem.n_dofs = problem.n_cps * length(problem.vars)

    # Check that each aux var has an existing associated var
    for aux_var in problem.aux_vars
        exist = false
        for var in problem.vars
            if var.sym == aux_var.ass_var
                exist = true
            end
        end
        if ~exist
            throw(ErrorException("Variable $(string(aux_var.ass_var)) does not exist!"))
        end
    end

    return
end