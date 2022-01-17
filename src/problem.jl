# abstract type Problem{T<:Real} end

# abstract type MechanicalProblem{T<:Real} <: Problem{T} end

abstract type AbstractProblem{T<:Real} end

mutable struct SteadyProblem{T<:Real} <: AbstractProblem{T}
    " The mesh of the problem"
    mesh::Mesh{T}

    " The order of the basis functions"
    order::Int64

    " Functions describing the stress residuals and jacobian"
    stress_residuals::Function
    stress_jacobian::Function

    " The total number of degrees of freedom"
    n_dofs::Int64

    " The collocation points coordinates"
    x::Vector{T}

    " The displacements"
    disp::Vector{T}

    " The stress"
    stress::Vector{T}
end

mutable struct TransientProblem{T<:Real} <: AbstractProblem{T}
    " The mesh of the problem"
    mesh::Mesh{T}

    " The order of the basis functions"
    order::Int64

    " Functions describing the stress residuals and jacobian"
    stress_residuals::Function
    stress_jacobian::Function

    " The total number of degrees of freedom"
    n_dofs::Int64

    " The collocation points coordinates"
    x::Vector{T}

    " The current time"
    time::T

    " The old time"
    time_old::T

    " The current time step"
    dt::T

    " The current time step"
    time_step::Int64

    " The displacements"
    disp::Vector{T}
    disp_old::Vector{T}

    " The stress"
    stress::Vector{T}
    stress_old::Vector{T}
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


function Problem(mesh::Mesh{T}, stress_res::Function, stress_jac::Function; order::Int64 = 0, transient::Bool = false) where {T<:Real}
    # Collocation point coordinates
    x = collocationPointsCoordinates(mesh, order)
    # Total number of dofs
    n_dofs = mesh.n_elems * (order + 1)

    # Displacements and stress
    disp = zeros(T, n_dofs)
    stress = zeros(T, n_dofs)

    if (transient)
        disp_old = zeros(T, n_dofs)
        stress_old = zeros(T, n_dofs)
        return TransientProblem{T}(mesh, order, stress_res, stress_jac, n_dofs, x, 0.0, 0.0, 0.0, 0, disp, disp_old, stress, stress_old)
    else
        return SteadyProblem{T}(mesh, order, stress_res, stress_jac, n_dofs, x, disp, stress)
    end
end

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
    for elem in problem.mesh.elems
        for i in 1:n_cp
            # Effective idx
            idx = (elem.id - 1) * n_cp + i

            # Displacements
            problem.disp_old[idx] = problem.disp[idx]

            # Stress
            problem.stress_old[idx] = problem.stress[idx]
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