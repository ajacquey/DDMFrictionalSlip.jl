abstract type Problem{T <: Real} end

struct MechanicalProblem{T <: Real} <: Problem{T}
    " The mesh of the problem"
    mesh::Mesh{T}

    " The order of the basis functions"
    order::Int64

    " Function describing the stress"
    stress_function::Function

    " The total number of degrees of freedom"
    n_dofs::Int64

    " The collocation points coordinates"
    x::Vector{T}
end

mutable struct HydroMechanicalProblem{T <: Real} <: Problem{T}
    " The mesh of the problem"
    mesh::Mesh{T}

    " The order of the basis functions"
    order::Int64

    " Functiom describing the fluid pressure"
    pressure_function::Function

    " The total number of degrees of freedom"
    n_dofs::Int64

    " The collocation points coordinates"
    x::Vector{T}

    " The displacements"
    disp::Vector{Vector{T}}
    disp_old::Vector{Vector{T}}

    " The stresses"
    stress::Vector{Vector{T}}
    stress_old::Vector{Vector{T}}

    " The fluid pressure"
    p::Vector{T}
    p_old::Vector{T}

    " Physical parameters"
    par::NamedTuple{}
end


function MechanicalProblem(mesh::Mesh{T}, stress_function::Function; order::Int64 = 0) where T  <: Real
    # Collocation point coordinates
    x = collocationPointsCoordinates(mesh, order)
    # Total number of dofs
    n_dofs = mesh.n_elems * (order + 1)

    return MechanicalProblem(mesh, order, stress_function, n_dofs, x)
end

function HydroMechanicalProblem(mesh::Mesh{T}, pressure_function::Function, param::NamedTuple{}; order::Int64 = 0) where T <: Real
    # Collocation point coordinates
    x = collocationPointsCoordinates(mesh, order)
    # Number of dofs per variable
    n_dofs_per_var = mesh.n_elems * (order + 1)
    # Total number of dofs
    n_dofs = 2 * n_dofs_per_var
    # Displacement Discontinuities
    disp = Vector{Vector{T}}(undef, 2)
    for i in 1:length(disp)
        disp[i] = zeros(T, n_dofs_per_var)
    end
    # Stresses
    stress = Vector{Vector{T}}(undef, 2)
    for i in 1:length(stress)
        stress[i] = zeros(T, n_dofs_per_var)
    end
    # Fluid pressure
    p = zeros(T, n_dofs_per_var)

    return HydroMechanicalProblem(mesh, order, pressure_function, n_dofs, x, disp, disp, stress, stress, p, p, param)
end

function reinit!(problem::HydroMechanicalProblem{T}) where T <: Real
    # Move current values as old
    for elem in problem.mesh.elems
        for i in 1:n_cp
            # Effective idx
            idx = (elem.id - 1) * n_cp + i

            # Displacements
            problem.disp_old[1][idx] = problem.disp[1][idx]
            problem.disp_old[2][idx] = problem.disp[2][idx]

            # Stress
            problem.stress_old[1][idx] = problem.stress[1][idx]
            problem.stress_old[2][idx] = problem.stress[2][idx]

            # Fluid pressure
            problem.p_old[idx] = problem.p[idx]
        end
    end
    return
end