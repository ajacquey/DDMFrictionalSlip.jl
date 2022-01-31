abstract type AbstractVariable{T<:Real} end

mutable struct Variable{T<:Real} <: AbstractVariable{T}
    " The variable ID"
    id::Int64

    " The variable symbol"
    sym::Symbol

    " The value of the variable"
    u::Vector{T}

    " The old value of the variable"
    u_old::Vector{T}

    # " The value of the time derivative of the variable"
    # u_dot::Vector{T}

    # " The value of the stress associated to this variable"
    # s::Vector{T}

    # " The old value of the stress associated to this variable"
    # s_old::Vector{T}

    # " The value of the time derivative of the stress associated variable"
    # s_dot::Vector{T}

    " The function describing the initial conditions"
    func_ic::Function
end

mutable struct AuxVariable{T<:Real} <: AbstractVariable{T}
    " The auxiliary variable ID"
    id::Int64

    " The auxiliary variable symbol"
    sym::Symbol

    " The associated variable"
    ass_var::Symbol

    " The value of the auxiliary variable"
    u::Vector{T}

    " The old value of the auxiliary variable"
    u_old::Vector{T}

    " The function to update the auxiliary variable and the residuals"
    func::Function

    " The function to update the jacobian"
    dfunc::Function

    # " The value of the time derivative of the variable"
    # u_dot::Vector{T}

    # " The value of the stress associated to this variable"
    # s::Vector{T}

    # " The old value of the stress associated to this variable"
    # s_old::Vector{T}

    " The function describing the initial conditions"
    func_ic::Function
end
