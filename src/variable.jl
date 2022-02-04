abstract type AbstractVariable{T<:Real} end

mutable struct Variable{T<:Real} <: AbstractVariable{T}
    " The variable ID"
    id::Int64

    " The variable symbol"
    sym::Symbol

    " The value of the variable"
    value::Vector{T}

    " The old value of the variable"
    value_old::Vector{T}

    " The function describing the initial conditions"
    func_ic::Function
end

mutable struct AuxVariable{T<:Real} <: AbstractVariable{T}
    " The auxiliary variable ID"
    id::Int64

    " The auxiliary variable symbol"
    sym::Symbol

    " The value of the auxiliary variable"
    value::Vector{T}

    " The old value of the auxiliary variable"
    value_old::Vector{T}

    " The function describing the initial conditions"
    func_ic::Function
end