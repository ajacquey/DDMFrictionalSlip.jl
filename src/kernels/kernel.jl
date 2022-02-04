abstract type AbstractKernel{T<:Real} end

function computeCpResidual(kernel::AbstractKernel{T}, t::T, x::T, cp::Integer) where {T<:Real}
    throw(MethodError(computeCpResidual, kernel, t, x, cp))
end

function computeCpJacobian(kernel::AbstractKernel{T}, v::Variable{T}, t::T, x::T, cp::Integer)::T where {T<:Real}
    return zero(T)
end