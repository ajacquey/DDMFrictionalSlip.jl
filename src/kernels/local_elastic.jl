struct LocalElasticKernel{T<:Real} <: AbstractKernel{T}
    " The variabe this kernel refers to"
    u::Variable{T}

    " The elastic modulus"
    μ::T

    " The fault thickness"
    h::T
end

function computeCpResidual(kernel::LocalElasticKernel{T}, t::T, x::T, cp::Integer)::T where {T<:Real}
    return kernel.μ / kernel.h * (kernel.u.value[cp] - kernel.u.value_old[cp])
end

function computeCpJacobian(kernel::LocalElasticKernel{T}, v::Variable{T}, t::T, x::T, cp::Integer)::T where {T<:Real}
    if v == kernel.u
        return kernel.μ / kernel.h
    else
        return zero(T)
    end
end