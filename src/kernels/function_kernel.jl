struct FunctionKernel{T<:Real} <: AbstractKernel{T}
    " The variabe this kernel refers to"
    u::Variable{T}

    " The function describing the kernel residuals"
    func::Function
end

function computeCpResidual(kernel::FunctionKernel{T}, t::T, x::T, cp::Integer)::T where {T<:Real}
    return kernel.func(t, x)
end

function computeCpJacobian(kernel::FunctionKernel{T}, v::Variable{T}, t::T, x::T, cp::Integer)::T where {T<:Real}
    return zero(T)
end