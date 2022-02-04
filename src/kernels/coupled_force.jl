struct CoupledForceKernel{T<:Real} <: AbstractKernel{T}
    " The variabe this kernel refers to"
    u::Variable{T}

    " The coupled auxiliary variable which rate describes the residuals"
    v::AuxVariable{T}
end

function computeCpResidual(kernel::CoupledForceKernel{T}, t::T, x::T, cp::Integer)::T where {T<:Real}
    return kernel.v.value[cp] - kernel.v.value_old[cp]
end

function computeCpJacobian(kernel::CoupledForceKernel{T}, v::Variable{T}, t::T, x::T, cp::Integer)::T where {T<:Real}
    return zero(T)
end