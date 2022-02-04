abstract type AbstractAuxKernel{T<:Real} end

function computeCpValue(aux_kernel::AbstractAuxKernel{T}, t::T, x::T, cp::Integer) where {T<:Real}
    throw(MethodError(computeCpValue, aux_kernel, t, x, cp))
end