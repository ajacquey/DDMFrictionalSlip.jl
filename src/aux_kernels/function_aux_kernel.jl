struct FunctionAuxKernel{T<:Real} <: AbstractAuxKernel{T}
    " The auxiliary variabe this auxiliary kernel refers to"
    u::AuxVariable{T}

    " The function to evaluate this auxiliary variable"
    func::Function

    " A symbol specifying when this auxiliary variable should be evaluated"
    execute_on::Symbol

    # # Constructor
    # FunctionAuxKernel{T}(u::AuxVariable{T}, func::Function; execute_on::Symbol=:non_linear) where {T<:Real} = new(u, func, execute_on)
end

function computeCpValue(aux_kernel::FunctionAuxKernel{T}, t::T, x::T, cp::Integer)::T where {T<:Real}
    return aux_kernel.func(t, x)
end