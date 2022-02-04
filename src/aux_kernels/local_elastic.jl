struct LocalElasticAuxKernel{T<:Real} <: AbstractAuxKernel{T}
    " The auxiliary variabe this auxiliary kernel refers to"
    u::AuxVariable{T}

    " The coupled variable describing the displacement discontinuity"
    v::Variable{T}

    " The elastic modulus"
    μ::T

    " The fault thickness"
    h::T

    " A symbol specifying when this auxiliary variable should be evaluated"
    execute_on::Symbol

    # Constructor
    # LocalElasticAuxKernel{T}(u::AuxVariable{T}, v::Variable{T}, μ::T, h::T; execute_on::Symbol=:non_linear) where {T<:Real} = new(u, v, μ , h, execute_on)
end

function computeCpValue(aux_kernel::LocalElasticAuxKernel, t::T, x::T, cp::Integer)::T where {T<:Real}
    return aux_kernel.u.value_old[cp] + aux_kernel.μ / aux_kernel.h * (aux_kernel.v.value[cp] - aux_kernel.v.value_old[cp])
end