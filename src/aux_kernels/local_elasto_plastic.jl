struct LocalElastoPlasticAuxKernel{T<:Real} <: AbstractAuxKernel{T}
    " The auxiliary variabe this auxiliary kernel refers to"
    u::AuxVariable{T}

    " The coupled variable describing the displacement discontinuity"
    v::Variable{T}

    " The shear elastic modulus"
    μ::T

    " The fault thickness"
    h::T

    " The friction coefficient"
    f::T

    " The normal effective stress auxiliary variable"
    σ::AuxVariable{T}

    " A symbol specifying when this auxiliary variable should be evaluated"
    execute_on::Symbol

    # Constructor
    # LocalElasticAuxKernel{T}(u::AuxVariable{T}, v::Variable{T}, μ::T, h::T; execute_on::Symbol=:non_linear) where {T<:Real} = new(u, v, μ , h, execute_on)
end

function computeCpValue(aux_kernel::LocalElastoPlasticAuxKernel, t::T, x::T, cp::Integer)::T where {T<:Real}
    # Trial stresses
    τ_tr = aux_kernel.u.value_old[cp] + aux_kernel.μ / aux_kernel.h * (aux_kernel.v.value[cp] - aux_kernel.v.value_old[cp])

    # Yield function
    if ((τ_tr - aux_kernel.f * aux_kernel.σ.value[cp]) >= 0.0)
        return aux_kernel.f * aux_kernel.σ.value[cp]
    else
        return τ_tr
    end
end