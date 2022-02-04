struct LocalElastoPlasticKernel{T<:Real} <: AbstractKernel{T}
    " The variabe this kernel refers to"
    u::Variable{T}

    " The coupled variable describing the opening"
    ϵ::Variable{T}

    " The normal elastic modulus"
    λ::T

    " The shear elastic modulus"
    μ::T

    " The fault thickness"
    h::T

    " The friction coefficient"
    f::T

    " The normal effective stress auxiliary variable"
    σ::AuxVariable{T}

    " The shear stress auxiliary variable"
    τ::AuxVariable{T}
end

function computeCpResidual(kernel::LocalElastoPlasticKernel{T}, t::T, x::T, cp::Integer)::T where {T<:Real}
    # Trial stresses
    τ_tr = kernel.τ.value_old[cp] + kernel.μ / kernel.h * (kernel.u.value[cp] - kernel.u.value_old[cp])
    σ_tr = kernel.σ.value_old[cp] + kernel.λ / kernel.h * (kernel.ϵ.value[cp] - kernel.ϵ.value_old[cp])

    # Yield function
    if ((τ_tr - kernel.f * σ_tr) >= 0.0)
        # println("Yield: ",  τ_tr - kernel.f * σ_tr)
        # println("Old shear stress: ", kernel.τ.value_old[cp])
        # println("New stress: ", kernel.f * σ_tr)
        # println("plastic!")
        return kernel.f * σ_tr - kernel.τ.value_old[cp]
    else
        return kernel.μ / kernel.h * (kernel.u.value[cp] - kernel.u.value_old[cp])
    end
end

function computeCpJacobian(kernel::LocalElastoPlasticKernel{T}, v::Variable{T}, t::T, x::T, cp::Integer)::T where {T<:Real}
    if (v == kernel.u || v == kernel.ϵ)
        # Trial stresses
        τ_tr = kernel.τ.value_old[cp] + kernel.μ / kernel.h * (kernel.u.value[cp] - kernel.u.value_old[cp])
        σ_tr = kernel.σ.value_old[cp] + kernel.λ / kernel.h * (kernel.ϵ.value[cp] - kernel.ϵ.value_old[cp])

        # Yield function
        if ((τ_tr - kernel.f * σ_tr) >= 0.0)
            if v == kernel.u
                return zero(T)
            elseif v == kernel.ϵ
                return kernel.f * kernel.λ / kernel.h
            end
        else
            if v == kernel.u
                return kernel.μ / kernel.h
            elseif v == kernel.ϵ
                return zero(T)
            end
        end
    else
        return zero(T)
    end
end