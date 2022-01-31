module TestTransientInjection

using DDMFrictionalSlip
using SpecialFunctions
using Interpolations
using Statistics
using Test

const λ = 1.0e+03::Float64
const μ = 1.0e+03::Float64
const h = 1.0e-03::Float64
const Δp = 0.4::Float64
const α = 0.04::Float64

@testset "Transient injection" begin
    @testset "Opening" begin
        function sigma_ic(x::Float64)
            return 1.0
        end

        function sigma(sigma_old::Float64, x::Float64, t::Float64, Δu::Float64)
            return sigma_old + λ / h * Δu
        end

        function dsigma(sigma_old::Float64, x::Float64, t::Float64, Δu::Float64)
            return λ / h
        end

        function p(p_old::Float64, x::Float64, t::Float64, Δu::Float64)
            return Δp * erfc(abs(x) / sqrt(α * t))
        end

        function dp(p_old::Float64, x::Float64, t::Float64, Δu::Float64)
            return 0.0
        end

        mesh = Mesh1D(-2.0, 2.0, 101)

        problem = Problem(mesh, order=2; transient=true, μ=λ)
        addVariable!(problem, :u)
        addAuxVariable!(problem, :p, :u, p, dp)
        addAuxVariable!(problem, :σ, :u, sigma, dsigma; func_ic = sigma_ic)

        solver = IterativeSolver(problem)

        time_seq = collect(range(0.5, stop=10.0, length=20))
        time_stepper = TimeSequence(time_seq; start_time=0.0, end_time=10.0)

        output = [DomainOutput("outputs/transient_opening"), MaximumOutput("outputs/transient_opening_max")]

        run!(problem, solver, time_stepper; log = false, outputs = output)
        # Analytical solution
        u_sol = -h / λ * Δp * erfc.(abs.(problem.x) / sqrt(α * problem.time))
        # Error
        err = mean(abs.(problem.vars[1].u - u_sol))
        # Error less than 1.0e-08
        @test err < 1.0e-08
    end
    # @testset "Uncoupled slip/opening" begin
    #     function sigma_ic(x::Float64)
    #         return 1.0
    #     end

    #     function tau_ic(x::Float64)
    #         return 0.4
    #     end

    #     function sigma(sigma_old::Float64, x::Float64, t::Float64, Δϵ::Float64)
    #         return sigma_old + λ / h * Δϵ
    #     end

    #     function dsigma(sigma_old::Float64, x::Float64, t::Float64, Δϵ::Float64)
    #         return λ / h
    #     end

    #     function tau(tau_old::Float64, x::Float64, t::Float64, Δδ::Float64)
    #         return tau_old + μ / h * Δδ
    #     end

    #     function dtau(tau_old::Float64, x::Float64, t::Float64, Δδ::Float64)
    #         return μ / h
    #     end

    #     function p(p_old::Float64, x::Float64, t::Float64, Δϵ::Float64)
    #         return Δp * erfc(abs(x) / sqrt(α * t))
    #     end

    #     function dp(p_old::Float64, x::Float64, t::Float64, Δϵ::Float64)
    #         return 0.0
    #     end

    #     mesh = Mesh1D(-2.0, 2.0, 101)

    #     problem = Problem(mesh, order=2; transient=true, μ=λ)
    #     addVariable!(problem, :ϵ)
    #     addVariable!(problem, :δ)
    #     addAuxVariable!(problem, :p, :ϵ, p, dp)
    #     addAuxVariable!(problem, :σ, :ϵ, sigma, dsigma; func_ic = sigma_ic)
    #     addAuxVariable!(problem, :τ, :δ, tau, dtau; func_ic = tau_ic)

    #     solver = IterativeSolver(problem)

    #     time_seq = collect(range(0.5, stop=10.0, length=20))
    #     time_stepper = TimeSequence(time_seq; start_time=0.0, end_time=10.0)

    #     output = [DomainOutput("outputs/transient_uncoupled"), MaximumOutput("outputs/transient_uncoupled_max")]

    #     run!(problem, solver, time_stepper; log = false, outputs = output)
    #     # Analytical solution
    #     u_sol = -h / λ * Δp * erfc.(abs.(problem.x) / sqrt(α * problem.time))
    #     # Error
    #     err = mean(abs.(problem.vars[1].u - u_sol))
    #     # Error less than 1.0e-08
    #     @test err < 1.0e-08
    # end
end
end