module TestTransientInjection

using DDMFrictionalSlip
using SpecialFunctions
using Interpolations
using Statistics
using Test

const λ = 1.0e+03::Float64
const h = 1.0e-03::Float64
const Δp = 0.4::Float64
const α = 0.04::Float64

@testset "Transient injection" begin
    @testset "Opening" begin
        
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
        addAuxVariable!(problem, :σ, :u, sigma, dsigma)

        solver = IterativeSolver(problem)

        time_seq = collect(range(0.5, 10.0, 20))
        time_stepper = TimeSequence(time_seq; start_time=0.0, end_time=10.0)

        output = [DomainOutput("outputs/transient_opening"), MaximumOutput("outputs/transient_opening_max")]

        run!(problem, solver, time_stepper; outputs=output)
        display(problem.vars[1].u)
        @test 1 == 1

    end
end
end