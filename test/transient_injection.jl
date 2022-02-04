module TestTransientInjection

using DDMFrictionalSlip
using SpecialFunctions
using Interpolations
using Statistics
using Test

const λ = 6.6667e+02::Float64
const μ = 6.6667e+02::Float64
const h = 1.0e-03::Float64
const f = 0.5::Float64
const Δp = 0.4::Float64
const α = 0.04::Float64

function sigma_ic(x::Vector{Float64})
    return ones(size(x))
end

function tau_ic(x::Vector{Float64})
    return 0.4 * ones(size(x))
end

function p_func(t::Float64, x::Float64)
    return Δp * erfc(abs(x) / sqrt(α * t))
end

@testset "Transient injection" begin
    @testset "Opening" begin
        mesh = Mesh1D(-2.0, 2.0, 101)

        problem = Problem(mesh, order = 2; transient = true, μ = λ)
        u = addVariable!(problem, :u)
        p = addAuxVariable!(problem, :p)
        σ = addAuxVariable!(problem, :σ; func_ic = sigma_ic)
        addKernel!(problem, LocalElasticKernel(u, λ, h))
        addKernel!(problem, CoupledForceKernel(u, p))
        addAuxKernel!(problem, FunctionAuxKernel(p, p_func, :time_step_begin))
        addAuxKernel!(problem, LocalElasticAuxKernel(σ, u, λ, h, :time_step_end))

        solver = IterativeSolver(problem)

        time_seq = collect(range(0.5, stop = 10.0, length = 20))
        time_stepper = TimeSequence(time_seq; start_time = 0.0, end_time = 10.0)

        output = [DomainOutput("outputs/transient_opening"), MaximumOutput("outputs/transient_opening_max")]

        run!(problem, solver, time_stepper; log = false, outputs = output)

        # Analytical solution
        u_sol = -h / λ * Δp * erfc.(abs.(problem.x) / sqrt(α * problem.time))
        # Error
        err = mean(abs.(u.value - u_sol))
        # Error less than 1.0e-08
        @test err < 1.0e-08
    end
    @testset "Uncoupled slip/opening" begin
        mesh = Mesh1D(-2.0, 2.0, 101)

        problem = Problem(mesh, order = 2; transient = true, μ = μ)
        ϵ = addVariable!(problem, :ϵ)
        δ = addVariable!(problem, :δ)
        p = addAuxVariable!(problem, :p)
        σ = addAuxVariable!(problem, :σ; func_ic = sigma_ic)
        τ = addAuxVariable!(problem, :τ; func_ic = tau_ic)
        addKernel!(problem, LocalElasticKernel(ϵ, λ, h))
        addKernel!(problem, CoupledForceKernel(ϵ, p))
        addKernel!(problem, LocalElasticKernel(δ, μ, h))
        addAuxKernel!(problem, FunctionAuxKernel(p, p_func, :time_step_begin))
        addAuxKernel!(problem, LocalElasticAuxKernel(σ, ϵ, λ, h, :time_step_end))
        addAuxKernel!(problem, LocalElasticAuxKernel(τ, δ, μ, h, :time_step_end))
        solver = IterativeSolver(problem)

        time_seq = collect(range(0.5, stop = 10.0, length = 20))
        time_stepper = TimeSequence(time_seq; start_time = 0.0, end_time = 10.0)

        output = [DomainOutput("outputs/transient_uncoupled"), MaximumOutput("outputs/transient_uncoupled_max")]

        run!(problem, solver, time_stepper; log = false, outputs = output)

        # Analytical solution
        ϵ_sol = -h / λ * Δp * erfc.(abs.(problem.x) / sqrt(α * problem.time))
        # Error
        err = mean(abs.(ϵ.value - ϵ_sol))
        # Error less than 1.0e-08
        @test err < 1.0e-08
    end
    @testset "Frictional slip" begin
        mesh = Mesh1D(-5.0, 5.0, 101)

        problem = Problem(mesh, order = 2; transient = true, μ = μ)
        ϵ = addVariable!(problem, :ϵ)
        δ = addVariable!(problem, :δ)
        p = addAuxVariable!(problem, :p)
        σ = addAuxVariable!(problem, :σ; func_ic = sigma_ic)
        τ = addAuxVariable!(problem, :τ; func_ic = tau_ic)
        addKernel!(problem, LocalElasticKernel(ϵ, λ, h))
        addKernel!(problem, CoupledForceKernel(ϵ, p))
        addKernel!(problem, LocalElastoPlasticKernel(δ, ϵ, λ, μ, h, f, σ, τ))
        addAuxKernel!(problem, FunctionAuxKernel(p, p_func, :time_step_begin))
        addAuxKernel!(problem, LocalElasticAuxKernel(σ, ϵ, λ, h, :time_step_end))
        addAuxKernel!(problem, LocalElastoPlasticAuxKernel(τ, δ, μ, h, f, σ, :time_step_end))
        solver = IterativeSolver(problem)

        time_seq = collect(range(0.5, stop = 10.0, length = 20))
        time_stepper = TimeSequence(time_seq; start_time = 0.0, end_time = 20.0)

        output = [DomainOutput("outputs/transient_slip"), MaximumOutput("outputs/transient_slip_max")]

        run!(problem, solver, time_stepper; log = false, outputs = output)

        # Analytical solution
        ϵ_sol = -h / λ * Δp * erfc.(abs.(problem.x) / sqrt(α * problem.time))
        # Error
        err = mean(abs.(ϵ.value - ϵ_sol))
        # Error less than 1.0e-08
        @test err < 1.0e-08
    end
end
end