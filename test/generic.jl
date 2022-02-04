module TestGeneric

using DDMFrictionalSlip
using Statistics
using Test

# Imposed stress
function stress(t::Float64, x::Float64)
    return 1.5 * (2.0 * x^2 - 1)
end

# Analytical solution
function δ_analytical(x::Float64)
    return (1.0 - x^2)^1.5
end

@testset "Generic example" begin
    @testset "PWC basis functions" begin
        mesh = Mesh1D(-1.0, 1.0, 96)

        problem = Problem(mesh; order = 0)
        solver = IterativeSolver(problem)
        u = addVariable!(problem, :u)
        addKernel!(problem, FunctionKernel(u, stress))
        # σ = addAuxVariable!(problem, :stress, :u, stress_res, stress_jac)
        run!(problem, solver; log = false)

        # Analytical solution
        δ_sol = δ_analytical.(problem.x)
        # Error
        err = mean(abs.(u.value - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end

    @testset "PWLC basis functions" begin
        mesh = Mesh1D(-1.0, 1.0, 48)

        problem = Problem(mesh; order = 1)
        solver = IterativeSolver(problem)
        u = addVariable!(problem, :u)
        addKernel!(problem, FunctionKernel(u, stress))
        # σ = addAuxVariable!(problem, :stress, :u, stress_res, stress_jac)
        run!(problem, solver, log = false)

        # Analytical solution
        δ_sol = δ_analytical.(problem.x)
        # Error
        err = mean(abs.(u.value - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end

    @testset "PWQC basis functions" begin
        mesh = Mesh1D(-1.0, 1.0, 32)

        problem = Problem(mesh; order = 2)
        solver = IterativeSolver(problem)
        u = addVariable!(problem, :u)
        addKernel!(problem, FunctionKernel(u, stress))
        # σ = addAuxVariable!(problem, :stress, :u, stress_res, stress_jac)
        run!(problem, solver; log = false)

        # Analytical solution
        δ_sol = δ_analytical.(problem.x)
        # Error
        err = mean(abs.(u.value - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end
end

end