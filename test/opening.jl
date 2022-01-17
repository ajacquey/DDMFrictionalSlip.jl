module TestOpening

using DDMFrictionalSlip
using Statistics
using Test

# Imposed stress
function stress_res(x::Float64, u::Float64)
    return -1.0
end

function stress_jac(x::Float64, u::Float64)
    return 0.0
end

# Analytical solution
function δ_analytical(x::Float64)
    return sqrt(1.0 - x^2)
end

@testset "Crack opening" begin
    @testset "PWC basis functions" begin
        mesh = Mesh1D(-1.0, 1.0, 96)

        problem = Problem(mesh, stress_res, stress_jac; order = 0)
        solver = IterativeSolver(problem)
        run!(problem, solver; log = false)

        # Analytical solution
        δ_sol = δ_analytical.(problem.x)
        # Error
        err = mean(abs.(problem.disp - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end

    @testset "PWLC basis functions" begin
        mesh = Mesh1D(-1.0, 1.0, 48)

        problem = Problem(mesh, stress_res, stress_jac; order = 1)
        solver = IterativeSolver(problem)
        run!(problem, solver; log = false)

        # Analytical solution
        δ_sol = δ_analytical.(problem.x)
        # Error
        err = mean(abs.(problem.disp - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end

    @testset "PWQC basis functions" begin
        mesh = Mesh1D(-1.0, 1.0, 32)

        problem = Problem(mesh, stress_res, stress_jac; order = 2)
        solver = IterativeSolver(problem)
        run!(problem, solver; log = false)

        # Analytical solution
        δ_sol = δ_analytical.(problem.x)
        # Error
        err = mean(abs.(problem.disp - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end
end

end