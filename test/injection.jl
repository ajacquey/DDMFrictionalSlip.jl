module TestInjection

using DDMFrictionalSlip
using SpecialFunctions
using Interpolations
using Statistics
using Test

include("injection_utils.jl")

@testset "Injection scaled" begin
    @testset "T = 0.1, PWC basis functions" begin
        T = 0.1

        # Imposed stress
        function stress_res(x::Float64, Δu::Float64)
            λ = lambda_analytical_gs(T, 500)
            return T - erfc(λ * abs(x))

        end

        function stress_jac(x::Float64, Δu::Float64)
            return 0.0
        end

        # Analytical solution
        function δ_analytical(x::Vector{Float64})
            x_a, δ_a, λ = injection_analytical_gs(T, length(x))
            itp = LinearInterpolation(x_a, δ_a) # create interpolation function
            return itp(x)
        end

        mesh = Mesh1D(-1.0, 1.0, 96)

        problem = Problem(mesh, stress_res, stress_jac; order=0)
        solver = IterativeSolver(problem)
        run!(problem, solver; log = false)
        # Analytical solution
        δ_sol = δ_analytical(problem.x)
        # Error
        err = mean(abs.(problem.disp - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end
    @testset "T = 0.1, PWLC basis functions" begin
        T = 0.1

        # Imposed stress
        function stress_res(x::Float64, Δu::Float64)
            λ = lambda_analytical_gs(T, 500)
            return T - erfc(λ * abs(x))

        end

        function stress_jac(x::Float64, Δu::Float64)
            return 0.0
        end

        # Analytical solution
        function δ_analytical(x::Vector{Float64})
            x_a, δ_a, λ = injection_analytical_gs(T, length(x))
            itp = LinearInterpolation(x_a, δ_a) # create interpolation function
            return itp(x)
        end

        mesh = Mesh1D(-1.0, 1.0, 96)

        problem = Problem(mesh, stress_res, stress_jac; order=1)
        solver = IterativeSolver(problem)
        run!(problem, solver; log = false)
        # Analytical solution
        δ_sol = δ_analytical(problem.x)
        # Error
        err = mean(abs.(problem.disp - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end
    @testset "T = 0.1, PWQC basis functions" begin
        T = 0.1

        # Imposed stress
        function stress_res(x::Float64, Δu::Float64)
            λ = lambda_analytical_gs(T, 500)
            return T - erfc(λ * abs(x))

        end

        function stress_jac(x::Float64, Δu::Float64)
            return 0.0
        end

        # Analytical solution
        function δ_analytical(x::Vector{Float64})
            x_a, δ_a, λ = injection_analytical_gs(T, length(x))
            itp = LinearInterpolation(x_a, δ_a) # create interpolation function
            return itp(x)
        end

        mesh = Mesh1D(-1.0, 1.0, 96)

        problem = Problem(mesh, stress_res, stress_jac; order=2)
        solver = IterativeSolver(problem)
        run!(problem, solver; log = false)
        # Analytical solution
        δ_sol = δ_analytical(problem.x)
        # Error
        err = mean(abs.(problem.disp - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end
    @testset "T = 0.5, PWC basis functions" begin
        T = 0.5

        # Imposed stress
        function stress_res(x::Float64, Δu::Float64)
            λ = lambda_analytical_gs(T, 500)
            return T - erfc(λ * abs(x))

        end

        function stress_jac(x::Float64, Δu::Float64)
            return 0.0
        end

        # Analytical solution
        function δ_analytical(x::Vector{Float64})
            x_a, δ_a, λ = injection_analytical_gs(T, length(x))
            itp = LinearInterpolation(x_a, δ_a) # create interpolation function
            return itp(x)
        end

        mesh = Mesh1D(-1.0, 1.0, 96)

        problem = Problem(mesh, stress_res, stress_jac; order=0)
        solver = IterativeSolver(problem)
        run!(problem, solver; log = false)
        # Analytical solution
        δ_sol = δ_analytical(problem.x)
        # Error
        err = mean(abs.(problem.disp - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end
    @testset "T = 0.5, PWLC basis functions" begin
        T = 0.5

        # Imposed stress
        function stress_res(x::Float64, Δu::Float64)
            λ = lambda_analytical_gs(T, 500)
            return T - erfc(λ * abs(x))

        end

        function stress_jac(x::Float64, Δu::Float64)
            return 0.0
        end

        # Analytical solution
        function δ_analytical(x::Vector{Float64})
            x_a, δ_a, λ = injection_analytical_gs(T, length(x))
            itp = LinearInterpolation(x_a, δ_a) # create interpolation function
            return itp(x)
        end

        mesh = Mesh1D(-1.0, 1.0, 96)

        problem = Problem(mesh, stress_res, stress_jac; order=1)
        solver = IterativeSolver(problem)
        run!(problem, solver; log = false)
        # Analytical solution
        δ_sol = δ_analytical(problem.x)
        # Error
        err = mean(abs.(problem.disp - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end
    @testset "T = 0.5, PWQC basis functions" begin
        T = 0.5

        # Imposed stress
        function stress_res(x::Float64, Δu::Float64)
            λ = lambda_analytical_gs(T, 500)
            return T - erfc(λ * abs(x))

        end

        function stress_jac(x::Float64, Δu::Float64)
            return 0.0
        end

        # Analytical solution
        function δ_analytical(x::Vector{Float64})
            x_a, δ_a, λ = injection_analytical_gs(T, length(x))
            itp = LinearInterpolation(x_a, δ_a) # create interpolation function
            return itp(x)
        end

        mesh = Mesh1D(-1.0, 1.0, 96)

        problem = Problem(mesh, stress_res, stress_jac; order=2)
        solver = IterativeSolver(problem)
        run!(problem, solver; log = false)
        # Analytical solution
        δ_sol = δ_analytical(problem.x)
        # Error
        err = mean(abs.(problem.disp - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end
    @testset "T = 0.9, PWC basis functions" begin
        T = 0.9

        # Imposed stress
        function stress_res(x::Float64, Δu::Float64)
            λ = lambda_analytical_gs(T, 500)
            return T - erfc(λ * abs(x))

        end

        function stress_jac(x::Float64, Δu::Float64)
            return 0.0
        end

        # Analytical solution
        function δ_analytical(x::Vector{Float64})
            x_a, δ_a, λ = injection_analytical_gs(T, length(x))
            itp = LinearInterpolation(x_a, δ_a) # create interpolation function
            return itp(x)
        end

        mesh = Mesh1D(-1.0, 1.0, 96)

        problem = Problem(mesh, stress_res, stress_jac; order=0)
        solver = IterativeSolver(problem)
        run!(problem, solver; log = false)
        # Analytical solution
        δ_sol = δ_analytical(problem.x)
        # Error
        err = mean(abs.(problem.disp - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end
    @testset "T = 0.9, PWLC basis functions" begin
        T = 0.9

        # Imposed stress
        function stress_res(x::Float64, Δu::Float64)
            λ = lambda_analytical_gs(T, 500)
            return T - erfc(λ * abs(x))

        end

        function stress_jac(x::Float64, Δu::Float64)
            return 0.0
        end

        # Analytical solution
        function δ_analytical(x::Vector{Float64})
            x_a, δ_a, λ = injection_analytical_gs(T, length(x))
            itp = LinearInterpolation(x_a, δ_a) # create interpolation function
            return itp(x)
        end

        mesh = Mesh1D(-1.0, 1.0, 96)

        problem = Problem(mesh, stress_res, stress_jac; order=1)
        solver = IterativeSolver(problem)
        run!(problem, solver; log = false)
        # Analytical solution
        δ_sol = δ_analytical(problem.x)
        # Error
        err = mean(abs.(problem.disp - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end
    @testset "T = 0.9, PWQC basis functions" begin
        T = 0.9

        # Imposed stress
        function stress_res(x::Float64, Δu::Float64)
            λ = lambda_analytical_gs(T, 500)
            return T - erfc(λ * abs(x))

        end

        function stress_jac(x::Float64, Δu::Float64)
            return 0.0
        end

        # Analytical solution
        function δ_analytical(x::Vector{Float64})
            x_a, δ_a, λ = injection_analytical_gs(T, length(x))
            itp = LinearInterpolation(x_a, δ_a) # create interpolation function
            return itp(x)
        end

        mesh = Mesh1D(-1.0, 1.0, 96)

        problem = Problem(mesh, stress_res, stress_jac; order=2)
        solver = IterativeSolver(problem)
        run!(problem, solver; log = false)
        # Analytical solution
        δ_sol = δ_analytical(problem.x)
        # Error
        err = mean(abs.(problem.disp - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end
end

end
#     @testset "PWC basis functions" begin
#         mesh = Mesh1D(-1.0, 1.0, 96)

#         problem = MechanicalProblem(mesh, stress; order=0)
#         solver = IterativeSolver(problem)
#         x, δ = run!(problem, solver)

#         # Analytical solution
#         δ_sol = δ_analytical(x)
#         # Error
#         err = mean(abs.(δ - δ_sol) ./ δ_sol)
#         # Error less than 2%
#         @test err < 0.02
#     end

#     @testset "PWLC basis functions" begin
#         mesh = Mesh1D(-1.0, 1.0, 48)

#         problem = MechanicalProblem(mesh, stress; order=1)
#         solver = IterativeSolver(problem)
#         x, δ = run!(problem, solver)

#         # Analytical solution
#         δ_sol = δ_analytical(x)
#         # Error
#         err = mean(abs.(δ - δ_sol) ./ δ_sol)
#         # Error less than 2%
#         @test err < 0.02
#     end

#     @testset "PWQC basis functions" begin
#         mesh = Mesh1D(-1.0, 1.0, 32)

#         problem = MechanicalProblem(mesh, stress; order=2)
#         solver = IterativeSolver(problem)
#         x, δ = run!(problem, solver)

#         # Analytical solution
#         δ_sol = δ_analytical(x)
#         # Error
#         err = mean(abs.(δ - δ_sol) ./ δ_sol)
#         # Error less than 2%
#         @test err < 0.02
#     end
# end

# end