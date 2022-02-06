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
        function stress(t::Float64, x::Float64)
            λ = lambda_analytical_gs(T, 500)
            return T - erfc(λ * abs(x))
        end

        # Analytical solution
        function δ_analytical(x::Vector{Float64})
            x_a, δ_a, λ = injection_analytical_gs(T, length(x))
            itp = LinearInterpolation(x_a, δ_a) # create interpolation function
            return itp(x)
        end

        mesh = Mesh1D(-1.0, 1.0, 96)

        problem = Problem(mesh; order = 0)
        u = addVariable!(problem, :u)
        addKernel!(problem, FunctionKernel(u, stress))
        run!(problem; log = false)

        # Analytical solution
        δ_sol = δ_analytical(problem.x)
        # Error
        err = mean(abs.(u.value - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end
    @testset "T = 0.1, PWLC basis functions" begin
        T = 0.1

        # Imposed stress
        function stress(t::Float64, x::Float64)
            λ = lambda_analytical_gs(T, 500)
            return T - erfc(λ * abs(x))
        end

        # Analytical solution
        function δ_analytical(x::Vector{Float64})
            x_a, δ_a, λ = injection_analytical_gs(T, length(x))
            itp = LinearInterpolation(x_a, δ_a) # create interpolation function
            return itp(x)
        end

        mesh = Mesh1D(-1.0, 1.0, 96)

        problem = Problem(mesh; order = 1)
        u = addVariable!(problem, :u)
        addKernel!(problem, FunctionKernel(u, stress))
        run!(problem; log = false)
        # Analytical solution
        δ_sol = δ_analytical(problem.x)
        # Error
        err = mean(abs.(u.value - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end
    @testset "T = 0.1, PWQC basis functions" begin
        T = 0.1

        # Imposed stress
        function stress(t::Float64, x::Float64)
            λ = lambda_analytical_gs(T, 500)
            return T - erfc(λ * abs(x))
        end

        # Analytical solution
        function δ_analytical(x::Vector{Float64})
            x_a, δ_a, λ = injection_analytical_gs(T, length(x))
            itp = LinearInterpolation(x_a, δ_a) # create interpolation function
            return itp(x)
        end

        mesh = Mesh1D(-1.0, 1.0, 96)

        problem = Problem(mesh; order = 2)
        u = addVariable!(problem, :u)
        addKernel!(problem, FunctionKernel(u, stress))
        run!(problem; log = false)
        # Analytical solution
        δ_sol = δ_analytical(problem.x)
        # Error
        err = mean(abs.(u.value - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end
    @testset "T = 0.5, PWC basis functions" begin
        T = 0.5

        # Imposed stress
        function stress(t::Float64, x::Float64)
            λ = lambda_analytical_gs(T, 500)
            return T - erfc(λ * abs(x))
        end

        # Analytical solution
        function δ_analytical(x::Vector{Float64})
            x_a, δ_a, λ = injection_analytical_gs(T, length(x))
            itp = LinearInterpolation(x_a, δ_a) # create interpolation function
            return itp(x)
        end

        mesh = Mesh1D(-1.0, 1.0, 96)

        problem = Problem(mesh; order = 0)
        u = addVariable!(problem, :u)
        addKernel!(problem, FunctionKernel(u, stress))
        run!(problem; log = false)
        # Analytical solution
        δ_sol = δ_analytical(problem.x)
        # Error
        err = mean(abs.(u.value - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end
    @testset "T = 0.5, PWLC basis functions" begin
        T = 0.5

        # Imposed stress
        function stress(t::Float64, x::Float64)
            λ = lambda_analytical_gs(T, 500)
            return T - erfc(λ * abs(x))
        end

        # Analytical solution
        function δ_analytical(x::Vector{Float64})
            x_a, δ_a, λ = injection_analytical_gs(T, length(x))
            itp = LinearInterpolation(x_a, δ_a) # create interpolation function
            return itp(x)
        end

        mesh = Mesh1D(-1.0, 1.0, 96)

        problem = Problem(mesh; order = 1)
        u = addVariable!(problem, :u)
        addKernel!(problem, FunctionKernel(u, stress))
        run!(problem; log = false)
        # Analytical solution
        δ_sol = δ_analytical(problem.x)
        # Error
        err = mean(abs.(u.value - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end
    @testset "T = 0.5, PWQC basis functions" begin
        T = 0.5

        # Imposed stress
        function stress(t::Float64, x::Float64)
            λ = lambda_analytical_gs(T, 500)
            return T - erfc(λ * abs(x))
        end

        # Analytical solution
        function δ_analytical(x::Vector{Float64})
            x_a, δ_a, λ = injection_analytical_gs(T, length(x))
            itp = LinearInterpolation(x_a, δ_a) # create interpolation function
            return itp(x)
        end

        mesh = Mesh1D(-1.0, 1.0, 96)

        problem = Problem(mesh; order = 2)
        u = addVariable!(problem, :u)
        addKernel!(problem, FunctionKernel(u, stress))
        run!(problem; log = false)
        # Analytical solution
        δ_sol = δ_analytical(problem.x)
        # Error
        err = mean(abs.(u.value - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end
    @testset "T = 0.9, PWC basis functions" begin
        T = 0.9

        # Imposed stress
        function stress(t::Float64, x::Float64)
            λ = lambda_analytical_gs(T, 500)
            return T - erfc(λ * abs(x))
        end

        # Analytical solution
        function δ_analytical(x::Vector{Float64})
            x_a, δ_a, λ = injection_analytical_gs(T, length(x))
            itp = LinearInterpolation(x_a, δ_a) # create interpolation function
            return itp(x)
        end

        mesh = Mesh1D(-1.0, 1.0, 96)

        problem = Problem(mesh; order = 0)
        u = addVariable!(problem, :u)
        addKernel!(problem, FunctionKernel(u, stress))
        run!(problem; log = false)
        # Analytical solution
        δ_sol = δ_analytical(problem.x)
        # Error
        err = mean(abs.(u.value - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end
    @testset "T = 0.9, PWLC basis functions" begin
        T = 0.9

        # Imposed stress
        function stress(t::Float64, x::Float64)
            λ = lambda_analytical_gs(T, 500)
            return T - erfc(λ * abs(x))
        end

        # Analytical solution
        function δ_analytical(x::Vector{Float64})
            x_a, δ_a, λ = injection_analytical_gs(T, length(x))
            itp = LinearInterpolation(x_a, δ_a) # create interpolation function
            return itp(x)
        end

        mesh = Mesh1D(-1.0, 1.0, 96)

        problem = Problem(mesh; order = 1)
        u = addVariable!(problem, :u)
        addKernel!(problem, FunctionKernel(u, stress))
        run!(problem; log = false)
        # Analytical solution
        δ_sol = δ_analytical(problem.x)
        # Error
        err = mean(abs.(u.value - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end
    @testset "T = 0.9, PWQC basis functions" begin
        T = 0.9

        # Imposed stress
        function stress(t::Float64, x::Float64)
            λ = lambda_analytical_gs(T, 500)
            return T - erfc(λ * abs(x))
        end

        # Analytical solution
        function δ_analytical(x::Vector{Float64})
            x_a, δ_a, λ = injection_analytical_gs(T, length(x))
            itp = LinearInterpolation(x_a, δ_a) # create interpolation function
            return itp(x)
        end

        mesh = Mesh1D(-1.0, 1.0, 96)

        problem = Problem(mesh; order = 2)
        u = addVariable!(problem, :u)
        addKernel!(problem, FunctionKernel(u, stress))
        run!(problem; log = false)
        # Analytical solution
        δ_sol = δ_analytical(problem.x)
        # Error
        err = mean(abs.(u.value - δ_sol) ./ δ_sol)
        # Error less than 2%
        @test err < 0.02
    end
end

end