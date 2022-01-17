abstract type AbstractOutput end

mutable struct DomainOutput <: AbstractOutput
    file::String
    initialized::Bool
    # p_func::Union{Function, nothing}
end

mutable struct MaximumOutput <: AbstractOutput
    file::String
    initialized::Bool
end

function DomainOutput(filename::String)
    file = checkFilename(filename)
    # return DomainOuput(file, p_func)
    return DomainOutput(file, false)
end

function MaximumOutput(filename::String)
    file = checkFilename(filename)
    return MaximumOutput(file, false)
end

function checkFilename(filename::String)
    # # Check if filename does not contain an extension
    # if (~endswith(filename, ".csv"))
    #     filename = string(filename, ".csv")
    # end
    # Check if filename contains a folder and if this folder exists
    if (contains(filename, '/'))
        output_dir = string()
        k_end = findlast(isequal('/'), filename)
        if (startswith(filename, '/')) # absolute path
            output_dir = filename[1:k_end]
        else # relative path
            output_dir = string(@__DIR__, "/", filename[1:k_end])
        end
        # Create folder if it doesn't exist
        if (!isdir(output_dir))
            mkpath(output_dir)
        end
        return string(output_dir, filename[k_end + 1:end])
    else
        return string(@__DIR__, filename)
    end
end

function initialize!(output::DomainOutput)
    output.file = checkFilename(output.file)
    output.initialized = true
    return
end

function initialize!(output::MaximumOutput)
    output.file = checkFilename(output.file)
    output.initialized = true
    open(string(output.file, ".csv"); write=true) do f
        write(f, "time,tau,delta,delta_dot\n")
    end
    return
end

function execute!(output::DomainOutput, problem::TransientProblem{T}) where {T<:Real}
    open(string(output.file, "_$(string(problem.time_step, pad=4)).csv"); write=true) do f
        # if isnothing(output.p_func)
        write(f, "x,tau,delta,delta_dot\n") # header
        writedlm(f, hcat(problem.x, problem.stress, problem.disp, (problem.disp .- problem.disp_old) / problem.dt), ',')
        # else
        #     write(f, "x,tau,p,delta,delta_dot\n") # header
        #     writedlm(f, hcat(problem.x, problem.stress, output.p_func.(problem.x, problem.t), problem.disp, (problem.disp .- problem.disp_old) / problem.dt), ',')
        # end
    end
    return
end

function execute!(output::MaximumOutput, problem::TransientProblem{T}) where {T<:Real}
    open(string(output.file, ".csv"); append=true) do f
        writedlm(f, [problem.time maximum(problem.stress) maximum(problem.disp) maximum(problem.disp .- problem.disp_old) / problem.dt ], ',')
    end
    return
end

function outputResults!(outputs::Vector{AbstractOutput}, problem::AbstractProblem)
    for output in outputs
        execute!(output, problem)
    end
    return 
end

function initializeOutputs!(outputs::Vector{AbstractOutput})
    for output in outputs
        initialize!(output)
    end
    return
end