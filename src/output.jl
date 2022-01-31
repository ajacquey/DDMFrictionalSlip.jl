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
            output_dir = string(pwd(), "/", filename[1:k_end])
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

function initialize!(output::DomainOutput, problem::AbstractProblem{T}) where {T<:Real}
    output.file = checkFilename(output.file)
    output.initialized = true
    return
end

function initialize!(output::MaximumOutput, problem::TransientProblem{T}) where {T<:Real}
    output.file = checkFilename(output.file)
    output.initialized = true
    open(string(output.file, ".csv"); write=true) do f
        header = "time"
        for aux_var in problem.aux_vars
            header = string(header, ",", string(aux_var.sym))
        end
        for var in problem.vars
            header = string(header, ",", string(var.sym))
        end
        header = string(header, "\n")
        write(f, header)
    end
    return
end

function execute!(output::DomainOutput, problem::TransientProblem{T}) where {T<:Real}
    open(string(output.file, "_$(string(problem.time_step, pad=4)).csv"); write=true) do f
        header = "x"
        data = problem.x
        for aux_var in problem.aux_vars
            header = string(header, ",", string(aux_var.sym))
            data = hcat(data, aux_var.u)
        end
        for var in problem.vars
            header = string(header, ",", string(var.sym))
            data = hcat(data, var.u)
        end
        header = string(header, "\n")
        ##### TODO #####
        # Add rate of var
        ################
        write(f, header) # write header
        writedlm(f, data, ',') # write data
    end
    return
end

function execute!(output::MaximumOutput, problem::TransientProblem{T}) where {T<:Real}
    open(string(output.file, ".csv"); append=true) do f
        data = problem.time
        for aux_var in problem.aux_vars
            data = hcat(data, maximum(aux_var.u))
        end
        for var in problem.vars
            data = hcat(data, maximum(var.u))
        end
        ##### TODO #####
        # Add rate of var
        ################
        writedlm(f, data, ',') # write data
    end
    return
end

function outputResults!(outputs::Vector{AbstractOutput}, problem::AbstractProblem{T}) where {T<:Real}
    for output in outputs
        execute!(output, problem)
    end
    return 
end

function initializeOutputs!(outputs::Vector{AbstractOutput}, problem::TransientProblem{T}) where {T<:Real}
    for output in outputs
        initialize!(output, problem)
    end
    return
end