abstract type AbstractOutput end

mutable struct DomainOutput <: AbstractOutput
    file::String
    initialized::Bool
end

mutable struct MaximumOutput <: AbstractOutput
    file::String
    initialized::Bool
end

function DomainOutput(filename::String)
    file = checkFilename(filename)
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
        return string(output_dir, filename[k_end+1:end])
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
    open(string(output.file, ".csv"); write = true) do f
        header = "time"
        # Vars and aux vars values
        for aux_var in problem.aux_vars
            header = string(header, ",", string(aux_var.sym))
        end
        for var in problem.vars
            header = string(header, ",", string(var.sym))
        end
        # Vars and aux vars rate values
        for aux_var in problem.aux_vars
            header = string(header, ",", string(aux_var.sym, "_dot"))
        end
        for var in problem.vars
            header = string(header, ",", string(var.sym, "_dot"))
        end
        header = string(header, "\n")
        write(f, header)
    end
    return
end

function execute!(output::DomainOutput, problem::TransientProblem{T}) where {T<:Real}
    open(string(output.file, "_$(string(problem.time_step, pad=4)).csv"); write = true) do f
        header = "x"
        data = problem.x
        # Vars and aux vars values
        for aux_var in problem.aux_vars
            header = string(header, ",", string(aux_var.sym))
            data = hcat(data, aux_var.value)
        end
        for var in problem.vars
            header = string(header, ",", string(var.sym))
            data = hcat(data, var.value)
        end
        # Vars and aux vars rate values
        for aux_var in problem.aux_vars
            header = string(header, ",", string(aux_var.sym, "_dot"))
            data = hcat(data, (aux_var.value - aux_var.value_old) / problem.dt)
        end
        for var in problem.vars
            header = string(header, ",", string(var.sym, "_dot"))
            data = hcat(data, (var.value - var.value_old) / problem.dt)
        end
        header = string(header, "\n")

        write(f, header) # write header
        writedlm(f, data, ',') # write data
    end
    return
end

function execute!(output::MaximumOutput, problem::TransientProblem{T}) where {T<:Real}
    open(string(output.file, ".csv"); append = true) do f
        data = problem.time
        # Vars and aux vars values
        for aux_var in problem.aux_vars
            data = hcat(data, maximum(aux_var.value))
        end
        for var in problem.vars
            data = hcat(data, maximum(var.value))
        end
        # Vars and aux vars values
        for aux_var in problem.aux_vars
            data = hcat(data, maximum((aux_var.value - aux_var.value_old) / problem.dt))
        end
        for var in problem.vars
            data = hcat(data, maximum((var.value - var.value_old) / problem.dt))
        end
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