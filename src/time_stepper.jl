abstract type TimeStepper{T<:Real} end

mutable struct TimeSequence{T<:Real} <: TimeStepper{T}
    # Start time
    start_time::T

    # End time
    end_time::T

    # Time sequence
    time_seq::Vector{T}

    # Tolerance
    tol::T
end

function TimeSequence(time_seq::Vector{T}; start_time::T = time_seq[1], end_time::T = time_seq[end], tol::T = 1.0e-08) where {T<:Real}
    if (~issorted(time_seq))
        throw(DomainError(time_seq, "Time sequence need to be sorted!"))
    end

    # Check that starting time is either first elem of lower
    if (start_time > time_seq[1])
        throw(DomainError(start_time, "Start time should be smaller or equal than first time in sequence!"))
    elseif (start_time == time_seq[1])
        pushfirst!(time_seq, 0.0)
    end

    # Check that end time is either last elem or bigger
    if (end_time < time_seq[end])
        throw(DomainError(end_time, "End time should be bigger or equal than last time in sequence!"))
    elseif (end_time > time_seq[end])
        push!(time_seq, end_time)
    end

    return TimeSequence(start_time, end_time, time_seq, tol)
end
