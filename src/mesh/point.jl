abstract type AbstractPoint end

struct Point{T<:Real} <: AbstractPoint
    x::T
    y::T
end

struct Node{T<:Real} <: AbstractPoint
    id::Int64
    x::T
    y::T
end

"""
Custom `show` function for `Point{T}` that prints some information.
"""
function Base.show(io::IO, P::Point{T} where {T<:Real})
    println("Point information:")
    println("  coordinates: ($(P.x), $(P.y))")
    println()
end

"""
Custom `show` function for `Node{T}` that prints some information.
"""
function Base.show(io::IO, N::Node{T} where {T<:Real})
    println("Node information:")
    println("  coordinates: ($(N.x), $(N.y))")
    println("           ID: $(N.id)")
    println()
end