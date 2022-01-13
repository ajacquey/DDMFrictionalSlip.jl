abstract type Elem{T<:Real} end

struct Edge1D{T<:Real} <: Elem{T}
    type::Symbol
    id::Int64
    centroid::Point{T}
    nodes::Vector{Node{T}}
end

"""
Custom `show` function for `Elem{T}` that prints some information.
"""
function Base.show(io::IO, E::Elem{T} where {T<:Real})
    println("Element information:")
    println("       type: $(E.type)")
    println("         ID: $(E.id)")
    println("  # ofnodes: $(length(E.nodes))")
    println()
end