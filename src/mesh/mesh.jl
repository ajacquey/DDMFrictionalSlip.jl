abstract type Mesh{T <: Real} end

struct Mesh1D{T <: Real} <: Mesh{T}
    dim::Int64
    n_elems::Int64
    n_nodes::Int64
    elems::Vector{Elem{T}}
    nodes::Vector{Node{T}}
end

"""
Custom `show` function for `Mesh{T}` that prints some information.
"""
function Base.show(io::IO, M::Mesh{T} where T <: Real)
    println("Mesh information:")
    println("  dimension: $(M.dim)")
	println("    n_elems: $(M.n_elems)")
    println("    n_nodes: $(M.n_nodes)")
end

function Mesh1D(x::Vector{T}) where T <: Real
    dim = 1
    n = length(x) - 1
    if (~issorted(x))
        throw(ErrorException("Array of coordinates need to be sorted!"))
    end
    elems = Vector{Elem{T}}(undef,n)
    elem_type = :Edge1D
    nodes = Vector{Node{T}}(undef, n + 1)
    # First point
    nodes[1] = Node{T}(1, x[1], 0.0)
    for i in 1:n
        centroid = Point{T}((x[i] + x[i+1]) / 2.0, 0.0)
        nodes[i+1] = Node{T}(i+1, x[i+1], 0.0)
        elems[i] = Edge1D{T}(elem_type, i, centroid, [nodes[i], nodes[i+1]])
    end
    return Mesh1D(dim, n, n + 1, elems, nodes)
end

function Mesh1D(xmin::T, xmax::T, n::Int64) where T <: Real
    if (xmin >= xmax)
        throw(ErrorException("xmax should be stricly bigger than xmin!"))
    end
    if (n <= 0)
        throw(ErrorException("Number of elements should be bigger than zero!"))
    end
    x = collect(range(xmin, xmax, length=n+1))
    return Mesh1D(x::Vector{T})
end