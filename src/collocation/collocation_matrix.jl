
function collocationMatrix(mesh::Mesh{T}, order::Int64; μ::T=1.0) where T <: Real
    # Check order
    @assert order > -1
    @assert order < 3

    # Declare matrix
    E = zeros(T, mesh.n_elems * (order + 1), mesh.n_elems * (order + 1))

    collocationMatrix!(E, mesh, order; μ = μ)
    return E
end

function collocationMatrix!(E::AbstractMatrix{T}, mesh::Mesh{T}, order::Int64; μ::T=1.0) where T <: Real
    # Check size of matrix
    @assert size(E, 1) == (mesh.n_elems * (order + 1))
    @assert size(E, 2) == (mesh.n_elems * (order + 1))

    # Check order
    @assert order > -1
    @assert order < 3

    # Number of collocation points per element
    n_cp = order + 1

    # Start building matrix element-wise
    Threads.@threads for elem_i in mesh.elems
        # Index range - i
        idx_i = (elem_i.id - 1) * n_cp + 1: elem_i.id * n_cp
        for elem_j in mesh.elems 
            # Index range - j
            idx_j = (elem_j.id - 1) * n_cp + 1: elem_j.id * n_cp
            # Compute local collocation matrix
            localCollocationMatrix!(view(E, idx_i, idx_j), elem_i, elem_j, order, μ)
        end
    end

    return
end

function localCollocationMatrix!(E::AbstractMatrix{T}, elem_i::Elem, elem_j::Elem, order::Int64, μ::T) where T <: Real
    
    if order == 0
        localPWCMatrix!(E, elem_i, elem_j, μ)
    elseif order == 1
        localPWLCMatrix!(E, elem_i, elem_j, μ)
    elseif order == 2
        localPWQCMatrix!(E, elem_i, elem_j, μ)
    else
        throw(DomainError(order, "Unsupported order for the basis functions!"))
    end
    
    return
end