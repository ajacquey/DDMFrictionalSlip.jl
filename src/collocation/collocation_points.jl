function collocationPointsCoordinates(mesh::Mesh{T}, order::Int64) where {T<:Real}
    # Number of collocation points
    n_cp = order + 1
    # Declare results
    c_points = zeros(T, mesh.n_elems * n_cp)


    Threads.@threads for elem in mesh.elems
        c_points[(elem.id-1)*n_cp+1:elem.id*n_cp] .= elemCollocationPointsCoordinates(elem, order)
    end

    return c_points
end

function elemCollocationPointsCoordinates(elem::Elem{T}, order::Int64) where {T<:Real}
    if elem.type == :Edge1D
        return elemCollocationPoints1DCoordinates(elem, order)
    else
        throw(DomainError(elem.type, "Unsupported element type!"))
    end
end

function elemCollocationPoints(elem::Elem{T}, order::Int64) where {T<:Real}
    if elem.type == :Edge1D
        return elemCollocationPoints1D(elem, order)
    else
        throw(DomainError(elem.type, "Unsupported element type!"))
    end
end

function elemCollocationPoints1DCoordinates(elem::Elem{T}, order::Int64) where {T<:Real}
    # Number of collocation points
    n_cp = order + 1
    # Declare vector for c_points
    c_points = zeros(T, n_cp)

    # Element half length
    aj = elementHalfLength(elem)

    if order == 0
        for i in 1:n_cp
            c_points[i] = elem.centroid.x
        end
    elseif order == 1
        for i in 1:n_cp
            c_points[i] = elem.centroid.x + (-1.0)^i / sqrt(2.0) * aj
        end
    elseif order == 2
        for i in 1:n_cp
            c_points[i] = elem.centroid.x + (i - 2) * sqrt(3.0) / 2.0 * aj
        end
    else
        throw(DomainError(order, "Unsupported order for the basis functions!"))
    end
    return c_points
end

function elemCollocationPoints1D(elem::Elem{T}, order::Int64) where {T<:Real}
    # Number of collocation points
    n_cp = order + 1
    # Declare vector for c_points
    c_points = Vector{Point{T}}(undef, n_cp)

    # Element half length
    aj = elementHalfLength(elem)

    if order == 0
        for i in 1:n_cp
            c_points[i] = Point(elem.centroid.x, 0.0)
        end
    elseif order == 1
        for i in 1:n_cp
            c_points[i] = Point(elem.centroid.x + (-1.0)^i / sqrt(2.0) * aj, 0.0)
        end
    elseif order == 2
        for i in 1:n_cp
            c_points[i] = Point(elem.centroid.x + (i - 2) * sqrt(3.0) / 2.0 * aj, 0.0)
        end
    else
        throw(DomainError(order, "Unsupported order for the basis functions!"))
    end
    return c_points
end