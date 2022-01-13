function localPWCMatrix!(E::AbstractMatrix{T}, elem_i::Elem{T}, elem_j::Elem{T}, μ::T) where T <: Real
    n_cp = 1
    # Element j half length
    aj = elementHalfLength(elem_j)
    # Collocation points element i
    c_points = elemCollocationPoints(elem_i, 0)
    
    for i in 1:n_cp
        uij = c_points[i].x - elem_j.centroid.x
        for j in 1:n_cp
            E[i,j] = μ * 2.0 * aj / (π * (uij^2 - aj^2))
        end
    end

    return
end

function localPWLCMatrix!(E::AbstractMatrix{T}, elem_i::Elem{T}, elem_j::Elem{T}, μ::T) where T <: Real
    n_cp = 2
    # Element j half length
    aj = elementHalfLength(elem_j)
    # Collocation points element i
    c_points = elemCollocationPoints(elem_i, 1)
    
    for i in 1:n_cp
        uij = c_points[i].x - elem_j.centroid.x
        for j in 1:n_cp
            E[i,j] = μ * aj / (π * (uij^2 - aj^2))
            E[i,j] += μ * 2.0 * (j - 1.5) * sqrt(2.0) / π * (uij / (uij^2 - aj^2) + 1.0 / (2.0 * aj) * log(abs(uij - aj) / abs(uij + aj)))
        end
    end

    return
end

function localPWQCMatrix!(E::AbstractMatrix{T}, elem_i::Elem{T}, elem_j::Elem{T}, μ::T) where T <: Real
    n_cp = 3
    # Element j half length
    aj = elementHalfLength(elem_j)
    # Collocation points element i
    c_points = elemCollocationPoints(elem_i, 2)
    
    for i in 1:n_cp
        uij = c_points[i].x - elem_j.centroid.x
        for j in 1:n_cp
            E[i,j] = μ * 2.0 * (abs(j - 2) - 1.0 / 3.0) *  aj / (π * (uij^2 - aj^2))
            E[i,j] += μ * (j - 2) * 2.0 / (π * sqrt(3.0)) * (uij / (uij^2 - aj^2) + 1.0 / (2.0 * aj) * log(abs(uij - aj) / abs(uij + aj)))
            E[i,j] += μ * 4.0 / π * (abs(j - 2) - 2.0 / 3.0) * (uij / aj^2 * log(abs(uij - aj) / abs(uij + aj)) + 2.0 / aj)
        end
    end

    return
end


