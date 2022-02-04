function elementHalfLength(elem::Elem{T})::T where {T<:Real}
    if elem.type == :Edge1D
        return abs(elem.nodes[1].x - elem.nodes[2].x) / 2.0
    else
        throw(DomainError(elem.type, "Unsupported element type!"))
    end
end
