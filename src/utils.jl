# """
# See tests for performance comparison. This is faster for sparse matrices, but Tullio is faster for dense.
# """
@inbounds function diagmul!(result, A, B)
    for (ind, r, c) in zip(eachindex(result), eachrow(A), eachcol(B))
        result[ind] = dot(r, c)
    end
    nothing
end

# """
# Need to receive transpose of A
# """
@inbounds function diagmul_t!(result, At, B)
    for (ind, r, c) in zip(eachindex(result), eachcol(At), eachcol(B))
        result[ind] = dot(r, c)
    end
    nothing
end

# "for dense matrices"
@inbounds function diagmul_d!(result, A, B)
    @tullio result[i] = A[i,j]*B[j,i]
    nothing
end

function alphafunction(t, t_s, τs::Float64)
    if t < t_s @warn "OPAAAA" end
    δt = (t-t_s)/τs
    return δt * exp(1 - δt)
end