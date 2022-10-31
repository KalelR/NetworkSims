# """
# See tests for performance comparison. This allocates but is nevertheless faster for small and big N!
# """
@inbounds function diagmul!(result, A, B)
    result .= dot.(eachrow(A), eachcol(B))
    nothing
end


function alphafunction(t, t_s, τs::Float64)
    δt = (t-t_s)/τs
    return δt * exp(1 - δt)
end