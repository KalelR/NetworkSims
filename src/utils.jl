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
    δt = (t-t_s)/τs
    return δt * exp(1 - δt)
end

logrange(x1, x2; length) = (10^y for y in range(log10(x1), log10(x2), length=length))

function is_fp(att)
    isapprox(att[end], att[end-1], atol=1e-5)
end

function _rescale_to_01(vec::AbstractVector)
    mini, maxi = extrema(vec)
    return map(f -> (f .- mini) ./ (maxi .- mini), vec)
end

"""
Convert a linear index idx into a 2d one (i, j). Linear follows (in 0-notation) idx = numcols*i + j; so i = div(idx, numcols)
"""
function linear_to_2d_idx(idx, numrows, numcols)
    idxrow = div( (idx-1), numcols) + 1
    idxcol = ((idx-1) - (idxrow - 1) * numcols) + 1
    return idxrow, idxcol
end 

function num_rows_cols(numpanels, maxcols)
    if numpanels < 1 numpanels = 1 end
    numcols = numpanels > maxcols ? maxcols : numpanels
    numrows = ceil(Int, numpanels / numcols)
    return numrows, numcols
end