function name_dir(pvals, ptypes; kwargs...)
    types=collect(values(ptypes))
    uniqueSPs = unique(types[findall(x->x[1:2]=="SP", types)])
    dirpath = ""
    sortstrvec!(uniqueSPs, 3)
    for SP in uniqueSPs
        params_type_SP = findall(x->values(x)==SP, ptypes)
        vals_type_SP = filter(x->x[1] in params_type_SP, pvals)
        SP_num = parse(Int, SP[3:end])
        dir = savename(vals_type_SP; kwargs...)
        if dirpath == "" dirpath = dir
        else dirpath = "$(dirpath)/$dir" end
    end
return dirpath
end

function sortstrvec!(v, idx)
    nums = [parse(Int, el[idx:end]) for el in v]
    idxs_sort = sortperm(nums)
    v .= v[idxs_sort]
    nothing
end

function name_file(savevariable, pvals, ptypes; kwargs...)
    cp = findfirst(x->values(x)=="CP", ptypes)
    vals_type_cp = filter(x->x[1] == cp, pvals)
    savename(savevariable, vals_type_cp; kwargs...)
end

"""
SPs go into directory names, with lowest number higher on top of the dir hierarchy; cp goes into variable name; any other type is naturally ignored in the name (but they are of course available for makesims in pvals)
"""
function name_result(savevariable, pvals, ptypes; fileformat="jld2", equals="_", connector="-", sort=true, ignores=(), digits=20, allowedtypes=(Real, String, Symbol, Vector), kwargs...)
    dirpath = name_dir(pvals, ptypes; equals, connector, sort, ignores, digits, allowedtypes, kwargs...);
    filename = name_file(savevariable, pvals, ptypes;  equals, connector, sort, ignores, digits, allowedtypes, kwargs...)
    fullfilename = "$(datadir())/sims/results/$(dirpath)/$(filename)"
    # fullfilename = "$(datadir())/sims/results/$(dirpath)/$(filename).$(fileformat)"
end

function plotname(plottitle, _pvals, ptypes; fileformat="png", equals="_", connector="-", sort=true, ignores=(), digits=20, allowedtypes=(Real, String, Symbol, Vector), padleft=4, padright=4, kwargs...)
    pvals = deepcopy(_pvals)
    cp = findfirst(x->values(x)=="CP", ptypes)
    pvals[cp] = pvals[cp] isa String ? pvals[cp] : number_to_padded_string(pvals[cp], padleft, padright)
    dirpath = name_dir(pvals, ptypes; equals, connector, sort, ignores, digits, allowedtypes, kwargs...);
    filename = name_file(plottitle, pvals, ptypes;  equals, connector, sort, ignores, digits, allowedtypes, kwargs...)
    fullfilename = "$(plotsdir())/sims/results/$(dirpath)/$(filename).$(fileformat)"
end

using Formatting
function number_to_padded_string(num::Number, padleft=4, padright=4)
    # num_s = string(num) 
    num_s = sprintf1("%.$(padright)f", num) 
    idx = findfirst(".", num_s)[1] #idx of the decimal point
    right = num_s[idx+1:end]
    left = num_s[1:idx-1] 
    left_s = lpad(left, 4, "0")
    right_s = rpad(right, 4, "0")
    num_s_pad = left_s * right_s
end

function number_to_padded_string(num::AbstractArray, padleft=4, padright=4)
    number_to_padded_string(num[1], padleft, padright)
end