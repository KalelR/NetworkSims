function saver(res, pvals, ptypes, savemode; equals="_", connector="-", sort=true, ignores=(), digits=20, kwargs...)
    if savemode == "allseparated"
        for (res_type, res_val) in res
            resfilename = name_result(res_type, pvals, ptypes; equals, connector, sort, ignores, digits, kwargs...)
            safesave(resfilename, Dict(res_type=>res_val); kwargs...)
        end
    elseif savemode == "alltogether"
            resfilename = name_result("all", pvals, ptypes; equals, connector, sort, ignores, digits)
            println(resfilename)
            safesave(resfilename, res; kwargs...)
    else
        resfilename = name_result("alltogether", pvals, ptypes; fileformat = "jld2", equals, connector, sort, ignores, digits)
        @warn("Specified savemode $savemode is unavailable, saving res into an emergency file $resfilename")
        safesave(resfilename, res; kwargs...)
    end
    return nothing
end

"""
Duplicate code, but can't be bother to improve haha. Just make sure code is the same!
"""
function savenames(resnames, pvals, ptypes, savemode; equals="_", connector="-", sort=true, ignores=(), digits=20, kwargs...)
    resfilenames = []
    if savemode == "allseparated"
        for res_type in resnames
            resfilename = name_result(res_type, pvals, ptypes; equals, connector, sort, ignores, digits)
            push!(resfilenames, resfilename)
        end
    elseif savemode == "alltogether"
            resfilename = name_result("all", pvals, ptypes; equals, connector, sort, ignores, digits)
            push!(resfilenames, resfilename)
    else
        resfilename = name_result("alltogether", pvals, ptypes; fileformat = "jld2", equals, connector, sort, ignores, digits)
        push!(resfilenames, resfilename)
    end
    return resfilenames
end
