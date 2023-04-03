function parse_inputs(args; printdictname=true)
    inputsdict = parse_args(args);

    dfull = getdict_all(inputsdict["dictname"], inputsdict["parameterlist_user"]; printdictname);
    pvals_all, ptypes = separatedicts(dfull);
    verify_inputs(pvals_all, ptypes);

    savemode = get_savemode(args);

    danalysis = getdict_all(inputsdict["analysisdictname"], inputsdict["analysisparameterlist_user"]; printdictname=false);
    # pvals_all, ptypes = separatedicts(danalysis);
    # verify_inputs(pvals_all, ptypes);
    # analysisnames = get_analysis(args);

    return pvals_all, ptypes, savemode, danalysis
end

function get_savemode(args)
    inputsdict = parse_args(args)
    return inputsdict["savemode"]
end

function get_analysis(args)
    inputsdict = parse_args(args)
    return inputsdict["analysismode"]
end

function parse_args(args)
    if length(args) == 0 @error("No parameter input given to run.jl! Pay attention, Kalel!") end
    dictname = ""; parameterlist_user = ""; savemode = ""; analysismode = "";
    if length(args) >= 1 dictname = args[1] end
    if length(args) >= 2 parameterlist_user = args[2] end
    if length(args) >= 3 savemode = args[3] end
    if length(args) >= 4 analysisdictname = args[4] end
    if length(args) >= 5 analysisparameterlist_user = args[5] end
    # return Dict("dictname"=>dictname, "parameterlist_user"=>parameterlist_user, "savemode"=>savemode, "analysisdictname"=>analysisdictname, "analy")
    return @strdict dictname parameterlist_user savemode analysisdictname analysisparameterlist_user
end

function getdict_all(args::Vector{String}; printdname=true)
    inputsdict = parse_args(args);
    return getdict_all(inputsdict["dictname"], inputsdict["parameterlist_user"]; printdname);
end

function getdict_all(dname, parameterlist_user=""; printdictname=true)
    if printdictname println("Running code for dict $dname in $parameterlist_user, with args $ARGS. Good luck!") end
    if parameterlist_user != ""
        include("$(datadir())/sims/inputs/$parameterlist_user");
        dfull = Base.invokelatest(getdict, dname)
        return dfull
    end
    parameterlists = ["parameters-example-devinterface.jl"]
    for parameterlist in parameterlists
        include("$(datadir())/sims/inputs/$parameterlist");
        dfull = Base.invokelatest(getdict, dname)
        return dfull
    end
    @error("No parameter found.")
end

""" From a Dict whose values are tuples, separate into two dicts whose values are each element of the tuple"""
function separatedicts(d)
    numvals = length(first(values(d)));
    ds = [Dict() for i=1:numvals]
    for (key, val) in d
        for i=1:numvals
            ds[i][key] = val[i]
        end
    end
    return ds
end

function verify_inputs(valsd, typesd)
    @assert sum("CP" .== values(typesd)) == 1 "No CP or more than 1 CP in inputs. Inspect thy inputs further, fool."
    nothing
end

"""
Check if a run was already made. Return true is it was already made, meaning that _all_ files that
would be saved already exist. Informs main() whether to run sim or continue.
"""
function check_existing_run(resnames, pvals, ptypes, savemode)
    resfilenames = savenames(resnames, pvals, ptypes, savemode)
    if all(isfile.(resfilenames)) return true end
    return false
end
