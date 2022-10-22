"""
analysisdict is a dictionary mapping the analysis name to a dictionary containing possible parameters for the analysis
"""
function analyse_solution(sol, analysisdict, pvals, odeprob)
    res = Dict();
    if "sol" in keys(analysisdict) #could include here the saving dt
        res["sol"] = sol
    end
    if "spiketimes" in keys(analysisdict)
        res["spiketimes"] = odeprob.p.params_coup.coup_type.spiketime_detection.spiketimes
    end
    return res
end