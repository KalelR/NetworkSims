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


function spiketimes(sol, pvals)
    @unpack Vth = pvals
    N = size(sol, 2)
    sts = [Float64[] for i=1:N]
    for i=1:N
        Vs = view(sol, 1, i, :)
        dist_to_th = Vs .- Vth #spike when sign changes from -1 to +1; so when diff = 1 - (-1) = 2
        signs = sign.(dist_to_th)
        idxs_spikes = findall(x->x==2, diff(signs))
        sts_i = sol.t[idxs_spikes]
        sts[i] = sts_i
    end
    return sts
end