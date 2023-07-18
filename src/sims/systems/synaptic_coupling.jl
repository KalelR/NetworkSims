



mutable struct SpikeTimeDetection
    spiketimes :: Matrix{Union{Float64, Nothing}} #{N, numts}
    current_spikes_idxs :: Vector{Int64}
    Vth :: Float64
    numstoredspikes :: Int64
end

function reinit_synaptic!(odeproblem)
    # p.spiketimedetection.spiketimes = Matrix{Union{Nothing, Float64}}(undef, (N, numstoredspiketimes))
    p = odeproblem.p
    fill!(p.spiketimedetection.spiketimes, nothing) 
    fill!(p.spiketimedetection.current_spikes_idxs, 1)
    fill!(p.Icoup, 0.0)
    fill!(p.ΔV, 0.0)
    fill!(p.AΔV, 0.0)
    fill!(p.Ct, 0.0)
    fill!(p.Csum, 0.0)
    nothing
end





"""
Current that each unit is *sending*, not receiving. This is different than what I used to
do. 
On my notes, C[i,j] = current that j sends to i. Here, I'm actually calculating
transpose(C); In this case, C[i,j] is the current that i sends to j Conversely, also
C[i,j] is the current that j receives from i

Objects *receiver_types: vector of size N, each element is a dictionary mapping each
connection type (e.g. :I, :E) to the neighborhood of i belonging to that type. This helps
to quickly assign the connection sum values to each connection type. Goes into the C=Ct
matrix.
"""
function syn_sum!(p, t)
    C = view(p.Ct, :, :); spiketimes = p.spiketimedetection.spiketimes; Csum = p.Csum; τs_vals = p.τs_vals; receivers_types = p.receivers_types
    spiketimes_rows = eachrow(spiketimes)
    @inbounds for (i, spiketimes_i) in enumerate(spiketimes_rows)
        fill!(Csum, 0.0)
        if all(isempty.(values(receivers_types[i]))) continue end 
        for t_s in spiketimes_i
            Csum[1] += alphafunction(t, t_s, τs_vals[:E])
            Csum[2] += alphafunction(t, t_s, τs_vals[:I])
        end
        C[i, receivers_types[i][:E]] .= Csum[1]; #transpose of actual C matrix; element [i,j] contains current that i sends to j
        C[i, receivers_types[i][:I]] .= Csum[2]; #transpose of actual C matrix; element [i,j] contains current that i sends to j
    end
    return nothing
end


function nextindex(idx, maxidx)
    if idx == maxidx return 1
    else return idx+1 end
    nothing
end

function condition_spike!(out, u, t, integrator)
    Vth = integrator.p.spiketimedetection.Vth
    N = integrator.p.N
    Vs =  view(u, 1:N)
    @. out = Vs - Vth
    nothing
end

function affect_spike!(integrator, idx)
    spiketimes = integrator.p.spiketimedetection.spiketimes
    current_spikes_idxs = integrator.p.spiketimedetection.current_spikes_idxs
    numstoredspikes = integrator.p.spiketimedetection.numstoredspikes
    # @show integrator.callback_cache
    # @show integrator.callback_cache.vector_event_idxs
    # spiketimes[idx, current_spikes_idxs[idx]] = integrator.t
    # current_spikes_idxs[idx] = nextindex(current_spikes_idxs[idx], numstoredspikes)
    # for idx_spiked::Int64 in findall(x->x==true, integrator.callback_cache.vector_event_idxs) #idx of all neurons that spiked
    for idx_spiked in findall(x->x==true, integrator.callback_cache.vector_event_idxs) #idx of all neurons that spiked
        # println("t = $(integrator.t), idx spiked = $idx_spiked, current idx = $(current_spikes_idxs[idx_spiked]), next =$(nextindex(current_spikes_idxs[idx_spiked], numstoredspikes))")
        spiketimes[idx_spiked, current_spikes_idxs[idx_spiked]] = integrator.t
        current_spikes_idxs[idx_spiked] = nextindex(current_spikes_idxs[idx_spiked], numstoredspikes)
    end
end

function alphafunction(t, t_s, p)
    δt = (t-t_s)/p.τs
    return δt * exp(1 - δt)
end

function alphafunction(t, t_s::Float64, τs::Float64)
    δt = (t-t_s)/τs
    return δt * exp(1 - δt)
end


function alphafunction(t, t_s::Nothing, τs::Float64)
    return 0.0
end

