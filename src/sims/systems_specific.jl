mutable struct SpikeTimeDetection{A, B, C}
    spiketimes :: Matrix{A} #{N, numts}
    current_st_idx :: Vector{B}
    Vth :: Float64
    numstoredspikes :: C
end

mutable struct SynType <: CoupType
    coup_f_sum! :: Function
    coup_f :: Function
    conductance :: Function #for coupling, eg conductance_alpha
    spiketime_detection :: SpikeTimeDetection
end

mutable struct ParamsSynAlpha <: ParamsCoupTypeSp #specific coupling parameter for each connection
    τs :: Float64
    gmax :: Float64
    E :: Float64
end

mutable struct ParamsFHN <: ParamsUnitDynamics
    a :: Float64
    b :: Float64
    c :: Float64
    I :: Float64
end

function fitzhugh_nagumo_rule!(du, u, p, t)
    V, w = u
    @unpack a, b, c, I = p
    du[1] = V*(a-V)*(V-1) -w + I
    du[2] = b*V - c*w
end

function coup_f_sum_syn!(i, du, u, p_coup, t)
    sum = 0.0;
    for (neighbor, params_sp) in p_coup.adjl[i]
        sum += p_coup.coup_type.coup_f(i, neighbor, u, params_sp, p_coup, t)
    end
    du[1] += sum #could of course generalize to add coupling terms to other dimensions by making sum a vector
end

"""
p_conn = params for this specific connection
p_conns = general parameters for all connections
"""
function coup_f_syn(post_unit, pre_unit, u, p_sp, p_coup, t)
    t_s = 0.0;
    # t_s = p_coup.ts[pre_unit, :];
    @unpack gmax, E = p_sp
    sum = 0.0
    for t_s in p_coup.coup_type.spiketime_detection.spiketimes[pre_unit, :]
        sum += gmax * p_coup.coup_type.conductance(t, t_s, p_sp) * (u[1, post_unit] - E)
    end
    return sum
end

function alphafunction(t, t_s, p)
    δt = (t-t_s)/p.τs
    return δt * exp(1 - δt)
end

function condition_spike(out, u, t, integrator)
    for i = 1:size(u,2)
        V = u[1, i];
        out[i] = V - integrator.p.params_coup.coup_type.spiketime_detection.Vth
    end
end

function affect_spike!(integrator, idx)
    # push!(integrator.p.spiketimes[idx], integrator.t) #TODO: what if multiple neurons spike at the same time-step??
    st_struct = integrator.p.params_coup.coup_type.spiketime_detection
    st_struct.spiketimes[idx, st_struct.current_st_idx[idx]] = integrator.t #TODO: what if multiple neurons spike at the same time-step??
    st_struct.current_st_idx[idx] = nextindex(st_struct.current_st_idx[idx], st_struct.numstoredspikes)
end

function nextindex(idx, maxidx)
    if idx == maxidx return 1
    else return idx+1 end
    nothing
end
