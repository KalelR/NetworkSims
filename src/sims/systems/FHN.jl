using FunctionWrappers
import FunctionWrappers: FunctionWrapper

mutable struct SpikeTimeDetection
    spiketimes :: Matrix{Float64} #{N, numts}
    current_spikes_idxs :: Vector{Int64}
    Vth :: Float64
    numstoredspikes :: Int64
end

@inbounds function f!(dVs, V, w, p, t)
    @unpack a, I = p
    @. dVs = V * (a - V) * (V - 1) - w + I
    nothing
end

@inbounds function g!(dws, V, w, p, t)
    @unpack b, c, d = p
    @. dws = d*(b*V - c*w)
    nothing
end

# """
# Current that each unit is *sending*, not receiving. This is different than what I used to do.
# On my notes, C[i,j] = current that j sends to i. Here, I'm actually
# calculating transpose(C);
# In this case, C[i,j] is the current that i sends to j
# Conversely, also  C[i,j] is the current that j receives from i
# """
@inbounds function syn_sum!(p, t)
    C = view(p.Ct, :, :)
    spiketimes = p.spiketimedetection.spiketimes
    conductance = p.conductance
    Csum = p.Csum
    τs_vals = p.τs_vals
    alltypes = keys(τs_vals)
    spiketimes_rows = eachrow(spiketimes)
    for (i, spiketimes_i) in enumerate(spiketimes_rows)
        fill!(Csum, 0.0)
        receivers = p.adjl[i]::Vector{Int64}
        # if isempty(receivers) view(C, i, :) .= 0.0; continue end #could change this to look at receiver types but oh well i have adjl already...
        if isempty(receivers) continue end #could change this to look at receiver types but oh well i have adjl already...
        for t_s in spiketimes_i
            # Csum[1] += conductance(t, t_s, τs_vals[:E])
            # Csum[2] += conductance(t, t_s, τs_vals[:I])
            Csum[1] += alphafunction(t, t_s, τs_vals[:E])
            Csum[2] += alphafunction(t, t_s, τs_vals[:I])
        end
        C[i, p.receivers_types[i][:E]] .= Csum[1]; #transpose of actual C matrix; element [i,j] contains current that i sends to j
        C[i, p.receivers_types[i][:I]] .= Csum[2]; #transpose of actual C matrix; element [i,j] contains current that i sends to j
    end
    return nothing
end


# ---------------------------------------------------------------------------- #
#                               Dense formulation                              #
# ---------------------------------------------------------------------------- #
@inbounds function fitzhugh_nagumo_mat_rule!(du, u, p, t)
    @unpack a, b, c, I, A, E, ΔV, Icoup, AΔV, Ct = p
    dVs = view(du, 1, :)
    dws = view(du, 2, :)
    Vs =  view(u, 1, :)
    ws =  view(u, 2, :)

    @. ΔV = E - Vs #driving force term, 1xN
    @. AΔV = A * ΔV #gmax*(E-V)*Aij, NxN
    syn_sum!(p, t) #Ct = sum of conductances, 1xN
    diagmul_d!(Icoup, AΔV, Ct) #Icoup, (allocates)
    f!(dVs, Vs, ws, p, t)
    @. dVs = dVs + Icoup
    g!(dws, Vs, ws, p, t)
    nothing
end

mutable struct ParametersFHN
    a :: Vector{Float64}
    b :: Vector{Float64}
    c :: Vector{Float64}
    d :: Vector{Float64}
    I :: Vector{Float64}
    A :: Matrix{Float64}
    adjl ::Vector{Vector{Int64}}
    E :: Matrix{Float64}
    τs_vals :: Dict{Symbol, Float64}
    receivers_types :: Vector{Dict{Symbol, Vector{Int64}}}
    conductance :: FunctionWrapper{Float64, Tuple{Float64, Float64, Float64}}
    spiketimedetection :: SpikeTimeDetection
    Ct :: Matrix{Float64}
    Csum :: Vector{Float64}
    Icoup :: Vector{Float64}
    ΔV :: Matrix{Float64}
    AΔV :: Matrix{Float64}
end

"""
Aij = weight of j->i
connsl[i] = [(unit, type), (unit2, type2), ...] of receivers of conns from node i
"""
function ParametersFHN(pvals, args...)
    @unpack a, b, c, d, I, numstoredspiketimes, Vth = pvals
    #topology
    A, adjl, E, receivertypes, τs_vals = get_coupling_parameters(pvals, args...)
    conductance = alphafunction
    N = length(adjl)
    #preallocated
    spiketimes = zeros(N, numstoredspiketimes);
    Icoup = zeros(Float64, N)
    ΔV = zeros(Float64, (N, N))
    AΔV = zeros(Float64, size(A))
    Ct = zeros(Float64, (N,N))

    Csum = zeros(Float64, 2) #one for Exc, another for Inh; in general length(unique(τsmat))
    current_spikes_idxs = ones(Int64, N)
    spiketimedetection = SpikeTimeDetection(spiketimes, current_spikes_idxs, Vth, numstoredspiketimes)

    return ParametersFHN(
    a, b, c, d, I, #units
    A, adjl, E, τs_vals, receivertypes, conductance, #coupling
    spiketimedetection, Ct, Csum, Icoup, ΔV, AΔV #prealloc
    )
end

# ---------------------------------------------------------------------------- #
#                              Sparse formulation                              #
# ---------------------------------------------------------------------------- #
mutable struct ParametersFHN_s
    a :: Vector{Float64}
    b :: Vector{Float64}
    c :: Vector{Float64}
    d :: Vector{Float64}
    I :: Vector{Float64}
    At :: SparseMatrixCSC{Float64, Int64} #sparse
    adjl ::Vector{Vector{Int64}}
    Et :: Matrix{Float64} #dense
    τs_vals :: Dict{Symbol, Float64}
    receivers_types :: Vector{Dict{Symbol, Vector{Int64}}}
    conductance :: FunctionWrapper{Float64, Tuple{Float64, Float64, Float64}}
    spiketimedetection :: SpikeTimeDetection
    Ct :: SparseMatrixCSC{Float64, Int64} #sparse
    Csum :: Vector{Float64}
    Icoup :: Vector{Float64}
    ΔVt :: Matrix{Float64} #dense
    AΔVt :: SparseMatrixCSC{Float64, Int64} #sparse
end

"""
Aij = weight of j->i
connsl[i] = [(unit, type), (unit2, type2), ...] of receivers of conns from node i
"""
function ParametersFHN_s(pvals, args...)
    @unpack a, b, c, d, I, numstoredspiketimes, Vth = pvals
    #topology
    A, adjl, E, receivertypes, τs_vals = get_coupling_parameters(pvals, args...)
    conductance = alphafunction
    N = length(adjl)
    #preallocated
    spiketimes = zeros(N, numstoredspiketimes);
    Icoup = zeros(Float64, N)
    ΔVt = zeros(Float64, (N, N))
    At = sparse(A');
    AΔVt = deepcopy(At) #sparse like A
    Ct = deepcopy(transpose(At)) #I think it has a form like A
    Et = E'

    Csum = zeros(Float64, 2) #one for Exc, another for Inh; in general length(unique(τsmat))
    current_spikes_idxs = ones(Int64, N)
    spiketimedetection = SpikeTimeDetection(spiketimes, current_spikes_idxs, Vth, numstoredspiketimes)

    return ParametersFHN_s(
    a, b, c, d, I, #units
    At, adjl, Et, τs_vals, receivertypes, conductance, #coupling
    spiketimedetection, Ct, Csum, Icoup, ΔVt, AΔVt #prealloc
    )
end

@inbounds function fitzhugh_nagumo_mat_rule_sparse!(du, u, p, t)
    @unpack a, b, c, I, At, Et, ΔVt, Icoup, AΔVt, Ct = p
    dVs = view(du, 1, :)
    dws = view(du, 2, :)
    Vs =  view(u, 1, :)
    Vst =  view(u, 1, :)'
    ws =  view(u, 2, :)

    @. ΔVt = Et - Vst #driving force term, 1xN
    @. AΔVt = ΔVt * At #gmax*(E-V)*Aij, NxN
    syn_sum!(p, t) #Ct = sum of conductances, 1xN
    diagmul_t!(Icoup, AΔVt, Ct) #Icoup, (allocates)

    f!(dVs, Vs, ws, p, t)
    @. dVs = dVs + Icoup
    g!(dws, Vs, ws, p, t)
    nothing
end






# ---------------------------------------------------------------------------- #
#                           Constant tau formulation                           #
                # (C is a vector, tested and works fine and is fast!)
# ---------------------------------------------------------------------------- #
function ParametersFHN_consttau(pvals)
    @unpack a, b, c, I, τs, numstoredspiketimes = pvals
    #topology
    A, adjl, E = get_coupling_parameters(pvals)
    conductance = alphafunction
    N = length(adjl)

    #preallocated
    spiketimes = zeros(N, numstoredspiketimes);
    C = zeros(Float64, N)
    Icoup = zeros(Float64, N)
    ΔV = zeros(Float64, (N, N))
    AΔV = zeros(Float64, size(A))

    return ParametersFHNconsttau(
    a, b, c, I, #units
    τs, A, adjl, E, conductance, #coupling
    spiketimes, C, Icoup, ΔV, AΔV #prealloc
    )
end

mutable struct ParametersFHNconsttau
    a :: Vector{Float64}
    b :: Vector{Float64}
    c :: Vector{Float64}
    I :: Vector{Float64}
    τs :: Float64
    A :: Matrix{Float64}
    adjl ::Vector{Vector{Int64}}
    E :: Matrix{Float64}
    conductance :: FunctionWrapper{Float64, Tuple{Float64, Float64, ParametersFHNconsttau}}
    spiketimes :: Matrix{Float64}
    C :: Vector{Float64}
    Icoup :: Vector{Float64}
    ΔV :: Matrix{Float64}
    AΔV :: Matrix{Float64}
end
@inbounds function fitzhugh_nagumo_mat_consttau_rule!(du, u, p, t)
    @unpack a, b, c, I, A, E, spiketimes, C, ΔV, Icoup, AΔV = p
    dVs = view(du, 1, :)
    dws = view(du, 2, :)
    Vs = view(u, 1, :)
    ws = view(u, 2, :)

    @. ΔV = E - Vs #1.1 driving force term, 1:N
    @. AΔV = A * ΔV #gmax*(E-V)*Aij
    syn_sum!(C, p, t) #2 g conductance term, 1:N
    mul!(Icoup, A, C) #3 Icoup = sum_j g(t,ts)*gmax*(E-v)*Aij coup term

    f!(dVs, Vs, ws, p, t)
    @. dVs = dVs + Icoup
    g!(dws, Vs, ws, p, t)
    nothing
end


@inbounds function syn_sum_consttau!(C, p, t)
    spiketimes = p.spiketimes
    conductance = p.conductance
    for i in eachindex(C)
        receivers = p.adjl[i]
        if length(receivers) == 0 C[i] = 0.0; continue end
        sum = 0.0
        for t_s in view(spiketimes, i, :)
            sum = sum + conductance(t, t_s, p)
        end
        C[i] = sum;
    end
    return nothing
end

function condition_spike!(out, u, t, integrator)
    Vth = integrator.p.spiketimedetection.Vth
    for i = 1:size(u,2)
        V = u[1, i];
        out[i] = V - Vth
    end
end

function affect_spike!(integrator, idx)
    spiketimes = integrator.p.spiketimedetection.spiketimes
    current_spikes_idxs = integrator.p.spiketimedetection.current_spikes_idxs
    numstoredspikes = integrator.p.spiketimedetection.numstoredspikes
    # @show integrator.callback_cache
    # @show integrator.callback_cache.vector_event_idxs
    # spiketimes[idx, current_spikes_idxs[idx]] = integrator.t
    # current_spikes_idxs[idx] = nextindex(current_spikes_idxs[idx], numstoredspikes)
    for idx_spiked in findall(x->x==true, integrator.callback_cache.vector_event_idxs)
        # println("idx spiked = $idx_spiked, current idx = $(current_spikes_idxs[idx_spiked]), next =$(nextindex(current_spikes_idxs[idx_spiked], numstoredspikes))")
        spiketimes[idx_spiked, current_spikes_idxs[idx_spiked]] = integrator.t
        current_spikes_idxs[idx_spiked] = nextindex(current_spikes_idxs[idx_spiked], numstoredspikes)
    end
end

function nextindex(idx, maxidx)
    if idx == maxidx return 1
    else return idx+1 end
    nothing
end