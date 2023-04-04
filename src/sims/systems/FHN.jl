using FunctionWrappers
import FunctionWrappers: FunctionWrapper

include("synaptic_coupling.jl")

mutable struct ParametersFHN
    f:: Function 
    g:: Function 
    a :: Vector{Float64}
    b :: Vector{Float64}
    c :: Vector{Float64}
    d :: Vector{Float64}
    I :: Vector{Float64}
    A :: Matrix{Float64}
    N :: Int64
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
    
    f = f_FHN! 
    g = g_FHN!
    
    #topology
    A, adjl, E, receivertypes, τs_vals = get_coupling_parameters(pvals, args...)
    conductance = alphafunction
    N = length(adjl)
    #preallocated
    spiketimes = Matrix{Union{Nothing, Float64}}(undef, (N, numstoredspiketimes))
    Icoup = zeros(Float64, N)
    ΔV = zeros(Float64, (N, N))
    AΔV = zeros(Float64, size(A))
    Ct = zeros(Float64, (N,N))
    Csum = zeros(Float64, 2) #one for Exc, another for Inh; in general length(unique(τsmat))
    current_spikes_idxs = ones(Int64, N)
    spiketimedetection = SpikeTimeDetection(spiketimes, current_spikes_idxs, Vth, numstoredspiketimes)

    return ParametersFHN(
    f, g, 
    a, b, c, d, I, #units
    A, N, adjl, E, τs_vals, receivertypes, conductance, #coupling
    spiketimedetection, Ct, Csum, Icoup, ΔV, AΔV #prealloc
    )
end

@inbounds function f_FHN!(dVs, V, w, p, t)
    @unpack a, I = p
    @. dVs = V * (a - V) * (V - 1) - w + I
    nothing
end

@inbounds function g_FHN!(dws, V, w, p, t)
    @unpack b, c, d = p
    @. dws = d*(b*V - c*w)
    nothing
end

@inbounds function bidimensional_synaptic_coupling_rule!(du, u, p, t)
    @unpack a, b, c, I, A, E, ΔV, Icoup, AΔV, Ct, N, f, g = p
    dVs = view(du, 1:N)
    dws = view(du, N+1:2N)
    Vs =  view(u, 1:N)
    ws =  view(u, N+1:2N)

    #compute synaptic coupling
    @. ΔV = E - Vs #driving force term, 1xN
    @. AΔV = A * ΔV #gmax*(E-V)*Aij, NxN
    syn_sum!(p, t) #Ct = sum of conductances, 1xN
    diagmul_d!(Icoup, AΔV, Ct) #Icoup, (allocates)
    
    f(dVs, Vs, ws, p, t)
    @. dVs = dVs + Icoup
    g(dws, Vs, ws, p, t)
    nothing
end



