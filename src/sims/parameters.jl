
function topologyname(pvals)
    @unpack topm = pvals
    if topm == "random"
        return "$(datadir())/sims/inputs/topology/random/adjl-N_$N-k_$k-topseed_$topseed.jld2"
    elseif topm == "121"
        @unpack τs, gmax_inh, gmax_exc, E_inh, E_exc = pvals
        return "$(datadir())/sims/inputs/topology/121/adjl-N_3-121-τs_$τs-gmaxinh_$gmax_inh-gmaxexc_$gmax_exc-Einh_$E_inh-Eexc_$E_exc.jld2"
    elseif topm == "masterslaves"
        @unpack τs, gmax_inh, gmax_exc, E_inh, E_exc = pvals
        return "$(datadir())/sims/inputs/topology/N_3/adjl-N_3-masterslaves-τs_$τs-gmaxinh_$gmax_inh-gmaxexc_$gmax_exc-Einh_$E_inh-Eexc_$E_exc.jld2"
    end
end

function get_topology(pvals)
    topfilename = topologyname(pvals) #
    return load(topfilename)["adjl"]
end

function get_initialconditions(pvals)
    @unpack ictype = pvals
    if ictype == "uniform"
        return get_initialconditions_uniform(pvals)
    end
    @error("No initial condition type could be found for pvals = $pvals.")
end

function get_initialconditions_uniform(pvals)
    @unpack N, icmin, icmax, icseed, ictype = pvals
    icfilename = "$(datadir())/sims/inputs/ics/N_$N-ictype_$ictype/ics-N_$N-type_$ictype-min_$icmin-max_$icmax-seedic_$icseed.jld2"
    if isfile(icfilename)
        return load(icfilename)["ics"]
    else
        @info("File $icfilename with initial conditions was not found. Generating new initial conditions and saving them.")
        numeqs = get_num_eqs(pvals)
        u0s = rand(MersenneTwister(icseed), Uniform(icmin, icmax), (numeqs, N))
        safesave(icfilename, Dict("ics"=>u0s))
        return u0s
    end
end

function get_num_eqs(pvals)
    @unpack unitm = pvals
    if unitm == "FHN"
        return 2
    elseif unitm == "izhikevich"
        return 2
    else
        @warn "incorrect unit model"
    end
end

# ---------------------------------------------------------------------------- #
#                                   TOPOLOGY                                   #
# ---------------------------------------------------------------------------- #
# I'll start with just two types of connetions: E and I; I can later implement other types eg E2, I2, ...
#1. generate adjl and save; one seed here
#2. assign types to the connections in adjl and save as connsl; another seed here

#adjl: pure connections
#connsl: connections + assignment of type; currently this is just an intermediate step
#A and E: matrices with the values of each type inside

function adjmat_to_adjlist(mat)
    N = size(mat)[1]
    v_adjl = [ Int64[] for i=1:N]
    for i=1:N
        v_inneighbors = findall(x->x==1, mat[i,:])
        v_adjl[i] = v_inneighbors
    end
    return v_adjl
end

function connlist_to_adjmat(adjl)
    A = zeros(Int64, (length(adjl), length(adjl)))
    for (i, neighbors) in enumerate(adjl)
        A[i, neighbors] .= 1
        A[neighbors, i] .= 1
    end
    return A
end

function get_coupling_parameters(pvals)
    @unpack ϵ, gmax_exc, E_exc, gmax_inh, E_inh, τs_exc, τs_inh  = pvals
    adjl = get_adjl(pvals)
    connsl = get_connsl(pvals)
    type_to_params = Dict(:E=>(gmax_exc, E_exc, τs_exc), :I=>(gmax_inh, E_inh, τs_inh))
    return get_coupling_parameters(pvals, adjl, connsl, type_to_params)
end

function get_coupling_parameters(pvals, adjl, connsl, type_to_params; ϵ=1.0)
    A, E = connlist_to_adjmat(connsl, type_to_params)
    @unpack ϵ = pvals
    A = ϵ .* A
    alltypes = keys(type_to_params)
    receivertypes = receiver_types(connsl, alltypes)

    τs_vals = Dict(keys(type_to_params) .=> map(x->x[3], values(type_to_params)))

    return A, adjl, E, receivertypes, τs_vals
end

function receiver_types(i, connsl, alltypes)
    recvtypes = Dict(alltypes .=> [Int[] for i in eachindex(collect(alltypes))])
    for (rec, recvtype) in connsl[i]
        push!(recvtypes[recvtype], rec)
    end
    return recvtypes
end

function receiver_types(connsl, alltypes)
    recvtypes_allunits = [receiver_types(i, connsl, alltypes) for i in eachindex(connsl)]
end

function get_connsl(pvals)
    @unpack N, topm = pvals
    if topm == "ER"
        @unpack k, topseed, types, probabilities, EIseed = pvals
        connslname = "$(datadir())/sims/inputs/N_$N/connsl-$topm-N_$N-k_$k-seed_$topseed-types_$types-probabilities_$probabilities-EIseed_$EIseed.jld2"
        if isfile(connslname)
            return load(connslname)["connsl"]
        else
            @info("File $connslname was not found. Generating it and saving now.")
            adjl = get_adjl(pvals)
            connsl = connectionslist(adjl, types, probabilities; EIseed)
            safesave(connslname, Dict("connsl"=>connsl))
            return connsl
        end
    elseif topm == "121"
        connsl = [[(1, :E), (2, :I)], [], []];
        return connsl
    else
        error("No other topm found")
    end
end

"""
adjl[i] contains receivers of connections starting from node i
"""
function get_adjl(pvals)
    @unpack topm = pvals
    if topm == "ER"
        @unpack N, k, topseed = pvals
        filename = "$(datadir())/sims/inputs/N_$N/graph-$topm-N_$N-k_$k-seed_$topseed.jld2"
        if isfile(filename)
            g = load(filename)["g"]
            adjl = g.fadjlist;
        else
            @info("File $filename was not found. Generating it and saving now.")
            Ne = k*N
            g = erdos_renyi(N, Ne; is_directed=true, rng=MersenneTwister(topseed))
            safesave(filename, Dict("g"=>g))
            adjl = g.fadjlist; #out neighbors == receivers of each node
            return adjl
        end
    elseif topm == "121"
        adjl = [[1,2], [], []]
    else
        error("No other topm found")
    end
end

numberofconnections(adjl) = sum(length.(adjl))

@inbounds function connectionslist(adjl, types::Vector{T}, probabilities::Vector{Float64}; EIseed=1) where {T}
    N = length(adjl)
    connsl = [Tuple{Int, T}[] for i=1:N]
    numconns = numberofconnections(adjl)
    types = StatsBase.sample(MersenneTwister(EIseed), types, Weights(probabilities), numconns) #pregenerate, so I can use the rng easily
    cont = 1
    for (i, receivers) in enumerate(adjl)
        conns = Tuple{Int, T}[]
        for (j, rec) in enumerate(receivers)
            type = types[cont]; cont+=1;
            push!(conns, (rec, type))
        end
        connsl[i] =  conns
    end
    return connsl
end

using SparseArrays
"""
type_to_params: Dict(:type=>(weight, e), :type2=>(weight2, e2))
adjl[i]: connections i->nodes; contains nodes that receives connections from i
A[i,j]: sender of i (weight of connection j->i); so A[i, :] control the senders of i; A[:, i] = receivers of i
"""
function connlist_to_adjmat(connl, type_to_params)
    A = zeros(Float64, (length(connl), length(connl)))
    E = zeros(Float64, (length(connl), length(connl)))
    # Tau = zeros(Float64, (length(connl), length(connl)))
    @inbounds for (i, conn) in enumerate(connl)
        for (rec, type) in conn
            weight, e, τ = type_to_params[type]
            A[rec, i] = weight
            E[rec, i] = e
            # Tau[i, rec] = τ
        end
    end
    # return sparse(A), sparse(E)
    return A, E
end