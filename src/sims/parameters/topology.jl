# ---------------------------------------------------------------------------- #
#                                   TOPOLOGY                                   #
# ---------------------------------------------------------------------------- #

#adjl: pure connections
#connsl: connections + assignment of type; currently this is just an intermediate step
#A and E: matrices with the values of each type inside

"""
Wrapper around pvals
"""
function get_coupling_parameters(pvals)
    coup_type = get(pvals, "coup_type", "synaptic")
    if coup_type == "synaptic"
        @unpack ϵ, gmax_exc, E_exc, gmax_inh, E_inh, τs_exc, τs_inh  = pvals
        connsl = get_connsl(pvals)
        type_to_params = _get_type_to_params(ϵ, gmax_exc, gmax_inh, τs_exc, τs_inh, E_exc, E_inh)
        
        @show connsl
        return get_coupling_parameters(connsl, type_to_params)
    elseif coup_type == "diffusive"
        adjl = get_adjl(pvals)
        @unpack ϵ = pvals
        return ϵ, adjl
    else 
        @error "no coupling type $coup_type"
    end
end

function _get_type_to_params(ϵ, gmax_exc, gmax_inh, τs_exc, τs_inh, E_exc, E_inh)
    gmax_exc *= ϵ; gmax_inh *= ϵ; #multiply by global coup strength already (saves up computations later)
    gmax_exc /= (exp(1)*τs_exc); gmax_inh /= (exp(1)*τs_inh); # normalize alpha function so that strength ofinteractions does not depend on tau
    type_to_params = Dict(:E=>(gmax_exc, E_exc, τs_exc), :I=>(gmax_inh, E_inh, τs_inh))
end

"""
Receives 
* `type_to_params` which is a Dict containing the information of each type of connection; e.g. :E=>(gmaxE, E_E, τE), :I=>(gmaxI, ...)
"""
function get_coupling_parameters(connsl, type_to_params)
    A, E = connlist_to_adjmat(connsl, type_to_params)
    alltypes = keys(type_to_params)
    receivertypes = receiver_types(connsl, alltypes)
    τs_vals = Dict(keys(type_to_params) .=> map(x->x[3], values(type_to_params))) #convenience dict

    return A, E, receivertypes, τs_vals
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
        adjl = [[2,3], [], []]
    elseif topm == "bidirectional"
        adjl = [[2], [1]]
    elseif topm == "wattsstrogatz"
        @unpack N, k, rewiring_prob, topseed = pvals
        filename = "$(datadir())/sims/inputs/N_$N/graph-$topm-N_$N-k_$k-rewiringprob_$rewiring_prob-seed_$topseed.jld2"
        if isfile(filename)
            g = load(filename)["g"]
            adjl = g.fadjlist;
        else
            g = watts_strogatz(N, k, rewiring_prob; is_directed=false, rng=MersenneTwister(topseed))
            safesave(filename, Dict("g"=>g))
            adjl = g.fadjlist; #out neighbors == receivers of each node
            return adjl
        end
    else
        error("No other topm found")
    end
end

"""
connsl: vector of vectors; the i-th entry contains tuples (out-neighbor::int, type_of_connection::symb)  
"""
function get_connsl(pvals)
    @unpack N, topm = pvals
    if topm isa Vector 
        return topm 
    elseif topm == "ER"
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
        connsl = [[(2, :E), (3, :I)], [], []];
        return connsl
    elseif topm == "bidirectional"
        connsl = [[(2, :E)], [(1, :E)]];
        return connsl
    elseif topm == "wattsstrogatz"
            @unpack k, rewiring_prob, topseed, types, probabilities, EIseed = pvals
            connslname = "$(datadir())/sims/inputs/N_$N/connsl-$topm-N_$N-k_$k-rewiringprob_$rewiring_prob-seed_$topseed-types_$types-probabilities_$probabilities-EIseed_$EIseed"
            if isfile(connslname)
                return load(connslname)["connsl"]
            else
                @info("File $connslname was not found. Generating it and saving now.")
                adjl = get_adjl(pvals)
                connsl = connectionslist(adjl, types, probabilities; EIseed)
                safesave(connslname, Dict("connsl"=>connsl))
                return connsl
            end
    elseif topm == "modular"
        @unpack Nmod, intermodseed, k, topseed, types, probabilities, EIseed, modules_topologies, num_connections, connection_decider = pvals
        num_modules = length(modules_topologies)
        connslname = "$(datadir())/sims/inputs/Nmod_$Nmod-nummods_$num_modules/connsl-$topm-nummods_$num_modules-intermodseed_$intermodseed-Nmod_$Nmod-k_$k-seed_$topseed-types_$types-probabilities_$probabilities-EIseed_$EIseed.jld2"
        data, _ = produce_or_load(get_connsl_modular, pvals; filename=connslname, force=true, tag=false)
        connsl = data["connsl"]
        @show connsl
    else
        error("No other topm found")
    end
end

function get_connsl_modular(pvals)
    @unpack N, Nmod, intermodseed, k, topseed, types, probabilities, EIseed, modules_topologies, num_connections, connection_decider = pvals
    num_modules = length(modules_topologies)
    modules_params = [Dict("N"=>Nmod, "k"=>k, "topseed"=>topseed, "types"=>types, "probabilities"=>probabilities, "EIseed"=>EIseed, "topm"=>modules_topologies[i]) for i in 1:num_modules]
    A = build_adjacency_matrix_modules(modules_params, num_connections, connection_decider; intermodseed)
    connsl = adjmat_to_connsl(A)
    data = @strdict connsl
end
    

"""
connsl[i] contains the receiver nodes from node i
adjl[i, :] also does
"""
function adjmat_to_connsl(A::Matrix{T}) where {T}
    symbol_for_null_connection = T <: Symbol ? :O : 0
    (N,N) = size(A)
    connsl = Vector{Vector{Tuple{Int64, Symbol}}}(undef, N)
    for (idx_node, receivers) in enumerate(eachrow(A)) 
        existing_conn_idxs = findall(x->x!=symbol_for_null_connection, receivers)
        existing_conn_types = receivers[existing_conn_idxs]
        connsl[idx_node] = collect( zip(existing_conn_idxs, existing_conn_types) )
    end 
    return connsl
end

function connsl_to_adjmat(connsl)
    N = length(connsl)
    A = Matrix{Symbol}(undef, (N,N))
    fill!(A, :O)
    @inbounds for (i, conn) in enumerate(connsl)
        for (rec, type) in conn
            A[rec, i] = type
        end
    end
    return A
end

"""
Vector of size N; each index i contains information about the neighborhood of i, each element is (index_of_receiving_unit, type_of_the_connection)
Assigns connection to its type randomly
"""
function connectionslist(adjl, types::Vector{T}, probabilities::Vector{Float64}; EIseed=1) where {T}
    N = length(adjl)
    connsl = [Tuple{Int, T}[] for i=1:N]
    # connsl = Vector{Tuple{Int, T}}(undef, N)
    numconns = numberofconnections(adjl)
    types = StatsBase.sample(MersenneTwister(EIseed), types, Weights(probabilities), numconns) #pregenerate, so I can use the rng easily
    cont = 1
    @inbounds for (i, receivers) in enumerate(adjl)
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

"""
Receives the connections list and 
Returns a vector of size N containing the receiver types for each unit. The receiver types a Dictionary
"""
function receiver_types(connsl, alltypes)
    recvtypes_allunits = [receiver_types(i, connsl, alltypes) for i in eachindex(connsl)]
end

function receiver_types(i, connsl, alltypes)
    recvtypes = Dict(alltypes .=> [Int[] for i in eachindex(collect(alltypes))])
    for (rec, recvtype) in connsl[i]
        push!(recvtypes[recvtype], rec)
    end
    return recvtypes
end

function adjmat_to_adjlist(mat)
    N = size(mat)[1]
    v_adjl = [ Int64[] for i=1:N]
    for i=1:N
        v_inneighbors = findall(x->x==1, mat[i,:])
        v_adjl[i] = v_inneighbors
    end
    return v_adjl
end


function multiplyglobalconstant!(adjl::Vector{Vector{Any}}, ϵ)
    for adjl_unit in adjl
        for (neighbor, params_conn) in adjl_unit
            params_conn.gmax *= ϵ
        end
    end
end




numberofconnections(adjl) = sum(length.(adjl))

# function connlist_to_adjmat(adjl)
#     A = zeros(Int64, (length(adjl), length(adjl)))
#     for (i, neighbors) in enumerate(adjl)
#         A[i, neighbors] .= 1
#         A[neighbors, i] .= 1
#     end
#     return A
# end

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