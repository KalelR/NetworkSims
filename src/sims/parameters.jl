
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


function adjmat_to_adjlist(mat)
    N = size(mat)[1]
    v_adjl = [ Int64[] for i=1:N]
    for i=1:N
        v_inneighbors = findall(x->x==1, mat[i,:])
        v_adjl[i] = v_inneighbors
    end
    return v_adjl
end

function adjlist_to_adjmat(adjl)
    A = zeros(Int64, (length(adjl), length(adjl)))
    for (i, neighbors) in enumerate(adjl)
        A[i, neighbors] .= 1
        A[neighbors, i] .= 1
    end
    return A
end

# using Graphs
# N = 10
# k = 3
# Ne = k*N
# g = erdos_renyi(N, Ne)
# adjm = adjacency_matrix(g)
# adjl = adjmat_to_adjlist(adjm)

# τss = [1.3]
# gmaxs = [-2, +2]
# Es = [-3, +3]
# using StatsBase
# # a=sample(gmaxs, Weights([0.2, 0.8]), 1000)
# connl = Vector{Vector{Tuple{Int64, ParamsCoupTypeSp}}}(undef, N)
# for (i, neighbors) in enumerate(adjl)
#     conn = Vector{Tuple{Int64, ParamsCoupTypeSp}}(undef, length(neighbors))
#     for (j, neigh) in enumerate(neighbors)
#         gmax=sample(gmaxs, Weights([0.2, 0.8]))
#         τs=sample(τss, Weights([1.0]))
#         E=sample(Es, Weights([0.2, 0.8]))
#         push!(conn, (neigh, ParamsSynAlpha(τs, gmax, E)))
#         conn[j] = (neigh, ParamsSynAlpha(τs, gmax, E))
#     end
#     connl[i] = conn
# end
# connl