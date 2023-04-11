"""
analysisdict is a dictionary mapping the analysis name to a dictionary containing possible parameters for the analysis
"""
function analyse_solution(sol, analysisdict, pvals, odeprob)
    res = Dict();
    if "sol" in keys(analysisdict) #could include here the saving dt
        res["sol"] = sol
    end
    if "spiketimes" in keys(analysisdict)
        # res["spiketimes"] = odeprob.p.params_coup.coup_type.spiketime_detection.spiketimes
        res["spiketimes"] = odeprob.p.spiketimedetection.spiketimes
    end
    
    if "syncerrorstats" in keys(analysisdict)
        res["syncerrorstats"] = sync_error_stats(sol)
    end
    
    if "phasesyncstats" in keys(analysisdict)
        spikethreshold = 0.5
        @unpack Δt = pvals
        res["phasesyncstats"] = order_parameter_stats(sol, spikethreshold; Δt)
    end

    if "phasesync" in keys(analysisdict)
        spikethreshold = 0.5
        @unpack Δt = pvals
        res["phasesync"] = order_parameter(sol, spikethreshold; Δt)
    end
    
    return res
end


# function spiketimes(sol, pvals)
#     @unpack Vth = pvals
#     N = size(sol, 2)
#     sts = [Float64[] for i=1:N]
#     for i=1:N
#         Vs = view(sol, 1, i, :)
#         dist_to_th = Vs .- Vth #spike when sign changes from -1 to +1; so when diff = 1 - (-1) = 2
#         signs = sign.(dist_to_th)
#         idxs_spikes = findall(x->x==2, diff(signs))
#         sts_i = sol.t[idxs_spikes]
#         sts[i] = sts_i
#     end
#     return sts
# end

using Distances 

function sync_error_stats(sol)
    ts, syncerrors = sync_error(sol)
    E_mean = mean(syncerrors)
    E_std =  std(syncerrors)
    return Dict("mean"=>E_mean, "std"=>E_std)
end


function sync_error(sol)
    N = floor(Int, size(sol, 1)/2)
    xs = sol[1:N, :]
    ys = sol[N+1:2N, :]
    syncerrors = sync_error(xs, ys)
    return sol.t, syncerrors
end

"""
receives xs and ys as matrices of size N x T, returns the sync error as the average of the std over time along each variable separeteldy
"""
function sync_error(xs, ys)
    errors_x = std(xs, dims=1)[1, :]
    errors_y = std(ys, dims=1)[1, :]
    errors = (errors_x .+ errors_y) ./ 2 
end

function sync_error_pairwise(xs, ys)
    traj1 = [xs[:, 1] ys[:, 1]]
    traj2 = [xs[:, 2] ys[:, 2]]
    dists = colwise(Euclidean(), traj1', traj2')
end



# ---------------------------------------------------------------------------- #
#                       Phases and phase synchronization                       #
# ---------------------------------------------------------------------------- #
function spiketimes(sol, Vth)
    ts = sol.t 
    N = floor(Int, size(sol, 1)/2)
    xs = sol[1:N, :]
    return spiketimes(ts, xs, Vth)
end

"""
variables_all as a matrix of size NxT
"""
function spiketimes(ts, variables_all::Matrix, Vth)
    variables_all_vec = eachrow(variables_all)
    sts = [eltype(ts)[] for i in variables_all_vec]
    for (i, variables) in enumerate(variables_all_vec)
        dist_to_th = variables .- Vth #spike when sign changes from -1 to +1; so when diff = 1 - (-1) = 2
        signs = sign.(dist_to_th)
        idxs_spikes = findall(x->x==2, diff(signs))
        sts_i = ts[idxs_spikes]
        sts[i] = sts_i
    end
    return sts
end


function order_parameter(sol, spikethreshold; Δt=0)
    ts = sol.t 
    N = floor(Int, size(sol, 1)/2)
    xs = sol[1:N, :]
    return order_parameter(ts, xs, spikethreshold; Δt)
end

function order_parameter_stats(sol, spikethreshold; Δt=0)
    ts, Rs = order_parameter(sol, spikethreshold; Δt)
    Rmean = mean(Rs)
    Rstd = std(Rs)
    return Dict("mean"=>Rmean, "std"=>Rstd)
end

"""
Calculates R(t), given matrix of phases of Nxtimes
"""
function order_parameter(θs::Matrix{T}) where T
    mapslices(order_parameter, θs, dims=1)[1,:]
end

function order_parameter(θs::Vector{T}) where T
    if isempty(θs) return [-1] end
    θ_sum  = zero(ComplexF64)
    for θ in θs
        θ_sum += exp(im * θ)
    end
    return abs(θ_sum) / length(θs)
end


"""
Variables_all is a matrix of size NxT
"""
function order_parameter(ts, variables_all, spikethreshold; Δt=0)
    if length(ts) != size(variables_all, 2) 
        if length(ts) == size(variables_all, 1)
            @warn "Matrix variables_all has the wrong size, current: $(size(variables_all)), length of ts is $(length(ts)). Transposing it!" 
            variables_all = Matrix(variables_all')
        else
            @error "Matrix variables_all has the wrong size, current: $(size(variables_all)), length of ts is $(length(ts)). " 
        end
    end
    spike_times = spiketimes(ts, variables_all, spikethreshold)
    ts_phases, ϕs = linearphases_network(ts, spike_times; Δt)
    if isempty(ϕs) return [0.0], [-1.0] end
    Rs = order_parameter(ϕs)
    return ts_phases, Rs
end

convert_time_to_idx(t, Δt, t0=0) = round(Int, (t-t0)/Δt + 1, RoundNearest)
convert_idx_to_time(idx, Δt, t0=0) = t0 + (idx-1)*Δt

function find_first_and_last_events(event_times)
    first_events_idx = [event_times_unit[1] for event_times_unit in event_times]
    last_events_idx = [event_times_unit[end] for event_times_unit in event_times]
    
    initial_event_idx = maximum(first_events_idx)
    last_event_idx = minimum(last_events_idx)
    return initial_event_idx, last_event_idx 
end


"""
Since time is usually float, the function mainly deals with the corresponding time indices, which are
Int and thus more exact. This should avoid several issues Ive encountered with small floating point errors in
the time vector leading to differently sized arrays and stuff. 
This is by far the most optimized way Ive tried, though I imagine that there is a more clever way of doing the linearphases_unit! part. 
receives `t` just for t0 
event_times as a vector of vectors, with vector i containing events of unit i 
delta t is the (uniform!) time step of t 
Returns the times at which the phases were computed and the phases as an NxT matrix. T is the size of the time vector, which ranges from the 
latest first event to the earliest last event across units.
"""
function linearphases_network(t, event_times; Δt=0)
    if any(isempty.(event_times)) @warn "Empty event times!!"; return Float64[], Float64[] end
    if Δt == 0 Δt = t[2] - t[1] end
    t0 = t[1] 
    event_times_idx = map(ets_unit->convert_time_to_idx.(ets_unit, Δt, t0), event_times)
    initial_event_idx, last_event_idx = find_first_and_last_events(event_times_idx)
    idxs = initial_event_idx:last_event_idx
    N = length(event_times)
    ϕs = zeros(Float64, (N, length(idxs)))
    
    @inbounds for i in 1:N
        ϕs_unit = view(ϕs, i, :)
        event_times_idx_unit = event_times_idx[i]
        linearphases_unit!(ϕs_unit, event_times_idx_unit, initial_event_idx, last_event_idx)
    end
   
    ts = convert_idx_to_time.(idxs, Δt, t0)
    return ts, ϕs
end

"""
Preallocated ϕs vector
"""
function linearphases_unit!(ϕs, event_times_idx, initial_event_idx, last_event_idx)
    idx_start = findlast(x->x <= initial_event_idx, event_times_idx) #index of start indexing event_times_idxs vector; the time-index corresponding ere is event_times_idxs[idx_start]
    event_times_start_idx = event_times_idx[idx_start]
    idx_end = findfirst(x->x >= last_event_idx, event_times_idx) -1 #index of end indexing event_times_idxs vector; the time-index corresponding ere is event_times_idxs[idx_start]
    @inbounds for i in idx_start:idx_end
        idx_tk = event_times_idx[i]
        idx_tk2 = event_times_idx[i+1]
        first_idx_ts = idx_tk < initial_event_idx ?  initial_event_idx : idx_tk
        last_idx_ts = idx_tk2 < last_event_idx ? idx_tk2 : last_event_idx
        idxs_ts = first_idx_ts : last_idx_ts
        ϕs[idxs_ts .- initial_event_idx .+ 1] .= @. 2π* (idxs_ts - idx_tk) ./ (idx_tk2 - idx_tk) #+ 2π*(i-idx_start)
    end
    return ϕs
end
