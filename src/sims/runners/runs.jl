function run_spatiotemporalprofile(pvals_list, ptypes)
    for idx in eachindex(pvals_list)
        pvals = pvals_list[idx];
        fig, axs = run_spatiotemporalprofile_single(pvals, ptypes)
    end
end

function run_spatiotemporalprofile_single(pvals, ptypes)
    analysisdict = Dict("sol"=>Dict(), "spiketimes"=>Dict(), "phasesync"=>Dict())
    res, odeprob = makesims_FHN(pvals, analysisdict);
    sol = res["sol"]
    fig, axs = plot_spatiotemporalprofile_single(sol, res)
    @unpack ϵ, probabilities = pvals
    ax = axs[1]
    ax.title = "ϵ = $ϵ, prop of E = $(probabilities[1])"
    figtitle = "spatiotemporalprofile"
    plotsave(figtitle, fig, pvals, ptypes)
    return fig, axs
end


function plot_spatiotemporalprofile_single(sol, res)
    fig, axs = subplotgrid(1, 1; resolution=(1000, 600))
    N = floor(Int, size(sol, 1)/2)
    Vs = sol[1:N, :]'
    t = sol.t

    ax = axs[1]
    plot_spatiotemporal_profile!(t, Vs; fig, ax)
    
    for possible_quantifier in ["phasesync"] 
        if possible_quantifier ∈ keys(res)
            t_quant, quantifier = res[possible_quantifier]
            ax2 = Axis(fig[1, 1]; ytickcolor =:red, yticklabelcolor = :red, ylabelcolor=:red, yaxisposition = :right, ylabel="R(t)")
            hidespines!(ax2)
            hidexdecorations!(ax2)
            lines!(ax2, t_quant, quantifier; color=:red)
            ylims!(ax2, 0, 1)
        end
    end

    return fig, axs
end


function collect_results_transition_to_sync(pvals_list)
    res_all = [Dict() for _ in eachindex(pvals_list)]
    for (idx_iterate, idx) in enumerate(eachindex(pvals_list))
        pvals = pvals_list[idx];
        analysisdict = Dict("syncerrorstats"=>Dict(), "phasesyncstats"=>Dict())
        res, odeprob = makesims_FHN(pvals, analysisdict);
        res_all[idx_iterate]  = res
    end
    return res_all
end

function run_transition_to_sync(pvals_list, ptypes)
    res_all = collect_results_transition_to_sync(pvals_list)
    ϵs = [pvals["ϵ"] for pvals in pvals_list]
    fig, axs = plot_transition_to_sync(res_all, ϵs)
    figtitle = "transitiontosync"
    pvals = pvals_list[1]
    pvals["ϵ"] = "all"
    plotsave(figtitle, fig, pvals, ptypes)
    return fig, axs 
end 

function plot_transition_to_sync(res_all, ϵs)
    fig, axs = subplotgrid(2,1; ylabels=["R", "E"], xlabel="ϵ", sharex=true, resolution=(800, 400))
    
    for (idx, quantifier_string) in  enumerate(["phasesyncstats", "syncerrorstats"])
        mean_quant = [res[quantifier_string]["mean"] for res in res_all]
        std_quant =  [res[quantifier_string]["std"] for res in res_all]
        supband = mean_quant .+ std_quant
        infband = mean_quant .- std_quant
        ax = axs[idx]
        band!(ax, ϵs, infband, supband; color=(:black, 0.6))
        scatterlines!(ax, ϵs, mean_quant; color=:black, markersize=10)
        ax.xscale = log10
        hidespines!(ax, :t, :r)
    end

    return fig, axs 
end 

function run_spatiotemporalprofile_and_transition_to_sync(pvals_list, ptypes)
    res_all = [Dict() for _ in eachindex(pvals_list)]
    ϵs = [pvals["ϵ"] for pvals in pvals_list]
    for idx in eachindex(pvals_list)
        pvals = pvals_list[idx];
        analysisdict = Dict("sol"=>Dict(), "spiketimes"=>Dict(), "phasesync"=>Dict(), "phasesyncstats"=>Dict(), "syncerrorstats"=>Dict())
        res, odeprob = makesims_FHN(pvals, analysisdict);
        sol = res["sol"]
        fig, axs = plot_spatiotemporalprofile_single(sol, res)
        @unpack ϵ, probabilities = pvals
        axs[1].title = "ϵ = $ϵ, prop of E = $(probabilities[1])"
        figtitle = "spatiotemporalprofile"
        plotsave(figtitle, fig, pvals, ptypes)
        
        res_quantifiers = deepcopy(res); delete!(res_quantifiers, "sol")
        res_all[idx]  = res_quantifiers
    end
    fig, axs = plot_transition_to_sync(res_all, ϵs)
    figtitle = "transitiontosync"
    pvals = pvals_list[1]
    pvals["ϵ"] = "all"
    plotsave(figtitle, fig, pvals, ptypes)
    return fig, axs
end

function run_allparams_study_sync_coordinates(pvals_list, ptypes; kwargs...)#, Ttrs=[1000.0], param_vals=[0.2], ϵs=[0.0, 1e-5, 5e-5])
    for pvals in pvals_list
        fig, axs = run_study_sync_coordinates_transients(pvals, ptypes; kwargs...)
    end
end



"""
Plot attractors on sync manifold + time series 
*colorsmode = "fromatts", "couplingcurrent" or "currents_ratio"
"""
function run_study_sync_coordinates_transients(pvals, ptypes; psymb="I", colorsmode="fromatts", gridlength = 100, numics = 10, kwargs...)

    @unpack ttrans, tend, Δt, ϵ = pvals
    T = tend-ttrans; Ttr = ttrans
    ics = _get_ics(pvals)
    unique_ics = Dict(1:length(ics) .=> ics)

    angles_projection3 = [ [0.92, 0.47], [4.7, 0.01], [4.7, 1.57]]
    
    fig, axs = study_sync_coordinates(pvals; T, Ttr, Δt, u0s = unique_ics, projections=["3"], colorsmode, angles_projection3, coordinate_grid = (1.0, 0.6, -0.8, -0.3))
    param_val = pvals[psymb]
    supertitle(fig, "$(string(psymb)) = $param_val, ϵ = $ϵ, Ttr = $ttrans")
    
    
    _plot_uncoupled!(pvals, ics[1], Ttr+500, Δt, axs; Ttr)

    figtitle = "3dstructures_synccoordinates"
    plotsave(figtitle, fig, pvals, ptypes; padright=10)
    
    return fig, axs
end

function _get_ics(pvals)
    ics = get(pvals, "ics", nothing)
    if isnothing(ics) 
        @unpack N = pvals
        return rand(Float64, 2N)
    elseif ics isa Vector{Vector{Vector{<:Number}}}
        return ics 
    elseif ics isa Tuple 
        ic_grid = ics 
        sampler, = statespace_sampler(Random.MersenneTwister(1234); min_bounds = minimum.(ic_grid), max_bounds = maximum.(ic_grid))
        @unpack extra_ics, numics = pvals
        ics_rand = [sampler() for i in 1:(numics-length(extra_ics))]
        # @show extra_ics
        _ics = vcat(extra_ics, ics_rand)
        # @show _ics
        return _ics
        # ics = StateSpaceSet(_ics)
        return ics
    end 
end

"""
Plot attractors on sync manifold given initial conditions; parameters defined in pvals
u0s is a dictionary mapping the ic label (which could be an att label) to an ic (vector)
"""
function study_sync_coordinates(pvals; u0s=nothing, colorsmode="fromatts", kwargs...)
    
    ds = get_ds(pvals)
    atts_integ, ts = _integrate_u0_or_att(ds, u0s; pvals, kwargs...)
    atts_integ_sync = sync_coordinates(atts_integ)
    
    colors = _select_colors(colorsmode, ds, atts_integ, ts; kwargs...)
    fig, axs = plot_atts_syncmanifold(atts_integ_sync, atts_integ, ts; u0s, colors, kwargs...) 
    
    return fig, axs 
end

function _select_colors(colorsmode, ds, atts_integ, ts; resize_Icoup=false, kwargs...)
    if colorsmode == "couplingcurrent" || colorsmode == "currents_ratio"
        p = ds.integ.sol.prob.p
        Icoups_atts = coupling_current(atts_integ, p, ts)
        if colorsmode == "currents_ratio"
            Idyns_atts = dynamics_current(atts_integ, p, ts)
            colors = currentratio_as_color(Icoups_atts, Idyns_atts; resize_Icoup)
        else 
            colors = couplingcurrent_as_color(Icoups_atts; resize_Icoup)
        end
    else 
        colors = nothing 
    end
    return colors 
end

function _integrate_u0_or_att(ds, u0s=nothing; pvals_all=nothing, kwargs...)
    if isnothing(u0s) #get ics from attractors
        # fs_curves, atts_info, idxfirstlc, idxbeginchaotic = get_info_continuation_network(pvals_all)
        # atts = atts_info[param_idx]
        # atts_integ, ts = integrate_atts(ds, atts; kwargs...)
    else
        atts_integ, ts = integrate_u0s(ds, u0s; kwargs...)
    end
    
    return atts_integ, ts
end


