
using Attractors
include("$(DIR)/sims/systems/interface_dynamical_systems.jl");

function find_attractors(pvals; numics=20, gridlength=200, reduce_method="uniform", grid=nothing, gridrange=[(-100, 30), (-0.5, 1.0)],force_non_adaptive=false, mx_chk_loc_att=100, mx_chk_fnd_att=100,kwargs...)
    @unpack N, Δt = pvals
    ds = get_ds(pvals)

    if isnothing(grid)
        xg = range(gridrange[1]...; length = gridlength)
        yg = range(gridrange[2]...; length = gridlength)
        grid = (ntuple(x->xg, N)..., ntuple(x->yg, N)...)
    end
    
    if reduce_method == "uniform"
        reduced_grid = map(g -> range(minimum(g), maximum(g); length = numics), grid)
    else 
        x2 = grid[2][1]; x4  = grid[4][1]
        reduced_grid = (grid[1], x2:x2, grid[3], x4:x4)
    end

    # @show force_non_adaptive
    mapper = AttractorsViaRecurrences(ds, grid; Δt,
        sparse = true, force_non_adaptive, mx_chk_loc_att, mx_chk_fnd_att
    )
    
    sampler, = statespace_sampler(Random.MersenneTwister(1234);
    min_bounds = minimum.(grid), max_bounds = maximum.(grid)
)

    ic_on_sync_manifold = ones(Float64, 2N)
    ics_rand = [sampler() for i in 1:(numics-1)]
    ics = [ic_on_sync_manifold, ics_rand...]
    ics = StateSpaceSet(ics)

    # basins, approx_atts, fs = Attractors.basins_of_attraction_all(mapper, reduced_grid; show_progress=false)
    fs, labels = Attractors.basins_fractions(mapper, ics; show_progress=false)
    atts = extract_attractors(mapper)
    res = @strdict atts fs labels ds
    return res
end



function run_find_attractors(pvals; numics=2, gridrange=[(-1.0, 1.5), (-1.0, 1.0)], gridlength=100, kwargs...)
    res = find_attractors(pvals; numics, gridrange, kwargs...)
    return res
end

function run_find_attractors_integrated(pvals; kwargs...)
    res = run_find_attractors(pvals; kwargs...)
    atts = res["atts"]
    u0s_from_atts = Dict(k => Vector(att[end]) for (k, att) in atts)
    @unpack tend, ttrans, Δt = pvals
    ds = res["ds"]
    atts_integ, ts_integ = integrate_u0s(ds, u0s_from_atts; T=tend-ttrans, Ttr=ttrans, Δt)
    full_res = @strdict atts_integ ts_integ res
    return full_res
end







# function find_attractors_plot_spatiotemporal_profile(pvals; kwargs...)
#     atts_integ, ts_integ, res = run_find_attractors_integrated(pvals; kwargs...)
#     measures = measures_atts(atts_integ, ts_integ; pvals, spikethreshold=0.5)
#     fig, axs = spatiotemporal_profile_attractors(atts_integ, ts_integ; measures, kwargs...)
#     run_details = @strdict atts_integ ts_integ res measures
#     return fig, axs, run_details
# end


function find_attractors_plot_spatiotemporal_profile(pvals, ptypes; force=false, recurrences_kwargs, kwargs...)
    # atts_integ, ts_integ, res = run_find_attractors_integrated(pvals; kwargs...)
    filename = name_result("attractorsinfo-$(string(recurrences_kwargs))", pvals, ptypes)
    full_res, file = produce_or_load(pvals, x->run_find_attractors_integrated(x; recurrences_kwargs..., kwargs...); filename, verbose=true, force, tag=false);
    @unpack atts_integ, ts_integ, res = full_res

    measures = measures_atts(atts_integ, ts_integ; pvals, spikethreshold=0.5)
    fig, axs = spatiotemporal_profile_attractors(atts_integ, ts_integ; measures, kwargs...)
    run_details = @strdict atts_integ ts_integ res measures
    return fig, axs, run_details
end



function spatiotemporal_profile_attractors(atts, ts; measures=nothing, kwargs...)
    # @show atts 
    # @show length(atts)
    numrows, numcols = num_rows_cols(length(atts), 4)
    fig, axs = subplotgrid(numrows, numcols; resolution=(300*numcols, 200*numrows), sharex=true, sharey=true)
    for (idx_iterate, (k, att)) in enumerate(atts)
        idxrow, idxcol = linear_to_2d_idx(idx_iterate, numrows, numcols)
        ax = axs[idxrow, idxcol]
        N = div(size(att, 2), 2)
        xs = Matrix(att[:, 1:N])
        plot_spatiotemporal_profile!(ts, xs; fig, ax, idxrow, idxcol, kwargs...)
        if measures isa Dict
            measure = measures[k]
            ts_R = measure[:, 1]; Rs = measure[:, 2]
            hidey = idxcol == numcols ? false : true
            ax2 = plot_quantifier!(ts_R, Rs; fig, ax, idxrow, idxcol, hidey)
            ttrans = ts[1]; tend = ts[end]
            xlims!(ax, ttrans, tend)
            xlims!(ax2, ttrans, tend)
        end
    end
    return fig, axs
end



function find_attractors_plot_statespace_subspace(pvals; kwargs...)
    atts_integ, ts_integ, res = run_find_attractors_integrated(pvals; kwargs...)
    measures = measures_atts(atts_integ, ts_integ; pvals, spikethreshold=0.5)
    fig, axs = statespace_subspace_attractors(atts_integ, ts_integ; measures)
    return fig, axs, res
end

function statespace_subspace_attractors(atts, ts; measures=nothing)
    # @show atts 
    # @show length(atts)
    numrows, numcols = num_rows_cols(length(atts), 4)
    fig, axs = subplotgrid(numrows, numcols; resolution=(300*numcols, 200*numrows), sharex=true, sharey=true)
    for (idx_iterate, (k, att)) in enumerate(atts)
        idxrow, idxcol = linear_to_2d_idx(idx_iterate, numrows, numcols)
        ax = axs[idxrow, idxcol]
        N = div(size(att, 2), 2)
        xs = Matrix(att[:, 1:N])
        ys = Matrix(att[:, N+1:2N])
        plot_state_space_subspace!(xs, ys; fig, ax, idxrow, idxcol, idx_t_snapshot=length(ts))
        measure = measures[k]; R = round(mean(measure[:, 2]), sigdigits=2)
        ax.title = "R = $R"
    end
    return fig, axs
end


"""
assumes variables in form TxN
"""
function plot_state_space_subspace!(xs, ys; fig=nothing, ax=nothing, idxrow=1, idxcol=1, label="V(t)", idx_t_snapshot=nothing, num_max_atts=Inf)
    N = size(xs, 2)
    for (idx, i) in enumerate(1:N)
        if idx > num_max_atts break end
        lines!(ax, xs[:, i], ys[:, i])#; color=colors[i])
        if idx_t_snapshot isa Number 
            scatter!(ax, xs[idx_t_snapshot, i], ys[idx_t_snapshot, i]; marker=:circle)#; color=colors[i])
        end
    end
    ax.ylabel = "y"
    ax.xlabel = "x"
    return fig, ax
end

function time_series_attractors(atts, ts; measures=nothing)
    # @show atts 
    # @show length(atts)
    numrows, numcols = num_rows_cols(length(atts), 4)
    fig, axs = subplotgrid(numrows, numcols; resolution=(300*numcols, 200*numrows), sharex=true, sharey=true)
    for (idx_iterate, (k, att)) in enumerate(atts)
        idxrow, idxcol = linear_to_2d_idx(idx_iterate, numrows, numcols)
        ax = axs[idxrow, idxcol]
        N = div(size(att, 2), 2)
        xs = Matrix(att[:, 1:N])
        # ys = Matrix(att[:, N+1:2N])
        plot_time_series!(ts, xs; fig, ax, idxrow, idxcol, idx_t_snapshot=length(ts))
        measure = measures[k]; R = round(mean(measure[:, 2]), sigdigits=2)
        ax.title = "R = $R"
    end
    return fig, axs
end


"""
assumes variables in form TxN
"""
function plot_timeseries!(ts, xs; fig=nothing, ax=nothing, idxrow=1, idxcol=1, ylabel="x(t)")
    N = size(xs, 2)
    for i in 1:N 
        lines!(ax, ts, xs[:, i])#; color=colors[i])
    end
    ax.ylabel = ylabel
    ax.xlabel = "t"
    return fig, ax
end

# function count_number_attractors(pvals_list; kwargs...)
#     num_atts = zeros(Int64, length(pvals_list))
#     for (idx, pvals) in enumerate(pvals_list)
#         atts_integ, ts_integ, res = run_find_attractors_integrated(pvals; kwargs...)
#         num_atts[idx] = length(atts_integ)




function continuation_network(pvals_list, param_idx="ϵ")
    @unpack ttrans = pvals_list[1]
    gridrange = [(-1.0, 1.5), (-1.0, 1.0)]; gridlength = 100
    recurrences_kwargs = (mx_chk_lost = 1000, Ttr=ttrans, mx_chk_att=10, mx_chk_loc_att=500, mx_chk_fnd_att=300, sparse=true)
    
    fs_curves, atts_info = attractors_continuation(pvals_list, param_idx; num_ics_per_parameter=100, gridlength, recurrences_kwargs...) 
    # res = @strdict fs_curves atts_info
    # jldsave(filename; res)
    return fs_curves, atts_info 
end


function plot_continuation_network_with_info(pvals_all, param_idx)
    fs_curves, atts_info = continuation_network(pvals_all, param_idx)
    # saved = load(filename)
    # @unpack fs_curves, atts_info = saved["res"]

    params_select = [4.8, 4.85, 4.9, 5]
    params_idx_selected = [findmin( abs.(pvals_all[psymb] .- param_select))[2] for param_select in params_select]
    fig, axs = plot_transition_info(fs_curves, atts_info, pvals_all[psymb], psymb; params_idx_selected)
    fig
    DataInspector(fig)
    fig

    makiesave("$(plotsdir())/inapk/desync_transition/desynctransition_recurrences-ϵ_$ϵ-$param_range-Ttr_$Ttr.png", fig)
    return fig, axs
end

function get_info_continuation_network(pvals_all)
    @unpack ϵ = pvals_all
    param_range = extrema(pvals_all[:I])
    filename = "$(datadir())/inapk/desync_transition/desynctransition_recurrences-ϵ_$ϵ-$param_range-Ttr_$Ttr.jld2"
    # filename = "$(datadir())/inapk/desync_transition/desynctransition_recurrences.jld2"
    saved = load(filename)
    @unpack fs_curves, atts_info = saved["res"]
    idxfirstlc = findfirst(haskey.(atts_info, 4))
    idxbeginchaotic = findlast(haskey.(atts_info, 2))-1
    return fs_curves, atts_info, idxfirstlc, idxbeginchaotic 
end

find_param_idx(paramsall, param_select; psymb=:I) = findmin( abs.(paramsall .- param_select))[2]

function plot_integrated_atts(pvals_all; params_select=nothing, use_sync_coordinates=false)
    fs_curves, atts_info, idxfirstlc, idxbeginchaotic = get_info_continuation_network(pvals_all)

    if !isnothing(params_select)
        psymb = :I
        params_idx_selected = [findmin( abs.(pvals_all[psymb] .- param_select))[2] for param_select in params_select]
        idxs_to_plot = params_idx_selected
    else
        idxs_to_plot = [idxfirstlc, idxbeginchaotic]
    end
    # T = 100; Ttr = 9900; Δt = 0.005
    # T = 100; Ttr = 9900; Δt = 0.001
    # T = 100; Ttr = 99900; Δt = 0.005
    # T = 40; Ttr = 9900; Δt = 0.001
    T = 40; Ttr = 0; Δt = 0.005
    fig, axs = plot_integrated_atts(pvals_all, atts_info, idxs_to_plot, :I; T, Ttr, Δt, use_sync_coordinates)
    supertitle(fig, "T = $(T+Ttr), Ttr = $Ttr, Δt = $Δt")
    # fig
    filename = "$(plotsdir())/inapk/desync_transition/attractors/attractors-T_$T-Ttr_$Ttr-Δt_$Δt-idxs_$idxs_to_plot-sync_$use_sync_coordinates.png"
    # xlims!(axs[2], 1980, 2000)
    return fig, axs, filename
end


"""
receive u0s in a Dictionary (usually for attractors)
"""
function integrate_u0s(ds, u0s::Dict{A, B}; T=1000, Ttr=500, Δt=0.1, kwargs...) where {A,B}
    atts_integ = Dict{A, StateSpaceSet}()
    ts_integ = Float64[]
    for (k, u0) in u0s
        if u0 isa StateSpaceSet u0 = u0[end] end
        # @show k, u0
        reinit!(ds, u0)
        if current_parameters(ds) isa ParametersFHN reinit_synaptic!(ds) end 
        att_integ, ts = trajectory(ds, T, u0; Ttr, Δt) 
        atts_integ[k] = att_integ
        ts_integ = ts
    end
    return atts_integ, ts_integ
end


""" 
Receive u0s in a vector of dictionaries (used on results of continuation analysis)
pvals_list must be the same as given to the continuation!
"""
function integrate_u0s(pvals_list, u0s_all::AbstractArray{A}; kwargs...) where A
    atts_integ_all = Vector{A}(undef, length(u0s_all))
    ts_integ = Float64[]
    for (idx, u0s) in enumerate(u0s_all)
        d = pvals_list[idx]
        @unpack tend, ttrans, Δt = d
        Ttr = ttrans; T = tend - ttrans;
        # @show idx
        ds = get_ds(d; abstol=1e-8, reltol=1e-8)
        atts_integ, ts_integ = integrate_u0s(ds, u0s; T, Ttr, Δt, kwargs...)
        atts_integ_all[idx] = atts_integ
    end 
    return atts_integ_all, ts_integ 
end
    
function integrate_u0s_resindict(pvals_list, u0s_all::AbstractArray{A}; kwargs...) where A
    atts_integ_all, ts_integ = integrate_u0s(pvals_list, u0s_all; kwargs...)
    @strdict atts_integ_all ts_integ 
end

_get_ics(ics::Function, i, fixedics::AbstractArray) = i <= length(fixedics) ? fixedics[i] : ics()
Attractors._get_ic(ics::Function, i) = _get_ics(ics, i, [ones(Float64, 22)])
    
function attractors_continuation(pvals_list, pidx; num_ics_per_parameter=50, gridlength=200, gridrange=[(-100, 30), (-0.5, 1.0)], ic_gridrange=nothing, kwargs...)
    dfirst = pvals_list[1]
    ds = get_ds(dfirst)
    N = dfirst["N"]
    
    xg = range(gridrange[1]...; length = gridlength)
    yg = range(gridrange[2]...; length = gridlength)
    grid = (ntuple(x->xg, N)..., ntuple(x->yg, N)...)
    if isnothing(ic_gridrange) ic_gridrange = grid end

    mapper = AttractorsViaRecurrences(ds, grid; Δt=0.1, kwargs...)
    
    u0_sync_manifold = ones(Float64, 2N)
    fixed_ics = [u0_sync_manifold]

    # Attractors._get_ic(ics::Function, i) = _get_ics(ics, i, fixed_ics)

    sampler, = statespace_sampler(Random.MersenneTwister(1); min_bounds = minimum.(grid), max_bounds = maximum.(grid))
    # rsc = RecurrencesSeededContinuation(mapper; threshold = Inf, distance=StrictlyMinimumDistance())
    # rsc = RecurrencesSeededContinuation(mapper; threshold = Inf, distance=Hausdorff())
    rsc = RecurrencesSeededContinuation(mapper)
    prange = [pvals[string(pidx)] for pvals in pvals_list]
    @info "running continuation"
    fractions_curves, attractors_info = continuation(rsc, prange, pidx, sampler; show_progress = false, samples_per_parameter = num_ics_per_parameter)
    fullres = @strdict fractions_curves attractors_info
    return fullres
end

function plot_transition_info(fs_curves, atts_info, prange, psymb; params_idx_selected=[1,2,3], plot_timescale_info=false, pvals_all=nothing, quantifiers=nothing, max_num_atts=4, access_all=nothing)
    
    fig = Figure(resolution=(1500, 1200)); axs = []
    ax1 = Axis(fig[1, 1:length(params_idx_selected)]); push!(axs, ax1)
    ukeys = unique_keys(fs_curves)
    @show ukeys
    COLORS = get_colors(length(ukeys), :distinguishable_colors)
    colors = Dict(k => (to_color(COLORS[i]), 0.65) for (i, k) in enumerate(ukeys))
    # quantifiers = attractors_quantifiers(prange, atts_info)
    ax1.ylabel = "R"
    
    for (j, k) in enumerate(ukeys)
        scatterlines!(ax1, prange, quantifiers[j]; color=colors[k], linewidth=4, markersize=4)
    end
    ax2 = Axis(fig[2, 1:length(params_idx_selected)]); push!(axs, ax2)
    basins_fractions_plot!(ax2, fs_curves, prange; colors)
    ax2.ylabel = "Basin fractions"
    ax2.xlabel = "I"
    linkxaxes!(ax1, ax2)
    
    access_1 = access_all isa AbstractArray ? access_all[1] : [1, 3]
    access_2 = access_all isa AbstractArray ? access_all[2] : [2, 4]
    
    for (idx, param_idx) in enumerate(params_idx_selected)
        atts = atts_info[param_idx]
        ax = Axis(fig[3, idx])
        plot_attractors(atts; fig, ax, colors, access=access_1, markersize=10, num_max_atts=4); push!(axs, ax)
        ax.title = "$(string(psymb)) = $(prange[param_idx])"
        ax = Axis(fig[4, idx])
        plot_attractors(atts; fig, ax, colors, access=access_2, markersize=10, num_max_atts=4); push!(axs, ax)
    end 

    return fig, axs
end

function attractors_quantifiers(prange, atts_info::AbstractVector, ts; spikethreshold=0.5, Δt=0.1)
    ukeys = unique_keys(atts_info)
    measures = [zeros(length(prange)) for _ in ukeys]
    for i in eachindex(atts_info)
        for (j, k) in enumerate(ukeys)
            if !haskey(atts_info[i], k) 
                measure = NaN
            else
                att = atts_info[i][k]
                N = floor(Int, length(att[1,:])/2)
                xs = Matrix(Matrix(att[:, 1:N])')
                measure =  mean(order_parameter(ts, xs, spikethreshold; Δt)[2])
            end
            measures[j][i] = measure
        end
    end
    return measures 
end

function measures_atts(atts::Dict{A, B}, ts; pvals=nothing, spikethreshold=nothing) where {A,B}
    measures = Dict{A, Any}()
    for (k, att) in atts 
        N = div(size(att, 2), 2)
        xs = Matrix(Matrix(att[:, 1:N])')
        @unpack Δt = pvals
        tsR, Rs =  order_parameter(ts, xs, spikethreshold; Δt)
        measures[k] = hcat(tsR, Rs)
    end 
    return measures
end

function plot_attractors(attractors::Dict; fig=nothing, ax=nothing, access = [1,2], markersize = 12, colors=nothing, num_max_atts=Inf, kwargs...)
    if isnothing(fig) fig = Figure(); ax = Axis(fig[1,1]) end
    ukeys = keys(attractors)
    # @show attractors 
    # @show ukeys
    if isnothing(colors)
        COLORS = get_colors(length(ukeys), :distinguishable_colors)
        colors = Dict(k => (to_color(COLORS[i]), 0.75) for (i, k) in enumerate(ukeys))
    end
    for (idx, k) in enumerate(ukeys)
        if idx > num_max_atts break end
        att = attractors[k]
        # scatter!(ax, vec(attractors[k][:, access]); color = colors[k], label = "$k", markersize = markersize, kwargs...)
        # lines!(ax, vec(att[:, access]); color = colors[k], label = "$k", markersize = markersize, linewidth=0.5, kwargs...)
        N = div(size(att, 2), 2)
        xs = Matrix(att[:, 1:N])
        ys = Matrix(att[:, N+1:2N])
        plot_state_space_subspace!(xs, ys; fig, ax, idxrow, idxcol, idx_t_snapshot=length(ts), num_max_atts)

    end
    axislegend(ax)
    return fig
end


function plot_basins(basins, grid; fig=nothing, ax=nothing)
    xg, yg = grid
    if isnothing(fig) fig = Figure(); ax = Axis(fig[1,1]) end
    ids = sort!(unique(basins))
    cmap = generate_cmap(length(ids))
    hmap = heatmap!(ax, xg, yg, basins; colormap = cmap, colorrange = (ids[1] - 0.5, ids[end]+0.5))
    return fig, ax
end

function plot_timeseries(attractors::Dict, ts; fig=nothing, ax=nothing, access = 1, colors=nothing)
    if isnothing(fig) fig = Figure(); ax = Axis(fig[1,1]) end
    ukeys = keys(attractors)
    if isnothing(colors)
        colors = Dict(k => (to_color(COLORS[i]), 0.75) for (i, k) in enumerate(ukeys))
    end
    for k in ukeys
        lines!(ax, ts, vec(attractors[k][:, access]); color = colors[k],
        label = "$k")
    end
    axislegend(ax)
    return fig
end


function fractions_to_cumulative(fractions_curves, prange)
    ukeys = unique_keys(fractions_curves)
    bands = [zeros(length(prange)) for _ in ukeys]
    for i in eachindex(fractions_curves)
        for (j, k) in enumerate(ukeys)
            bands[j][i] = get(fractions_curves[i], k, 0)
        end
    end
    # transform to cumulative sum
    for j in 2:length(bands)
        bands[j] .+= bands[j-1]
    end
    return ukeys, bands
end

function basins_fractions_plot!(ax, fractions_curves, prange; colors=nothing,
        add_legend = false
    )
    ukeys, bands = fractions_to_cumulative(fractions_curves,prange)

    for (j, k) in enumerate(ukeys)
        if j == 1
            low, upp = [0.0 for _ in prange], bands[j]
        else
            low, upp = bands[j-1], bands[j]
        end
        color = isnothing(colors) ? Cycled(j) : colors[j]
        band!(ax, prange, low, upp; color, label = "$k")
    end
    ylims!(ax, 0, 1); xlims!(ax, minimum(prange), maximum(prange))
    add_legend && axislegend(ax; position = :lt)
    return
end


function fixed_point_analysis(pvals)
    ds = get_ds(pvals)
    
    xs = interval(-3.0, 3.0)
    ys = interval(-2.0, 2.0)
    box = IntervalBox()
    J = jacobian_FHN 
    method = IntervalRootFinding.Krawczyk
    
    fp, eigs, stable = fixedpoints(ds, box, J; method)
end