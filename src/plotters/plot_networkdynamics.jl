"""
assumes variables in form TxN
"""
function plot_spatiotemporal_profile!(t, variables; fig=nothing, ax=nothing, idxrow=1, idxcol=1, label="V(t)", plotcolorbar=true, kwargs...)
    # fig, ax = util_create_fig_ax(fig, ax, idxrow, idxcol)
    N = size(variables, 2)
    hm = heatmap!(ax, t, 1:N, variables)
    ax.ylabel = "Unit index"
    ax.xlabel = "t"
    if plotcolorbar Colorbar(fig[idxrow+1, idxcol], hm; label, tellwidth=false, vertical=false, flip_vertical_label=true) end
    return fig, ax
end



function plot_quantifier!(t, measure; fig=nothing, ax=nothing, idxrow=1, idxcol=1, ylabel="R(t)", hidey=false)
    ax = Axis(fig[idxrow, idxcol]; ytickcolor =:red, yticklabelcolor = :red, ylabelcolor=:red, yaxisposition = :right, ylabel)
    hidespines!(ax)
    hidexdecorations!(ax)
    if hidey hideydecorations!(ax) end
    lines!(ax, t, measure; color=:red)
    ylims!(ax, 0, 1)
    return ax 
end



function syncerror(ts, xs, ys; fig=nothing, ax=nothing)
    dists = syncerror_pairwise(xs, ys)
    if isnothing(fig) fig, axs = subplotgrid(1,1); ax = axs[1,1] end
    lines!(ax, ts, dists; color=:black)
    return fig, ax
end

function couplingcurrent_as_color(Icoups_atts::Dict{A, B}; resize_Icoup=false) where {A,B}
    colors = Dict{A, Any}()
    for (k, Icoups) in Icoups_atts
        if resize_Icoup 
            for (idx, col) in enumerate(eachcol(Icoups))
                Icoups[:, idx] = _rescale_to_01(collect(col))
            end
        end
        colors[k] = mapslices(x->norm(x[1:2]), Icoups, dims=2)[:, 1] 
    end
    
    return colors
end
    
function currentratio_as_color(Icoups_atts::Dict{A, B}, Idyns_atts; resize_Icoup=false) where {A,B}
    colors = Dict{A, Any}()
    for (k, Icoups) in Icoups_atts
        Idyns = Idyns_atts[k]
        if resize_Icoup 
            for (idx, col) in enumerate(eachcol(Icoups))
                Icoups[:, idx] = _rescale_to_01(collect(col))
            end
        end
        Icoup_norm = mapslices(x->norm(x[1:2]), Icoups, dims=2)[:, 1] 
        Idyn_norm = mapslices(x->norm(x[1:2]), Idyns, dims=2)[:, 1] 
        if any(Idyn_norm .== 0.0) 
            colors[k] = nothing 
        else
            colors[k] = log10.(Icoup_norm ./ Idyn_norm)
        end
        @show colors[k]
    end
    
    return colors
end


function plot_atts_syncmanifold(atts_syncmanifold, atts_pure, ts; u0s=nothing, projections=["1","2","3"], colors=nothing, angles_projection3 = [[1.2, 0.45]], coordinate_grid = (-5, 0.8, -80, -0.1), kwargs...)
    
    if isnothing(colors)
        ukeys = keys(atts_pure)
        colors = Dict(k => (to_color(COLORS[i]), 0.75) for (i, k) in enumerate(ukeys))
    end
    
    fig = Figure(resolution=(1800, 800)); axs = []
    idxcol = 1
    
    if "1" ∈ projections
        ax = Axis3(fig[1:3, idxcol]; azimuth=4.05, elevation=0.88); push!(axs, ax)
        for (k, att) in atts_syncmanifold
            if is_fp(att) 
                scatter!(ax, att[:, 1], att[:, 2], att[:, 3]; color=colors[k], markersize=20)
            else
                lines!(ax, att[:, 1], att[:, 2], att[:, 3]; color=colors[k])
            end
        end
        ax.xlabel = "x2-x1"; ax.ylabel = "y2-x1"; ax.zlabel="x2+x1"
        x_sm = zeros(10); y_sm = x_sm; z_sm = range(-65, 0, length=10)
        lines!(ax, x_sm, y_sm, z_sm; color=(:black, 0.4))
        idxcol += 1
    end
    
    if "2" ∈ projections
        ax = Axis3(fig[1:3, idxcol]; azimuth=4.05, elevation=0.88); push!(axs, ax)
        for (k, att) in atts_syncmanifold
            if is_fp(att)
                scatter!(ax, att[:, 1], att[:, 2], att[:, 4]; color=colors[k], markersize=20)
            else
                lines!(ax, att[:, 1], att[:, 2], att[:, 4]; color=colors[k])
            end
        end
        ax.xlabel = "x2-x1"; ax.ylabel = "y2-x1"; ax.zlabel="y2+y1"
        x_sm = zeros(10); y_sm = x_sm; z_sm = range(-0.2, 0.8, length=10)
        lines!(ax, x_sm, y_sm, z_sm; color=(:black, 0.4))
        idxcol += 1
    end
    
    if "3" ∈ projections
        for (idx, angles) in enumerate(angles_projection3)
            azimuth, elevation = angles
            ax = Axis3(fig[1:3, idxcol]; azimuth, elevation); push!(axs, ax)
            for (k, att) in atts_pure
                # if is_fp(att)
                    # s=scatter!(ax, att[:, 1], att[:, 3], att[:, 2]; color=colors[k], markersize=20)
                # else
                    l=lines!(ax, att[:, 1], att[:, 3], att[:, 2]; color=colors[k])
                # end
            end
            ax.xlabel = "x1"; ax.ylabel = "y1"; ax.zlabel="x2"
            # xmin = minimum(minimum([att[:, 1] for (k, att) in atts_pure]))
            # xmax = maximum(maximum([att[:, 1] for (k, att) in atts_pure]))
            # ymin = minimum(minimum([att[:, 3] for (k, att) in atts_pure]))
            # ymax = maximum(maximum([att[:, 3] for (k, att) in atts_pure]))
            xmax, ymax, xmin, ymin = coordinate_grid
            x = range(xmin, 1.0*xmax, length=10)
            y = range(ymin, 1.0*ymax, length=10)
            zfunc(x,y) = x 
            surface!(ax, x, y, zfunc; highclip=(:black, 0.4), colorrange=(-1e4, -1e4), transparency=true)

            # if !isnothing(u0s)
            #     # ukeys = keys(u0s)
            #     # colors = Dict(k => (to_color(COLORS[i]), 0.75) for (i, k) in enumerate(ukeys))
            #     for (k, u0) in u0s 
            #         @show u0
            #         scatter!(ax, u0[1], u0[3], u0[2]; color=colors[k], markersize=20)
            #         scatter!(ax, u0[1], u0[4], u0[2]; color=colors[k], markersize=20, marker=:rect)
            #     end
            # end
            # idxcol += 1
            # ax = Axis3(fig[1:3, idxcol]; azimuth=1.2, elevation=0.45); push!(axs, ax)
            # cbar = Colorbar(ax, l)
            idxcol += 1
        end

        
        colgap!(fig.layout, 0.01)
        
    end
            
        idxs = [[1,2], [3,4]] #idxs of variables to plot
        # idx_ic = 3 #red
        idx_ic = 1 
        for (idx_iterate, idx) in enumerate(idxs)
            ax = Axis(fig[3+idx_iterate,:]); push!(axs, ax)
            att = atts_pure[idx_ic] 
            lines!(ax, ts, att[:, idx[1]], color=colors[idx_ic])
            lines!(ax, ts, att[:, idx[2]], color=colors[idx_ic], linestyle=:dash)
        end
    
    return fig, axs 
end

function _plot_uncoupled!(_pvals, u0, T, Δt, axs; Ttr)
    pvals = deepcopy(_pvals)
    pvals["ϵ"] = 0.0
    ds = get_ds(pvals)
    att, ts = trajectory(ds, T, u0; Ttr, Δt)
    lines!(axs[1], att[:, 1], att[:, 3], att[:, 2]; color=:blue)
    lines!(axs[2], att[:, 1], att[:, 3], att[:, 2]; color=:blue)
    lines!(axs[3], att[:, 1], att[:, 3], att[:, 2]; color=:blue)
    
    
    ds = get_ds(pvals)
    att, ts = trajectory(ds, T, u0; Ttr, Δt)
    lines!(axs[4], ts, att[:, 1], color=(:blue, 0.75))
    axs[4].ylabel = "x"
    xlims!(axs[4], Ttr, T)
    lines!(axs[5], ts, att[:, 3], color=(:blue, 0.75))
    axs[5].ylabel = "y"; axs[5].xlabel = "t"
    xlims!(axs[5], Ttr, T)
end

"""
att is a Tx2N matrix
"""
function sync_coordinates(att::Matrix)
    diff_x1 = (att[:, 2] .- att[:, 1] ) ./ 2
    diff_x2 = (att[:, 4] .- att[:, 3]) ./ 2
    sum_x1 = (att[:, 2] .+ att[:, 1]) ./ 2
    sum_x2 = (att[:, 4] .+ att[:, 3]) ./ 2
    return hcat(diff_x1, diff_x2, sum_x1, sum_x2)
end



function sync_coordinates(atts::Dict{A,B}) where {A,B}
    atts_sync = Dict{A,B}()
    for (k, att) in atts 
        atts_sync[k] = sync_coordinates(att)
    end 
    return atts_sync 
end

function sync_coordinates(atts_all::Vector{Dict{A,B}}) where {A,B}
    map(atts->sync_coordinates(atts), atts_all)
end


function plot_time_series!(t, variables; fig=nothing, ax=nothing, idxrow=1, idxcol=1, ylabel="R(t)", color=:red)
    if isnothing(ax) ax = Axis(fig[idxrow, idxcol]) end
    hidespines!(ax, :r, :t)
    lines!(ax, t, variables; color)
    return ax 
end
    

function plot_voltages(sols::Vector; axs=nothing, fig=nothing)
    if isnothing(fig) fig = Figure() end
    for (idx, sol) in enumerate(sols)
        ax = 
        if isnothing(axs)
            idxrow = idx; idxcol = 1
            Axis(fig[idxrow, idxcol])
        else 
            axs[idx]
        end
        plot_voltages!(sol; ax)
    end
    return fig 
end

function plot_voltages!(sol::ODESolution; ax=nothing)
    N = size_network(sol)
    Vs = view(sol, 1:N, :)
    ts = sol.t
    colors = get_colors(size(Vs, 1), :distinguishable_colors )
    for (i, Vs_unit) in enumerate(eachrow(Vs))
        plot_time_series!(ts, Vs_unit; ax, color=colors[i])
    end
end


function plot_state_space(sols::Vector; fig=nothing, axs=nothing, idxrow=nothing, idxcol=nothing)
    if isnothing(fig) fig = Figure() end
    for (idx, sol) in enumerate(sols)
        ax = 
        if isnothing(axs)
            idxrow = idx; idxcol = 1
            Axis(fig[idxrow, idxcol])
        else 
            axs[idx]
        end
        plot_state_space!(sol; ax)
    end
    return fig 
end


function plot_state_space!(sol; ax=nothing)
    N = size_network(sol)
    Vs = view(sol, 1:N, :)
    ws = view(sol, N+1:2N, :)
    colors = get_colors(size(Vs, 1), :distinguishable_colors )
    for i in 1:N
        lines!(ax, Vs[i, :], ws[i,:], color=colors[i])
    end
end


function plot_voltages_state_space(sols)
    numrows = length(sols); numcols=  2
    fig, axs = subplotgrid(numrows, numcols)
    plot_voltages(sols; axs=axs[:, 1])
    plot_state_space(sols; axs=axs[:, 2])
    return fig
end