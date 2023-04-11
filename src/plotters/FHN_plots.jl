function plot_FHN(sol, ts, params; spiketimes = nothing, plotspikes=false, res=(800,600), haslegend=false, xlims=nothing, ttrans=nothing)
    if plotspikes
        if isnothing(spiketimes) spiketimes = params.spiketimedetection.spiketimes; end #for i=1:size(spiketimes, 1) filter!(x->x!=0.0, spiketimes[i, :]) end
        Vth = params.spiketimedetection.Vth
    end
    if !isnothing(ttrans) idx = findfirst(x->x>ttrans, ts) else idx = 1 end
    N = floor(Int, size(sol,1)/2)
    Vs = sol[1:N, :]
    ws = sol[N+1:end, :]
    fig, axs = subplotgrid(3, 1; sharex = false, sharey=true, resolution = res)
    for i=1:N
        ax = axs[1,1]; ax.ylabel = "V"; ax.xlabel="t";
        lines!(ax, ts[idx:end], Vs[i, idx:end], label="$i")
        if plotspikes
            spiketimes_i = spiketimes isa Matrix ? spiketimes[i, :] : spiketimes[i]
            filter!(x->!isnothing(x), spiketimes_i)
            spiketimes_i = Array{Float64}(spiketimes_i)
            # @show spiketimes_i 
            if isempty(spiketimes_i) continue end
            scatter!(ax, spiketimes_i, [Vth for i=1:length(spiketimes_i)])
        end
        ax = axs[2,1]; ax.ylabel = "w"; ax.xlabel="t";
        lines!(ax, ts[idx:end], ws[i, idx:end])
        ax = axs[3,1]; ax.ylabel = "w"; ax.xlabel="V";
        lines!(ax, Vs[i, idx:end], ws[i, idx:end])
    end
    if haslegend
        axislegend(axs[1,1])
    end
    if !isnothing(xlims)
        for i=1:length(xlims) xlims!(axs[i,1], xlims[i]) end
    end
    # if !isnothing(legend)
    # end
    return fig, axs
end