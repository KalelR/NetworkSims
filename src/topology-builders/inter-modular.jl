using Random, StatsBase

"""
Considers entry A[idx_source, idx_target].
#TODO: support non-directed networks
"""
function connect_modules(As::Vector{Matrix{A}}, num_connections::Matrix{Int}; kwargs...) where {A}
    sizes_modules = map(Ai->size(Ai, 1), As)
    Ntot = sum(sizes_modules)
    Afull = zeros(A, (Ntot, Ntot))
    
    #copy each module in As into A
    for (idx_module, Ai) in enumerate(As)
        copy_adj!(Ai, Afull, idx_module, sizes_modules)
    end
    
    connect_modules!(Afull, As, num_connections; kwargs...)

    
    return Afull
end

function copy_adj!(A_source, A_target, module_idx, sizes_modules)
    idxs = _idxs_module(A_source, module_idx, sizes_modules)
    A_target[idxs, idxs] .= A_source
    nothing 
end

function connect_modules!(A, As, num_connections::Matrix{Int}; net_seed = 1, kwargs...)
    sizes_modules = map(Ai->size(Ai, 1), As)
    rng = MersenneTwister(net_seed)
    number_connections_ij = Int64[]
    number_modules = length(As)
    number_pairs_modules = number_modules^2
    for idx_module_i in 1:length(As)
        Ai = As[idx_module_i]
        for idx_module_j in 1:length(As)
            if idx_module_i == idx_module_j continue end
            Aj = As[idx_module_j]
            number_connections_ij = num_connections[idx_module_i, idx_module_j]
            connect_modules_pairwise!(A, Ai, idx_module_i, Aj, idx_module_j, number_connections_ij, sizes_modules, rng; kwargs...)
        end
    end
    nothing
end

function connect_modules_pairwise!(A, Ai, i, Aj, j, num_connections, sizes_modules, rng; is_directed::Bool = false, kwargs...)
    idxs_all = _idxs_module(A, 1, sizes_modules)
    idxs_source_module = _idxs_module(Ai, i, sizes_modules)
    idxs_target_module = _idxs_module(Aj, j, sizes_modules)
    # idxs_others = setdiff(idxs_all, idxs_module)
    
    idxs_source, idxs_target = _decide_source_and_target_idxs(rng, num_connections, idxs_source_module, idxs_target_module, A)
    @show idxs_source
    @show idxs_target
    _assign_connections!(idxs_source, idxs_target, A)
    
    nothing
end

"""
Find the empty entries in the adjacency matrix A that are valid connections - this means excluding the on-diagonal entries and the intra-modular entries.
"""
function _decide_source_and_target_idxs(rng, num_connections_ij, idxs_source_module, idxs_target_module, A)
    idxs_offdiag = offdiag_idxs(A) #cartesian
    idxs_empty_entries_offdiag = idxs_offdiag[findall(x->A[x] == 0, idxs_offdiag)] 
    idxs_empty_entries_permissible = [idx for idx in idxs_empty_entries_offdiag if ( !(idx[2] ∈ idxs_source_module ) && (idx[1] ∈ idxs_source_module) && (idx[2] ∈ idxs_target_module) )] #removed off diag, now remove idxs belonging to the same (source) module #TODO: this is not correct, finding many repeated targets
    max_number_conns = length(idxs_empty_entries_permissible)
    
    if max_number_conns < num_connections_ij 
        @warn "maximum number of connections $max_number_conns is smaller than the number of requested connections $num_connections_ij - the matrix is too full. Filling it with all possible connections!"
    end
    num_connections = clamp(num_connections_ij, 0, max_number_conns)
    chosen_connections = sample(rng, idxs_empty_entries_permissible, num_connections; replace=false)
    idxs_source = [chosen_connection[1] for chosen_connection in chosen_connections]
    idxs_target = [chosen_connection[2] for chosen_connection in chosen_connections]
    
    return idxs_source, idxs_target
end

function _assign_connections!(idxs_s, idxs_t, A; is_directed = false)
    num_conns = length(idxs_s)
    for idx_connection in 1:num_conns
        idx_t = idxs_t[idx_connection]
        idx_s = idxs_s[idx_connection]
        A[idx_s, idx_t] = 1 #todo: check ordering (matters for directed nets!)
        if is_directed 
            A[idx_t, idx_s] = 1 
        end
    end 
    nothing 
end
            

function _idxs_module(A_module, module_idx, sizes_modules)
    (N_module, N_module) = size(A_module)
    sizes_modules_before = module_idx == 1 ? [0] : sizes_modules[1:module_idx-1]
    idx_offset = sum(sizes_modules_before) + 1
    idxs = idx_offset : idx_offset + N_module - 1
end

function offdiag_idxs(A::AbstractMatrix)
    [ι for ι in CartesianIndices(A) if ι[1] ≠ ι[2]]
end 

function plot_adjacency_matrix!(A, ax, fig; idxcol=1)
    (N, N) = size(A)
    sources = 1:N
    targets = 1:N
    heatmap!(ax, targets, sources, A; colormap=:grays)
    ax.yticks = 1:N 
    ax.xticks = 1:N 
    ax.xlabel = "sources"
    ax.ylabel = "targets"
    colsize!(fig.layout, idxcol, Aspect(1, 1.0))
    return nothing
end

function plot_adjacency_matrix(A, As)
    fig = Figure(); 
    for (i, Ai) in enumerate(As)
        ax = Axis(fig[1,i])
        plot_adjacency_matrix!(Ai, ax, fig; idxcol=i)
    end 
    ax = Axis(fig[2:3,1:end])
    plot_adjacency_matrix!(A, ax, fig; idxcol=1)
    
    return fig, ax 
end


using Test
using CairoMakie
@testset "2 modules" begin
    A1 = [0 1 0 1; 0 0 0 1; 0 1 0 0; 0 0 1 0]
    A2 = [0 0 1 1 0; 1 0 0 1 1; 1 1 0 0 0; 0 1 0 0 1; 1 0 1 1 0]
    As = [A1, A2]
    num_connections = [0 4; 4 0]
    # A3 = connect_modules(A1, A2, num_connections)
    A3 = connect_modules(As, num_connections)
    @test A3[1:4, 1:4] == A1 
    @test A3[5:end, 5:end] == A2

    l1=sum(A3[1:4, 5:end])
    l2=sum(A3[5:end, 1:4])
    @test l1+l2 == sum(num_connections)
    fig, axs = plot_adjacency_matrix(A3, As); fig
end

@testset "3 modules" begin
    A1 = [0 1 0 1; 0 0 0 1; 0 1 0 0; 0 0 1 0]
    A2 = [0 0 1 1 0; 1 0 0 1 1; 1 1 0 0 0; 0 1 0 0 1; 1 0 1 1 0]
    A3 = [0 1 0 0 0; 0 0 1 0 1; 0 1 0 0 1; 1 1 0 0 1; 1 1 1 1 0]
    As = [A1, A2, A3]
    num_connections = [0 4 2; 4 0 1; 0 0 0]
    
    Afull = connect_modules(As, num_connections)
    @test  Afull[1:4, 1:4] == A1 
    @test  Afull[5:9, 5:9] == A2
    @test  Afull[10:14, 10:14] == A3
    
    idxs_modules = [1:4, 5:9, 10:14]
    ls = zeros(Int64, (3,3))
    for i=1:3 
        for j=1:3
            if i == j continue  end
            ls[i,j]=sum(Afull[idxs_modules[i], idxs_modules[j]])
        end
    end
    @test ls == num_connections
    ls 
    num_connections
    

    fig, axs = plot_adjacency_matrix(Afull, As); fig
end