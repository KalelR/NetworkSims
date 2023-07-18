using Random, StatsBase

"""
can return matrix of symbols
`topology_of_modules`: index i is the `topm` for module i
`num_connections`: matrix of NxN

Currently implemented connection_deciders:
- _decide_source_and_targets_only_E_to_I 
- _decide_source_and_targets_uniform
"""
function build_adjacency_matrix_modules(modules_params, num_connections, connection_decider=_decide_source_and_targets_only_E_to_I; intermodseed=1)

    As = Vector{Matrix{Symbol}}(undef, length(modules_params))
    for i in eachindex(modules_params)
        module_params = modules_params[i]
        connsl = get_connsl(module_params)
        @show connsl
        A = connsl_to_adjmat(connsl)
        check_solitary_nodes(A)
        @show A
        As[i] = A
    end
    inter_modular_specifier = Dict("num_connections"=>num_connections) 
    
    Afull = connect_modules(As, inter_modular_specifier, connection_decider; net_seed=intermodseed)
    @show Afull
    
    return Afull
end

function check_solitary_nodes(A)
    (N,N) = size(A)
    for i = 1:N
        incoming = A[:, i]
        filter!(x->x!=:O, incoming)
        if length(incoming) == 0 
            @warn "node $i does not receive any input"
        end
        
        outgoing = A[i, :]
        filter!(x->x!=:O, outgoing)
        if length(outgoing) == 0 
            @warn "node $i does not send any input"
        end
    end
    nothing
end

    

"""
Considers entry A[idx_source, idx_target].
#TODO: support non-directed networks
"""
function connect_modules(As::Vector{Matrix{T}}, inter_modular_specifier::Dict, connections_decider=_decide_source_and_targets_uniform; kwargs...) where {T}
    @show As
    sizes_modules = map(Ai->size(Ai, 1), As)
    Ntot = sum(sizes_modules)
    Afull = Matrix{Symbol}(undef, (Ntot, Ntot)); fill!(Afull, :O)
    
    #copy each module in As into A
    for (idx_module, Ai) in enumerate(As)
        copy_adj!(Ai, Afull, idx_module, sizes_modules)
    end
    
    connect_modules!(Afull, As, inter_modular_specifier, connections_decider; kwargs...)

    
    return Afull
end

function copy_adj!(A_source, A_target, module_idx, sizes_modules)
    idxs = _idxs_module(A_source, module_idx, sizes_modules)
    A_target[idxs, idxs] .= A_source
    nothing 
end

function connect_modules!(A, As, inter_modular_specifier, connections_decider; net_seed = 1, kwargs...)
    sizes_modules = map(Ai->size(Ai, 1), As)
    rng = MersenneTwister(net_seed)
    for idx_module_i in 1:length(As)
        Ai = As[idx_module_i]
        for idx_module_j in 1:length(As)
            if idx_module_i == idx_module_j continue end
            Aj = As[idx_module_j]
            connect_modules_pairwise!(A, Ai, idx_module_i, Aj, idx_module_j, inter_modular_specifier, sizes_modules, rng, connections_decider; kwargs...)
        end
    end
    nothing
end

function connect_modules_pairwise!(A, Ai, i, Aj, j, inter_modular_specifier, sizes_modules, rng, connections_decider; is_directed::Bool = false, kwargs...)
    idxs_all = _idxs_module(A, 1, sizes_modules)
    idxs_source_module = _idxs_module(Ai, i, sizes_modules)
    idxs_target_module = _idxs_module(Aj, j, sizes_modules)
    
    idxs_source, idxs_target, connection_types = source_and_targets(rng, idxs_source_module, idxs_target_module, A, connections_decider, inter_modular_specifier, i, j; kwargs...)
    @show idxs_source
    @show idxs_target
    @show connection_types
    _assign_connections!(A, idxs_source, idxs_target, connection_types)
    
    nothing
end

function source_and_targets(rng, idxs_source_module, idxs_target_module, A, connections_decider::Function, inter_modular_specifier, idx_module_i, idx_module_j; kwargs...)
    idxs_permissible_connections = _available_connections(idxs_source_module, idxs_target_module, A)
    return connections_decider(rng, idxs_permissible_connections, inter_modular_specifier, idx_module_i, idx_module_j, A; kwargs...)
end

"""
Find the empty entries in the adjacency matrix A that are valid connections - this means excluding the on-diagonal entries and the intra-modular entries.
"""
function _available_connections(idxs_source_module, idxs_target_module, A)
    idxs_offdiag = offdiag_idxs(A) #cartesian
    idxs_empty_entries_offdiag = idxs_offdiag[findall(x->A[x] == :O, idxs_offdiag)] 
    idxs_empty_entries_permissible = [idx for idx in idxs_empty_entries_offdiag if ( !(idx[2] ∈ idxs_source_module ) && (idx[1] ∈ idxs_source_module) && (idx[2] ∈ idxs_target_module) )] 
    return idxs_empty_entries_permissible
end

function _decide_source_and_targets_uniform(rng, idxs_permissible_connections, inter_modular_specifier, idx_module_i, idx_module_j, A; kwargs...)
    num_connections_ij = inter_modular_specifier["num_connections"][idx_module_i, idx_module_j]
    
    max_number_conns = length(idxs_permissible_connections)
    if max_number_conns < num_connections_ij 
        @warn "maximum number of connections $max_number_conns is smaller than the number of requested connections $num_connections_ij - the matrix is too full. Filling it with all possible connections!"
    end
    num_connections = clamp(num_connections_ij, 0, max_number_conns)
    @show idxs_permissible_connections 
    @show num_connections
    chosen_connections = sample(rng, idxs_permissible_connections, num_connections; replace=false)
    idxs_source = [chosen_connection[1] for chosen_connection in chosen_connections]
    idxs_target = [chosen_connection[2] for chosen_connection in chosen_connections]
    connection_types = [:E for _ in eachindex(idxs_source)]
    
    return idxs_source, idxs_target, connection_types
end

function _decide_source_and_targets_only_E_to_I(rng, idxs_permissible_connections, inter_modular_specifier, idx_module_i, idx_module_j, A; kwargs...)
    node_types = identify_node_types(A) 
    idxs_connections = [idx for idx in idxs_permissible_connections if ( (node_types[idx[1]] == :E) && (node_types[idx[2]] == :I)   )] #excitatory to inh only!
    @show idxs_connections
    return _decide_source_and_targets_uniform(rng, idxs_connections, inter_modular_specifier, idx_module_i, idx_module_j; kwargs...)
end

function identify_node_types(A)
    (N,N)=size(A)
    node_types = Vector{Symbol}(undef, N)
    for (idx_node, row) in enumerate(eachrow(A))
        unique_conns = unique(row)
        filter!(x->x!=:O, unique_conns)
        if length(unique_conns) == 1 
            node_type = unique_conns[1]
        else 
            @error "More than one connection type in row for neuron $idx_node, with row $row"
        end
        
        node_types[idx_node] = node_type
    end
    return node_types
end



function _assign_connections!(A, idxs_s, idxs_t, connection_types=[:E for _ in eachindex(idxs_s)]; is_directed = false)
    num_conns = length(idxs_s)
    for idx_connection in 1:num_conns
        idx_t = idxs_t[idx_connection]
        idx_s = idxs_s[idx_connection]
        conn_type = connection_types[idx_connection]
        A[idx_s, idx_t] = conn_type #todo: check ordering (matters for directed nets!)
        if is_directed 
            A[idx_t, idx_s] = conn_type
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
    @show A
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

function plot_adjacency_matrix!(_A::Matrix{Symbol}, args...; kwargs...)
    A = convert_symbol_to_int(_A)
    plot_adjacency_matrix!(A, args...; kwargs...)
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



function count_connections(A)
    length(findall(x->x!=:O, A))
end

function convert_int_to_symbol(_A)
    A = replace(_A, 1=>:E, 0=>:O, -1=>:I)
    return convert.(Symbol, A)
end


function convert_symbol_to_int(_A)
    A = replace(_A, :E=>1, :I=>-1, :O=>0)
    return convert.(Int, A)
end