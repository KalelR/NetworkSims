using Random
function connect_modules(A1::Matrix{A}, A2::Matrix{A}, num_connections::Int64; kwargs...) where {A}
    (N1, N1) = size(A1)
    (N2, N2) = size(A2)
    Ntot = N1+N2
    A3 = zeros(A, (Ntot, Ntot))
    sizes_modules = [N1, N2]
    copy_adj!(A1, A3, 1, sizes_modules)
    copy_adj!(A2, A3, 2, sizes_modules)
    connect_modules!(A3, A1, A2, num_connections, sizes_modules; kwargs...)
    
    return A3
end

function copy_adj!(A_source, A_target, module_idx, sizes_modules)
    idxs = _idxs_module(A_source, module_idx, sizes_modules)
    A_target[idxs, idxs] .= A_source
    nothing 
end

"""
so far, only directed
"""
function connect_modules!(A, A1, A2, num_connections, sizes_modules; net_seed = 1, is_directed::Bool = false, kwargs...)
    idxs_all = _idxs_module(A, 1, sizes_modules)
    rng = MersenneTwister(net_seed)
    num_modules = 2
    for idx_module in 1:num_modules
        Ai = idx_module == 1 ? A1 : A2
        num_connections_ij = floor(Int64, num_connections/2) #TODO: generalize this
        idxs_module = _idxs_module(Ai, idx_module, sizes_modules)
        idxs_others = setdiff(idxs_all, idxs_module)
        
        idxs_source, idxs_target = _decide_source_and_target_idxs(rng, num_connections_ij, idxs_module, A)
        @show idxs_source
        @show idxs_target
        _assign_connections!(idxs_source, idxs_target, A)
    end
    
    nothing
end


"""
Find the empty entries in the adjacency matrix A that are valid connections - this means excluding the on-diagonal entries and the intra-modular entries.
"""
function _decide_source_and_target_idxs(rng, num_connections_ij, idxs_module, A)
    idxs_offdiag = offdiag_idxs(A) #cartesian
    idxs_empty_entries_offdiag = idxs_offdiag[findall(x->A[x] == 0, idxs_offdiag)] 
    idxs_empty_entries_permissible = [idx for idx in idxs_empty_entries_offdiag if ( !(idx[2] ∈ idxs_module ) && (idx[1] ∈ idxs_module))] #removed off diag, now remove idxs belonging to the same (source) module #TODO: this is not correct, finding many repeated targets
    max_number_conns = length(idxs_empty_entries_permissible)
    
    if max_number_conns < num_connections_ij 
        @warn "maximum number of connections $max_number_conns is smaller than the number of requested connections $num_connections_ij - the matrix is too full. Filling it with all possible connections!"
    end
    num_connections = clamp(num_connections_ij, 1, max_number_conns)
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


using Test
@testset "2 modules" begin
    A1 = [0 1 0 1; 0 0 0 1; 0 1 0 0; 0 0 1 0]
    A2 = [0 0 1 1 0; 1 0 0 1 1; 1 1 0 0 0; 0 1 0 0 1; 1 0 1 1 0]
    num_connections = 4
    A3 = connect_modules(A1, A2, num_connections)
    @test A3[1:4, 1:4] == A1 
    @test A3[5:end, 5:end] == A2

    l1=sum(A3[1:4, 5:end])
    l2=sum(A3[5:end, 1:4])
    @test l1+l2 == num_connections
end