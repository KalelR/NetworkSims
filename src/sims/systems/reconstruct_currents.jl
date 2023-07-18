
function _coupling_current(u, p, t)
    N = floor(Int, length(u)/2)
    Vs = view(u, 1:N)
    ns = view(u, N+1:2N)
    @unpack coup_params = p
    diffusive_coupling_all_prealloc!(u, t, Vs, ns, coup_params);
    Icoupx = get_tmp(coup_params.Icoupx, u)
    Icoupy = get_tmp(coup_params.Icoupy, u)
    return hcat(Icoupx..., Icoupy...)
end

function coupling_current(traj, p, ts)
    N = floor(Int, size(traj, 2)/2)
    Icoups = zeros(Float64, (length(ts), 2N))
    for idx_t in eachindex(ts)
        t = ts[idx_t]
        u = traj[idx_t, :]
        Icoup = _coupling_current(u, p, t)
        Icoups[idx_t, :] = Icoup
    end 
    return Icoups 
end

function coupling_current(atts::Dict{A, B}, p, ts) where {A,B}
    Icoups = Dict{A, Matrix{Float64}}()
    for (k, att) in atts 
        Icoups[k] = coupling_current(att, p, ts)
    end 
    return Icoups
end

function dynamics_current!(du, u, p, t)
    @unpack I, gl, El, gna, Ena, gk, Ek, τ, Vhm, km, Vhn, kn, coup_params = p
    C = 1
    N = floor(Int, length(u)/2)
    dVs = view(du, 1:N)
    dns = view(du, N+1:2N)
    Vs = view(u, 1:N)
    ns = view(u, N+1:2N)
    minf = get_tmp(coup_params.minf, u)
    ninf = get_tmp(coup_params.ninf, u)
    
    @. minf = activationfunc(Vs, Vhm, km)
    @. ninf = activationfunc(Vs, Vhn, kn)
    @. dVs = (I - gl*(Vs- El) - gna*minf*(Vs - Ena) - gk * ns * (Vs - Ek) ) / C
    @. dns = (ninf - ns) / τ  
    return nothing
end
    
function dynamics_current(traj, p, ts)
    N = floor(Int, size(traj, 2)/2)
    Icoups = zeros(Float64, (length(ts), 2N))
    d_current = zeros(Float64, size(traj))
    for idx_t in eachindex(ts)
        t = ts[idx_t]
        u = traj[idx_t, :]
        du = view(d_current,idx_t, :)
        dynamics_current!(du, u, p, t)
    end 
    return d_current 
end

function dynamics_current(atts::Dict{A, B}, p, ts) where {A,B}
    Ids = Dict{A, Matrix{Float64}}()
    for (k, att) in atts 
        Ids[k] = dynamics_current(att, p, ts)
    end 
    return Ids
end