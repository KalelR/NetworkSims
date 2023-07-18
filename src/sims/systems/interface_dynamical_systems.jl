using Attractors 

function get_ds_synaptic(pvals; abstol=1e-9, reltol=1e-9, kwargs...)
    u0 = get_initialconditions(pvals)
    p = ParametersFHN(pvals)
    ode_function = bidimensional_synaptic_coupling_rule!
    solver = get_solver(pvals)
    details = get_integration_details(pvals)
    diffeq = (alg = solver, abstol, reltol, details...)
    ds = CoupledODEs(ode_function, u0, p; diffeq)
end

function get_ds_diffusive(pvals; abstol=1e-9, reltol=1e-9, kwargs...)
    u0 = get_initialconditions(pvals)
    p = ParametersFHNDiffusive(pvals)
    ode_function = bidimensional_diffusive_coupling_rule!
    solver = get_solver(pvals)
    diffeq = (alg = solver, abstol, reltol)
    ds = CoupledODEs(ode_function, u0, p; diffeq)
end

function get_ds(pvals; kwargs...)
    coup_type = get(pvals, "coup_type", "synaptic")
    if coup_type == "synaptic"
        return get_ds_synaptic(pvals; kwargs...)
    elseif coup_type == "diffusive"
        return get_ds_diffusive(pvals; kwargs...)
    else
        @error "no found coup type"
    end
end


function makesims_FHN_ds(pvals, analysisdict)
    @unpack tend, ttrans, Δt = pvals
    T = tend-ttrans; Ttr = ttrans
    ds = get_ds(pvals)
    traj, ts = trajectory(ds, T; Ttr, Δt)
    sol = Matrix(transpose(Matrix(traj)))
    # res = analyse_solution(sol, analysisdict, pvals, odeprob)
    res = Dict("sol"=>sol, "ts"=>ts)
    odeprob = nothing
    return res, odeprob
end

function sync_coordinates(att::StateSpaceSet)
    StateSpaceSet(sync_coordinates(Matrix(att)); warn=false)
end

function reinit_synaptic!(ds::CoupledODEs)
    reinit_synaptic!(ds.integ.sol.prob)
end
function SciMLBase.reinit!(ds::Attractors.ContinuousTimeDynamicalSystem, u = initial_state(ds);
        p = current_parameters(ds), t0 = initial_time(ds)
    )
    if current_parameters(ds) isa ParametersFHN reinit_synaptic!(ds) end 
    isnothing(u) && return
    set_parameters!(ds, p)
    reinit!(ds.integ, u; reset_dt = true, t0)
    return ds
end