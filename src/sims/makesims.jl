function makesims_FHN(pvals, analysisdict)
    odeprob = construct_odeproblem_FHN(pvals)
    details = get_integration_details(pvals)
    solver = get_solver(pvals)
    @unpack ttrans, Δt, tend = pvals
    sol = solve(odeprob, solver, saveat=ttrans:Δt:tend; details...)
    res = analyse_solution(sol, analysisdict, pvals, odeprob)
    return res, odeprob
end

function construct_odeproblem_FHN(pvals)
    u0s = get_initialconditions(pvals)
    params = ParametersFHN(pvals)
    @unpack tend = pvals
    return ODEProblem(bidimensional_synaptic_coupling_rule!, u0s, (0, tend), params)
end

function get_integration_details(pvals; skipinterpolation=false)
    if haskey(pvals, "skipinterpolation") skipinterpolation = pvals["skipinterpolation"] end
    details = NamedTuple()
    @unpack unitm, N = pvals
    if unitm == "FHN"
        cb = VectorContinuousCallback(condition_spike!, affect_spike!, N; affect_neg! = nothing, skipinterpolation, save_positions=(false, false));
        
        shock_times = get(pvals, "shock_times", nothing)
        if !isnothing(shock_times) 
            shock_magnitudes = prepare_shock_magnitudes(pvals)
            affect!(integrator) = shock_affect!(integrator, shock_magnitudes)
            cb_shock, solve_kwargs_shock = prepare_shock_intervention(affect!, shock_times)
            cb = CallbackSet(cb_shock, cb) 
        else 
            solve_kwargs_shock = ()
        end
        
        if !skipinterpolation
            details = (callback=cb, abstol=1e-9, reltol=1e-9, solve_kwargs_shock...)
            # @info "normal, not skipping interpolation"
        else
            @unpack Δt = pvals
            # @info "skipping interpolation"
            details = (callback=cb, adaptive=false, dt=Δt, solve_kwargs_shock...)
        end
    end

    return details
end

function get_solver(pvals)
    @unpack solver = pvals
    if solver == "tsit5"
        return Tsit5()
    elseif solver == "vern9"
        return Vern9()
    elseif solver == "autotsi5trbdf2"
        solver = AutoTsit5(TRBDF2())
    else
        @error "No solver found, $solver"
    end
end

function shock_affect!(integrator, shock_magnitudes)
    integrator.u .+= shock_magnitudes
end 

function prepare_shock_intervention(affect!, shock_times)
    cb = PresetTimeCallback(shock_times, affect!)
    solve_kwargs = (tstops=shock_times, )
    # solve_kwargs = (callback=cb)
    return cb, solve_kwargs
end

function prepare_shock_magnitudes(pvals)
    shock_magnitudes = get(pvals, "shock_magnitudes", nothing)
    if isnothing(shock_magnitudes)
        shock_magnitude_info = get(pvals, "shock_magnitude_info", nothing)
        @unpack N = pvals
        shock_magnitudes = zeros(Float64, 2N); shock_magnitudes_V = view(shock_magnitudes, 1:N); shock_magnitudes_w = view(shock_magnitudes, N+1:2N)
        for (idx_neuron_str, (Vshock, wshock)) in shock_magnitude_info 
            idx_neuron = parse(Int64, idx_neuron_str)
            shock_magnitudes_V[idx_neuron] = Vshock
            shock_magnitudes_w[idx_neuron] = wshock
        end         
    end
    @show shock_magnitudes
    
    return shock_magnitudes
end


# 
# 
# 
# 
# 
# 
# #TODO: check if loading all parameters every time slows code down significantly!; also if reusing some allocated matrices makes sense!
# function makesims(pvals, analysisdict)
#     odeprob = construct_odeproblem(pvals)
#     details = get_integration_details(pvals)
#     solver = get_solver(pvals)
#     @unpack ttrans, Δt, tend = pvals
#     sol = solve(odeprob, solver, saveat=ttrans:Δt:tend; details...)
#     res = analyse_solution(sol, analysisdict, pvals, odeprob)
#     return res
# end
# 
# function construct_odeproblem(pvals)
#     u0s = get_initialconditions(pvals)
#     paramscoup = get_coupling(pvals)
#     paramsunits = get_paramsunits(pvals);
#     paramssystem = ParamsSystem(paramscoup, paramsunits)
#     @unpack tend = pvals
#     return ODEProblem(system!, u0s, (0, tend), paramssystem)
# end
# 
# function get_coupling(pvals)
#     adjl = get_topology(pvals)
#     coup_type = get_coupling_type(pvals)
#     @unpack ϵ = pvals
#     multiplyglobalconstant!(adjl, ϵ)
#     paramscoup = ParamsCoup(adjl, coup_type)
#     return paramscoup
# end
# 
# 
# function get_coupling_type(pvals)
#     @unpack couplingm = pvals;
#     if couplingm == "synapse"
#         @unpack synapsem = pvals;
#         coup_f_sum! = coup_f_sum_syn!
#         coup_f = coup_f_syn
#         if synapsem == "alpha"
#             conductance = alphafunction
#         end
#         @unpack numstoredspiketimes, N, Vth = pvals
#         spiketimes = zeros(N, numstoredspiketimes);
#         spiketime_detection = SpikeTimeDetection(spiketimes, ones(UInt16, N), Vth, numstoredspiketimes)
#         couptype = SynType(coup_f_sum!, coup_f, conductance, spiketime_detection)
#         return couptype
#     else
#         @error "No coupling mode found corresponding to $couplingm"
#     end
# end
# 
# function get_paramsunits(pvals)
#     @unpack N, unitm = pvals
#     if unitm == "FHN"
#         @unpack a, b, c, I = pvals
#         unit_f! = fitzhugh_nagumo_rule!
#         params_units_dyn = Vector{ParamsFHN}(undef, N)
#         params_units_dyn .= ParamsFHN.(a, b, c, I)
#         unit_f! = fitzhugh_nagumo_rule!
#         params_units = ParamsUnit.(unit_f!, params_units_dyn)
#     else
#         @error "Unavailable model."
#     end
#     return params_units
# end



