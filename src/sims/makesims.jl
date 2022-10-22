
#TODO: check if loading all parameters every time slows code down significantly!; also if reusing some allocated matrices makes sense!
function makesims(pvals, analysisdict)
    odeprob = construct_odeproblem(pvals)
    details = get_integration_details(pvals)
    solver = get_solver(pvals)
    @unpack ttrans, Δt, tend = pvals
    sol = solve(odeprob, solver, saveat=ttrans:Δt:tend; details...)
    res = analyse_solution(sol, analysisdict, pvals, odeprob)
    return res
end

function construct_odeproblem(pvals)
    u0s = get_initialconditions(pvals)
    paramscoup = get_coupling(pvals)
    println(paramscoup)
    paramsunits = get_paramsunits(pvals);
    paramssystem = ParamsSystem(paramscoup, paramsunits)
    @unpack tend = pvals
    return ODEProblem(system!, u0s, (0, tend), paramssystem)
end

function get_coupling(pvals)
    adjl = get_topology(pvals)
    coup_type = get_coupling_type(pvals)
    @unpack ϵ = pvals
    multiplyglobalconstant!(adjl, ϵ)
    paramscoup = ParamsCoup(adjl, coup_type)
    return paramscoup
end

function multiplyglobalconstant!(adjl::Vector{Vector{Any}}, ϵ)
    for adjl_unit in adjl
        for (neighbor, params_conn) in adjl_unit
            params_conn.gmax *= ϵ
        end
    end
end

function get_coupling_type(pvals)
    @unpack couplingm = pvals;
    if couplingm == "synapse"
        @unpack synapsem = pvals;
        coup_f_sum! = coup_f_sum_syn!
        coup_f = coup_f_syn
        if synapsem == "alpha"
            conductance = alphafunction
        end
        @unpack numstoredspiketimes, N, Vth = pvals
        spiketimes = zeros(N, numstoredspiketimes);
        spiketime_detection = SpikeTimeDetection(spiketimes, ones(UInt16, N), Vth, numstoredspiketimes)
        couptype = SynType(coup_f_sum!, coup_f, conductance, spiketime_detection)
        return couptype
    else
        @error "No coupling mode found corresponding to $couplingm"
    end
end

function get_paramsunits(pvals)
    @unpack N, unitm = pvals
    if unitm == "FHN"
        @unpack a, b, c, I = pvals
        unit_f! = fitzhugh_nagumo_rule!
        params_units_dyn = Vector{ParamsFHN}(undef, N)
        params_units_dyn .= ParamsFHN.(a, b, c, I)
        unit_f! = fitzhugh_nagumo_rule!
        params_units = ParamsUnit.(unit_f!, params_units_dyn)
    else
        @error "Unavailable model."
    end
    return params_units
end

function get_integration_details(pvals)
    details = NamedTuple()
    @unpack unitm, N = pvals
    if unitm == "FHN"
        cb = VectorContinuousCallback(condition_spike,affect_spike!, N; affect_neg! = nothing);
        details = (callback=cb, )
    end

    return details
end

function get_solver(pvals)
    @unpack solver = pvals
    if solver == "tsit5"
        return Tsit5()
    elseif solver == "vern9"
        return Vern9()
    else
        @error "No solver found, $solver"
    end
end