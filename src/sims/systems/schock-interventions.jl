function shock_affect!(integrator)
    shocks = integrator.p.shocks 
    integrator.u .+= shocks 
end 

function prepare_shock_intervention(affect!, shock_times)
    cb = PresetTimeCallback(shock_times, affect!)
    solve_kwargs = (callback=cb, tstops=shock_times)
end