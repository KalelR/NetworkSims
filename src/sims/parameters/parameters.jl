
include("topology.jl")

function get_initialconditions(pvals)
    # ics = get(pvals, "ics", nothing)
    # if ics isa Vector return ics end 
    
    @unpack ictype = pvals
    if ictype == "uniform"
        return get_initialconditions_uniform(pvals)
    elseif ictype isa AbstractArray
        return ictype
    end
    @error("No initial condition type $ictype could be found for pvals = $pvals.")
end

function get_initialconditions_uniform(pvals)
    @unpack N, icmin, icmax, icseed, ictype = pvals
    icfilename = "$(datadir())/sims/inputs/ics/N_$N-ictype_$ictype/ics-N_$N-type_$ictype-min_$icmin-max_$icmax-seedic_$icseed.jld2"
    if isfile(icfilename)
        return load(icfilename)["ics"]
    else
        @info("File $icfilename with initial conditions was not found. Generating new initial conditions and saving them.")
        numeqs = get_num_eqs(pvals)
        # u0s = rand(MersenneTwister(icseed), Uniform(icmin, icmax), (numeqs, N))
        u0s = rand(MersenneTwister(icseed), Uniform(icmin, icmax), (numeqs*N))
        safesave(icfilename, Dict("ics"=>u0s))
        return u0s
    end
end

function get_num_eqs(pvals)
    @unpack unitm = pvals
    if unitm == "FHN"
        return 2
    elseif unitm == "izhikevich"
        return 2
    else
        @warn "incorrect unit model"
    end
end



