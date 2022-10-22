module NetworkSim

using DiffEqCallbacks, OrdinaryDiffEq, LinearAlgebra, DelimitedFiles, Distributions, Random, CSVFiles, DrWatson
export parse_inputs, check_existing_run, makesims, saver, savenames
export SpikeTimeDetection, SynType, CoupType, ParamsSynAlpha, ParamsCoupTypeSp, ParamsFHN, ParamsUnitDynamics, ParamsSystem

include("input-interface.jl");
include("output-interface.jl");
include("savers.jl");
include("sims/systems_generic.jl");
include("sims/systems_specific.jl");
include("sims/parameters.jl");
include("sims/analysis.jl")
include("sims/makesims.jl");

end
