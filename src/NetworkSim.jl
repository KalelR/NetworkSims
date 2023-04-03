module NetworkSim

using DiffEqCallbacks, OrdinaryDiffEq, LinearAlgebra, DelimitedFiles,
    Distributions, Random, CSVFiles, DrWatson, Graphs, StatsBase,
    SparseArrays
using Tullio
export parse_inputs, check_existing_run, makesims, saver, savenames
export SpikeTimeDetection, SynType, CoupType, ParamsSynAlpha, ParamsCoupTypeSp, ParamsFHN, ParamsUnitDynamics, ParamsSystem, ParamsCoup, coup_f_sum_syn!, coup_f_syn, alphafunction, fitzhugh_nagumo_rule!, ParamsUnit
export condition_spike, affect_spike!, system!
export ParametersFHN, makesims_FHN
export separatedicts, plotname, spiketimes

include("input-output/input-interface.jl");
include("input-output/output-interface.jl");
include("input-output/savers.jl");
include("sims/systems_generic.jl");
include("sims/systems/FHN.jl");
include("utils.jl");
include("sims/parameters.jl");
include("sims/analysis.jl")
include("sims/makesims.jl");

end
