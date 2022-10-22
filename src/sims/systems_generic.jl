abstract type ParamsCoupSp end
abstract type ParamsCoupTypeSp end
abstract type CoupType end
abstract type ParamsUnitDynamics end
abstract type ParamsUnitDynamics end

#ParamsCoup: general coupling parameters;paramscouptype: type of coupling made (eg synaptic); paramscouptypesp: parameter values for that coupling type
mutable struct ParamsCoup #general coupling parameters
    adjl :: Vector{Vector{Tuple{Int64, ParamsCoupTypeSp}}} #adjl[i] contains vector with each element a (neighbor, params_coup_sp of connection neighbor->i) of unit i
    coup_type :: CoupType
end

mutable struct ParamsUnit
    unit_f! :: Function
    params_model :: ParamsUnitDynamics
end

mutable struct ParamsSystem
    params_coup :: ParamsCoup
    params_units :: Vector{ParamsUnit}
end

"""
Generic system: unit_f + coup_f
#TODO: check if this is more efficient than directly implementing system!
"""
function system!(du, u, p, t)
    for i = 1:size(du, 2)
        du_view = view(du, 1:2, i)
        params_unit = p.params_units[i]
        params_unit.unit_f!(du_view, u, params_unit.params_model, t)
        p.params_coup.coup_type.coup_f_sum!(i, du_view, u, p.params_coup, t)
    end
    return nothing
end
