module GridapFokkerPlanck

using Gridap
using GridapGmsh
using Configurations
using ArgCheck
using Glob

export FokkerPlanckStationary,
    FokkerPlanckTransient,
    to_stationary_config,
    ThetaMethodIntegrator,
    FESpacesSpecs,
    generate,
    AutonomousFP,
    probability_current,
    autonomous_operator,
    crank_nicholson_integrator,
    create_FE_spaces,
    save_stationary_solution,
    save_transient_solution,
    normalized_IC,
    load_config,
    sanitize,
    solveFP,
    _get_script_name


include("AutonomousFP.jl")
include("utils.jl")
include("configs.jl")


end
