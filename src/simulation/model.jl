
"""
    StoichiometryModel

TK2 model parameters

# Fields:
- `N::Float64`: Typical population size
- `D::Float64`: diffusion constant (VK expansion parameter)
- `kappa::Float64`: death rate bias
- `beta::Float64`: birth rate bias
"""
@option "Stoich" struct StoichiometryModel
    N::Int = 100
    D::Float64 = 1e-3
    kappa::Float64 = 0.5
    beta::Float64 = 0.0
end

ParamPair = Pair{Symbol,Float64}

"""
    initial(model::StoichiometryModel)

Get the default initial system configuration as a list of params with their corresponding population sizes. Default to 50/50 split.
"""
function initial(params::StoichiometryModel)::Vector{ParamPair}
    N = params.N
    x1 = Int(floor(N / 2))
    return [:X₁ => x1, :X₂ => N - x1]
end
