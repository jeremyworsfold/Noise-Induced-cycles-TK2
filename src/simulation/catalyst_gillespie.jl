mutable struct ModelCatalyst
    rs::ReactionSystem
    pars::Vector{ParamPair}
end


function ModelCatalyst(params::StoichiometryModel)::ModelCatalyst
    @unpack N, kappa, beta, D = params
    k1 = D * (1 + kappa)
    k2 = D * (1 - kappa)
    b1 = D * (1 + beta) * N / 2
    b2 = D * (1 - beta) * N / 2

    rs = @reaction_network begin
        @species X₂(t) X₁(t)
        k1, X₁ --> ∅
        k2, X₂ --> ∅
        b1, ∅ --> X₁
        b2, ∅ --> X₂
        r, X₁ + X₂ --> 2X₂
        r, X₂ + X₁ --> 2X₁
    end

    p = [:r => 1 / N, :b1 => b1, :b2 => b2, :k1 => k1, :k2 => k2]

    return ModelCatalyst(rs, p)
end


function do_sim(model::ModelCatalyst, T_end, u0)
    @unpack rs, pars = model
    jprob = JumpProblem(rs, DiscreteProblem(rs, u0, (0.0, T_end), pars), Direct())
    return solve(jprob, SSAStepper())
end
