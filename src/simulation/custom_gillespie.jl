
mutable struct ModelCustom
    X::MVector{2,Int}
    t::Float64
    stoich::SMatrix{2,6,Int}
    rates::MVector{6,Float64}
    pars::SVector{5,ParamPair}
end

function ModelCustom(params::StoichiometryModel)::ModelCustom
    @unpack N, kappa, beta, D = params
    k1 = D * (1 + kappa)
    k2 = D * (1 - kappa)
    b1 = D * (1 + beta) * N / 2
    b2 = D * (1 - beta) * N / 2

    stoich = SMatrix{2,6,Int}(hcat([1, 0], [0, 1], [-1, 0], [0, -1], [1, -1], [-1, 1]))

    rates = MVector{6,Float64}(zeros(Float64, 6))
    rates[1] = b1
    rates[2] = b2

    u0 = initial(params)

    return ModelCustom(
        [x[2] for x in u0],
        0.0,
        stoich,
        rates,
        [:r => 1 / N, :b1 => b1, :b2 => b2, :k1 => k1, :k2 => k2],
    )
end


function ensemble_custom_sim(conf::SimParams, nsims::Int)
    @unpack N, D, kappa, beta = conf.params
    ϕ2 = 1 / (1 - kappa)
    # last bin dumps all events out of range
    bins = 0:1:Int(ceil(N * ϕ2 + 5 * sqrt(N)))
    Nmax = length(bins) - 1

    output = Vector{OutputHists}(undef, nsims)

    prog = Progress(nsims)
    @threads for i = 1:nsims
        output[i] = OutputHists(Nmax, N)
        model = ModelCustom(conf.params)
        do_sim!(output[i], model, conf.T_end, show_prog = i == 1)
        next!(prog)
    end

    h = hist_mean(output, N)
    data = results_dict(h, bins ./ N)
    return data, conf
end


function ensemble_custom_sim(scriptname::String, nsims::Int)
    conf = load(SimParams, scriptname)
    ensemble_custom_sim(conf, nsims)
end


function do_sim!(output::OutputHists, model::ModelCustom, t_end::Real; show_prog = false)
    prog = Progress(Int(ceil(t_end)), desc = "simulation on thread 1", enabled = show_prog)
    t = 0.0
    while t < t_end
        t += gill_step!(model, output)
        update!(prog, Int(round(t)))
    end
    normalize!(output, t)
end

function gill_step!(model::ModelCustom, output::OutputHists)
    update_rates!(model)
    s = sum(model.rates)
    Δτ = -log(rand()) / s
    update_output!(output, model.X, Δτ)

    idx = rand(Distributions.Categorical(model.rates ./ s))
    @. model.X = model.X + model.stoich[:, idx]
    return Δτ
end


function update_rates!(model::ModelCustom)
    @unpack X, pars, rates = model
    r = pars[1][2]
    k1 = pars[4][2]
    k2 = pars[5][2]

    rates[3] = k1 * X[1]
    rates[4] = k2 * X[2]
    rates[5] = r * X[1] * X[2]
    rates[6] = rates[5]
end
