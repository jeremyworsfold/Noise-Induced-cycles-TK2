if isdefined(@__MODULE__, :LanguageServer)
    include("../src/TK2.jl")
    using .TK2
else
    using TK2
end

include("plotting/plotting.jl")
theme = include("plotting/theme.jl")
theme isa Attributes && set_theme!(theme)

################################################################
# Simulation

SCRIPTNAME = "05-relative-prop"
SIMULATE = false
NUM_SIMS = 100
SAVE = true
SAVEFIG = true

expo = -4:0.5:2
lambdas = 10.0 .^ expo

base_config = load(SimParams, SCRIPTNAME)
N = base_config.params.N
κ = base_config.params.kappa
β = base_config.params.beta
DIFFS = 2 .* lambdas ./ N
configs = [@set base_config.T_end = 4e3 / D for D in DIFFS]
configs = [@set conf.params.D = D for (conf, D) in zip(configs, DIFFS)]

y = collect(range(-1, 1, N + 3))

mean_y = []
stacked_prob = []
if SIMULATE
    for conf in configs
        dump(conf)
        data, conf = ensemble_custom_sim(conf, NUM_SIMS)
        append!(mean_y, sum(data["y"] .* y) / N)
        append!(stacked_prob, [data["y"]])
    end
    if SAVE
        save_res_and_conf(
            SCRIPTNAME,
            Dict("lambda" => lambdas, "y" => mean_y, "prob" => stacked_prob),
            configs[1],
        )
    end
else
    res, conf = load_res_and_conf(SCRIPTNAME)
    append!(mean_y, res["y"])
    append!(stacked_prob, res["prob"])
    lambdas = res["lambda"]
end


################################################################
# Theory

kappa = 0:0.001:1
beta = -1:0.001:1

prob_tele(β) = (1 + β) / 2
prob_plus(κ, β) =
    κ == 0 ? prob_tele(β) :
    (κ - 1 + (1 - κ)^((1 - β) / 2) * (1 + κ)^((1 + β) / 2)) / (2 * κ)
expectation_y(κ, β) = 2 * prob_plus(κ, β) - 1
y_star(κ, β) = (β - κ) / (1 - β * κ)

y_diff = reshape(
    [expectation_y(k, b) - y_star(k, b) for k in kappa for b in beta],
    (length(beta), length(kappa)),
)
y_exp = reshape(
    [expectation_y(k, b) for k in kappa for b in beta],
    (length(beta), length(kappa)),
)


################################################################
# Plotting

pwr10_label(a::AbstractString) = rich("10", superscript(a))
xticksfmt = LogTicks(WilkinsonTicks(4, k_min=3))

f = fig_in_cm(10, 17.2)

##### Top Left
begin
    ax1 = Axis(
        f[1, 1],
        xlabel=L"\lambda",
        ylabel=L"\langle y\rangle",
        xscale=log10,
        yticks=[-0.3, -0.2, -0.1],
        xticks=xticksfmt,
    )
    hlines!(ax1, expectation_y(κ, β), color=COLORS[1], label=L"\langle y\rangle_{\textrm{PDMP}}", linewidth=1.5, linestyle=:dot)
    hlines!(
        ax1,
        y_star(κ, β),
        color=COLORS[2],
        linestyle=:dash,
        label=L"y^*",
        linewidth=1.5,
    )
    expo_ = -4:0.02:2
    lambdas_ = 10.0 .^ expo_
    lines!(ax1, lambdas_, (@. (lambdas_ * y_star(κ, β) + 0.5 * expectation_y(κ, β)) / (0.5 + lambdas_)), color=COLORS[3])
    scatter!(ax1, lambdas, mean_y, color=(COLORS[4], 0.5), strokewidth=0.5, strokecolor=:black, markersize=5)
    axislegend(ax1, position=:rt, padding=3, rowgap=-1, patchsize=(12, 12))
    xlims!(ax1, extrema(lambdas)...)
end

##### Bottom Left
begin
    ax2 = Axis(f[2, 1], xlabel=L"y", ylabel=L"P^*(y)", yscale=log10)
    local lbl = pwr10_label.(["-4", "-2", "0", "2"])
    local lnstyl = [:solid, :dash, :dot, :dashdot]

    for (i, prob) in enumerate(stacked_prob[1:4:end])
        last_index = length(prob)
        # smooth out bumps due to precision errors in binning
        vals = [
            (j > 2 && j < last_index) ? (prob[j-1] + 2 * prob[j] + prob[j+1]) / 4 : prob[j] for j in eachindex(prob)
        ]
        lines!(ax2, y, prob, label=lbl[i], linestyle=lnstyl[i])
    end
    xlims!(ax2, -1, 1)
    ylims!(ax2, high=1e3)
    axislegend(
        ax2,
        titlegap=-1,
        padding=0,
        rowgap=2,
        colgap=5,
        patchsize=(10, 10),
        nbanks=2,
    )
end

##### Top Right
begin
    axmain1 = Axis(
        f[1, 2],
        xlabel=L"\kappa",
        ylabel=L"\beta",
        yticks=[-1, 0, 1],
        xticks=xticksfmt,
    )
    local cf =
        contourf!(axmain1, kappa, beta, y_exp', colormap=:PRGn, levels=range(-1, 1, 10))
    scatter!(axmain1, κ, β, marker=:star5, markersize=10, color=:white, strokewidth=0.7, strokecolor=:black,)
    limits!(axmain1, 0, 0.9, -1, 1)
    Colorbar(f[1, 3], cf, label=L"\langle y\rangle_{\textrm{PDMP}}")
end

##### Bottom Right

begin
    axmain2 = Axis(
        f[2, 2],
        xlabel=L"\kappa",
        ylabel=L"\beta",
        yticks=[-1, 0, 1],
        xticks=xticksfmt,
    )
    local cf =
        contourf!(axmain2, kappa, beta, y_diff', colormap=Reverse(:grays), levels=15)
    scatter!(axmain2, κ, β, marker=:star5, markersize=10, color=:white, strokewidth=0.7, strokecolor=:black,)
    Colorbar(f[2, 3], cf, label=L"\langle y\rangle_{\textrm{PDMP}} - y^*")
    linkaxes!(axmain1, axmain2)
end
f

if SAVEFIG
    save_fig(SCRIPTNAME, f)
end
f