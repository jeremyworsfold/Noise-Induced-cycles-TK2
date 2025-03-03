if isdefined(@__MODULE__, :LanguageServer)
    include("../src/TK2.jl")
    using .TK2
else
    using TK2
end

include("plotting/plotting.jl")
theme = include("plotting/theme.jl")
theme isa Attributes && set_theme!(theme)
CairoMakie.activate!()

SIMULATE = true
NUM_SIMS = 100
SAVE = true
SAVEFIG = true

#################################################

scriptnames = ["04-LNA-sym", "04-LNA-asym"]
f1 = fig_in_cm(9, 8.6)
axs = [
    Axis(f1[1, 1], ylabel = L"P^\ast_\pm(s)", xticklabelsvisible = false),
    Axis(f1[2, 1], xlabel = L"s", ylabel = L"P^\ast_\pm(s)"),
]
linkxaxes!(axs...)
rowgap!(f1.layout, Relative(0.02))

for (ax, sname) in zip(axs, scriptnames)
    if SIMULATE
        data, conf = ensemble_custom_sim(sname, NUM_SIMS)
        if SAVE
            save_res_and_conf(sname, data, conf)
        end
    else
        data, conf = load_res_and_conf(sname)
    end

    #################################################
    @unpack N, kappa, beta = conf.params
    approx = from_params(PDMPAsymFirstOrder, kappa, beta)

    ϕ1, ϕ2 = TK2.bounds(approx)
    x = collect(range(ϕ1 - 3 / sqrt(N), ϕ2 + 4 / sqrt(N), 2000))

    Πp, Πm, ϕ = theoretical_pdfs(approx)

    ΠpLNA, ΠmLNA = theoretical_LNA_pdfs(approx, x, N)
    Π0LNA = ΠpLNA .+ ΠmLNA

    manual_hist!(ax, data["edges"], data["x1"]; color = COLORS[1])
    manual_hist!(ax, data["edges"], data["x2"]; color = COLORS[2])
    plot_theory!(ax, x, ΠpLNA, var = :x1)
    plot_theory!(ax, x, ΠmLNA, var = :x2)
    xlims!(ax, extrema(x))
    ylims!(ax, 0, max(maximum(ΠpLNA), maximum(ΠmLNA)) * 1.05)
end

f1

if SAVEFIG save_fig("asym-sym-LNA-comp", f1) end
