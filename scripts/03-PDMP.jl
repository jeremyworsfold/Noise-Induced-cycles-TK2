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

SCRIPTNAME = "03-PDMP"
SIMULATE = false
NUM_SIMS = 100
SAVE = true
SAVEFIG = true

if SIMULATE
    data, conf = TK2.ensemble_custom_sim(SCRIPTNAME, NUM_SIMS)
    if SAVE
        save_res_and_conf(SCRIPTNAME, data, conf)
    end
else
    data, conf = load_res_and_conf(SCRIPTNAME)
end

#################################################
# Theory

@unpack N, kappa, D = conf.params
approx = from_params(PDMPStationarySol, kappa)
ϕ1, ϕ2 = TK2.bounds(approx)
x = collect(range(ϕ1 - 3 / sqrt(N), ϕ2 + 5 / sqrt(N), 2000))

Πp, Πm, ϕ = theoretical_pdfs(approx)
p_plus = probability_fraction(approx, ϕ)
ΠpLNA, ΠmLNA = theoretical_LNA_pdfs(approx, x, N)
Π0LNA = ΠpLNA .+ ΠmLNA

#################################################
# Plotting

# STATIONARY DISTRIBUTION FOR X+,X- WITH LNA 
f1 = fig1()
ax = Axis(f1[1, 1], xlabel = L"\phi", ylabel = L"\Pi^\ast_\pm(\phi)")
manual_density!(ax, data["edges"], data["x1"], fillcolor = COLORS[1], label=L"+")
manual_density!(ax, data["edges"], data["x2"], fillcolor = COLORS[2], label=L"-")
plot_PDMP!(ax, ϕ, Πp, ϕ1, var = :x1, label=L"+")
plot_PDMP!(ax, ϕ, Πm, ϕ2, var = :x2, label=L"-")
axislegend(ax, 
    position = :ct,
    padding = 0,
    rowgap = 2,
    patchsize = (10, 10), patchlabelgap=5, merge=true
)
xlims!(ax, extrema(x))
ylims!(ax, 0, maximum(Π0LNA) * 1.05)
f1


# STATIONARY PROB OF X+ OR X- GIVEN A PARTICULAR s
indices = (ϕ2 + 0.15 .> data["edges"][1:end-1] .> ϕ1 - 0.15)
xvals = data["edges"][1:end-1][indices]
frac = data["x1"][indices] ./ data["s"][indices]

f2 = fig1()
ax = Axis(f2[1, 1], xlabel = L"\phi", ylabel = L"\Pi^*(\pm|\phi)")
inset_ax = inset_axis(f2[1, 1], xlabel = L"\phi", ylabel = L"\Pi^*_0(\phi)")

fill_between!(ax, xvals, 0, frac, color = (COLORS[1], 0.5))
fill_between!(ax, xvals, frac, 1, color = (COLORS[2], 0.5))
lines!(ax, ϕ, p_plus, color = :black, linestyle = :dash)
xlims!(ax, TK2.bounds(approx) .+ (-0.1, 0.1))
ylims!(ax, 0, 1)

manual_hist!(inset_ax, data["edges"], data["s"], color = COLORS[4])
plot_PDMP!(inset_ax, ϕ, Πp .+ Πm, [ϕ1, ϕ2], var = :s, linewidth = 0.6)
xlims!(ax, TK2.bounds(approx) .+ (-0.1, 0.1))
limits!(inset_ax, ϕ1 - 0.2, ϕ2 + 0.3, 0, maximum(Π0LNA) * 1.1)
f2

if SAVEFIG
    save_fig_and_conf(SCRIPTNAME, f1, conf)
    save_fig_and_conf("$SCRIPTNAME-prop", f2, conf)
end
