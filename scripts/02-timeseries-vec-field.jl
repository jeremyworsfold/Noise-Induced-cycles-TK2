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

dy(s, y, β, κ) = (β - y) / s - κ * (1 - y^2)
ds(s, y, β, κ) = 1 - s * (1 + κ * y)
s_nullcine(y, β, κ) = 1 / (1 + y * κ)

SCRIPTNAME = "02-timeseries"
SIMULATE = true
NUM_SIMS = 100
SAVE = true
SAVEFIG = true

kappas = [0.0, 0.5]
ts_Tend = 1e4

base_config = load(SimParams, SCRIPTNAME)
configs = [@set base_config.params.kappa=κ for κ in kappas]
# configs = [SimParams(1e7, StoichiometryModel(200, 0.0004, κ, 0)) for κ in kappas]
fnames = ["$SCRIPTNAME-sym", "$SCRIPTNAME-asym"]

# Simulate for stationary state
data_ss = []
for (conf, fname) in zip(configs, fnames)
    if SIMULATE
        data, conf = ensemble_custom_sim(conf, NUM_SIMS)
        if SAVE save_res_and_conf(fname, data, conf) end
    else
        data, conf = load_res_and_conf(fname)
    end
    append!(data_ss, [data])
end

# Simulate shorter timeseries
ts_Tend = 1e4
configs_ts = [@set conf T_end=ts_Tend for conf in configs]
data_ts = []
for conf in configs
    @info conf.params.kappa, conf.params.beta
    @unpack N, D, kappa = conf.params
    model = TK2.ModelCatalyst(conf.params)

    sol = TK2.do_sim(model, ts_Tend, TK2.initial(conf.params))
    d, t = get_variables(sol, N)

    append!(data_ts, [Dict("t" => t, "x1" => d[:x1].vals, "x2" => d[:x2].vals)])
end

#########################################################
# PLOTTING 

CairoMakie.activate!()
f = fig_in_cm(13, 10.2)
ts_grid1 = f[1, 1] = GridLayout()
ts_grid2 = f[2, 1] = GridLayout()
axtop1 = Axis(ts_grid1[1, 1], height = 0.9 * PIXEL_SCALE)
axtop2 = Axis(ts_grid2[1, 1], height = 0.9 * PIXEL_SCALE)
ts_axes = [
    Axis(
        ts_grid1[2, 1],
        yaxisposition = :right,
        yticklabelsvisible = false,
        xticks = Int.(floor.([0, ts_Tend / 2, ts_Tend])),
        xlabel = L"t",
        xtickformat = values -> [power10fmt(v) for v in values],
    ),
    Axis(
        ts_grid2[2, 1],
        yaxisposition = :right,
        yticklabelsvisible = false,
        xticks = [0, ts_Tend / 2, ts_Tend],
        xlabel = L"t",
        xtickformat = values -> [power10fmt(v) for v in values],
    ),
]

hidedecorations!(axtop1)
hidespines!(axtop1)
rowgap!(ts_grid1, 1, Relative(0.03))
hidedecorations!(axtop2)
hidespines!(axtop2)
rowgap!(ts_grid2, 1, Relative(0.03))

text!(
    axtop1,
    Point3f(0.5, 0.5, 0),
    text = "Noise-induced bistability",
    space = :relative,
    align = (:center, :center),
)
text!(
    axtop2,
    Point3f(0.5, 0.5, 0),
    text = "Noise-induced cycles",
    space = :relative,
    align = (:center, :center),
)
linkxaxes!(ts_axes...)


bounds = [(0.5, 1.5), (1 / (1 + 0.5) - 0.2, 1 / (1 - 0.5) + 0.2)]

for (i, (conf, ts_ax, dat_ts, dat_ss, bs)) in
    enumerate(zip(configs, ts_axes, data_ts, data_ss, bounds))
    # Axes
    gr = f[i, 2] = GridLayout()
    axmain = Axis(gr[2, 1], xlabel = L"y", ylabel = L"s", xticks = [-1, 0, 1])
    axtop = Axis(gr[1, 1], height = 0.8 * PIXEL_SCALE)
    axright = Axis(gr[2, 2], width = 0.8 * PIXEL_SCALE)
    linkyaxes!(axmain, axright)
    linkxaxes!(axmain, axtop)

    # Calculate vector fields
    β = conf.params.beta
    κ = conf.params.kappa

    y_star = (κ - β) / (β * κ - 1)
    s_star = (β * κ - 1) / (κ^2 - 1)

    N = 10
    svals = LinRange(bs..., N)
    yvals = LinRange(-1, 1, N)
    dy_vals = [dy(s, y, β, κ) for y in yvals, s in svals]
    ds_vals = [ds(s, y, β, κ) for y in yvals, s in svals]


    # Vector field
    band!(axmain, yvals, 0.0, s_nullcine.(yvals, β, κ), color = (:black, 0.2))
    lines!(axmain, yvals, s_nullcine.(yvals, β, κ), color = (:black, 0.4))
    arrows!(axmain, yvals, svals, dy_vals, ds_vals, normalize = false)
    scatter!(
        axmain,
        y_star,
        s_star,
        markersize = 10
    )

    # y marginal histogram
    manual_density!(
        axtop,
        range(-1, 1, length(dat_ss["y"]) - 1),
        dat_ss["y"][1:end-2],
        color = (COLORS[3], 0.5),
    )
    lines!(
        axtop,
        range(-1, 1, length(dat_ss["y"]) - 2),
        dat_ss["y"][1:end-2],
        color = COLORS[3],
    )
    ylims!(axtop, 0, 1.2)
    xlims!(axtop, -1, 1)

    # s marginal histogram
    x = dat_ss["edges"][1:end-1]
    lower = Point2f.(dat_ss["x1"], 0.0)
    upper = Point2f.(dat_ss["x1"], x)
    band!(axright, lower, upper, color = (COLORS[1], 0.5))
    lines!(axright, upper, label = "s", color = COLORS[1])
    lower = Point2f.(dat_ss["x2"], 0.0)
    upper = Point2f.(dat_ss["x2"], x)
    band!(axright, lower, upper, color = (COLORS[2], 0.5))
    lines!(axright, upper, label = "s", color = COLORS[2])
    ylims!(axright, bs...)
    xlims!(
        axright,
        low = 0.0,
        high = maximum(vcat(dat_ss["x1"][10:end], dat_ss["x2"][10:end])),
    )

    # formatting
    hidedecorations!(axtop)
    hidespines!(axtop, :t, :r, :l)
    hidedecorations!(axright)
    hidespines!(axright, :t, :r, :b)
    colgap!(gr, 1, Relative(0.03))
    rowgap!(gr, 1, Relative(0.03))


    t = dat_ts["t"]
    x1 = TimeSol(dat_ts["x1"] ./ conf.params.N, L"x_+")
    x2 = TimeSol(dat_ts["x2"] ./ conf.params.N, L"x_-")
    stacked!(ts_ax, [x1, x2], t, linecolor = RGBf(0.45, 0.45, 0.45))
    ts_ax.tellheight = true
    xlims!(ts_ax, 0, ts_Tend)
    linkyaxes!(ts_ax, axmain)
end

rowgap!(f.layout, 1, Relative(0.02))
colgap!(f.layout, 1, Relative(0.02))
f

if SAVEFIG
    save_fig_and_conf("vec_flow_with_timeseries-2", f, configs[2])
end




###### 3D histograms


using GLMakie
GLMakie.activate!()

figs = [fig_in_cm(12, 10) for i = 1:2]
yticks = [[0.5, 1.0], [1.0, 2.0]]

for (f, conf, data, ytcks) in zip(figs, configs, data_ss, yticks)
    ax = Axis3(
        f[1, 1],
        xlabel = L"y",
        ylabel = L"s",
        xgridvisible = false,
        ygridvisible = false,
        zgridvisible = false,
        xticks = [-1, 0.0, 1.0],
        yticks = ytcks,
    )

    edges_y = collect(range(-1, 1, size(data["x_prime"])[1]))
    colors = reshape(
        repeat(tanh.(edges_y), inner = length(data["edges"])),
        (length(data["edges"]), length(edges_y)),
    )

    srfc = surface!(
        ax,
        edges_y,
        data["edges"],
        data["x_prime"],
        color = colors',
        colormap = [(COLORS[2], 0.5), :white, (COLORS[1], 0.5)],
        shading = NoShading,
    )

    hidespines!(ax)
    hidezdecorations!(ax)
end

figs[2]


if SAVEFIG
    save_fig_and_conf("3d-sym-bistability", figs[1], conf, fmt = "png", px_per_unit = 4)
    save_fig_and_conf("3d-asym-cycles", figs[2], conf, fmt = "png", px_per_unit = 4)
end
