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
SIMULATE_TS = false
NUM_SIMS = 100
SAVE = false
SAVEFIG = true

KAPPAS = [0.0, 0.5]
ts_Tend = 2e3

base_config = load(SimParams, SCRIPTNAME)
configs = [@set base_config.params.kappa = κ for κ in KAPPAS]
fnames = ["$SCRIPTNAME-sym", "$SCRIPTNAME-asym"]

# Simulate for stationary state
data_ss = []
for (conf, fname) in zip(configs, fnames)
    if SIMULATE
        data, conf = ensemble_custom_sim(conf, NUM_SIMS)
        if SAVE
            save_res_and_conf(fname, data, conf)
        end
    else
        data, conf = load_res_and_conf(fname)
    end
    append!(data_ss, [data])
end

data_ts = []
if SIMULATE_TS
    # Simulate shorter timeseries
    for (conf, fname) in zip(configs, fnames)
        @info conf.params.kappa, conf.params.beta
        @unpack N, D, kappa = conf.params
        model = TK2.ModelCatalyst(conf.params)

        sol = TK2.do_sim(model, ts_Tend, TK2.initial(conf.params))
        d, t = get_variables(sol, N)

        dat = Dict("t" => t, "x1" => d[:x1].vals, "x2" => d[:x2].vals)
        append!(data_ts, [dat])
        if SAVE
            save_res_and_conf("$fname-ts", dat, conf)
        end
    end
else
    data_ts = [load_res_and_conf("$fname-ts")[1] for fname in fnames]
end


#########################################################
# PLOTTING 

CairoMakie.activate!()
f = fig_in_cm(13, 12.2)
# ts_grid1 = f[1, 1] = GridLayout()
# ts_grid2 = f[2, 1] = GridLayout()
# axtop1 = Axis(ts_grid1[1, 1], height=PIXEL_SCALE)
# axtop2 = Axis(ts_grid2[1, 1], height=PIXEL_SCALE)
# axtop1 = Axis(f[1, 1], height=0.9 * PIXEL_SCALE)
# axtop2 = Axis(f[2, 1], height=0.9 * PIXEL_SCALE)
ts_axes = [
    Axis(
        f[1, 1],
        # height=PIXEL_SCALE,
        # yaxisposition=:right,
        # yticklabelsvisible=false,
        xticks=LogTicks(WilkinsonTicks(3, k_min=2)),
        xlabel=L"t",
        ylabel=L"s",
    ),
    Axis(
        f[2, 1],
        # height=0.9 * PIXEL_SCALE,
        # yaxisposition=:right,
        # yticklabelsvisible=false,
        xticks=LogTicks(WilkinsonTicks(3, k_min=2)),
        xlabel=L"t",
        ylabel=L"s",
    ),
]

# hidedecorations!(axtop1)
# hidespines!(axtop1)
# rowgap!(ts_grid1, 1, Relative(0.03))
# hidedecorations!(axtop2)
# hidespines!(axtop2)
# rowgap!(ts_grid2, 1, Relative(0.03))

# text!(
#     axtop1,
#     Point3f(0.5, 0.5, 0),
#     text="Fluctuating Population",
#     space=:relative,
#     align=(:center, :center),
# )
# text!(
#     axtop2,
#     Point3f(0.5, 0.5, 0),
#     text="Population growth and decay",
#     space=:relative,
#     align=(:center, :center),
# )
linkxaxes!(ts_axes...)


bounds = [(0.5, 1.5), (1 / (1 + 0.5) - 0.2, 1 / (1 - 0.5) + 0.2)]

for (i, (conf, ts_ax, dat_ts, dat_ss, bs)) in
    enumerate(zip(configs, ts_axes, data_ts, data_ss, bounds))
    # Axes
    # gr = f[i, 2] = GridLayout()
    # axmain = Axis(gr[1, 1], xlabel=L"y", ylabel=L"s", xticks=[-1, 0, 1])
    # axtop = Axis(gr[1, 1], height=0.8 * PIXEL_SCALE)
    # axright = Axis(gr[2, 2], width=0.8 * PIXEL_SCALE)
    axmain = Axis(f[i, 2], xticks=[-1, 0, 1], yticklabelsvisible=false)
    axright1 = Axis(f[i, 3], width=0.8 * PIXEL_SCALE)
    axright2 = Axis(f[i, 4], width=0.8 * PIXEL_SCALE)
    # axright1 = Axis(f[i, 4], ylabel=L"s", xticks=[-1, 0, 1], width=0.8 * PIXEL_SCALE)
    # hideydecorations!(axmain)
    linkyaxes!(axmain, axright1, axright2)
    # linkxaxes!(axmain, axtop)

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

    func(pt) = Point2f(dy(pt[2], pt[1], β, κ), ds(pt[2], pt[1], β, κ))
    # Vector field
    band!(axmain, yvals, 0.0, s_nullcine.(yvals, β, κ), color=(:black, 0.2))
    # lines!(axmain, yvals, s_nullcine.(yvals, β, κ), color = (:black, 0.4))
    # arrows!(axmain, yvals, svals, dy_vals, ds_vals, normalize = true)
    streamplot!(axmain, func, yvals, svals, arrow_size=6, color=x -> 1, colormap=["#333", "#444"], density=0.45, linewidth=0.7, colorscale=Makie.pseudolog10, quality=100, gridsize=(35, 35))
    scatter!(axmain, y_star, s_star, color=COLORS[4], markersize=10, strokewidth=0.7, strokecolor=:black)

    # y marginal histogram
    # manual_density!(
    #     axtop,
    #     range(-1, 1, length(dat_ss["y"]) - 1),
    #     dat_ss["y"][1:end-2],
    #     color=(COLORS[3], 0.5),
    # )
    # lines!(
    #     axtop,
    #     range(-1, 1, length(dat_ss["y"]) - 2),
    #     dat_ss["y"][1:end-2],
    #     color=COLORS[3],
    # )
    # ylims!(axtop, 0, 1.2)
    xlims!(axmain, -1, 1)

    # s marginal histogram
    x = dat_ss["edges"][1:end-1]
    x1 = vec(sum(dat_ss["x_prime"][(Int(length(dat_ss["y"]) ÷ 2)):end, :], dims=1))
    x2 = vec(sum(dat_ss["x_prime"][1:(Int(length(dat_ss["y"]) ÷ 2)), :], dims=1))
    lower = Point2f.(dat_ss["x1"], 0.0)
    upper = Point2f.(dat_ss["x1"], x)
    lower = Point2f.(x1, 0.0)
    upper = Point2f.(x1, x)
    band!(axright1, lower, upper, color=(COLORS[1], 0.5))
    lines!(axright1, upper, label="s", color=COLORS[1])
    ylims!(axright1, bs...)
    lower = Point2f.(dat_ss["x2"], 0.0)
    upper = Point2f.(dat_ss["x2"], x)
    lower = Point2f.(x2, 0.0)
    upper = Point2f.(x2, x)
    band!(axright2, lower, upper, color=(COLORS[2], 0.5))
    lines!(axright2, upper, label="s", color=COLORS[2])
    ylims!(axright2, bs...)
    xlims!(
        axright1,
        low=0.0,
        # high=maximum(vcat(dat_ss["x1"][10:end], dat_ss["x2"][10:end])) / 2,
    )
    xlims!(
        axright2,
        low=0.0,
        # high=maximum(vcat(dat_ss["x1"][10:end], dat_ss["x2"][10:end])) / 2,
    )

    # formatting
    # hidedecorations!(axtop)
    # hidespines!(axtop, :t, :r, :l)
    hidedecorations!(axright1)
    hidespines!(axright1, :t, :r, :b)
    hidedecorations!(axright2)
    hidespines!(axright2, :t, :r, :b)
    # colgap!(gr, 1, Relative(0.03))
    # rowgap!(gr, 1, Relative(0.03))


    t = dat_ts["t"]
    x1 = TimeSol(dat_ts["x1"] ./ conf.params.N, L"x_+")
    x2 = TimeSol(dat_ts["x2"] ./ conf.params.N, L"x_-")
    stacked!(ts_ax, [x1, x2], t, linecolor=RGBf(0.45, 0.45, 0.45))
    # ts_ax.tellheight = true
    xlims!(ts_ax, 0, ts_Tend)
    ylims!(ts_ax, bs...)
    linkyaxes!(ts_ax, axmain)
end

f

if SAVEFIG
    save_fig_and_conf("$SCRIPTNAME-ts-flows", f, configs[2])
end




###### 3D histograms


using GLMakie
GLMakie.activate!()

figs = [fig_in_cm(12, 10) for i = 1:2]
yticks = [[0.5, 1.0], [1.0, 2.0]]

for (f, conf, data, ytcks) in zip(figs, configs, data_ss, yticks)
    ax = Axis3(
        f[1, 1],
        xlabel=L"y",
        ylabel=L"s",
        xgridvisible=false,
        ygridvisible=false,
        zgridvisible=false,
        xticks=[-1, 0.0, 1.0],
        yticks=ytcks,
    )

    edges_y = collect(range(-1, 1, size(data["x_prime"])[1]))
    colors = reshape(
        repeat(tanh.(edges_y), inner=length(data["edges"])),
        (length(data["edges"]), length(edges_y)),
    )

    srfc = surface!(
        ax,
        edges_y,
        data["edges"],
        data["x_prime"],
        color=colors',
        colormap=[(COLORS[2], 0.5), :white, (COLORS[1], 0.5)],
        shading=NoShading,
    )

    hidespines!(ax)
    hidezdecorations!(ax)
end

figs[1]


if SAVEFIG
    save_fig_and_conf("$SCRIPTNAME-sym", figs[1], configs[1], fmt="png", px_per_unit=4)
    save_fig_and_conf("$SCRIPTNAME-asym", figs[2], configs[2], fmt="png", px_per_unit=4)
end