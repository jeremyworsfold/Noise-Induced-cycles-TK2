using CairoMakie
using LaTeXStrings
using ColorSchemes
using Parameters
using Configurations
using Format
using Parameters
using Setfield


COLORS = colorschemes[:seaborn_muted].colors

COLOR_DICT = Dict(:x1 => COLORS[1], :x2 => COLORS[2], :s => COLORS[4])
LABEL_P_DICT = Dict(:x1 => L"\Pi^*_+(s)", :x2 => L"\Pi^*_-(s)", :s => L"\Pi^*(s)")
LABEL_X_DICT = Dict(:x1 => L"+", :x2 => L"-", :s => L"x_1+x_2")

PIXEL_SCALE = 28.3465

"""
    fig_in_cm(height::Real, width::Real; padding, args)

Create a figure with correct resolution for adding to documents

# Arguments:
- `height`: height in cm of figure
- `width`: width in cm of figure
- `padding`: padding to ensure axes fit
"""
fig_in_cm(height::Real, width::Real; padding = (2.5, 5, 2.5, 5), args...) =
    Figure(size = (width * PIXEL_SCALE, height * PIXEL_SCALE), figure_padding = padding)

fig1(; args...) = fig_in_cm(6, 8.6; args...)

fig2(; args...) = fig_in_cm(12, 8.6; args...)

function inset_axis(f::Union{GridLayout,GridPosition}; args...)
    ax = Axis(
        f,
        width = Relative(0.4),
        height = Relative(0.4),
        halign = 0.95,
        valign = 0.95,
        backgroundcolor = :white;
        args...,
    )
    translate!(ax.scene, 0, 0, 10)
    translate!(ax.elements[:background], 0, 0, 9)
    return ax
end

##############################################################################

function plot_theory!(ax, x, y; var = :x1, args...)
    return lines!(ax, x, y; color = COLOR_DICT[var], label = LABEL_P_DICT[var], args...)
end

function plot_PDMP!(ax, x, y, bound; var = :x1, args...)
    vlines!(ax, bound, color = COLOR_DICT[var], linestyle = :dash; args...)
    return lines!(
        ax,
        x,
        y;
        color = COLOR_DICT[var],
        label = LABEL_P_DICT[var],
        linestyle = :dash,
        args...,
    )
end

function add_rect_legend(scene; vars = (:x1, :x2), title = nothing, args...)
    group_color = [
        MarkerElement(marker = :rect, markersize = 10, color = (COLOR_DICT[v], 0.5)) for
        v in vars
    ]
    if !isnothing(title)
        return Legend(scene, group_color, [LABEL_X_DICT[v] for v in vars], title; args...)
    end
    return Legend(scene, group_color, [LABEL_X_DICT[v] for v in vars]; args...)
end

function save_fig(name::AbstractString, fig::Figure; path = "figures", fmt = "pdf", args...)
    CairoMakie.save("$path/$name.$fmt", fig; pt_per_unit = 1, args...)
end

function save_fig_and_conf(
    name::AbstractString,
    fig::Figure,
    conf::SimParams;
    path = "figures",
    fmt = "pdf",
    args...,
)
    save_fig(name, fig; path = path, pt_per_unit = 1, fmt = fmt, args...)
    to_toml("$path/$name.toml", conf)
end


##############################################################################

function hist_pts(edges, heights)
    @assert length(edges) == length(heights) + 1
    xs = zeros(Float64, 2 * length(edges))
    ys = similar(xs)
    xs[1] = edges[1]
    for (i, y) in enumerate(heights)
        idx = 2 * i
        xs[idx:idx+1] = edges[i:i+1]
        ys[idx:idx+1] = [y, y]
    end
    xs[end] = edges[end]
    return xs, ys
end


"""
    manual_hist!(ax::Axis, edges, heights, color = COLORS[1]; args...)

Plot a histogram with either an outline or a fill color. 
Does not show internal lines like the standard.

# Arguments:
- `ax`: axis for plotting
- `edges`: edges used to plot vertical lines
- `heights`: heights of the histogram
- `fillcolor`: color of histogram with α=0.5. Default is `COLORS[1]`
"""
function manual_hist!(ax::Axis, edges, heights; color = COLORS[1], args...)
    xs, ys = hist_pts(edges, heights)
    fill_between!(ax, xs, 0, ys, color = (color, 0.5); args...)
end


function manual_stephist!(ax::Axis, edges, heights; args...)
    xs, ys = hist_pts(edges, heights)
    lines!(ax, xs, ys; args...)
end

function manual_density!(ax::Axis, edges, heights; fillcolor = COLORS[1], args...)
    x = edges[1:end-1]
    return fill_between!(ax, x, 0, heights, color = (fillcolor, 0.5); args...)
end

function stacked_hist!(ax::Axis, edges, heights; colors = COLORS, args...)
    last = zeros(2 * length(edges))
    for (c, h) in zip(colors, heights)
        xs, ys = hist_pts(edges, h)
        stacked_h = last .+ ys
        fill_between!(ax, xs, last, stacked_h, color = (c, 0.5); args...)
        last .= stacked_h
    end
end



"""
    stacked!(ax::Axis, xs::Vector{TimeSol}, t; colors = COLORS)

Plot stacked timeseries for a vector of `TimeSol` in order.

# Arguments:
- `ax`: axis for plotting
- `xs`: vector of values over time
- `t`: time points
- `colors`: colors used for plotting the fill and lines
"""
function stacked!(ax::Axis, xs::Vector{TimeSol}, t; colors = COLORS, linecolor = nothing)
    L = length(xs)
    last = zeros(length(t))
    i = 1
    while i ≤ L
        x, i = iterate(xs, i)
        y = last .+ x.vals
        fill_between!(ax, t, last, y, color = (colors[i-1], 0.5), overdraw = false)
        last .= y
    end
    i = 1
    last .= 0.0
    while i ≤ L
        x, i = iterate(xs, i)
        y = last .+ x.vals
        lc = isnothing(linecolor) ? colors[i-1] : linecolor
        lines!(ax, t, y, color = lc, linewidth = 0.3)
        last .= y
    end
    xlims!(ax, 0.0, maximum(t))
    ylims!(ax, low = 0.0)
end
