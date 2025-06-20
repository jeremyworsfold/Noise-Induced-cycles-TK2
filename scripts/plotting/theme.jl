Attributes(
    size = (8.6 * 28.3465, 6 * 28.3465),
    fontsize = 9,
    fonts = (; regular = "OPTITimes-Roman.otf",),
    palette = (
        color = colorschemes[:seaborn_colorblind].colors,
        marker = [:rect, :circle, :utriangle, :diamond],
        linestyle = [:dash, :dot, :dashdot, :dashdotdot],
    ),
    linewidth = 1,
    colormap = :PRGn,
    strokewidth=0.1,
    Scatter = (
        color = COLORS[4],
        markersize = 7,
        marker = :circle,
        # strokewidth = 0.7,
        # strokecolor = :black,
    ),
    grid = false,
    Axis = (
        xgridvisible = false,
        ygridvisible = false,
        spinewidth = 0.5,
        xtickwidth = 0.5,
        ytickwidth = 0.5,
    ),
    Legend = (
        framevisible = false,
        tellheight = false,
        tellwidth = false,
        halign = :right,
        valign = :top,
        padding = 5,
        rowgap = -8,
        patchlabelgap = 2,
    ),
    Arrows = (arrowsize = 1, lengthscale = 0.08, arrowcolor = (:black,0.0), normalize = true),
)
