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


f = fig_in_cm(2, 5.46)

Colorbar(f[1, 1], limits=(4e-3, 4e1), colormap=Reverse(:grays),
    vertical=false, scale=log10, spinewidth=0.3, size=8)

f
save_fig("colorbar", f, fmt="pdf")