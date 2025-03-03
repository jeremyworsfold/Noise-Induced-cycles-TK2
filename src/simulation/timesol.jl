
replace_nan(v) = map(x -> isnan(x) ? zero(x) : x, v)

replace_inf(v) = map(x -> isinf(x) ? one(x) : x, v)

replace_out_range(v) = map(x -> x > 2 || x < -2 ? zero(x) : x, v)

"""
    TimeSol

Timeseries solution of a variable at discrete (not necessarily equispaced) time points.

# Fields:
- `vals::Vector{Float64}`: timeseries of the variable
- `label::Union{LaTeXString, String}`: string representation of variable
"""
struct TimeSol
    vals::Vector{Float64}
    label::Union{LaTeXString,String}
end

function get_variables(sol, Omega; vars = (:x1, :x2))
    d = Dict{Symbol,TimeSol}()
    if :x1 ∈ vars
        d[:x1] = TimeSol(sol[:X₁], L"X_1")
    end
    if :x2 ∈ vars
        d[:x2] = TimeSol(sol[:X₂], L"X_2")
    end
    if :s ∈ vars
        d[:s] = TimeSol((sol[:X₁] + sol[:X₂]), L"s")
    end
    if :z ∈ vars
        d[:z] = TimeSol((sol[:X₁] - sol[:X₂]), L"z")
    end
    if :y ∈ vars
        d[:y] = TimeSol(
            map(
                replace_out_range,
                map(replace_nan, (sol[:X₁] - sol[:X₂]) ./ (sol[:X₁] + sol[:X₂])),
            ),
            L"y",
        )
    end

    return d, sol.t
end
