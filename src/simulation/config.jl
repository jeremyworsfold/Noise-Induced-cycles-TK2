"""
    SimParams

Collection of parameters to generate and solve a particular reaction network

# Fields:
- `T_end::Float64` : Final time of simulation
- `model::StoichiometryModel`: reaction network model
"""
@option struct SimParams
    T_end::Float64 = 1e3
    params::StoichiometryModel
end

function save_res_and_conf(
    name::AbstractString,
    dict::Dict,
    conf::SimParams;
    path = "results",
)
    JLD.save("$path/$name.jld", dict)
    save(conf, name, path = path)
end


function load_res_and_conf(name::AbstractString, path = "results"; show = true)
    data = JLD.load("$path/$name.jld")
    conf = load(SimParams, name, path = path, show = show)
    return data, conf
end


"""
    load(t::Type{SimParams}, name::AbstractString; path = "parameters", show = true)

Load simulation parameters from a `.toml` file
"""
function load(t::Type{SimParams}, name::AbstractString; path = "parameters", show = true)
    conf = from_toml(t, "$path/$name.toml")
    if show == true
        dump(conf)
    end
    return conf
end

"""
    save(conf::SimParams, name; dir = results/raw)

Save simulation parameters to a `.toml` file
"""
function save(conf::SimParams, name::AbstractString; path = "results")
    open("$path/$name.toml", "w") do io
        to_toml(io, conf)
    end
end
