@option "CartesianModel" struct CartesianSpecs
    limits::Vector{Vector{Float64}}
    N_pts::Vector{Int}
end

function generate(p::CartesianSpecs)
    limits = tuple(reduce(vcat, p.limits)...)
    return CartesianDiscreteModel(limits, tuple(p.N_pts...))
end

function sanitize(p::CartesianSpecs)
    @argcheck length(p.limits) == length(p.N_pts)
    @argcheck all([length(a) == 2 for a in p.limits])
    @argcheck all([a[2] - a[1] > 0 for a in p.limits])
end

@option "MeshFromFile" struct MeshFromFile
    directory::String = "examples/models"
    filename::String
end

function generate(p::MeshFromFile)
    dir = p.directory
    name, ext = splitext(p.filename)
    if ext == ".msh"
        return GmshDiscreteModel("$dir/$name.msh")
    elseif ext == ".json" || ext == ".jld"
        return DiscreteModelFromFile("$dir/$name$ext")
    end
end

function sanitize(p::MeshFromFile)
    _, ext = splitext(p.filename)
    @argcheck ext in [".msh", ".json", ".jld"]
    @argcheck ispath("$(p.directory)/$(p.filename)")
end

#################################################################

@option "ThetaMethod" struct ThetaMethodIntegrator
    delta_t::Float64
    theta::Float64 = 0.5
    T_final::Float64 = 1.0
end

function sanitize(p::ThetaMethodIntegrator)
    @argcheck p.delta_t > 0.0
    @argcheck 0.0 < p.theta < 1.0
    @argcheck p.T_final > p.delta_t
end

generate(p::ThetaMethodIntegrator) = crank_nicholson_integrator(p.delta_t; theta=p.theta)

#################################################################
@option "fe_specs" struct FESpacesSpecs
    polynomial_order::Int = 1
    integration_degree::Int = 2
    conformity::String = "H1"
    precision::String = "Float64"
end
_precision_dict = Dict("float64" => Float64, "float32" => Float32, "float16" => Float16)

function generate(model::DiscreteModel, p::FESpacesSpecs; lagrange_normalization=false)
    type = lowercase(p.precision)
    return create_FE_spaces(
        model,
        integration_degree=p.integration_degree,
        polynomial_order=p.polynomial_order,
        conformity=Symbol(p.conformity),
        precision=_precision_dict[type],
        lagrange_normalization=lagrange_normalization
    )
end

function sanitize(p::FESpacesSpecs)
    @argcheck p.polynomial_order > 0
    @argcheck p.integration_degree ≥ p.polynomial_order
    allowed_precisions = keys(_precision_dict)
    @argcheck lowercase(p.precision) in allowed_precisions
end

#################################################################

@option "UniformIC" struct UniformIC end

generate(p::UniformIC) = x -> 1

function sanitize(p::UniformIC) end

@option "GaussianIC" struct GaussianIC
    x0::Vector{Float64}
    sigma::Float64
end

function generate(p::GaussianIC)
    σ2 = p.sigma^2
    a = VectorValue(p.x0...)
    return x -> exp(-inner((x - a), (x - a)) / (2 * σ2))
end

function sanitize(p::GaussianIC)
    @argcheck p.sigma > 0.0
end

#################################################################

@option struct FokkerPlanckStationary
    model::Union{CartesianSpecs,MeshFromFile}
    fe_specs::FESpacesSpecs
end


function sanitize(fp::FokkerPlanckStationary; verbose=true)
    if verbose
        @info "Checking configurations of stationary solver for correct inputs"
    end
    sanitize(fp.model)
    sanitize(fp.fe_specs)
    if verbose
        @info "Checking complete!"
        dump(fp)
    end
end


function solveFP(fokker_planck::AutonomousFP, config::FokkerPlanckStationary)
    model = generate(config.model)
    Ω, dΩ, U, V = generate(model, config.fe_specs, lagrange_normalization=true)

    operator = ss_operator(fokker_planck, dΩ, U, V)
    u, λ = solve(operator)
    return Dict(:u => u, :lambda => λ, :omega => Ω, :domega => dΩ)
end

#################################################################


@option struct FokkerPlanckTransient
    model::Union{CartesianSpecs,MeshFromFile}
    fe_specs::FESpacesSpecs
    scheme::ThetaMethodIntegrator
    initial_cond::Union{UniformIC,GaussianIC}
end

function sanitize(fp::FokkerPlanckTransient; verbose=true)
    if verbose
        @info "Checking configurations of transient solver for correct inputs"
    end
    sanitize(fp.model)
    sanitize(fp.fe_specs)
    sanitize(fp.scheme)
    sanitize(fp.initial_cond)
    if verbose
        @info "Checking complete!"
        dump(fp)
    end
end

function solveFP(fokker_planck::AutonomousFP, config::FokkerPlanckTransient)
    model = generate(config.model)
    Ω, dΩ, U, V = generate(model, config.fe_specs, lagrange_normalization=false)
    integrator = generate(config.scheme)
    init_cond = generate(config.initial_cond)
    operator = transient_operator(fokker_planck, dΩ, U, V)
    # u0 = interpolate_everywhere(init_cond, U)
    u0 = normalized_IC(init_cond, U, dΩ)

    sol = solve(integrator, operator, 0.0, config.scheme.T_final, u0)

    return Dict(:u => sol, :u0 => u0, :omega => Ω, :domega => dΩ)
end

#################################################################


function load_config(
    p::Union{Type{FokkerPlanckStationary},Type{FokkerPlanckTransient}},
    fname::AbstractString;
    dir="examples/parameters",
    check_pars::Bool=true,
    verbose=true,
)
    name, _ = splitext(fname)
    pth = isnothing(dir) || dir == "" ? pth = name : "$dir/$name"

    conf = from_toml(p, "$(pth).toml")
    if check_pars
        sanitize(conf, verbose=verbose)
    end
    return conf
end

to_stationary_config(p::FokkerPlanckTransient) = FokkerPlanckStationary(p.model, p.fe_specs)