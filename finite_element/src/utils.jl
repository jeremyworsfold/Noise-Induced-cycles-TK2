_get_script_name(file) = splitext(splitpath(file)[end])[1]

function crank_nicholson_integrator(Δt::Float64; theta::Float64=0.5)
    return ThetaMethod(LUSolver(), Δt, theta)
end

function normalized_IC(init_cond::Function, U::FESpace, dΩ::Measure)
    u0_unnormalized = interpolate_everywhere(init_cond, U)
    A = 1 / sum(∫(u0_unnormalized)dΩ)
    return interpolate_everywhere(x -> A * init_cond(x), U)
end

function create_FE_spaces(
    model::DiscreteModel;
    integration_degree::Int=2,
    polynomial_order::Int=1,
    conformity=:H1,
    precision=Float64,
    lagrange_normalization=false
)
    Ω = Triangulation(model)
    dΩ = Measure(Ω, integration_degree)
    refFE = ReferenceFE(lagrangian, precision, polynomial_order)
    V_p = TestFESpace(model, refFE; conformity=conformity)
    U_p = TrialFESpace(V_p)

    if !lagrange_normalization
        return Ω, dΩ, U_p, V_p
    end

    V_λ = ConstantFESpace(model)
    U_λ = TrialFESpace(V_λ)
    U = MultiFieldFESpace([U_p, U_λ])
    V = MultiFieldFESpace([V_p, V_λ])
    return Ω, dΩ, U, V
end
