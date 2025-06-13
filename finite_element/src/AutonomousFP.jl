
struct AutonomousFP{M<:Function,S<:Function}
    μ::M
    Σ::S
end

@inline J(fp::AutonomousFP, u) = fp.μ * u - fp.Σ ⋅ ∇(u)

function residuals(fp::AutonomousFP, dΩ::Measure)
    f(t) = x -> 0  # probability mass is conserved to no forcing term
    stiffness(t, u, v) = ∫(∇(v) ⋅ (-J(fp, u)))dΩ
    mass(t, dtu, v) = ∫(v ⋅ dtu)dΩ
    res(t, v) = ∫(v ⋅ f(t))dΩ
    return stiffness, mass, res
end


function transient_operator(fp::AutonomousFP, dΩ::Measure, U::FESpace, V::FESpace)
    s, m, r = residuals(fp, dΩ)
    return Gridap.TransientLinearFEOperator((s, m), r, U, V, constant_forms=(true, true))
end

function ss_bilinear(fp::AutonomousFP, dΩ::Measure, mean_val)
    a((p, λ), (q, vλ)) = ∫(J(fp, p) ⋅ ∇(q) + (λ * q + vλ * p))dΩ
    l((q, vλ)) = ∫(mean_val * vλ)dΩ
    return a, l
end

function ss_operator(fp::AutonomousFP, dΩ::Measure, U::FESpace, V::FESpace)::AffineFEOperator
    mean_val = 1 / sum(∫(1)dΩ)
    a, l = ss_bilinear(fp, dΩ, mean_val)
    return AffineFEOperator(a, l, U, V)
end


function save_transient_solution(fname::String, u, u0, fp::AutonomousFP, Ω::Triangulation; dir="examples/results")
    directory = isnothing(dir) || dir == "" ? fname : "$dir/$fname"
    if !ispath(directory)
        @warn "Path to directory $directory not found, creating new directory."
        mkpath("$directory")
    else
        rm.(glob("$directory/*"), force=true)
        @info "Path to directory $directory found, removing existing files."
    end
    pth = "$directory/$fname"
    createpvd(pth) do pvd
        pvd[0.0] = createvtk(
            Ω,
            "$(pth)_0.0.vtu",
            cellfields=["u" => u0, "j" => J(fp, u0)],
        )
        for (t, uᵢ) in u
            t_approx = round(t, sigdigits=3)
            pvd[t] = createvtk(
                Ω,
                "$(pth)_$(t_approx).vtu",
                cellfields=["u" => uᵢ, "j" => J(fp, uᵢ)],
            )
        end
    end
end


function save_stationary_solution(fname::String, p, fp::AutonomousFP, Ω::Triangulation; dir="examples/results")
    directory = isnothing(dir) || dir == "" ? "." : dir
    if !ispath(directory)
        @warn "Path to directory $directory not found, creating new directory."
        mkpath("$directory")
    end
    pth = "$directory/$fname"
    writevtk(Ω,
        "$(pth).vtu",
        cellfields=["u" => p, "j" => J(fp, p)]
    )
end
