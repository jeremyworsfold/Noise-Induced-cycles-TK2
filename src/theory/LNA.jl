
EPS = 5e-15

LNA_fluctuation(ξ, ϕ) = exp(-ξ^2 / (2ϕ)) / sqrt(2π * ϕ)

function LNA_convolution_p(p::PDMPStationarySol, x::Real, Omega::Real)
    domain = bounds(p) .+ (EPS, -EPS)
    function f(ϕ, p)
        p_, Omega = p
        ξ = sqrt(Omega) * (x - ϕ)
        return LNA_fluctuation(ξ, ϕ) * pdf_p(p_, ϕ)
    end
    return solve(IntegralProblem(f, domain, (p, Omega)), QuadGKJL()).u
end


function LNA_convolution_m(p::PDMPStationarySol, x::Real, Omega::Real)
    domain = bounds(p) .+ (EPS, -EPS)
    function f(ϕ, p)
        p_, Omega = p
        ξ = sqrt(Omega) * (x - ϕ)
        return LNA_fluctuation(ξ, ϕ) * pdf_m(p_, ϕ)
    end
    return solve(IntegralProblem(f, domain, (p, Omega)), QuadGKJL()).u
end


"""
    theoretical_LNA_stationary_profiles(p::PDMPStationarySol, x::AbstractVector, Omega::Real)

Perform linear noise approximation on PDMP stationary state.

# Arguments:
- `p`: PDMP approximation parameters
- `x`: positions to be evaluated
- `Omega`: system size (typical population)
"""
function theoretical_LNA_pdfs(p::PDMPStationarySol, x::AbstractVector, Omega::Real)
    ΠpLNA = [LNA_convolution_p(p, x_i, Omega) for x_i in x]
    ΠmLNA = [LNA_convolution_m(p, x_i, Omega) for x_i in x]

    problem = SampledIntegralProblem(Π1LNA .+ Π2LNA, x)
    method = TrapezoidalRule()
    sLNA = solve(problem, method)
    ΠpLNA ./= sLNA
    ΠmLNA ./= sLNA

    return vec(ΠpLNA), vec(ΠmLNA)
end
