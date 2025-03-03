"""
    PDMPStationarySol

Second order PDMP approximation parameters

# Fields:
- `κ::Float64`: bias in the degradation rates
- `β::Float64`: boas in the influx rates
- `k1::Float64`: decay rate in state 1
- `k2::Float64`: decay rate in state 2
"""
struct PDMPStationarySol
    κ::Float64
    β::Float64
    k1::Float64
    k2::Float64
end

bounds(p::PDMPStationarySol) = (1 / p.k1, 1 / p.k2)

const_norm(p::PDMPStationarySol) =
    cos(π * p.β / 2) * (1 + p.κ) * ((1 - p.κ) / (1 + p.κ))^((1 - p.β) / 2) / (2π * p.κ)

probability_fraction(p::PDMPStationarySol, x) = @. (1 + x * (p.κ - 1)) / (2 * p.κ * x)

function pdf_p(p::PDMPStationarySol, x)
    @unpack β, k1, k2 = p
    expo = -(1 + β) / 2
    return (k1 * x - 1)^expo * (1 - k2 * x)^(-expo) * const_norm(p) / x
end

function pdf_m(p::PDMPStationarySol, x)
    @unpack β, k1, k2 = p
    expo = (1 - β) / 2
    return (k1 * x - 1)^expo * (1 - k2 * x)^(-expo) * const_norm(p) / x
end

"""Create stationary distribution params. Beta defaults to zero (symmetric switching)."""
from_params(t::Type{PDMPStationarySol}, κ, β) = PDMPStationarySol(κ, β, 1 + κ, 1 - κ)

from_params(t::Type{PDMPStationarySol}, κ) = PDMPStationarySol(κ, 0, 1 + κ, 1 - κ)


"""
    theoretical_pdfs(approx::PDMPStationarySol; npts::Int = 2000, ϵ = 1.0e-9)

Calculate stationary states for a PDMP.

# Arguments:
- `approx`: PDMP approximation
- `npts`: spatial resolution for evaluation
- `ϵ`: small reduction in bounds to avoid singularities at the boundaries.
"""
function theoretical_pdfs(approx::PDMPStationarySol; npts::Int = 2000, ϵ = 1e-9)
    ϕ1, ϕ2 = bounds(approx)
    ϕ = collect(range(ϕ1 + ϵ, ϕ2 - ϵ, npts))
    Πp = [pdf_p(approx, ϕᵢ) for ϕᵢ in ϕ]
    Πm = [pdf_m(approx, ϕᵢ) for ϕᵢ in ϕ]
    return Πp, Πm, ϕ
end
