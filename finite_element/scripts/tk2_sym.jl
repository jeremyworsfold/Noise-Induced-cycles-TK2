using GridapFokkerPlanck
using Gridap

FNAME = _get_script_name(@__FILE__)

Vol = 100000
λ = 0.5
κ = 0.0
β = 0.0

drift1(x) = 1 - x[1] * (1 + κ * x[2])
drift2(x) = -x[2] / x[1] - κ * (1 - x[2]^2)

B11(x) = (1 + x[1] * (1 + κ * x[2])) / Vol
B12(x) = (κ * x[1] * (1 - x[2]^2) - x[2] / x[1]) / Vol
B21(x) = B12(x)
B22(x) = (((1 - κ * x[2]) * (1 - x[2]^2) + (1 + x[2]^2)) / x[1]^2 + Vol / λ * (1 - x[2]^2)) / Vol

effective_drift1(x) = (1 - 1 / x[1] + κ) / Vol
effective_drift2(x) = (κ + (x[2] - κ * (1 + (x[1]^2 - 3) * x[2]^2)) / x[1]^2 - 2 * Vol * x[2] / λ) / Vol

struct μ <: Function end
(μf::μ)(x) = VectorValue(
    drift1(x) - effective_drift1(x),
    drift2(x) - effective_drift2(x)
)

struct Σ <: Function end
(Σf::Σ)(x) = TensorValue{2,2,Float64}(
    B11(x),
    B12(x),
    B21(x),
    B22(x)
)


fokker_planck = AutonomousFP{μ,Σ}(μ(), Σ())

#################################################################

config = load_config(FokkerPlanckStationary, FNAME)
sol = solveFP(fokker_planck, config)


struct ρs_ <: Function end
(ρf::ρs_)(x) = x[1]
struct ρy_ <: Function end
(ρf::ρy_)(x) = x[2]
ρfs = ρs_()
ρfy = ρy_()

dΩ = sol[:domega]
u = sol[:u]

fe_s = sum(∫(u * ρfs)dΩ)
fe_y = sum(∫(u * ρfy)dΩ)

av_s = (1 - β * κ) / (1 - κ^2)
av_y = iszero(κ) ? β : (k2^(b2) * k1^(b1) - 1) / κ
@info "Theoretical <s, y>:" av_s av_y
@info "FE solution <s, y>:" fe_s fe_y

save_stationary_solution("$FNAME", sol[:u], fokker_planck, sol[:omega])
