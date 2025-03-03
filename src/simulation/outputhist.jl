
mutable struct OutputHists
    x::Matrix{Float64}
    x_prime::Matrix{Float64}
    Nmax::Int
    Ny::Int
end

function OutputHists(Nmax::Int, Nx::Int)
    Ny = Nx + 2
    return OutputHists(
        zeros(Float64, (Nmax + 1, Nmax + 1)),
        zeros(Float64, (Ny + 1, Nmax + 1)),
        Nmax,
        Ny,
    )
end

function update_output!(output, X, dt)
    @unpack x, x_prime, Nmax, Ny = output

    x[min(X[1], Nmax)+1, min(X[2], Nmax)+1] += dt
    Y = (X[1] - X[2]) / (X[1] + X[2])
    idx = isnan(Y) ? Ny + 2 : Int(floor(Ny * (Y + 1) / 2)) + 1
    x_prime[idx, min(X[1] + X[2], Nmax)+1] += dt
end

function hist_mean(output::Vector{OutputHists}, N)
    x = reshape(mean(stack([o.x for o in output]), dims = 3), size(output[1].x)) .* N
    x_prime =
        reshape(
            mean(stack([o.x_prime for o in output]), dims = 3),
            size(output[1].x_prime),
        ) .* N
    return OutputHists(x, x_prime, output[1].Nmax, output[1].Ny)
end

function results_dict(h::OutputHists, edges)
    data = Dict(
        "edges" => vec(edges),
        "x" => h.x[1:end-1, 1:end-1],
        "x_prime" => h.x_prime[1:end, 1:end-1],
        "x1" => vec(sum(h.x[1:end-1, :], dims = 2)),
        "x2" => vec(sum(h.x[:, 1:end-1], dims = 1)),
        "s" => vec(sum(h.x_prime[:, 1:end-1], dims = 1)),
        "y" => vec(sum(h.x_prime[1:end, :], dims = 2)),
    )
    return data
end

function normalize!(o::OutputHists, t)
    @. o.x = o.x / t
    @. o.x_prime = o.x_prime / t
end
