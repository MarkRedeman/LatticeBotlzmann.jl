
density(q::Quadrature, f::Array{Float64,N}) where {N} = sum(f, dims = N)
density(q::Quadrature, f::Array{Float64,1}) = sum(f)

function momentum(q::Quadrature, f::Array{Float64,3})
    j = Array{Float64}(undef, size(f, 1), size(f, 2), dimension(q))

    @inbounds for x = 1:size(f, 1), y = 1:size(f, 2)
        j[x, y, :] = momentum(q, f[x, y, :])
    end

    return j
end
function momentum(q::Quadrature, f::Array{Float64,1})::Array{Float64,1}
    j = zeros(dimension(q))
    momentum!(q, f, j)
    return j
end

function momentum!(q::Quadrature, f::Array{Float64,1}, j::Array{Float64,1})
    @inbounds for d = 1:dimension(q)
        j[d] = 0.0
        for idx = 1:length(f)
            j[d] += f[idx] * q.abscissae[d, idx]
        end
    end
    return
end

function velocity!(q::Quadrature, f::Array{Float64,1}, ρ::Float64, u::Array{Float64,1})
    @inbounds for d = 1:dimension(q)
        u[d] = 0.0
        for idx = 1:length(f)
            u[d] += f[idx] * q.abscissae[d, idx]
        end
        u[d] / ρ
    end
    return
end

function pressure(
    q::Quadrature,
    f::Array{Float64},
    ρ::Float64,
    u::Array{Float64,1},
)::Float64
    D = dimension(q)
    E = 0.0
    @inbounds for idx = 1:length(f)
        E += f[idx] * (q.abscissae[1, idx]^2 + q.abscissae[2, idx]^2)
    end

    p = q.speed_of_sound_squared * (E - ρ * (u[1]^2 + u[2]^2)) / D
    return p
end
function momentum_flux(q::Quadrature, f::Array{Float64}, ρ::Float64, u::Array{Float64,1})
    D = dimension(q)
    P = zeros(D, D)
    @inbounds for x_idx = 1:D, y_idx = 1:D
        P[x_idx, y_idx] = sum(
            # Pressure tensor
            f[f_idx] *
                (q.abscissae[x_idx, f_idx] - u[x_idx]) *
                (q.abscissae[y_idx, f_idx] - u[y_idx])

            # Stress tensor
            # f[f_idx] * (q.abscissae[x_idx, f_idx]) * (q.abscissae[y_idx, f_idx])
            for f_idx = 1:length(f)
        ) #- ρ * u[x_idx] * u[y_idx]
    end

    return q.speed_of_sound_squared * P #- I(2) * pressure(q, f, ρ, u)

    E = 0.0
    @inbounds for idx = 1:length(f)
        E += f[idx] * (q.abscissae[1, idx]^2 + q.abscissae[2, idx]^2)
    end

    p = q.speed_of_sound_squared * (E - ρ * (u[1]^2 + u[2]^2)) / D
    return p
end

# total energy density
function total_energy(q::Quadrature, f)
    return sum(
        [
            f[x, y, f_idx] * (q.abscissae[:, f_idx]' * q.abscissae[:, f_idx])
            for x = 1:size(f, 1), y = 1:size(f, 2), f_idx = 1:size(f, 3)
        ],
        dims = 3,
    )
end

# internal energy density
function internal_energy(q::Quadrature, f, ρ, u)
    return total_energy(q, f) - kinetic_energy(q, f, ρ, u)
end

# kinetic energy density
function kinetic_energy(q::Quadrature, f, ρ, u)
    return [
        ρ[x, y] * (u[x, y, 1] .^ 2 + u[x, y, 2] .^ 2)
        for x = 1:size(f, 1), y = 1:size(f, 2)
    ]
end

function dimension(q::Quadrature)::Int
    return size(q.abscissae, 1)
end

function temperature(q::Quadrature, f, ρ, u)
    return pressure(q, f, ρ[:, :, 1], u) ./ ρ
    # return internal_energy(q, f, ρ, u) * (2 / dimension(q))
end
function pressure(
    q::Quadrature,
    f::Array{Float64,3},
    ρ::Array{Float64,2},
    u::Array{Float64,3},
)
    return [
        pressure(
            q,
            f[x_idx, y_idx, :],
            ρ[x_idx, y_idx, 1],
            u[x_idx, y_idx, :],
        ) ./ ρ[x_idx, y_idx]
        for x_idx = 1:size(f, 1), y_idx = 1:size(f, 2)
    ]
end
# function temperature(q::Quadrature, f::Array{Float64, 3}, ρ, u)

#     return [
#         pressure(q, f[x_idx, y_idx, :], ρ[x_idx, y_idx, 1], u[x_idx, y_idx, :]) ./ ρ for x_idx = 1:size(f, 1), y_idx = 1:size(f, 2)
#     ]
#     # return internal_energy(q, f, ρ, u) * (2 / dimension(q))
# end

function temperature(q::Quadrature, f::Vector{Float64}, ρ::Float64, u::Vector{Float64})
    return pressure(q, f, ρ, u) ./ ρ
end
