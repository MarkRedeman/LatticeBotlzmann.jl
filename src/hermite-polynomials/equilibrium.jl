function equilibrium(
    q::Quadrature,
    ρ,#::Array{Float64, 2},
    u,#::Array{Float64, 3},
    T,#::Array{Float64, 2}
)
    fs = zeros(size(ρ, 1), size(ρ, 2), length(q.weights))
    f = zeros(length(q.weights))
    for x_idx = 1:size(ρ, 1), y_idx = 1:size(ρ, 2)
        equilibrium!(q, ρ[x_idx, y_idx], u[x_idx, y_idx, :], T[x_idx, y_idx], f)
        fs[x_idx, y_idx, :] .= f
    end
    return fs
end

function equilibrium(q::Quadrature, ρ::Float64, u::Array{Float64,1}, T::Float64)
    u_squared = sum(u .^ 2)
    f = zeros(length(q.weights))
    hermite_based_equilibrium!(q, ρ, u, T, f)
    return f

    for idx = 1:length(q.weights)
        cs = q.abscissae[1, idx] .* u[1] .+ q.abscissae[2, idx] .* u[2]

        f[idx] = _equilibrium(
            q,
            ρ,
            q.weights[idx],
            cs,
            u_squared,
            T,
            q.abscissae[1, idx]^2 + q.abscissae[2, idx]^2,
        )
    end
    return f
end

function equilibrium!(q::Quadrature, ρ::Float64, u::Array{Float64,1}, T::Float64, f)
    u_squared = 0.0
    for d = 1:dimension(q)
        u_squared += u[d] .^ 2
    end

    for idx = 1:length(q.weights)
        u_dot_xi = q.abscissae[1, idx] .* u[1] .+ q.abscissae[2, idx] .* u[2]

        f[idx] = _equilibrium(
            q,
            ρ,
            q.weights[idx],
            u_dot_xi,
            u_squared,
            T,
            q.abscissae[1, idx]^2 + q.abscissae[2, idx]^2,
        )
    end

    return
end

# Truncated upto order 2
function _equilibrium(q::Quadrature, ρ, weight, u_dot_xi, u_squared, T, xi_squared)
    cs = q.speed_of_sound_squared
    # a = 0.0
    a_H_0 = 1.0
    a_H_1 = cs * u_dot_xi

    # H_2_temperature = cs * ( T .- 1) .* (xi_squared * cs - dimension(q))
    H_2_temperature = 0.0

    a_H_2 = cs^2 * (u_dot_xi .* u_dot_xi) .+ H_2_temperature .+ -cs * u_squared
    return ρ .* weight .* (a_H_0 .+ a_H_1 .+ (1 / 2) * a_H_2)
end
