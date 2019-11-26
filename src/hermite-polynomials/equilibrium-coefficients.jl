# Use F.2 to generate them

function hermite_based_equilibrium(q, f)
    ρ = sum(f)
    u = velocity(q, f)
    T = temperature(q, f, ρ, u)

    hermite_based_equilibrium(q, ρ, u, T)
end
function hermite_based_equilibrium(q, ρ, u, T)
    @warn "HUH"
    N = div(lbm.order(q), 2)
    cs = 1 / q.speed_of_sound_squared
    D = dimension(q)

    # Compute (get!?) the hermite polynomials for this quadrature
    Hs = [[hermite(Val{n}, q.abscissae[:, i], q) for i = 1:length(q.weights)] for n = 1:N]

    as = [equilibrium_coefficient(Val{n}, q::Quadrature, ρ, u, T) for n = 1:N]

    return [
        q.weights[f_idx] * (ρ + sum([
            sum(as[n] .* Hs[n][f_idx]) / (factorial(n) * cs^n) for n = 1:N
        ]))
        for f_idx = 1:length(q.weights)
    ]
end


function hermite_based_equilibrium!(q::Q, ρ, u, T, f) where {Q<:Quadrature}
    N = div(lbm.order(q), 2)

    # NOTE: we skip H_0 since its one (also because of Julia's 1 based index....)
    #Hs = [[hermite(n, q.abscissae[:, i], q) for i = 1:length(q.weights)] for n = 1:N]

    cs = 1 / q.speed_of_sound_squared
    D = dimension(q)
    δ(α, β) = α == β ? 1 : 0

    a0 = ρ
    H1 = [hermite(Val{1}, q.abscissae[:, i], q) for i = 1:length(q.weights)]
    a1 = u

    Hs = [[hermite(Val{n}, q.abscissae[:, i], q) for i = 1:length(q.weights)] for n = 1:N]
    a_eq = [equilibrium_coefficient(Val{n}, q, ρ, u, T) for n = 1:N]


    @inbounds for f_idx = 1:length(f)
        f[f_idx] =
            q.weights[f_idx] * (
                ρ + sum([dot(a_eq[n], Hs[n][f_idx]) / (factorial(n) * cs^n) for n = 1:N])
            )
    end
    return
end

# equilibrium(Val{1}, ρ, u, T)

function equilibrium_coefficient(::Type{Val{0}}, q::Quadrature, ρ, u, T)
    return ρ
end

function equilibrium_coefficient(::Type{Val{1}}, q::Quadrature, ρ, u, T)
    return ρ * u
end

function equilibrium_coefficient(::Type{Val{2}}, q::Quadrature, ρ, u, T)
    cs = 1 / q.speed_of_sound_squared
    return ρ * (u * u' + cs * (T - 1) * I(2))
end
function equilibrium_coefficient(::Type{Val{3}}, q::Quadrature, ρ, u, T)
    cs = 1 / q.speed_of_sound_squared
    D = dimension(q)
    return ρ * [
        u[a] *
        u[b] *
        u[c] + cs * (T - 1) * (u[a] * δ(b, c) + u[b] * δ(a, c) + u[c] * δ(a, b))
        for a = 1:D, b = 1:D, c = 1:D
    ]
end
function equilibrium_coefficient(::Type{Val{4}}, q::Quadrature, ρ, u, T)
    cs = 1 / q.speed_of_sound_squared
    D = dimension(q)
    return ρ * [
        u[a] * u[b] * u[c] * u[d] +
        cs *
        (T - 1) *
        (
            u[a] * u[b] * δ(c, d) +
            u[a] * u[c] * δ(b, d) +
            u[a] * u[d] * δ(b, d) +

            u[b] * u[c] * δ(a, d) +
            u[b] * u[d] * δ(a, d) +
            u[c] * u[d] * δ(a, b)
        ) +
        cs^2 * (T - 1)^2 * (δ(a, b) * δ(c, d) + δ(a, c) * δ(b, d) + δ(a, d) * δ(b, c))
        for a = 1:D, b = 1:D, c = 1:D, d = 1:D
    ]
end


# function temperature(q::Quadrature, f::Vector{Float64}, ρ::Float64, u::Vector{Float64})
#     return pressure(q, f, ρ, u) ./ ρ


#     Hs = [[hermite(n, q.abscissae[:, i], q) for i = 1:length(q.weights)] for n = 1:2]
#     ρ_f = sum(f)
#     a_f = [sum([f[idx] * Hs[n][idx] for idx = 1:length(q.weights)]) for n = 1:2]

#     D = dimension(q)
#     P = Array{Float64}(undef, D, D)
#     u_f = a_f[1] / ρ_f
#     P_f = q.speed_of_sound_squared * a_f[2] - ρ_f * (u_f * u_f' - I(2))
#     T_f = tr(P_f) / (D * ρ_f)

#     # No need to allocate P
#     T_2 = (
#         q.speed_of_sound_squared * (a_f[2][1, 1] + a_f[2][2, 2])
#         - ρ_f * (u_f[1]^2 + u_f[2]^2 - D)
#     ) / (D * ρ_f)

#     @show T T_f T_2

#     return T_f
# end

"""
Compute the temperature from hermite coefficients
"""
function temperature(
    q::Quadrature,
    f::Vector{Float64},
    a_0::Float64,
    a_1::Vector{Float64},
    a_2::Matrix{Float64},
)
    D = dimension(q)
    P = Array{Float64}(undef, D, D)
    ρ = a_0
    u = a_1 / ρ
    P = q.speed_of_sound_squared * a_2 - ρ * (u * u' - I(2))

    T = tr(P) / (D * ρ)

    return T
end
