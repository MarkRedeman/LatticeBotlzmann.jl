density(q::Quadrature, fs::P) where {N <: Int, P <: AbstractArray{<:Real, N}} =
    sum(fs, dims = N)
density(q::Quadrature, f::P) where {P <: AbstractVector{<:Real}} = sum(f)
@inbounds density(q::D2Q9, f::P) where {P <: AbstractVector{<:Real}} = f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9]


function velocity!(
    q::Quadrature,
    f::P,
    ρ::T,
    u::VT,
) where {T <: Real, VT <: AbstractVector{T}, P <: AbstractVector{T}}
    @inbounds for d in 1:dimension(q)
        u[d] = zero(T)
        for idx in 1:length(q.weights)
            u[d] += f[idx] * q.abscissae[d, idx]
        end
        u[d] /= ρ
    end
    return
end
@inbounds function velocity!(
    q::D2Q9,
    f::P,
    ρ::T,
    u::VT,
) where {T <: Real, VT <: AbstractVector{T}, P <: AbstractVector{T}}
    u[1] = (f[5] + f[6] + f[7] - (f[2] + f[3] + f[4])) / ρ
    u[2] = (f[2]  + f[8] + f[9] - (f[4] + f[5] + f[6])) / ρ
    return
end

function pressure(
    q::Quadrature,
    f::P,
    ρ::T,
    u::VT,
) where {T <: Real, VT <: AbstractVector{T}, P <: AbstractVector{T}}
    a_2 =
        sum(f[idx] * hermite(Val{2}, q.abscissae[:, idx], q) for idx in 1:length(q.weights))
    D = dimension(q)

    p = (tr(a_2) - ρ * (u[1]^2 + u[2]^2 - D)) / D
    return p
end

function momentum_flux(
    q::Quadrature,
    f::Population,
    ρ::T,
    u::VT,
) where {T <: Real, VT <: AbstractVector{T}, Population <: AbstractVector{T}}
    D = dimension(q)
    P = zeros(D, D)
    @inbounds for x_idx in 1:D, y_idx in 1:D
        P[x_idx, y_idx] = sum(
            # Pressure tensor
            f[f_idx] *
            (q.abscissae[x_idx, f_idx] - u[x_idx]) *
            (q.abscissae[y_idx, f_idx] - u[y_idx])

            # Stress tensor
            # f[f_idx] * (q.abscissae[x_idx, f_idx]) * (q.abscissae[y_idx, f_idx])
            for f_idx in 1:length(f)
        ) #- ρ * u[x_idx] * u[y_idx]
    end

    return q.speed_of_sound_squared * P #- I * pressure(q, f, ρ, u)

    E = 0.0
    @inbounds for idx in 1:length(f)
        E += f[idx] * (q.abscissae[1, idx]^2 + q.abscissae[2, idx]^2)
    end

    p = q.speed_of_sound_squared * (E - ρ * (u[1]^2 + u[2]^2)) / D
    return p
end

function temperature(
    q::Quadrature,
    f::P,
    ρ::T,
    u::VT,
) where {T <: Real, VT <: AbstractVector{T}, P <: AbstractVector{T}}
    return pressure(q, f, ρ, u) ./ ρ
end

"""
Computes the deviatoric tensor σ

τ is the relaxation time such that ν = cs^2 τ
"""
function deviatoric_tensor(
    q::Quadrature,
    τ::T,
    f::P,
    ρ::T,
    u::VT,
) where {T <: Real, VT <: AbstractVector{T}, P <: AbstractVector{T}}
    D = dimension(q)

    a_bar_2 = sum([
        f[idx] * hermite(Val{2}, q.abscissae[:, idx], q) for idx in 1:length(q.weights)
    ])
    a_eq_2 = equilibrium_coefficient(Val{2}, q, ρ, u, 1.0)
    σ = (a_bar_2 - a_eq_2) / (1 + 1 / (2 * τ))
    return σ - I * tr(σ) / D
end

"""
Compute the temperature from hermite coefficients
"""
function temperature(
    q::Q,
    f::VT,
    a_0::VT,
    a_1::VT,
    a_2::MT,
) where {Q <: Quadrature, T <: Real, VT <: AbstractVector{T}, MT <: AbstractMatrix{T}}
    D = dimension(q)
    P = Array{T}(undef, D, D)
    ρ = a_0
    u = a_1 / ρ
    P = q.speed_of_sound_squared * a_2 - ρ * (u * u' - I(2))

    return tr(P) / (D * ρ)
end
