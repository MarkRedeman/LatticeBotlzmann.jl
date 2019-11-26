struct TRT{Force} <: CollisionModel
    τ_symmetric
    τ_asymmetric
    force::Force

    feq::Array{Float64,1}
    f::Array{Float64,1}
    feq_symmetric::Array{Float64,1}
    feq_asymmetric::Array{Float64,1}
    f_symmetric::Array{Float64,1}
    f_asymmetric::Array{Float64,1}
end

function CollisionModel(cm::Type{<:TRT}, q::Quadrature, problem::FluidFlowProblem)
    τ = q.speed_of_sound_squared * lattice_viscosity(problem) + 0.5
    Λ = 1 / 6
    τ_ = 0.5 + Λ / (τ - 0.5)
    @show τ_

    if has_external_force(problem)
        force_field = (x_idx, y_idx, t) -> lattice_force(problem, x_idx, y_idx, t)
        return TRT(
            τ,
            τ_,
            force_field,
            similar(q.weights),
            similar(q.weights),
            similar(q.weights),
            similar(q.weights),
            similar(q.weights),
            similar(q.weights),
        )
        return TRT(τ, τ_, force_field)
    else
        return TRT(
            τ,
            τ_,
            nothing,
            similar(q.weights),
            similar(q.weights),
            similar(q.weights),
            similar(q.weights),
            similar(q.weights),
            similar(q.weights),
        )
    end
end

function collide!(
    collision_model::TRT{Force},
    q::Quadrature,
    f_in,
    f_out;
    time = 0.0,
) where {Force}
    τ_s = collision_model.τ_symmetric
    τ_a = collision_model.τ_asymmetric

    cm = collision_model
    feq = cm.feq
    f = cm.f
    feq_symmetric = cm.feq_symmetric
    feq_asymmetric = cm.feq_asymmetric
    f_symmetric = cm.f_symmetric
    f_asymmetric = cm.f_asymmetric

    u = zeros(dimension(q))
    if !(Force <: Nothing)
        F = zeros(dimension(q))
    end

    nx, ny, nf = size(f_in)
    @inbounds for x = 1:nx, y = 1:ny
        @inbounds for f_idx = 1:nf
            f[f_idx] = f_in[x, y, f_idx]
        end

        ρ = density(q, f)

        # Momentum
        velocity!(q, f, ρ, u)

        # Temperature
        # T = temperature(q, f, ρ, u)
        T = 1.0

        if Force <: Nothing
            equilibrium!(q, ρ, u, T, feq)
        else
            F .= collision_model.force(x, y, time)

            equilibrium!(q, ρ, u + τ_s * F, T, feq)
        end

        @inbounds for f_idx = 1:nf
            opposite_idx = opposite(q, f_idx)
            feq_symmetric[f_idx] = 0.5 * (feq[f_idx] + feq[opposite_idx])
            feq_asymmetric[f_idx] = 0.5 * (feq[f_idx] - feq[opposite_idx])
            f_symmetric[f_idx] = 0.5 * (f[f_idx] + f[opposite_idx])
            f_asymmetric[f_idx] = 0.5 * (f[f_idx] - f[opposite_idx])
        end

        @inbounds for f_idx = 1:nf
            f_out[x, y, f_idx] =
                f[f_idx] + (
                    -(1 / τ_s) * (f_symmetric[f_idx] - feq_symmetric[f_idx]) -
                    (1 / τ_a) * (f_asymmetric[f_idx] - feq_asymmetric[f_idx])
                )
        end
    end
    return
end
