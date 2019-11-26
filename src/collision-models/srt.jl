struct SRT{Force} <: CollisionModel
    τ::Float64
    force::Force
end
SRT(τ) = SRT(τ, nothing)

function CollisionModel(cm::Type{<:SRT}, q::Quadrature, problem::FluidFlowProblem)
    τ = q.speed_of_sound_squared * lattice_viscosity(problem) + 0.5

    if has_external_force(problem)
        force_field = (x_idx, y_idx, t) -> lattice_force(problem, x_idx, y_idx, t)
        return SRT(τ, force_field)
    else
        return SRT(τ, nothing)
    end
end

function collide!(
    collision_model::SRT{Force},
    q::Quadrature,
    f_in,
    f_out;
    time = 0.0,
) where {Force}
    τ = collision_model.τ

    feq = Array{Float64}(undef, size(f_in, 3))
    f = Array{Float64}(undef, size(f_in, 3))
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

            equilibrium!(q, ρ, u + τ * F, T, feq)
        end

        @inbounds for f_idx = 1:nf
            f_out[x, y, f_idx] = (1 - 1 / τ) * f[f_idx] + (1 / τ) * feq[f_idx]
        end
    end
    return
end
