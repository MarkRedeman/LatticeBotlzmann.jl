abstract type StopCriteria end
struct NoStoppingCriteria <: StopCriteria end
mutable struct MeanVelocityStoppingCriteria <: StopCriteria
    old_mean_velocity::Float64
    tolerance::Float64
    problem::FluidFlowProblem
end
StopCriteria(problem::FluidFlowProblem) = NoStoppingCriteria()
StopCriteria(problem::PoiseuilleFlow) = MeanVelocityStoppingCriteria(0.0, 1e-12, problem)
StopCriteria(problem::CouetteFlow) = MeanVelocityStoppingCriteria(0.0, 1e-12, problem)
StopCriteria(problem::LidDrivenCavityFlow) =
    MeanVelocityStoppingCriteria(0.0, 1e-5, problem)
StopCriteria(problem::DecayingShearFlow) =
    problem.static ? MeanVelocityStoppingCriteria(0.0, 1e-8, problem) : NoStoppingCriteria()

should_stop!(::StopCriteria, q, f_in) = false
function should_stop!(stop_criteria::MeanVelocityStoppingCriteria, q, f_in)
    f = Array{Float64}(undef, size(f_in, 3))
    u = zeros(dimension(q))
    Nx = size(f_in, 1)
    Ny = size(f_in, 2)

    divide_by = 0
    u_mean = 0.0

    @inbounds for x_idx = 1:Nx, y_idx = 1:Ny
        divide_by += 1

        # Calculated
        @inbounds for f_idx = 1:size(f_in, 3)
            f[f_idx] = f_in[x_idx, y_idx, f_idx]
        end
        ρ = density(q, f)
        velocity!(q, f, ρ, u)

        u_mean += u[1]
    end

    u_mean /= divide_by

    converged = abs(u_mean / stop_criteria.old_mean_velocity - 1)

    if (converged < stop_criteria.tolerance)
        return true
    end

    if (isnan(u_mean))
        @warn "nan in velocity profile"
        return true
    end

    stop_criteria.old_mean_velocity = u_mean

    return false
end
