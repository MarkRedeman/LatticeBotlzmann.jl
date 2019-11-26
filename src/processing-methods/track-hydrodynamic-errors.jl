struct TrackHydrodynamicErrors{T} <: ProcessingMethod
    problem::FluidFlowProblem
    should_process::Bool
    n_steps::Int64
    stop_criteria::StopCriteria
    df::T
end
TrackHydrodynamicErrors(problem, should_process, n_steps) = TrackHydrodynamicErrors(
    problem,
    should_process,
    n_steps,
    StopCriteria(problem),
    Vector{
        NamedTuple{
            (
                :timestep,
                :error_ρ,
                :error_u,
                :error_p,
                :error_σ_xx,
                :error_σ_xy,
                :error_σ_yy,
                :error_σ_yx,
            ),
            Tuple{Int64,Float64,Float64,Float64,
                Float64,Float64,Float64,Float64},
        },
    }(),
)


function next!(process_method::TrackHydrodynamicErrors, q, f_in, t::Int64)
    if (!process_method.should_process)
        if (t != process_method.n_steps)
            return false
        end
    end

    problem = process_method.problem
    nx, ny, nf = size(f_in)
    x_range, y_range = range(problem)

    f = Array{Float64}(undef, size(f_in, 3))
    u = zeros(dimension(q))
    expected_u = zeros(dimension(q))

    error_ρ = 0.0
    error_u = 0.0
    error_p = 0.0
    error_σ_xx = 0.0
    error_σ_xy = 0.0
    error_σ_yy = 0.0
    error_σ_yx = 0.0

    time = t * delta_t(problem)
    Δ = Float64(y_range.step) * Float64(x_range.step)
    @inbounds for x_idx = 1:nx, y_idx = 1:ny
        @inbounds for f_idx = 1:size(f_in, 3)
            f[f_idx] = f_in[x_idx, y_idx, f_idx]
        end

        ρ = density(q, f)
        velocity!(q, f, ρ, u)
        p = pressure(q, f, ρ, u)

        # Adding the forcing term moves the optimal tau for poiseuille flows
        # F = cm.force(x_idx, y_idx, 0.0)
        # u += cm.τ * F

        ρ = dimensionless_density(problem, ρ)
        u = dimensionless_velocity(problem, u)
        p = dimensionless_pressure(q, problem, p)

        # Analytical
        x = x_range[x_idx]
        y = y_range[y_idx]

        # Compute statistics
        expected_ρ = lbm.density(q, problem, x, y, time)
        expected_u = lbm.velocity(problem, x, y, time)

        expected_p = lbm.pressure(q, problem, x, y, time)
        expected_ϵ = (dimension(q) / 2) * expected_p / expected_ρ
        expected_T = expected_p / expected_ρ

        error_ρ += Δ * (ρ - expected_ρ)^2
        error_p += Δ * (p - expected_p)^2
        error_u += Δ * ((u[1] - expected_u[1])^2 + (u[2] - expected_u[2])^2)
    end

    push!(
        process_method.df,
        (
            timestep = t,
            error_ρ = sqrt(error_ρ),
            error_u = sqrt(error_u),
            error_p = sqrt(error_p),
            error_σ_xx = error_σ_xx,
            error_σ_xy = error_σ_xy,
            error_σ_yy = error_σ_yy,
            error_σ_yx = error_σ_yx,
        ),
    )

    if mod(t, 100) == 0
        if (should_stop!(process_method.stop_criteria, q, f_in))
            @info "Stopping after $t steps out of $process_method.n_steps"
            return true
        end
    end

    if mod(t, 1) == 0
        should_visualize = false
        if (process_method.should_process)
            if t == process_method.n_steps
                should_visualize = true
            end

            if mod(t, max(10, round(Int, process_method.n_steps / 5))) == 0
                should_visualize = true
            end
        end


        if (should_visualize)
            Δt = delta_t(process_method.problem)
            visualize(process_method.problem, q, f_in, time, process_method.df)
        end
    end

    return false
end
