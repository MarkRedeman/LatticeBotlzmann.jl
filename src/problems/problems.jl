abstract type FluidFlowProblem end
abstract type SteadyStateProblem <: FluidFlowProblem end
abstract type TimeDependantProblem <: FluidFlowProblem end

abstract type InitialValueProblem <: FluidFlowProblem end
abstract type DoubleDistributionProblem <: FluidFlowProblem end

has_external_force(::FluidFlowProblem) = false
function acceleration(::FluidFlowProblem, x::Float64, y::Float64, timestep::Float64 = 0.0)
    return zeros(2, 2)
end

function initial_equilibrium(
    quadrature::Quadrature,
    problem::FluidFlowProblem,
    x::Float64,
    y::Float64,
)
    return equilibrium(
        quadrature,
        lattice_density(quadrature, problem, x, y),
        lattice_velocity(quadrature, problem, x, y),
        lattice_temperature(quadrature, problem, x, y),
    )
end

function initial_condition(q::Quadrature, problem::FluidFlowProblem, x::Float64, y::Float64)
    initial_equilibrium(q, problem, x, y)
end

function range(problem::FluidFlowProblem)
    Δx = problem.domain_size[1] / problem.NX
    Δy = problem.domain_size[2] / problem.NY

    x_range = range(Δx / 2, problem.domain_size[1] - Δx / 2, length = problem.NX)
    y_range = range(Δy / 2, problem.domain_size[1] - Δy / 2, length = problem.NY)

    return x_range, y_range
end

function initialize(quadrature::Quadrature, problem::FluidFlowProblem, cm = SRT)
    f = Array{Float64}(undef, problem.NX, problem.NY, length(quadrature.weights))

    x_range, y_range = range(problem)
    for x_idx = 1:problem.NX, y_idx = 1:problem.NY
        f[x_idx, y_idx, :] =
            initial_equilibrium(quadrature, problem, x_range[x_idx], y_range[y_idx])
    end

    return f
end

boundary_conditions(problem::FluidFlowProblem) = BoundaryCondition[]

function force(problem::FluidFlowProblem, x_idx::Int64, y_idx::Int64, time::Float64 = 0.0)
    x_range, y_range = range(problem)

    x = x_range[x_idx]
    y = y_range[y_idx]

    return force(problem, x, y, time)
end

function is_steady_state(problem::FluidFlowProblem)
    return problem.static
end
function is_time_dependant(problem::FluidFlowProblem)
    return !problem.static
end

# Dimensionless
function viscosity(problem) #::FluidFlowProblem)
    return problem.ν * delta_x(problem)^2 / delta_t(problem)
end

function heat_diffusion(problem) #::FluidFlowProblem)
    return problem.κ * delta_x(problem)^2 / delta_t(problem)
end

function delta_t(problem::FluidFlowProblem)
    return delta_x(problem) * problem.u_max
end

function delta_x(problem::FluidFlowProblem)
    return problem.domain_size[1] * (1 / problem.NX)
end
function reynolds(problem::FluidFlowProblem)
    return problem.NY * problem.u_max / problem.ν
end

lattice_viscosity(problem) = problem.ν #::FluidFlowProblem)
lattice_density(q, problem::FluidFlowProblem, x, y, t = 0.0) = density(q, problem, x, y, t)
lattice_velocity(q, problem::FluidFlowProblem, x, y, t = 0.0) =
    problem.u_max * velocity(problem, x, y, t)
lattice_pressure(q, problem::FluidFlowProblem, x, y, t = 0.0) =
    problem.u_max^2 * pressure(q, problem, x, y, t)
lattice_force(problem::FluidFlowProblem, x, y, t = 0.0) =
    problem.u_max * delta_t(problem) * force(problem, x, y, t)
lattice_temperature(q, problem::FluidFlowProblem, x, y, t = 0.0) =
    pressure(q, problem, x, y) / density(q, problem, x, y)

dimensionless_viscosity(problem) = problem.ν * delta_x(problem)^2 / delta_t(problem)
dimensionless_density(problem::FluidFlowProblem, ρ) = ρ
dimensionless_velocity(problem::FluidFlowProblem, u) = u / problem.u_max
dimensionless_pressure(q, problem::FluidFlowProblem, p) = p# (p - 1.0) / (q.speed_of_sound_squared * problem.u_max^2 )
dimensionless_temperature(q, problem::FluidFlowProblem, T) = T #* q.speed_of_sound_squared
dimensionless_force(problem::FluidFlowProblem, F) = F / (problem.u_max * delta_t(problem))

include("taylor-green-vortex-decay.jl")
include("decaying-shear-flow.jl")
include("poiseuille.jl")
include("couette-flow.jl")
include("lid-driven-cavity.jl")
include("linear-hydrodynamics-modes.jl")

# error(::Val{:density}, node, solution) = density(node) - density(solution)
