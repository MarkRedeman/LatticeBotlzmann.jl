struct PoiseuilleFlow <: FluidFlowProblem
    rho_0::Float64
    u_max::Float64
    ν::Float64
    NX::Int64
    NY::Int64
    k::Float64
    domain_size::Tuple{Float64,Float64}
    G::Float64
end

function PoiseuilleFlow(
    ν = 1.0 / 6.0,
    scale = 2,
    NX = 5 * scale + 0,
    NY = NX,
    domain_size = (1.0, 1.0);
    static = true,
)
    # Since our velocity profile is constant over x it is sufficient to only take
    # 3 nodes in the x direction
    NX = 3
    u_max = 0.1 / scale

    return PoiseuilleFlow(1.0, u_max, ν, NX, NY, 1.0, domain_size, 1.0)
end

function density(
    q::Quadrature,
    problem::PoiseuilleFlow,
    x::Float64,
    y::Float64,
    timestep::Float64 = 0.0,
)
    return 1.0
end

function pressure(
    q::Quadrature,
    problem::PoiseuilleFlow,
    x::Float64,
    y::Float64,
    timestep::Float64 = 0.0,
)
    return 1.0
end

function velocity(problem::PoiseuilleFlow, x::Float64, y::Float64, timestep::Float64 = 0.0)
    return [
        y * (problem.domain_size[2] - y) * (problem.G / 2)
        0.0
    ]
end
function force(problem::PoiseuilleFlow, x::Float64, y::Float64, time::Float64 = 0.0)
    return [
        viscosity(problem) * problem.G
        0.0
    ]
end

function delta_x(problem::PoiseuilleFlow)
    return problem.domain_size[2] / problem.NY
end

has_external_force(problem::PoiseuilleFlow) = true

"""
Apply a bottom and top wall bounce back boundary condition
"""
boundary_conditions(problem::PoiseuilleFlow) = [
    BounceBack(North(), 1:problem.NX, 1:problem.NY),
    BounceBack(South(), 1:problem.NX, 1:problem.NY),
]
