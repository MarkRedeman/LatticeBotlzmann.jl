struct LidDrivenCavityFlow <: FluidFlowProblem
    rho_0::Float64
    u_max::Float64
    ν::Float64
    NX::Int64
    NY::Int64
    domain_size::Tuple{Float64,Float64}
end

function LidDrivenCavityFlow(
    ν = 1.0 / 6.0,
    scale = 2,
    NX = 16 * scale,
    NY = NX,
    domain_size = (1.0, 1.0),
)
    # # Since our velocity profile is constant over x it is sufficient to only take
    # # 3 nodes in the x direction
    # NX = 5
    u_max = 0.01 / scale

    return LidDrivenCavityFlow(1.0, u_max, ν, NX, NY, domain_size)
end

function density(
    q::Quadrature,
    problem::LidDrivenCavityFlow,
    x::Float64,
    y::Float64,
    timestep::Float64 = 0.0,
)
    return 1.0
end

function pressure(
    q::Quadrature,
    problem::LidDrivenCavityFlow,
    x::Float64,
    y::Float64,
    timestep::Float64 = 0.0,
)
    return 1.0
end

function velocity(
    problem::LidDrivenCavityFlow,
    x::Float64,
    y::Float64,
    timestep::Float64 = 0.0,
)
    return [
        0.0
        0.0
    ]
end

"""
Apply a moving wall to the top of the domain and bounce back to all other
sides of the domain
"""
boundary_conditions(problem::LidDrivenCavityFlow) = [
    BounceBack(East(), 1:problem.NX, 1:problem.NY),
    BounceBack(South(), 1:problem.NX, 1:problem.NY),
    BounceBack(West(), 1:problem.NX, 1:problem.NY),
    MovingWall(North(), 1:problem.NX, 1:problem.NY, [problem.u_max, 0]),
]
