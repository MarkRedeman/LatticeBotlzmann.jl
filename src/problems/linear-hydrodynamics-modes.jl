# This problem attempts to reproduce the results from the article
# A General Multiple-Relaxation-Time Boltzmann Collision Model
# by Xiaowen Shan and Hudong Chen

# struct LinearHydrodynamicsModes <: FluidFlowProblem
#     rho_0::Float64
#     u_max::Float64
#     ν::Float64
#     NX::Int64
#     NY::Int64
#     domain_size::Tuple{Float64, Float64}
# end

"""
Models thermal diffusion given by,

ρ = ρ_0 + ρ̃ sin(2pi y / Ly)
u = 0
θ = ρ_0 θ_0 / ρ

It is expected that the temperature θ decays at the rate of exp(- κ t)
"""
struct LinearizedThermalDiffusion <: FluidFlowProblem
    ρ_0::Float64
    θ_0::Float64

    ρ̃::Float64
    ũ::Float64

    ν::Float64
    κ::Float64

    domain_size::Tuple{Float64,Float64}

    u_max::Float64
    NX::Int64
    NY::Int64
end
function LinearizedThermalDiffusion(ν, κ, scale, NY = 4 * scale)
    u_max = 0.01 / scale

    return LinearizedThermalDiffusion(
        1.0,
        1.0,
        0.001,
        0.001,

        ν,
        κ,

        (2pi, 2pi),

        u_max,
        # 3,
        NY,
        NY,
    )
end

function delta_x(problem::LinearizedThermalDiffusion)
    return problem.domain_size[2] / problem.NY
end


function density(
    q::Quadrature,
    problem::LinearizedThermalDiffusion,
    x::Float64,
    y::Float64,
    timestep::Float64 = 0.0,
)
    ρ_0 = problem.ρ_0
    ρ̃ = problem.ρ̃

    return ρ_0 + ρ̃ * sin(y) * exp(-heat_diffusion(problem) * timestep)
end

function pressure(
    q::Quadrature,
    problem::LinearizedThermalDiffusion,
    x::Float64,
    y::Float64,
    timestep::Float64 = 0.0,
)
    ρ_0 = problem.ρ_0
    θ_0 = problem.θ_0

    # The constant pressure will enforce a temperature given by
    # θ = ρ_0 * θ_0 / ρ
    cs = q.speed_of_sound_squared
    return ρ_0 * θ_0
end

function velocity(
    problem::LinearizedThermalDiffusion,
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
Models the lienarized transverse shear wave,

ρ = ρ_0
u = [ ũ sin(2pi y / Ly), 0 ]
θ = θ_0

It is expected that the temperature θ decays at the rate of exp(- κ t)
"""
struct LinearizedTransverseShearWave <: FluidFlowProblem
    ρ_0::Float64
    θ_0::Float64

    ρ̃::Float64
    ũ::Float64

    ν::Float64
    κ::Float64

    domain_size::Tuple{Float64,Float64}

    u_max::Float64
    NX::Int64
    NY::Int64
end
function LinearizedTransverseShearWave(ν, κ, scale, NY = 8 * scale)
    u_max = 0.01 / scale

    return LinearizedTransverseShearWave(
        1.0,
        1.0,
        0.001,
        0.001,

        ν,
        κ,

        (2pi, 2pi),

        u_max,
        NY,
        NY,
    )
end

function delta_x(problem::LinearizedTransverseShearWave)
    return problem.domain_size[2] / problem.NY
end


function density(
    q::Quadrature,
    problem::LinearizedTransverseShearWave,
    x::Float64,
    y::Float64,
    timestep::Float64 = 0.0,
)
    return problem.ρ_0
end

function pressure(
    q::Quadrature,
    problem::LinearizedTransverseShearWave,
    x::Float64,
    y::Float64,
    timestep::Float64 = 0.0,
)
    return problem.θ_0
end

function velocity(
    problem::LinearizedTransverseShearWave,
    x::Float64,
    y::Float64,
    timestep::Float64 = 0.0,
)
    return [
        problem.ũ * sin(y) * exp(-viscosity(problem) * timestep)
        0.0
    ]
end
