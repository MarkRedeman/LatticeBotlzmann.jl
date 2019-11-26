struct TaylorGreenVortexExample <: lbm.FluidFlowProblem
    rho_0::Float64
    u_max::Float64
    ν::Float64
    NX::Int64
    NY::Int64
    k_x::Float64
    k_y::Float64
    domain_size::Tuple{Float64,Float64}
    static::Bool
end
function TaylorGreenVortexExample(
    ν = 1.0 / 6.0,
    scale = 2,
    NX = 16 * scale,
    NY = NX,
    domain_size = (2pi, 2pi);
    static = true,
)
    u_max = 0.02 / scale
    Re = NX * u_max / ν
    @show Re
    p = TaylorGreenVortexExample(
        1.0,
        u_max,
        ν,
        NX,
        NY,
        domain_size[1] / NX,
        domain_size[2] / NY,
        domain_size,
        static,
    )

    if (p.u_max > sqrt(2 / 3) * delta_x(p) / delta_t(p))
        @warn p.u_max, sqrt(2 / 3) * delta_x(p) / delta_t(p)
    end
    p
end

function density(
    q::Quadrature,
    tgv::TaylorGreenVortexExample,
    x::Float64,
    y::Float64,
    timestep::Float64 = 0.0,
)
    return pressure(q, tgv, x, y, timestep)

    # If not athermal
    return 1.0
end

function pressure(
    q::Quadrature,
    tgv::TaylorGreenVortexExample,
    x::Float64,
    y::Float64,
    timestep::Float64 = 0.0,
)
    # return 1.0
    P = -0.25 *
        tgv.rho_0 *
        ((tgv.k_y / tgv.k_x) * cos(2.0 * x) + (tgv.k_x / tgv.k_y) * cos(2.0 * y)) *
        decay(tgv, x, y, timestep)^2

    return 1.0 + q.speed_of_sound_squared * tgv.u_max^2 * P
end
function pressure_tensor(
    q::Quadrature,
    problem::TaylorGreenVortexExample,
    x::Float64,
    y::Float64,
    time::Float64 = 0.0,
)
# q.speed_of_sound_squared * ()
    return pressure(q, problem, x, y, time) * [1.0 0.0; 0.0 1.0] +
           deviatoric_tensor(q, problem, x, y, time)
    return pressure(q, problem, x, y, time) * [1.0 0.0; 0.0 1.0]
    # - q.speed_of_sound_squared [
    #     0.0 1.0
    #     1.0 0.0
    # ]
end
function deviatoric_tensor(
    q::Quadrature,
    problem::TaylorGreenVortexExample,
    x::Float64,
    y::Float64,
    time::Float64 = 0.0,
)
    a = acceleration(problem, x, y, time)

    ν = problem.ν
    # ν = viscosity(problem)

    return ν * [
        2 * a[1, 1] a[1, 2] + a[2, 1]
        a[1, 2] + a[2, 1] 2 * a[2, 2]
    ]
    return q.speed_of_sound_squared^2 *
           ν *
           [
               2 * a[1, 1] a[1, 2] + a[2, 1]
               a[1, 2] + a[2, 1] 2 * a[2, 2]
           ]
    return q.speed_of_sound_squared *
           ν *
           [
               2 * a[1, 1] a[1, 2] + a[2, 1]
               a[1, 2] + a[2, 1] 2 * a[2, 2]
           ]


    a = 1.0
    A = 1.0
    b = 1.0
    B = -1.0

    σ_xx = -2 * a * A * sin(a * x) * cos(b * y)
    σ_yx = (a * B + b * A) * cos(a * x) * cos(b * y)
    σ_xy = (a * B + b * A) * cos(a * x) * cos(b * y)
    σ_yy = -2 * b * B * sin(a * x) * cos(b * y)


    u_x = -a * A * sin(x) * sin(y)
    u_y = b * A * cos(x) * cos(y)
    v_x = a * B * cos(x) * cos(y)
    v_y = -b * B * sin(x) * sin(y)

    # @show [σ_xx σ_yx; σ_xy σ_yy]
    σ_xx = 2 * u_x
    σ_yx = v_x + u_y
    σ_xy = u_y + v_x
    σ_yy = 2 * v_y
    # @show [σ_xx σ_yx; σ_xy σ_yy]

    F = decay(problem, x, y, time)
    return F * q.speed_of_sound_squared * problem.ν * [σ_xx σ_yx; σ_xy σ_yy]
end

function velocity(
    tgv::TaylorGreenVortexExample,
    x::Float64,
    y::Float64,
    timestep::Float64 = 0.0,
)
    return decay(tgv, x, y, timestep) * [
        -sqrt(tgv.k_y / tgv.k_x) * cos(x) * sin(y),
        sqrt(tgv.k_x / tgv.k_y) * sin(x) * cos(y),
    ]
end
function acceleration(
    tgv::TaylorGreenVortexExample,
    x::Float64,
    y::Float64,
    timestep::Float64 = 0.0,
)
    a = 1.0
    A = 1.0
    b = 1.0
    B = -1.0
    u_x = -sqrt(tgv.k_y / tgv.k_x) * sin(x) * sin(y)
    u_y = sqrt(tgv.k_y / tgv.k_x) * cos(x) * cos(y)
    v_x = sqrt(tgv.k_x / tgv.k_y) * cos(x) * cos(y)
    v_y = sqrt(tgv.k_x / tgv.k_y) * sin(x) * sin(y)

    u_x = -a * A * sin(x) * sin(y)
    u_y = b * A * cos(x) * cos(y)
    v_x = a * B * cos(x) * cos(y)
    v_y = -b * B * sin(x) * sin(y)

    return decay(tgv, x, y, timestep) * [u_x v_x; u_y v_y]

    return decay(tgv, x, y, timestep) * [
        -u_max * sqrt(tgv.k_y / tgv.k_x) * cos(x) * sin(y),
        u_max * sqrt(tgv.k_x / tgv.k_y) * sin(x) * cos(y),
    ]
end
function decay(tgv::TaylorGreenVortexExample, x::Float64, y::Float64, timestep::Float64)
    return tgv.static ? 1.0 : exp(-2.0 * timestep * viscosity(tgv))
end

function force(tgv::TaylorGreenVortexExample, x::Float64, y::Float64, time::Float64 = 0.0)
    return tgv.static ? 2 * viscosity(tgv) * velocity(tgv, x, y, 0.0) : [0.0 0.0]
end

has_external_force(problem::TaylorGreenVortexExample) = problem.static
