@testset "Poiseuille Problem" begin #for q in quadratures
    q = D2Q9()

    scale = 1
    ν = 1.0 / 6.0
    problem = PoiseuilleFlow(ν, scale, static = true)
    f_in = lbm.initialize(q, problem)
    f_out = copy(f_in)
    collision_operator = CollisionModel(SRT, q, problem)

    feq = Array{Float64}(undef, size(f_in, 3))
    f = Array{Float64}(undef, size(f_in, 3))
    u = zeros(dimension(q))


    for x_idx in 2, y_idx = 1:size(f_in, 2)
        f = f_in[x_idx, y_idx, :]

        ρ = lbm.density(q, f)
        lbm.velocity!(q, f, ρ, u)
        T = lbm.temperature(q, f, ρ, u)
        # @show y_idx, u, f
    end

    # @show size(f_in)
    t = 0
    Δt = lbm.delta_t(problem)
    o = size(f_in, 2)
    p = size(f_in, 2) - 1
    # @show f_in[4, 1, :]
    # @show f_in[4, 2, :]
    # @show f_in[4, o, :]
    # @show f_in[4, p, :]
    # @show f_in[2, 1, :] .- f_in[2, o, :]
    # @show f_in
    # @show f_in[4, 2, :] .- f_in[4, p, :]
    # lbm.collide!(collision_operator, q, f_in, f_out, time = t * Δt, problem = problem)
    lbm.collide!(collision_operator, q, time = t * Δt, f_old = f_in, f_new = f_out)
    # f_out = copy(f_in)

    # @show "after collision"
    # @show f_out[2, 1, :] .- f_out[2, o, :]
    # @show f_out
    # @show f_out[:, :, 2]
    # @show f_out[:, :, 4]
    # @show f_out[:, :, 6]
    # @show f_out[:, :, 8]
    # @show f_in[4, 1, :]
    # @show f_in[4, 2, :]
    # @show f_in[4, o, :]
    # @show f_in[4, p, :]
    # @show f_in[4, 2, :] .- f_in[4, p, :]

    # for x_idx = 4, y_idx = 1 : size(f_in, 2)
    #     @show f_in[3, y_idx, :]
    #     @show f_in[4, y_idx, :]
    #     @show f_in[5, y_idx, :]
    # end

    # lbm.apply_boundary_conditions!(q, problem, f_new = f_in, f_old = f_out, time = t * Δt)
    # lbm.apply_boundary_conditions!(q, problem, f_out, f_in, time = t * Δt)
    # @show "after boundary"
    # @show f_out[2, 1, :] .- f_out[2, o, :]
    # @show f_out

        # stream!(quadrature, f_out, f_in)
    lbm.stream!(q, f_new = f_in, f_old = f_out)
    # lbm.stream!(q, f_out, f_in)
    # @show "after stream"
    # @show f_in[2, 1, :]
    # @show f_in[2, o, :]
    # @show f_in[2, 1, :] .- f_in[2, o, :]
    # @show f_in

    # @show f_in[:, :, 2]
    # @show f_in[:, :, 4]
    # @show f_in[:, :, 6]
    # @show f_in[:, :, 8]

    for x_idx in 2, y_idx = 1:size(f_in, 2)
        f = f_in[x_idx, y_idx, :]

        ρ = lbm.density(q, f)
        lbm.velocity!(q, f, ρ, u)
        T = temperature(q, f, ρ, u)
        # @show y_idx, u, f
    end

    # @show f_in[4, 1, :]
    # @show f_in[4, 2, :]
    # @show f_in[4, o, :]
    # @show f_in[4, p, :]
    # @show f_in[4, 2, :] .- f_in[4, p, :]
    # @show "First step"
    for x_idx in 2, y_idx = 1:size(f_in, 2)
        f = f_in[x_idx, y_idx, :]

        ρ = lbm.density(q, f)
        lbm.velocity!(q, f, ρ, u)
        T = temperature(q, f, ρ, u)

        # @show y_idx, u, f
    end

    # @show "Simulate"
    # result = lbm.siumlate(problem, q, base = 200, should_process=false)
    # @show "DONE"

    # TODO: when an a thermal model is considered, the density fluctuates
    # TODO: when a thermal model is considerd we have constant density
end

# struct LatticeNode{Q, T}
#     distributions::Array{Q.d}
#     ρ::Float64
#     u::Vector{Float64}
#     T::Float64


# end
