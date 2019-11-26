@testset "Symmetry requirements of lattice $q" for q in quadratures
    w = q.weights
    ξs = q.abscissae
    D = 2
    Q = length(w)
    δ(α, β) = α == β ? 1.0 : 0.0

    # If we didn't manage to choose correct weights all other tests will likely fail
    @test sum(q.weights) ≈ 1.0

    @testset "Symmetry $n" for n = 1:div(lbm.order(q), 2)
        H = sum([
            q.weights[f_idx] * hermite(n, q.abscissae[:, f_idx], q)
            for f_idx = 1:length(q.weights)
        ])

        @test all(isapprox.(H, 0.0, atol = 1e-15))
    end
end

@testset "moments of the equilibrium distribution of the $q lattice" for q in quadratures
    N = 1
    ρ = fill(1.0, N, N)
    u = fill(0.001, N, N, lbm.dimension(q))
    T = fill(1.0, N, N)

    @testset "equilibira" begin
        f = lbm.equilibrium(q, ρ, u, T)

        @test all(lbm.density(q, f) .≈ ρ)

        @test all(lbm.momentum(q, f) ./ ρ .≈ u)
        @test all(isapprox.(lbm.temperature(q, f, ρ, u), T, atol = 1e-5))
        # @code_warntype lbm.equilibrium(q, ρ, u, T)
        @inferred lbm.equilibrium(q, ρ, u, T)
    end

    @testset "temperature" begin
        N = 1
        ρ = fill(1.0, N, N)
        u = fill(0.0, N, N, lbm.dimension(q))
        T = fill(1.0, N, N)
        f = lbm.equilibrium(q, ρ, u, T)
        # @show lbm.temperature(q, f, ρ, u) * q.speed_of_sound_squared
        @test all(lbm.temperature(q, f, ρ, u) .≈ T)
    end

    @testset "moments of equilibira in a single point" begin
        ρ = 1.0
        u = [0.1, 0.1]
        T = 1.0
        f = lbm.equilibrium(q, ρ, u, T)

        @test lbm.density(q, f) ≈ ρ
        @test all(lbm.momentum(q, f) ./ ρ .≈ u)

        @inferred lbm.density(q, f)
        @inferred lbm.momentum(q, f)

        j = [0.0, 0.0]
        lbm.momentum!(q, f, j)
        @test all(j .≈ u)
    end

    @testset "SRT collision operator" begin
        N = 1
        ρ = fill(1.0, N, N)
        u = fill(0.01, N, N, lbm.dimension(q))
        T = fill(1.05, N, N)

        f = lbm.equilibrium(q, ρ, u, T)
        f_out = copy(f)
        # @show q.speed_of_sound_squared
        # @show 1 / q.speed_of_sound_squared
        # @show ρ_1 = lbm.density(q, f)
        # @show j_1 = lbm.momentum(q, f)[:]
        # @show lbm.temperature(q, f, ρ_1, j_1 ./ ρ_1)

        # @show "COLLIDE"
        # When τ = 1.0 we immediatly get the equilibrium distribution
        lbm.collide!(SRT(1.0), q, f_old = f, f_new = f_out, time = 0.0)
        # @show ρ = lbm.density(q, f)
        # @show j = lbm.momentum(q, f)
        # @show lbm.temperature(q, f, ρ, j ./ ρ)

        f_ = f_out[1, 1, :]
        u_ = [0.0, 0.0]
        ρ_ = lbm.density(q, f_)
        lbm.velocity!(q, f_, ρ_, u_)
        T_ = lbm.temperature(q, f_, ρ_, u_)
        @test all(isapprox(f, f_out, atol = 1e-4))

        @inferred lbm.equilibrium(q, ρ, u, T)
        @inferred lbm.equilibrium(q, ρ, u, 1.0)
        # @code_warntype lbm.equilibrium(q, ρ, u, 1.0);
        # @code_warntype collide(collision_model, q, f)
        @inferred lbm.collide!(SRT(1.0), q, f_old = f, f_new = f_out, time = 0.0)
    end

    @testset "Stream and collide" begin
        N = 1
        ρ = fill(1.0, N, N)
        u = fill(0.1, N, N, lbm.dimension(q))
        T = fill(1.0, N, N)

        f = lbm.equilibrium(q, ρ, u, T)

        f_inn = copy(f)
        f_out = copy(f)
        f_fff = f
        τ = 1.0
        for t = 0:3
            lbm.collide!(SRT(τ), q, f_old = f_inn, f_new = f_out, time = 0.0)
            lbm.stream!(q, f_new = f_inn, f_old = f_out)
        end

        @test all(isapprox.(f, f_inn, atol = 1e-5))

        # @test all(isapprox.(H, 0.0))
    end
end

@testset "Opposite directions of $q" for q in quadratures
    for idx = 1:length(q.weights)
        opposite_idx = opposite(q, idx)

        @test all(q.abscissae[:, idx] .+ q.abscissae[:, opposite_idx] .== 0)
    end
end

# There was a regression where the stream function did not correctly set a
# period boundary when using a multispeed quadrature that had a speed exceeding
# the amount of distributions in the x or y direction
@testset "Streaming quadratures with multispeeds" for q in quadratures
    q = D2Q13()
    N = 1
    ρ = fill(1.0, N, N)
    u = fill(0.1, N, N, lbm.dimension(q))
    T = fill(1.0, N, N)

    feq = lbm.equilibrium(q, 1.0, [0.0, 0.0], 1.0)
    f = Array{Float64}(undef, 1, 1, length(q.weights))
    f[1, 1, :] = feq

    f_new = copy(f)

    lbm.stream!(q, f_new = f_new, f_old = f)

    @test all(isapprox.(f, f_new))
end
