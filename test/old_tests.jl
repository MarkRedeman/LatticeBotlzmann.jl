@testset "2 dimensional lattice boltzmann" begin
    return

    const tol = 1e-8

    @test sum(LBM.weights) ≈ 1.0
    @test LBM.density(LBM.weights) ≈ 1.0
    @test LBM.velocity(LBM.weights) ≈ [0.0, 0.0]

    # When constructing an equilibrium function, it should have
    # conserve the density and velocity that were used constructing it
    @testset "moments of the equilibrium distribution" begin
        # Wrappers to easily calculate density and velocity of an equilibrium
        equilibrium(ρ, u) = LBM.equilibrium(ρ, u)
        density(ρ, u) = LBM.density(equilibrium(ρ, u))
        velocity(ρ, u) = LBM.velocity(equilibrium(ρ, u))

        pressure(ρ, u) = LBM.pressure(equilibrium(ρ, u))

        @test density(1.0, [0.0, 0.0]) ≈ 1.0
        @test velocity(1.0, [0.0, 0.0]) ≈ [0.0, 0.0]

        # Todo: analytically check the numerical error bounds
        for u ∈ logspace(2, -9, 12)
            for v ∈ [[0, 0], [u, 0.0], [0.0, u], [u, u]]
                @test density(1.0, v) ≈ 1.0
                @test velocity(1.0, v) ≈ v
            end
        end
    end

    @testset "single relaxation bgk collision" begin
        @testset "equilibrium distributions should remain in equilibrium" begin
            f = LBM.equilibrium(1.0, [0.0, 0.0])

            @test LBM.collide(f, 1.0) ≈ f
        end

        # Create a distribution that is perturbed
        f = LBM.equilibrium(1.0, [0.0, 0.0])
        f[7] *= 3
        f[:] = f[:] / sum(f)
        @test LBM.velocity(f) ≈ [1.0 / 5.5, 0.0]
        @test f ≉ LBM.equilibrium(f)
        @test copy(f) ≉ LBM.collide(f, 1.1)

        # Check that we can collide multiple distributions using broadcastign
        fs = [LBM.equilibrium(1.0, [0.0, 0.0]) for x = 1:10, y = 1:10]
        @test isapprox(norm((fs.-LBM.collide.(fs, 1.0))[:]), 0, atol = 1e-14)
    end

    @testset "streaming" begin
        f = LBM.equilibrium(1.0, [0.0, 0.0])

        @test LBM.stream(
            [[1, 2, 3] [4, 5, 6] [7, 8, 9]],
            [0, 0],
        ) == [[1, 2, 3] [4, 5, 6] [7, 8, 9]]
        @test LBM.stream(
            [[1, 2, 3] [4, 5, 6] [7, 8, 9]],
            [0, 1],
        ) == [[7, 8, 9] [1, 2, 3] [4, 5, 6]]
        @test LBM.stream(
            [[1, 2, 3] [4, 5, 6] [7, 8, 9]],
            [1, 0],
        ) == [[3, 1, 2] [6, 4, 5] [9, 7, 8]]
        @test LBM.stream(
            [[1, 2, 3] [4, 5, 6] [7, 8, 9]],
            [1, 1],
        ) == [[9, 7, 8] [3, 1, 2] [6, 4, 5]]


        # Check that if we have a 2x2 square with periodic boundaries then
        # streaming twice should give the same result as no streaming
        for size = 2:10
            fs = [LBM.equilibrium(1.0, [Float64(x), Float64(y)]) for x = 1:size, y = 1:size]
            const original = deepcopy(fs)

            for idx = 1:size
                fs = LBM.stream(fs)
            end

            @test fs == original
        end
    end

    @testset "a small simulation" begin
        const Re = 10.0
        const U = 0.05
        const lx, ly = 10, 10
        const ν = U * lx / Re
        const ω = 1.0 / (3 * ν + 1 / 2.0)

        ρ, u = 1.0, [0.0, 0.0]
        fs = [LBM.equilibrium(ρ, u) for x = 1:lx, y = 1:ly]

        const original_density = sum(LBM.density.(fs))
        const original_velocity = sum(LBM.velocity.(fs))

        @test original_density ≈ lx * ly
        @test original_velocity ≈ [0.0, 0.0]

        # 100x100 1000
        # 52.668638 seconds (539.06 M allocations: 35.709 GB, 27.51% gc time)
        @time @testset "Naive streaming" begin
            @inbounds for t = 1:1000
                fs .= LBM.collide.(fs, ω)
                fs = LBM.stream(fs)
            end

            @test sum(LBM.density.(fs)) ≈ original_density
            @test isapprox(sum(LBM.velocity.(fs)), original_velocity, atol = 1e-16)
        end

        # # 100x100 1000
        # # 49.988281 seconds (500.07 M allocations: 33.235 GB, 21.75% gc time)
        # @time @testset "Efficient memory streaming" begin
        #     @inbounds for t = 1:1000
        #         fs .= LBM.collide_and_swap.(fs, ω)
        #         fs = LBM.stream_by_swapping(fs)
        #     end

        #     @test sum(LBM.density.(fs)) ≈ original_density
        #     @test isapprox(sum(LBM.velocity.(fs)), original_velocity, atol=1e-16)
        # end

        # # 100x100 1000
        # # 44.297459 seconds (500.06 M allocations: 33.235 GB, 22.48% gc time)
        # @time @testset "Efficient looping and memory streaming" begin
        #     @inbounds for t = 1:1000
        #         fs = LBM.collide_and_stream(fs, ω)
        #     end

        #     @test sum(LBM.density.(fs)) ≈ original_density
        #     @test isapprox(sum(LBM.velocity.(fs)), original_velocity, atol=1e-16)
        # end

        fs_as_array = zeros(9, lx, ly)
        for x = 1:lx, y = 1:ly
            fs_as_array[:, x, y] = fs[x, y]
        end
        # about 120% faster
        @inbounds for t = 1:150
            for x = 1:lx, y = 1:ly
                fs_as_array[:, x, y] = LBM.collide(fs_as_array[:, x, y]', ω)
            end

            # STREAMING STEP
            for idx = 1:9
                fs_as_array[idx, :, :] =
                    LBM.stream(fs_as_array[idx, :, :], LBM.abscissae[:, idx])
            end
        end
        @test sum(fs_as_array) ≈ original_density
    end

    @testset "the lattice is applicable for efficient streaming" begin
        # For details see
        # [Technical report: How to implement your DdQq dynamics with only q variables per node (instead of 2q)]
        # http://optilb.com/openlb/wp-content/uploads/2011/12/olb-tr1.pdf

        abscissae = LBM.abscissae

        for idx = 2:5
            @test abscissae[1, idx] < 0 || (abscissae[1, idx] == 0 && abscissae[2, idx] < 0)
        end

        # f = LBM.weights
        # u = LBM.velocity(f)
        # sum(f .* abscissae, 2) / density(f)
    end


    # @testset "plotting functionallity" begin
    #     pressure = (x, y) -> LBM.pressure(fs[x, y])
    # end

    # @testset "initial conditions" begin
    #     initial_density(x, y) = 1
    #     initial_velocity(x, y) = [sin(x), cos(y)]

    #     for x = linspace(0, 2π, 10), y = linspace(0, 2π, 10)
    #         f = LBM.initial_distributions(density=initial_density, velocity=initial_velocity)
    #     end
    # end
end
