using Test

@test 1 == 1

using lbm

quadratures = [D2Q4(), D2Q5(), D2Q9(), D2Q13(), D2Q17(), D2Q21(), D2Q37()]

include("problems/taylor-green-vortex.jl")
include("problems/poiseuille.jl")

@testset "LBM.jl" begin

    @testset "Multinomial hermite coefficients" begin

        ξ = [1.0 0.0]
        @test hermite(0, ξ) ≈ 1.0
        @test all(hermite(1, ξ) .≈ [1.0 0.0])
        @test all(hermite(2, ξ) .≈ [0.0 0.0; 0.0 -1.0])

        H_3 = hermite(3, ξ)
        @test all(H_3[:, :, 1] .≈ [-2.0 0.0; 0.0 -1.0])
        @test all(H_3[:, :, 2] .≈ [0.0 -1.0; -1.0 0.0])

        H_4 = hermite(4, ξ)

    # Check that our hermite polynomials obey the recurrance rule in 1 dimension
        ξ = 1.0
        for n = 3:4
            H = hermite(n, [ξ])[1]

            @test H ≈ ξ * hermite(n - 1, [ξ])[1] .- (n - 1) * hermite(n - 2, [ξ])[1]
        end

    # Check H4
    # H4 = [
    #     ξ[α_4] .* ξ[α_3] .* ξ[α_2] .* ξ[α_1] .- (
    #         # 6 terms
    #         (ξ[α_1] .* ξ[α_2]) .* δ(α_3, α_4) .+
    #         (ξ[α_1] .* ξ[α_3]) .* δ(α_2, α_4) .+
    #         (ξ[α_1] .* ξ[α_4]) .* δ(α_3, α_2) .+
    #         (ξ[α_2] .* ξ[α_3]) .* δ(α_1, α_4) .+
    #         (ξ[α_2] .* ξ[α_4]) .* δ(α_1, α_3) .+
    #         (ξ[α_3] .* ξ[α_4]) .* δ(α_1, α_2)
    #     ) .+ (δ(α_1, α_2) * δ(α_3, α_4) + δ(α_1, α_3) * δ(α_2, α_4) + δ(α_1, α_4) * δ(α_2, α_3))
    #     for α_1 = 1:D, α_2 = 1:D, α_3 = 1:D, α_4 = 1:D
    # ]
    end

    include("quadrature.jl")



    logspace(start, stop, length) = exp10.(range(start, stop = stop, length = length))
    @testset "moments of the equilibrium distribution of the $q lattice" for q in quadratures
        N = 1
        density_field = fill(1.0, N, N)
        velocity_field = fill(0.0, N, N, lbm.dimension(q))
        temperature_field = fill(1.0, N, N)

        f = lbm.equilibrium(q, density_field, velocity_field, temperature_field)

        @test lbm.density(q, f) ≈ density_field
        @test isapprox(lbm.momentum(q, f), velocity_field, atol = 1e-16)

        equilibrium(ρ, u) = begin
            v = fill(0.0, 1, 1, lbm.dimension(q))
            for d = 1:lbm.dimension(q)
                v[1, 1, d] = u[d]
            end

            lbm.equilibrium(q, fill(ρ, 1, 1), v, 1.0)
        end
        density(ρ, u) = lbm.density(q, equilibrium(ρ, u))
        momentum(ρ, u) = lbm.momentum(q, equilibrium(ρ, u))

        pressure(ρ, u) = LBM.pressure(equilibrium(ρ, u))

        @test density(1.0, [0.0, 0.0]) ≈ fill(1.0, 1, 1)
        @test isapprox(
            momentum(1.0, [0.0, 0.0]),
            fill(0.0, 1, 1, lbm.dimension(q)),
            atol = 1e-16,
        )

    # Todo: analytically check the numerical error bounds

        for u ∈ logspace(2, -9, 12)
            for v ∈ [[0, 0], [u, 0.0], [0.0, u], [u, u]]
                @test isapprox(density(1.0, v)[1, 1], 1.0, rtol = 1e-10, atol = 1e-8)
                for d = 1:lbm.dimension(q)
                    @test isapprox(
                        momentum(1.0, v)[1, 1, d],
                        v[d],
                        atol = 1e-11,
                        rtol = 1e-6,
                    )
                end
            end
        end
    end
end
