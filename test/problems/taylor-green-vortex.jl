import LinearAlgebra: I, tr

@testset "Taylor Green Vortex Problem" begin #for q in quadratures
    q = D2Q9()

    scale = 1
    ν = 1.0 / 6.0
    problem = TaylorGreenVortexExample(ν, scale)
    f = lbm.initialize(q, problem)

    ρ = lbm.density(q, f)
    j = lbm.momentum(q, f)
    T = lbm.temperature(q, f, ρ, j ./ ρ)

    # TODO: when an a thermal model is considered, the density fluctuates
    # TODO: when a thermal model is considerd we have constant density
    @testset "Hermite Pressure component" begin
        q = D2Q9()
        D = dimension(q)
        N = div(lbm.order(q), 2)
        Hs = [
            [hermite(Val{n}, q.abscissae[:, i], q) for i = 1:length(q.weights)] for n = 1:N
        ]

        scale = 1
        ν = 1.0 / 6.0
        problem = DecayingShearFlow(ν, scale, static = false)
        problem = TaylorGreenVortexExample(ν, scale, static = true)

        x = 0.30
        y = 0.30
        f = lbm.initial_equilibrium(q, problem, x, y)

        ρ = sum(f)

        τ = problem.ν
        τ = q.speed_of_sound_squared * lbm.lattice_viscosity(problem) + 0.5

        a_f = [sum([f[idx] * Hs[n][idx] for idx = 1:length(q.weights)]) for n = 1:N]

        @show a_f[2]
        P = a_f[2] - (a_f[1] * a_f[1]') / ρ - ρ * I(2)
        P = P * (1 - 1 / (2 * τ))

        σ_lb = (P - I(2) * tr(P) / (D)) / (problem.u_max)
        σ = lbm.deviatoric_tensor(q, problem, x, y, 0.0)

        @show σ σ_lb
        @show σ - σ_lb
        @show σ ./ σ_lb
        @show σ_lb ./ σ
        # return


        # @inbounds for f_idx = 1 : nf
        #     f_out[x_idx, y_idx, f_idx] = q.weights[f_idx] * (
        #         ρ +
        #         cs * ρ * sum(u .* Hs[1][f_idx]) +
        #         sum([
        #             cs^n * dot(a_coll[n], Hs[n][f_idx]) / (factorial(n))
        #             for n = 2:N
        #         ])
        #     )
        # end
    end
    @testset "Pressure component" begin
        scale = 1
        ν = 2.0 / 6.0
        problem = TaylorGreenVortexExample(ν, scale, static = false)
        problem = DecayingShearFlow(ν, scale, static = false)

        x = 0.30
        y = 0.30
        f = lbm.initial_equilibrium(q, problem, x, y)

        ρ = lbm.density(q, f)

        cs = 1 / q.speed_of_sound_squared
        # cs = q.speed_of_sound_squared
        d_u = problem.u_max * lbm.acceleration(problem, x, y, 0.0)
        τ = problem.ν
        @show d_u
        τ = q.speed_of_sound_squared * lbm.lattice_viscosity(problem) + 0.5
        @show τ
        for f_idx = 1:length(f)
            f[f_idx] +=
                -(q.weights[f_idx] * ρ * τ / cs) *
                sum(lbm.hermite(Val{2}, q.abscissae[:, f_idx], q) .* d_u)
        end
        @show 1 / delta_t(problem)
        @show 1 / lbm.delta_x(problem)
        @show 1 / viscosity(problem)


        j = lbm.momentum(q, f)
        u = j ./ ρ
        T = lbm.temperature(q, f, ρ, u)
        @show lbm.momentum_flux(q, f, ρ, u)
        @show lbm.pressure(q, f, ρ, u)

        @show -(
            q.speed_of_sound_squared / problem.ν
        ) * (lbm.momentum_flux(q, f, ρ, u) - I(2) * lbm.pressure(q, f, ρ, u)) ./
              problem.u_max
        @show ρ * u * u'

        @show lbm.deviatoric_tensor(q, problem, x, y, 0.0)
        # @show lbm.pressure_tensor(q, problem, x, y, 0.0)
        @show problem.u_max
        @show lbm.pressure_tensor(q, problem, x, y, 0.0) -
              lbm.deviatoric_tensor(q, problem, x, y, 0.0)
        @show lbm.pressure(q, f, ρ, u) - lbm.pressure(q, problem, x, y, 0.0)


        @warn "Check deviatoric stress"
        σ_lb = lbm.momentum_flux(q, f, ρ, u) - I(2) * lbm.pressure(q, f, ρ, u)
        σ = lbm.deviatoric_tensor(q, problem, x, y, 0.0)
        @show σ_lb σ
        @show (1 - 1 / (2 * τ))
        @show 0.5 * (1 / cs) * (1 / problem.u_max) * σ_lb / (1 - 1 / (2 * τ))
        # @show (1 / (problem.ν^2 * problem.u_max)) * σ_lb / (1 + 1 / (2 * problem.ν))
        # @show (q.speed_of_sound_squared / (problem.ν * problem.u_max) ) * (
        #     σ_lb
        # )
        # @show σ

        # @show (1 / (problem.ν^2 * problem.u_max)) * σ_lb / (1 + 1 / (2 * problem.ν)) + σ/2
        # @show (q.speed_of_sound_squared / (problem.ν * problem.u_max) ) * (
        #     σ_lb
        # ) + σ
        # return

    end
end

# struct LatticeNode{Q, T}
#     distributions::Array{Q.d}
#     ρ::Float64
#     u::Vector{Float64}
#     T::Float64


# end
