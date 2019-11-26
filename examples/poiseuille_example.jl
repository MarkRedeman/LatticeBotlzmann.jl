using DataFrames
using lbm
using Plots


# relaxation time (BGK model) (tau=sqrt(3/16)+0.5 gives exact solution)
# tau=sqrt(3/16)+0.5;
# nu=(2*tau-1)/6;

let
    # ν = (1.4383 - 0.5) / q.speed_of_sound_squared
    q = D2Q9()
    s = []
    for q in lbm.Quadratures
        stats = DataFrame([Float64[], Float64[]], [:τ, :error_u])
        for τ = 0.5:0.01:2.0
            ν = (τ - 0.5) / q.speed_of_sound_squared

            problem = DecayingShearFlow(ν, 1, static = true)
        # problem = PoiseuilleFlow(ν, 2, static = true)

            result = lbm.siumlate(
                problem,
                q,
                t_end = 1.0,
                should_process = false,
                collision_model = SRT,
            )
            if (!isnan(result.processing_method.df[end].error_u))
                push!(stats, [τ, result.processing_method.df[end].error_u])
            end
        end

    # plot(stats.τ, stats.error_u, yscale=:log10)
    # gui()
        push!(s, stats)
    end

    # p = plot()
    for stats in s_srt
        @show stats.τ[argmin(stats.error_u)]
        plot!(stats.τ, stats.error_u, yscale = :log10, linestyle = :dot)
    end
    gui()

# stats.τ[argmin(stats.error_u)] = 1.438
# 1.438
# julia> @show stats.τ[argmin(stats.error_u)]
# stats.τ[argmin(stats.error_u)] = 1.76
# 1.76
end
# @show result[2]
    # return result


νs = (0.0:2.0:6.0) ./ 6.0
    # νs = [1.0 / 6.0]
scales = [1, 2, 4, 6, 8, 10, 12]#, 8, 10, 12]
scales = [1, 2, 4] #, 8, 16]#, 8, 10, 12]

    # νs = (0.0:0.5:4.0) ./ 6.0
    # scales = [1, 2, 4]
    # νs = [1.0 / 6.0]
for ν in νs, scale in scales
    continue

        # example = TaylorGreenVortexExample(ν, scale, static = false)
    example = DecayingShearFlow(ν, scale, static = true)
        # example = TaylorGreenVortexExample(τ, scale, static = true)

    result = lbm.siumlate(example, quadrature, base = 20)
        # @show result[2]
    push!(stats, [ν, scale, result[2]])
end
    # return stats

    # end
    # end

s = stats
using Plots, LaTeXStrings
nu_scale_error = Array{Float64}([s.nu s.scale map(i -> i.error_u[end], s.stats)])
plot()
for nu_idx = 1:length(νs)
    nu = nu_scale_error[nu_idx*length(scales), 1]
    nu_round = round(nu, digits = 2)

    plot!(
        map(i -> 16 * i, nu_scale_error[(nu_idx-1)*length(scales).+(1:length(scales)), 2]),
        nu_scale_error[(nu_idx-1)*length(scales).+(1:length(scales)), 3],
        label = "nu = $nu_round",
    )
end
plot!(x -> 10.0 * x^(-2), label = "y(x) = 10x^(-2)", linestyle = :dash)
    # plot!(x -> 10.0 * x^(-3), label="y(x) = 10x^(-3)", linestyle=:dash)
    # plot!(x -> 10.0 * x^(-4), label="y(x) = 10x^(-4)", linestyle=:dash)
    # plot!(x -> 10.0 * x^(-5), label="y(x) = 10x^(-5)", linestyle=:dash)
plot!(scale = :log10, legend = :bottomleft, legendfontsize = 5)
gui()

@show -log.(nu_scale_error[2:end, 3] ./ nu_scale_error[1:(end-1), 3]) ./
      log.(nu_scale_error[2:end, 2] ./ nu_scale_error[1:(end-1), 2])


    # nu_idx = 0
    # plot(nu_scale_error[nu_idx * length(scales) .+ (1:length(scales)), 2], nu_scale_error[nu_idx * length(scales) .+ (1:length(scales)), 3])

# Some thoughs: hide the storage of the distributions f inside of an interface
# so that we can do:
# for x in xs, y in ys, z in zs
#  f = fs[x, y ,z]::Vector
#
# the interface could hide it by returning
# return @view f_internal[x, y, z, :]


# Idea: introduce an Initial Value Problem
# problem = FluidFlowProblem
# solution = solve(problem, LBM(Lattice, CollisionModel))
# LBM(Lattice, CollisionModel) can be a solution method
