function test_convergence(problem::Problem, scales) where {Problem<:FluidFlowProblem} end

function analyze_convergence(q::Quadrature, p, viscosity::Float64, N = 2, t_end = 2pi)
    stats = Vector{NamedTuple{(:nu, :scale, :u_error),Tuple{Float64,Int,Float64}}}()

    for scale = 0:N
        problem = p(2^scale, viscosity)
        result = siumlate(problem, q, should_process = false, t_end = t_end)

        push!(
            stats,
            (
                nu = viscosity,
                scale = 16 * 2^scale,
                u_error = result.processing_method.df[end].u_error,
            ),
        )
    end

    plot_convergence(stats, viscosity)
    gui()

    # @show -log.(
    #     stats.u_error[2:end] ./ stats.u_error[1:(end - 1)]
    # ) ./ log.(
    #     stats.nu[2:end] ./ stats.nu[1:(end - 1)]
    # )

    return stats
end

function plot_convergence(stats, viscosity)
    nu_round = round(viscosity, digits = 2)

    p = plot()
    plot!(p, getfield.(stats, :scale), getfield.(stats, :u_error), label = "nu = $nu_round")
    plot!(
        p,
        x -> stats[1].u_error * 10.0 * x^(-2),
        label = "y(x) = 10x^(-2)",
        linestyle = :dash,
    )
    plot!(p, scale = :log10, legend = :bottomleft, legendfontsize = 5)
    return p
end
