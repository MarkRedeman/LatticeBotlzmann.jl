abstract type InitializationStrategy2 end
InitializationStrategy(problem) = FromEquilibrium2()

struct ZeroInitialCondition2 <: InitializationStrategy2 end

struct FromEquilibrium2 <: InitializationStrategy2 end

struct WithNonEquilibrium2 <: InitializationStrategy2 end

struct MeiEtAl2 <: InitializationStrategy2 end

# function initial_equilibrium(quadrature::Quadrature, problem::FluidFlowProblem, x::Float64, y::Float64)
#     return equilibrium(
#         quadrature,
#         lattice_density(quadrature, problem, x, y),
#         lattice_velocity(quadrature, problem, x, y),
#         lattice_temperature(quadrature, problem, x, y)
#     )
# end

# function initial_condition(q::Quadrature, problem::FluidFlowProblem, x::Float64, y::Float64)
#     initial_equilibrium(q, problem, x, y)
# end

# function initialize(quadrature::Quadrature, problem::FluidFlowProblem, cm = SRT)
#     f = Array{Float64}(undef, problem.NX, problem.NY, length(quadrature.weights))

#     x_range, y_range = range(problem)
#     for x_idx in 1:problem.NX, y_idx in 1:problem.NY
#         f[x_idx, y_idx, :] = initial_equilibrium(
#             quadrature,
#             problem,
#             x_range[x_idx],
#             y_range[y_idx]
#         )
#     end

#     return f
# end
