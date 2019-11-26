abstract type BoundaryCondition end

"""
Apply a group of boundary conditions after streaming
"""
function apply!(
    boundary_conditions::Vector{<:BoundaryCondition},
    q::Quadrature,
    f_new,
    f_old;
    time = 0.0,
)
    for boundary_condition in boundary_conditions
        apply!(boundary_condition, q, f_new, f_old, time = time)
    end
end

function apply!(bc::BoundaryCondition, q::Quadrature, f_new, f_old; time = 0.0)
    apply!(bc::BoundaryCondition, q::Quadrature, f_new, f_old)
end

"""
Often we need special rules when a boundary is placed in a specific direction

Currently we only support a few 2D boundary conditions
In the future we could add additional directions to support specific corners
"""
abstract type Direction end
struct North <: Direction end
struct East <: Direction end
struct South <: Direction end
struct West <: Direction end

include("boundary-conditions/bounce-back.jl")
include("boundary-conditions/moving-wall.jl")
