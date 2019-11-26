# Halfway Bounce Back rules
struct BounceBack{D<:Direction,Ints} <: BoundaryCondition
    direction::D
    xs::Ints
    ys::Ints
end

function apply!(bc::BounceBack{<:North}, q::Quadrature, f_new, f_old)
    nx, ny, nf = size(f_new)
    for f_idx = 1:size(f_new, 3)
        # Bounce back the unkown populations
        opposite_f_idx = opposite(q, f_idx)

        for y_idx in bc.ys
            if y_idx + q.abscissae[2, opposite_f_idx] <= ny
                continue
            end
            for x_idx in bc.xs
                f_new[x_idx, y_idx, f_idx] = f_old[x_idx, y_idx, opposite_f_idx]
            end
        end
    end
end

function apply!(bc::BounceBack{<:South}, q::Quadrature, f_new, f_old)
    nx, ny, nf = size(f_new)
    for f_idx = 1:size(f_new, 3)
        opposite_f_idx = opposite(q, f_idx)

        for y_idx in bc.ys
            if y_idx + q.abscissae[2, opposite_f_idx] > 0
                continue
            end
            for x_idx in bc.xs
                f_new[x_idx, y_idx, f_idx] = f_old[x_idx, y_idx, opposite_f_idx]
            end
        end
    end
end

function apply!(bc::BounceBack{<:East}, q::Quadrature, f_new, f_old)
    nx, ny, nf = size(f_new)
    for f_idx = 1:size(f_new, 3)
        # Bounce back the unkown populations
        opposite_f_idx = opposite(q, f_idx)

        for x_idx in bc.xs
            if x_idx + q.abscissae[1, opposite_f_idx] <= nx
                continue
            end
            for y_idx in bc.ys
                f_new[x_idx, y_idx, f_idx] = f_old[x_idx, y_idx, opposite_f_idx]
            end
        end
    end
end
function apply!(bc::BounceBack{<:West}, q::Quadrature, f_new, f_old)
    nx, ny, nf = size(f_new)
    for f_idx = 1:size(f_new, 3)
        # Bounce back the unkown populations
        opposite_f_idx = opposite(q, f_idx)

        for x_idx in bc.xs
            if x_idx + q.abscissae[1, opposite_f_idx] > 0
                continue
            end
            for y_idx in bc.ys
                f_new[x_idx, y_idx, f_idx] = f_old[x_idx, y_idx, opposite_f_idx]
            end
        end
    end
end
