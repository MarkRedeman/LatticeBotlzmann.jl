"""
Show all quadratures currently implemented in lbm.jl
"""
function plot_quadratures(color_each_group_separatly = false)
    plot(map(q -> plot_quadrature(q, color_each_group_separatly), Quadratures)...)
end

"""
Show the arrows of a quadrature, useful for debugging
"""
function plot_quadrature(q::Quadrature, color_each_group_separatly = false)
    unique_weights = unique(q.weights)
    colors = distinguishable_colors(length(unique_weights))

    limit = max(q.abscissae[1, :]...) + 1

    p = plot(
        title = string(q),
        legend = nothing,
        framestyle = :grid,
        tickfontcolor = :gray,
        xlims = (-limit, limit),
        ylims = (-limit, limit),
        gridstyle = :dash,
    )

    for f_idx = 1:length(q.weights)
        color_idx = findfirst(isequal(q.weights[f_idx]), unique_weights)

        if q.abscissae[:, f_idx] != [0, 0]
            quiver!(
                p,
                [0.0],
                [0.0],
                quiver = ([q.abscissae[1, f_idx]], [q.abscissae[2, f_idx]]),
                color = color_each_group_separatly ? colors[color_idx] : :black,
            )
        end
    end

    for f_idx = 1:length(q.weights)
        if q.abscissae[:, f_idx] != [0, 0]
            continue
        end
        color_idx = findfirst(isequal(q.weights[f_idx]), unique_weights)
        scatter!(
            p,
            [(q.abscissae[1, f_idx], q.abscissae[2, f_idx])],
            # color=colors[color_idx],
            color = color_each_group_separatly ? colors[color_idx] : :gray,
            foreground_color_border = :black,
            markersize = 2,
        )
    end

    return p
end
