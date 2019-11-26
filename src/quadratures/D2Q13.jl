struct D2Q13 <: Quadrature
    abscissae::Array{Int64,2}
    weights::Array{Float64,1}
    speed_of_sound_squared::Float64

end
function D2Q13()
    # High Knudsen Number Thermal Flows with the D2Q13 Lattice Boltzmann Model
    # P. Lopez & Y. Bayazitoglu
    w_0 = 3 / 8
    w_1 = 1 / 12
    w_2 = 1 / 16
    w_3 = 1 / 96
    cs = 2

    return D2Q13(
        [   #A #B          #C          #D
            0 1 -1 0 0 1 -1 1 -1 2 -2 0 0
            0 0 0 1 -1 1 -1 -1 1 0 0 -2 2
        ],
        [w_0, w_1, w_1, w_1, w_1, w_2, w_2, w_2, w_2, w_3, w_3, w_3, w_3],
        cs,
    )
end

order(q::D2Q13) = 5
Base.show(io::IO, q::D2Q13) = show(io, "D2Q13")
