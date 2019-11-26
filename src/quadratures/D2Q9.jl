# weights(lattice::D2Q9) = [4/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9]'
# abscissaes(lattice::D2Q9) = [0  0; -1  1; -1  0; -1 -1; 0 -1; 1 -1; 1  0; 1  1; 0  1]'
# opposite(lattice::D2Q9) = [1 2 3 4 5 6 7 8 9]
# opposite(lattice::D2Q9, idx) = idx + (9 - 1) / 2
# dimension(lattice::D2Q9) = 2
# order(lattice::D2Q9) = [1, 7, 9, 3, 5, 8, 2, 4, 6]

# D2Q9
# const original_order = [1, 7, 9, 3, 5, 8, 2, 4, 6]

struct D2Q9 <: Quadrature
    abscissae::Array{Int64,2}
    weights::Array{Float64,1}
    speed_of_sound_squared::Float64
end
D2Q9() = D2Q9(
    [0 -1 -1 -1 0 1 1 1 0 0 1 0 -1 -1 -1 0 1 1],
    [4 / 9, 1 / 36, 1 / 9, 1 / 36, 1 / 9, 1 / 36, 1 / 9, 1 / 36, 1 / 9],
    3.0,
)

order(q::D2Q9) = 5
function opposite(q::D2Q9, idx::Int64)
    if idx == 1
        return 1
    end
    if idx <= 5
        return idx + 4
    end
    return idx - 4
end

Base.show(io::IO, q::D2Q9) = show(io, "D2Q9")
