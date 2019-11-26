# D2Q21
# See lectures in LBM

# Same as D2Q17, but vertical and horizontal vectors are extended by one

struct D2Q21 <: Quadrature
    abscissae::Array{Int64,2}
    weights::Array{Float64,1}
    speed_of_sound_squared::Float64
end
function D2Q21()
    cs = sqrt(2 / 3)
    w_0 = 91 / 324
    w_1 = 1 / 12
    w_2 = 2 / 27
    w_3 = 7 / 360
    w_4 = 1 / 432
    w_5 = 0 / 36
    w_6 = 1 / 1620

    # https://strathprints.strath.ac.uk/29195/1/Zhang_YH_Pure_Accuracy_analysis_of_high_order_lattice_Boltzmann_models_for_rarefied_gas_flows_21_Jan_2011.pdf
    w_0 = 91 / 324
    w_1 = 1 / 12
    w_2 = 2 / 27
    w_3 = 7 / 360
    w_4 = 1 / 432
    w_5 = 1 / 1620
    w_6 = 0
    cs = 3 / 2

    return D2Q21(
        [   #A #B          #C          #D          #E          #F
            0 1 -1 0 0 1 -1 1 -1 2 -2 0 0 2 -2 2 -2 3 -3 0 0 3 -3 -3 3
            0 0 0 1 -1 1 -1 -1 1 0 0 2 -2 2 -2 -2 2 0 0 3 -3 3 -3 3 -3
        ],
        [
            w_0,
            w_1,
            w_1,
            w_1,
            w_1,
            w_2,
            w_2,
            w_2,
            w_2,
            w_3,
            w_3,
            w_3,
            w_3,
            w_4,
            w_4,
            w_4,
            w_4,
            w_5,
            w_5,
            w_5,
            w_5,
            w_6,
            w_6,
            w_6,
            w_6,
        ],
        cs,
    )
end

order(q::D2Q21) = 7
Base.show(io::IO, q::D2Q21) = show(io, "D2Q21")
