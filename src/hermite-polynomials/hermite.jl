# functions to determine hermite coefficients of a function

hermite(::Type{Val{0}}) = 1.0
hermite(::Type{Val{0}}, ξ) = hermite(Val{0})
hermite(::Type{Val{1}}, ξ) = ξ .* hermite(Val{0})

hermite(::Type{Val{0}}, ξ, q) = hermite(Val{0})
hermite(::Type{Val{1}}, ξ, q) = ξ .* hermite(Val{0})

hermite(N::Int, ξ) = hermite(Val{N}, ξ)
hermite(N::Int, ξ, q::Quadrature) = hermite(Val{N}, ξ, q)

function hermite(::Type{Val{2}}, ξ)
    D = length(ξ)
    return [ξ[β] * ξ[α] - δ(α, β) for α = 1:D, β = 1:D]
end
function hermite(::Type{Val{3}}, ξ)
    D = length(ξ)
    return [
        ξ[γ] * ξ[β] * ξ[α] - (ξ[α] * δ(β, γ) + ξ[β] * δ(α, γ) + ξ[γ] * δ(α, β))

        for α = 1:D, β = 1:D, γ = 1:D
    ]
end

function hermite(::Type{Val{4}}, ξ)
    D = length(ξ)
    return [
        ξ[ϵ] *
        ξ[γ] *
        ξ[β] *
        ξ[α] - (ξ[α] * ξ[β] * δ(γ, ϵ) +
                ξ[α] * ξ[γ] * δ(β, ϵ) +
                ξ[α] * ξ[ϵ] * δ(β, γ) +

                ξ[β] * ξ[γ] * δ(α, ϵ) +
                ξ[β] * ξ[ϵ] * δ(α, γ) +
                ξ[γ] * ξ[ϵ] * δ(α, β)) + (δ(α, β) * δ(γ, ϵ) +
                                          δ(α, γ) * δ(β, ϵ) +
                                          δ(α, ϵ) * δ(β, γ))

        for α = 1:D, β = 1:D, γ = 1:D, ϵ = 1:D
    ]
end

function hermite(::Type{Val{2}}, ξ, q)
    D = length(ξ)
    cs = 1 / q.speed_of_sound_squared
    return [ξ[β] * ξ[α] - cs * δ(α, β) for α = 1:D, β = 1:D]
end
function hermite(::Type{Val{3}}, ξ, q)
    D = length(ξ)
    cs = 1 / q.speed_of_sound_squared
    return [
        ξ[γ] * ξ[β] * ξ[α] - cs * (ξ[α] * δ(β, γ) + ξ[β] * δ(α, γ) + ξ[γ] * δ(α, β))

        for α = 1:D, β = 1:D, γ = 1:D
    ]
end

function hermite(::Type{Val{4}}, ξ, q)
    D = length(ξ)
    cs = 1 / q.speed_of_sound_squared
    return [
        ξ[ϵ] *
        ξ[γ] *
        ξ[β] *
        ξ[α] - cs * (ξ[α] * ξ[β] * δ(γ, ϵ) +
                     ξ[α] * ξ[γ] * δ(β, ϵ) +
                     ξ[α] * ξ[ϵ] * δ(β, γ) +

                     ξ[β] * ξ[γ] * δ(α, ϵ) +
                     ξ[β] * ξ[ϵ] * δ(α, γ) +
                     ξ[γ] * ξ[ϵ] * δ(α, β)) + cs^2 * (δ(α, β) * δ(γ, ϵ) +
                                                      δ(α, γ) * δ(β, ϵ) +
                                                      δ(α, ϵ) * δ(β, γ))

        for α = 1:D, β = 1:D, γ = 1:D, ϵ = 1:D
    ]
end
