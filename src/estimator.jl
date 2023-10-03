mutable struct Estimator
    const H::OpString
    const legs_first::Matrix{Leg}
    const legs_last::Matrix{Leg}
    const ψ0::Matrix{Bool}
    ######################## numbers
    const Λ::Int64
    n::Int64
    β::Float64
    ξ::Float64
    ######################## Sz buffers
    Sz1::Float64
    Sz2::Float64
    ######################## observable buffers
    const ψk::Matrix{ComplexF64}
    const Sk::Matrix{Float64}
    const ψT::Vector{ComplexF64}
    const ST::Vector{Float64}
    const ψiψj::Matrix{Bool}
end
function Estimator(ψ0::Matrix{Bool}, Λ::Int64)
    return Estimator(
        [Op(i) for i ∈ 1:Λ],
        null_legs(ψ0),
        null_legs(ψ0),
        deepcopy(ψ0),
        #
        Λ, 0, 1.0, 1.0,
        #
        0.0, 0.0,
        # 
        zeros(ComplexF64, size(ψ0)),
        zeros(Float64, size(ψ0)),
        #
        zeros(ComplexF64, Λ),
        zeros(Float64, Λ),
        #
        zeros(Bool, length(ψ0), length(ψ0))
    )
end

# function calψT!(X::Estimator, i::Int64)
#     ψcur = ψ0[i]
#     l0 = X.legs_last[i]
#     if l0 == null_leg
#         return false
#     end
#     X.ψT .= 1im
#     cur = l0
#     Λ = X.Λ
#     p0 = cur.p
#     Λs = X.Λτ
#     Λ0 = X.Λ
#     while true
#         cur = cur.next
#         if cur == l0
#             return true
#         end

#     end
#     l_cur ≠ 

# end


