mutable struct Estimator
    const H::OpString
    const legs_first::Matrix{Leg}
    const legs_last::Matrix{Leg}
    const ψ0::Matrix{Bool}
    const ψt::Matrix{Bool}
    ######################## numbers
    n::Int64
    β::Float64
    ξ::Float64
    μ::Float64
    ######################## Sz buffers
    Sz1::Float64
    Sz2::Float64
    ######################## observable buffers
    const ψk::Matrix{ComplexF64}
    const Sk::Matrix{Float64}
end
function Estimator(ψ0::Matrix{Bool}, Λ0::Int64)
    return Estimator(
        [Op(i) for i ∈ 1:Λ0],
        null_legs(ψ0),
        null_legs(ψ0),
        deepcopy(ψ0),
        deepcopy(ψ0),
        #
        0, 0.0, 0.0, 0.0,
        #
        0.0, 0.0,
        # 
        zeros(ComplexF64, size(ψ0)),
        zeros(Float64, size(ψ0)),
    )
end

