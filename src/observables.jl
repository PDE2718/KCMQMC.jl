include("accumulator.jl")

Base.@kwdef struct Obs
    t::Int64
    L::Tuple{Int64,Int64}
    β::Float64
    ξ::Float64
    μ::Float64
    # N::Int64 = prod(L)
    # nwltrace::Vector{Int64} = Int64[]
    E::Accum{Float64} = Accum(0.0)
    ρ::Accum{Float64} = Accum(0.0)
    nwl1::Accum{Int64} = Accum(0)
    nwl2::Accum{Int64} = Accum(0)
    Sz1::Accum{Float64} = Accum(0.)
    Sz2::Accum{Float64} = Accum(0.)
    # ρleaf::Accum{Float64} = Accum(0.)
    # ρflip::Accum{Float64} = Accum(0.)
    Sk::Accum{Matrix{Float64}} = Accum(zeros(Float64, L), 0)
    ψ̄::Accum{Matrix{Float64}} = Accum(zeros(Float64, L), 0)
    ρtrace::Vector{Float64} = Float64[]
    ψsnapshots::Vector{BitMatrix} = BitMatrix[]
end

Obs(X::Estimator, t::Int64=0) = Obs(t,
    L=size(X.ψ0), β=X.β, ξ=X.ξ, μ=X.μ
)