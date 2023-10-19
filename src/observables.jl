include("accumulator.jl")

@with_kw_noshow struct Obs
    L::Tuple{Int64,Int64}
    β::Float64
    ξ::Float64
    μ::Float64
    N::Int64 = prod(L)
    ψsnapshots::Vector{BitMatrix} = BitMatrix[]
    nwltrace::Vector{Int64} = Int64[]
    ρtrace::Vector{Float64} = Float64[]
    E::Accum{Float64} = Accum(0.0)
    ρ::Accum{Float64} = Accum(0.0)
    nwl1::Accum{Int64} = Accum(0)
    nwl2::Accum{Int64} = Accum(0)
    Sz1::Accum{Float64} = Accum(0.)
    Sz2::Accum{Float64} = Accum(0.)
    ρleaf::Accum{Float64} = Accum(0.)
    ρflip::Accum{Float64} = Accum(0.)
    Sk::Accum{Matrix{Float64}} = Accum(zeros(Float64, L), 0)
end

Obs(X::Estimator) = Obs(
    L=size(X.ψ0), β=X.β, ξ=X.ξ, μ=X.μ
)

###################################test
# using LinearAlgebra

# ψt = (rand(Bool, 48, 48))
# ψiψj = rand(Bool,length(ψt), length(ψt))
# using BenchmarkTools

# function calψij(ψiψj, ψt)
#     for i ∈ eachindex(ψt)
#         for j ∈ eachindex(ψt)
#             ψiψj[i,j] = ψt[i]*ψt[j]
#         end
#     end
# end
# @btime calψij(ψiψj,ψt)

# @btime kron!(ψiψj, ψt, ψt)


