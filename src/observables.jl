include("accumulator.jl")

@with_kw_noshow struct Obs
    L::Tuple{Int64,Int64}
    β::Float64
    ψsnapshots::Vector{BitMatrix} = BitMatrix[]
    N::Int64 = prod(L)
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
    ST::Accum{Matrix{Float64}}
    ψiψj::Accum{Matrix{Float64}} = Accum(zeros(Float64, N, N), 0)
end

function Obs(X::Estimator)
    return Obs(
        L = size(X.ψ0), β = X.β,
        ST = Accum(zeros(Float64, X.Λ)
        )
    )
end

# @with_kw_noshow mutable struct Obs
#     para::Union{Tuple,NamedTuple}
#     ψ0::AbstractArray{Bool}
#     ntrace::Vector{Float64} = Float64[]
#     Nwl::Float64 = 0.0
#     Nwl2::Float64 = 0.0
#     Nσx::Float64 = 0.0
#     E::Float64 = 0.0
#     ρ::Float64 = 0.0
#     ρ2::Float64 = 0.0
#     ρflip::Float64 = 0.0
#     ρleaf::Float64 = 0.0
#     Sk::Matrix{Float64} = zeros(para[1])
#     ninj::Matrix{Float64} = zeros(para[1])
# end



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


