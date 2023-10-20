# number of flipable site / flipable up state
function obs_leaves(ψ)::Tuple{Float64,Float64}
    nflip = nleaf = cnt = 0
    q = false
    L = size(ψ)
    for (r, x) ∈ enumerate(ψ)
        cnt = 0
        for rp in udlr(r, L)
            cnt += ψ[rp]
        end
        q = (cnt == 1)
        nflip += q
        nleaf += q && x
    end
    return (nflip, nleaf) ./ length(ψ)
end

function track_trace!(X::Estimator, O::Obs)
    push!(O.ρtrace, X.Sz1 / O.N)
    push!(O.nwltrace, X.n)
end
function onestep_measure!(X::Estimator, O::Obs)
    ψ0, ψk, Sk, β, ξ, μ = X.ψ0, X.ψk, X.Sk, X.β, X.ξ, X.μ
    @assert O.L == size(ψ0)
    N = length(ψ0)
    nwl1::Int64 = X.n ; push!(O.nwl1, nwl1) ; push!(O.nwltrace, nwl1)
    E = - nwl1 / β / N + μshift(μ) ; push!(O.E,E)
    nwl2::Int64 = nwl1^2 ; push!(O.nwl2, nwl2)
    Sz1 = X.Sz1 ; push!(O.Sz1, Sz1) ; push!(O.Sz2, X.Sz2)
    ρ = Sz1 / N ; push!(O.ρ, ρ) ; push!(O.ρtrace, ρ)
    ρflip, ρleaf = obs_leaves(ψ0)
    push!(O.ρleaf, ρleaf) ; push!(O.ρflip, ρflip)

    ψk .= ψ0 ; fft!(ψk) ; map!(abs2, Sk, ψk) ; push!(O.Sk, Sk)
end

function snapshot_ψ!(X::Estimator, O::Obs)
    push!(O.ψsnapshots, X.ψ0 .== true)
end

obs_ninj(Sk::Matrix{Float64})::Matrix{Float64} = real.(bfft(Sk)) ./ (length(Sk)^2)

function obs_cicj(Sk::Matrix{Float64})::Matrix{Float64}
    S̃ = deepcopy(Sk)
    S̃[1] = 0
    return real.(bfft(S̃)) ./ (length(S̃)^2)
end
# Auxiliary functions for visualization
function pbcfill(A)
    B = fftshift(A)
    if iseven(size(B, 1))
        B = vcat(B, B[1:1, 1:end])
    end
    if iseven(size(B, 2))
        B = hcat(B, B[1:end, 1:1])
    end
    return B
end
function pbcgrid(A)
    @assert size(A) |> allequal
    n = size(A, 1)
    @assert isodd(n)
    return -(n ÷ 2):(n÷2)
end


# function track_ψiψj!(X::Estimator, O::Obs)
#     ψt = view(X.ψ0,:)
#     kron!(view(X.ψiψj,:), ψt, ψt)
#     push!(O.ψiψj, X.ψiψj)
# end


# # energy
# obs_E(nwl, β, N) = -nwl / β / N + 1

# # density
# obs_n̄(ψ) = mean(ψ)

# # Structural factor
# obs_Sk(ψ)::Matrix{Float64} = abs2.(fft(ψ))

# # correlation function
# obs_ninj(ψ::AbstractArray{Bool})::Matrix{Float64} = real.(bfft(obs_Sk(ψ))) ./ (length(ψ)^2)