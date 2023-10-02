# using Plots, Measures, LaTeXStrings
using Revise, KCMQMC, ProgressMeter
using Plots
# Ob1 = onesimu(L, β, randsector(L, 1); tth=10000, ts=20000)
# using Plots
# plot(Ob1.ρtrace)

# Ob1 = onesimu((4, 4), 10.0, randsector((4, 4), 1); tth=10000, ts=20000)

L = (48, 48)
β = 48.0
# sec = randsector(L, 30)
tth = 10000
ts = 30000
Λ = ceil(Int, prod(L) * β * 3)
X = Estimator(zeros(Bool, L), Λ)
# for c ∈ sec
#     X.ψ0[c] = true
# end
X.ψ0 .= false
# X.ψ0[1] = true
for i ∈ 1:2:L[1], j ∈ 1:2:L[2]
    if ~(isodd(i ÷ 2) ⊻ isodd(j ÷ 2))
        X.ψ0[i, j] = true
    end
end
heatmap(X.ψ0)

Ob = Obs(X)
@showprogress for t ∈ 1:tth
    X.β = annealing(t, β, tth)
    sweep!(X)
    only_track_trace!(X, Ob)
    if t % 100 == 0
        snapshot_ψ!(X, Ob)
    end
end
using Plots
heatmap(X.ψ0)
X.β = β
Ob.ρtrace |> plot
using Statistics

@showprogress for t ∈ 1:ts
    sweep!(X)
    onestep_measure!(X, Ob)
    if t % 100 == 0
        snapshot_ψ!(X, Ob)
        track_ψiψj!(X, Ob)
    end
end
ρ̄ = mean(Ob.ρ)
mean(mean.(Ob.ψsnapshots[end÷2:end]))
mean(Ob.E)

ST = mean(Ob.ST)
# ST[1] = 0
using FFTW
Sτ = real.(bfft(ST)) ./ length(ST)^2
τgrid = LinRange(0.0, β, length(Sτ) + 1)
Sτ .-= ρ̄^2
push!(Sτ, Sτ[1])
plot(τgrid, Sτ, ylims=(0, 0.25))
heatmap(Ob.ψsnapshots[end])
Skbar = mean(Ob.Sk)
# Skbar[1] = 0
using FFTW
Cij = real.(bfft(Skbar)) ./ prod(L)^2 .- ρ̄^2
heatmap()
plot(1:(size(Cij,1)-1), abs.(Cij[1,2:end]), yscale = :log10)