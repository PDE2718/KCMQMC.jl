# using Plots, Measures, LaTeXStrings
using Revise, KCMQMC, Statistics, FFTW, ProgressMeter
using Plots, Measures, LaTeXStrings
# Ob1 = onesimu(L, β, randsector(L, 1); tth=10000, ts=20000)
# using Plots
# plot(Ob1.ρtrace)

# Ob1 = onesimu((4, 4), 10.0, randsector((4, 4), 1); tth=10000, ts=20000)

L = (32, 32)
β = 20.0
ξ = 2/β
# sec = randsector(L, 30)
tth = 5000
ts = 20000
Λ = ceil(Int, prod(L) * β * 3)
X = Estimator(zeros(Bool, L), Λ)
X.ξ = ξ
# for c ∈ sec
#     X.ψ0[c] = true
# end
# X.ψ0 .= rand(L) .< 0.3
X.ψ0[1] = true
# for i ∈ 1:2:L[1], j ∈ 1:2:L[2]
#     if ~(isodd(i ÷ 2) ⊻ isodd(j ÷ 2))
#         X.ψ0[i, j] = true
#     end
# end
heatmap(X.ψ0)

Ob = Obs(X)
X.β = β
@showprogress for t ∈ 1:tth
    # X.β = annealing(t, β, tth)
    X.β = (t/tth) * β
    X.ξ = 1/X.β
    sweep!(X)
    only_track_trace!(X, Ob)
    if t % 100 == 0
        snapshot_ψ!(X, Ob)
    end
end
heatmap(X.ψ0, aspect_ratio=:equal, size = [300,300], colorbar=false)
X.β = β
X.ξ = ξ
Ob.ρtrace |> plot
X.H

@showprogress for t ∈ 1:ts
    sweep!(X)
    onestep_measure!(X, Ob)
    if t % 100 == 0
        snapshot_ψ!(X, Ob)
        track_ψiψj!(X, Ob)
    end
end
ρ̄ = mean(Ob.ρ)
# mean(mean.(Ob.ψsnapshots[end÷2:end]))
mean(Ob.E)

ST = mean(Ob.ST)
# ST[1] = 0
Sτ = real.(bfft(ST)) ./ length(ST)^2
# Sτ .-= ρ̄^2
τgrid = LinRange(0.0, β, length(Sτ) + 1)
push!(Sτ, Sτ[1])
plot(τgrid, Sτ, ylims=(0, 0.25))

heatmap(Ob.ψsnapshots[end])
Skbar = mean(Ob.Sk)
Skbar[1] = 0
heatmap(pbcfill(Skbar), aspect_ratio=:equal)
Cij = real.(bfft(Skbar)) ./ prod(L)^2 .- ρ̄^2
heatmap(pbcfill(Cij))
plot(1:(size(Cij,1)-1), abs.(Cij[1,2:end]), yscale = :log10)

# ani = @animate for ψ ∈ Ob.ψsnapshots
#     heatmap(ψ, aspect_ratio=:equal, size=[300, 300], colorbar=false)
# end
# ψ_ani = gif(ani, fps = 24)
# Ob.ψsnapshots
ψ̄t = mean(Ob.ψsnapshots[end-50:end])
heatmap(ψ̄t, aspect_ratio=:equal, size=[300, 300], colorbar=true)


