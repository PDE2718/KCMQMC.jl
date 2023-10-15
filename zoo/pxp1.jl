# using Plots, Measures, LaTeXStrings
using Revise, KCMQMC, Statistics, FFTW, ProgressMeter
using Plots, Measures, LaTeXStrings
# Ob1 = onesimu(L, β, randsector(L, 1); tth=10000, ts=20000)
# using Plots
# plot(Ob1.ρtrace)

# Ob1 = onesimu((4, 4), 10.0, randsector((4, 4), 1); tth=10000, ts=20000)

L = (16, 16)
β = 16.0
ξmax = 0.2
ξmin = 0.2/β
# sec = randsector(L, 30)
tth = 5000
ts = 40000
Λ = ceil(Int, prod(L) * β * 2)
X = Estimator(zeros(Bool, L), Λ)
X.ξ = ξmin
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
    X.β = annealing(t, β, tth)
    # X.β = (t/tth) * β
    X.ξ = (1-t/tth) * ξmax + ξmin
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

###
ST = mean(Ob.ST)
# ST[1] = 0
Sτ = real.(bfft(ST)) ./ length(ST)^2
Sτ .-= ρ̄^2
τgrid = LinRange(0.0, β, length(Sτ) + 1)
push!(Sτ, Sτ[1])
# plot(τgrid, Sτ, ylims=(0, 0.25))
plot(τgrid, Sτ)

heatmap(Ob.ψsnapshots[end])
Skbar = mean(Ob.Sk)
# Skbar[1] = 0
heatmap(pbcfill(Skbar), aspect_ratio=:equal)
Cij = real.(bfft(Skbar)) ./ prod(L)^2 .- ρ̄^2
heatmap(pbcfill(Cij))
plot(1:(size(Cij,1)-1), abs.(Cij[1,2:end]), yscale = :log10, xscale = :log10)

# ani = @animate for ψ ∈ Ob.ψsnapshots
#     heatmap(ψ, aspect_ratio=:equal, size=[300, 300], colorbar=false)
# end
# ψ_ani = gif(ani, fps = 24)
# Ob.ψsnapshots
Ob.ψsnapshots
tslice = 106
ψ̄t = mean(Ob.ψsnapshots[tslice:tslice+200])
heatmap(ψ̄t, aspect_ratio=:equal, size=[500, 500], colorbar=true)


# X.ξ = 0.0001
ψss = typeof(X.ψ0)[]
for t ∈ 1:1000
    sweep!(X)
end
@showprogress for t ∈ 1:1000
    sweep!(X)
    push!(ψss, deepcopy(X.ψ0))
    # onestep_measure!(X, Ob)
    # if t % 100 == 0
    #     snapshot_ψ!(X, Ob)
    #     track_ψiψj!(X, Ob)
    # end
end

t_mean = 40
ij00 = (1,1)
mag = 1
ani = @animate for i ∈ 1:2:1000-t_mean
    heatmap(mean(ψss[i:i+t_mean])[ij00[1]:mag:end, ij00[2]:mag:end],
    color = :RdBu_9,
    aspect_ratio=:equal, size=[350, 400], colorbar=false, margins = 5mm)
    annotate!((1.0, 1.0), text("t = $(i)", 11, :right, :bottom, :green))
end
gif(ani, fps = 20)