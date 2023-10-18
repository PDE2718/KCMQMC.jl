# using Plots, Measures, LaTeXStrings
using Revise, KCMQMC, Statistics, FFTW, ProgressMeter
using Plots, Measures, LaTeXStrings
# Ob1 = onesimu(L, β, randsector(L, 1); tth=10000, ts=20000)
# using Plots
# plot(Ob1.ρtrace)

# Ob1 = onesimu((4, 4), 10.0, randsector((4, 4), 1); tth=10000, ts=20000)
L = (16, 16)
β = 24.0
# sec = randsector(L, 30)
tth = 20_000
ts = 20_000
Λ = ceil(Int, prod(L) * β * 10.0)
X = Estimator(zeros(Bool, L), Λ)
X.β = β
X.μ = -3.0
X.ξ = 0.005

Ob = Obs(X)
# μgrid = vcat([fill(μ, 1000) for μ ∈ -3.0:0.1:2.5]...)
# μgrid = vcat(μ1,μ2)
tmax = 100_000
tgrid = 1:tmax
μgrid = LinRange(-3,2.5,tmax)
μgrids = repeat([μgrid,reverse(μgrid)],3)
for μgrid ∈ μgrids
    @showprogress for μ ∈ μgrid
        X.μ = μ
        sweep!(X)
        only_track_trace!(X, Ob)
    end
end
ρgrids = [Ob.ρtrace[(i*tmax).+tgrid] for i ∈ 0:5]

Etrace = -Ob.nwltrace./(β*prod(L)) .+ μshift.(μgrid)
binn(x,Lbin=50) = mean(x[i:Lbin:end] for i ∈ 1:Lbin)
plot(μgrid |> binn, Etrace |> binn)
plot(μgrid |> binn, Ob.nwltrace |> binn)
plot(binn.(μgrids), binn.(ρgrids))

plot([μgrid[1:tmax], μgrid[tmax+1:end]],
    [Ob.ρtrace[1:tmax],Ob.ρtrace[tmax+1:end]],
    # ylims = (0.3,0.35)
)
# heatmap(X.ψ0, aspect_ratio=:equal, size = [300,300], colorbar=false)
plt1 = plot(Ob.ρtrace)
plt2 = plot(μgrid)
plt = plot(plt1,plt2)

mean(Ob.ρtrace)

@showprogress for t ∈ 1:ts
    sweep!(X)
    onestep_measure!(X, Ob)
    # if t % 100 == 0
    #     snapshot_ψ!(X, Ob)
    #     track_ψiψj!(X, Ob)
    # end
end
# Ob.ρtrace |> plot
# mean(mean.(Ob.ψsnapshots[end÷2:end]))
mean(Ob.E)
ρ̄ = mean(Ob.ρ)

X.n
# heatmap(Ob.ψsnapshots[end])
Skbar = mean(Ob.Sk)
Skbar[1] = 0
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
    onestep_measure!(X, Ob)
    if t % 100 == 0
        snapshot_ψ!(X, Ob)
        # track_ψiψj!(X, Ob)
    end
end

t_mean = 10
ij00 = (1,1)
mag = 1
ani = @animate for i ∈ 1:2:1000-t_mean
    heatmap(mean(ψss[i:i+t_mean])[ij00[1]:mag:end, ij00[2]:mag:end],
    color = :RdBu_9,
    aspect_ratio=:equal, size=[350, 400], colorbar=false, margins = 5mm)
    annotate!((1.0, 1.0), text("t = $(i)", 11, :right, :bottom, :green))
end
gif(ani, fps = 20)