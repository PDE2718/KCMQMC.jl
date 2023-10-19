# using Plots, Measures, LaTeXStrings
using Revise, KCMQMC, Statistics, FFTW, ProgressMeter
using Plots, Measures, LaTeXStrings

L = (48,48)
β = 48.0
μ = -1.2
ξ = 0.005
tth = 2000
tbin = 5_000
nbin = 1
Λ = ceil(Int, prod(L) * β * 1)
X = Estimator(zeros(Bool, L), Λ)
# for i ∈ eachindex(X.ψ0)
#     X.ψ0[i] = metro(0.2)
# end
X.β = β ; X.μ = μ ; X.ξ = ξ
# X.ψ0 .= rand(Bool, size(X.ψ0))

Ob_t = Obs(X)
Ob_m = [Obs(X) for i ∈ 1:nbin]
@showprogress for t ∈ 1:2tth
    sweep!(X)
    check_increment!(X)
    track_trace!(X, Ob_t)
end
X.n / length(X.H)
# heatmap(X.ψ0)
plot(Ob_t.ρtrace,
    # ylims = (0.25,0.35)
)
ylims!(0.3,0.4)
for n ∈ 1:nbin
    Obcurr = Ob_m[n]
    @showprogress for t ∈ 1:tbin
        sweep!(X)
        onestep_measure!(X, Obcurr)
    end
end
mean(Ob_m[1].ρ)
Ob_m[1].ρtrace
plot(Ob_m[1].ρtrace)

Sk = [mean(o.Sk) for o ∈ Ob_m]
Skbar = mean(Sk)
Cij = obs_cicj(Skbar)

heatmap(pbcfill(Cij))
Cijxy = Cij[1,:] + Cij[:,1]
using LinearAlgebra
CijD = Cij[diagind(Cij)]

rg = 1:L[1]-1
ydata = Cijxy[rg.+1]
plot(rg, abs.(ydata),
    marker = (3),
    color = (x->x>0 ? 1 : 2).(ydata),
    xlims = (1,24), ylims = (1e-4,0.05),
    xscale = :log,
    yscale = :log,
)
plot!(rg, x->0.05 * x^-1.1)
Ob_m[10].ψsnapshots
# heatmap(Sk[1])
# Skbar = mean(Sk)
# Skbar[1] = 0

heatmap(pbcfill(Skbar), clims = (0,2000))

@showprogress for t ∈ 1:1000
    sweep!(X)
    snapshot_ψ!(X, Ob_t)
end

t_mean = 4
ani = @animate for i ∈ 1:2:1000-t_mean
    heatmap(mean(Ob_t.ψsnapshots[i:i+t_mean]),
        # color=:RdBu_9,
        aspect_ratio=:equal, size=[350, 400], colorbar=false, margins=5mm)
    annotate!((1.0, 1.0), text("t = $(i)", 11, :right, :bottom, :green))
end
gif(ani, fps=20)

XXX = [1.]
function modXXX()
    XXX = [5]
end
modXXX()
XXX