function onesimu(L::Tuple{Int,Int},
    β::Float64, ξ::Float64, μ::Float64,
    tth::Int, tbin::Int, nbin::Int, sector
    )
    # 
    Λ0 = ceil(Int, prod(L) * β)
    X = Estimator(zeros(Bool, L), Λ0)
    for s ∈ sector
        X.ψ0[s] = true
    end
    X.β = β
    X.μ = μ
    X.ξ = ξ

    Ob_t = Obs(X)
    M = [Obs(X) for _ ∈ 1:nbin]

    println("thermalizing...")
    for t ∈ 1:2tth

        X.β, X.ξ = linear_annealing(t, β, ξ, tth)
        sweep!(X)
        check_increment!(X)
        track_trace!(X, Ob_t)
    end
    println("thermalized!")
    for (i, obcurr) ∈ enumerate(M)
        println("simulating bin $(i)/$(nbin)")
        for t ∈ 1:tbin
            sweep!(X)
            onestep_measure!(X, obcurr)
        end
        snapshot_ψ!(X, Ob_t)
    end
    println("finished!")
    return Ob_t, M
end