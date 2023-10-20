function onesimu(L::Tuple{Int64,Int64}, β::Float64, sector; tth = 10000, ts = 20000)
    Λ = ceil(Int, prod(L) * β * 3)
    X = Estimator(zeros(Bool,L), Λ)
    for c ∈ sector
        X.ψ0[c] = true
    end
    Ob = Obs(X)
    for t ∈ 1:tth
        X.β = annealing(t, β, tth)
        sweep!(X)
        only_track_trace!(X,Ob)
        if t%100 == 0
            snapshot_ψ!(X,Ob)
        end
    end
    X.β = β
    for t ∈ 1:ts
        sweep!(X)
        onestep_measure!(X,Ob)
        if t%100 == 0
            snapshot_ψ!(X,Ob)
            track_ψiψj!(X,Ob)
        end
    end
    return Ob
end
