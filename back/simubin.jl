using Parameters

@with_kw_noshow mutable struct Obs
    para::Union{Tuple,NamedTuple}
    ψ0::AbstractArray{Bool}
    ntrace::Vector{Float64} = Float64[]
    E::Float64 = 0.0
    Nwl::Float64 = 0.0
    Nwl2::Float64 = 0.0
    ρ::Float64 = 0.0
    ρ2::Float64 = 0.0
    ρflip::Float64 = 0.0
    ρleaf::Float64 = 0.0
    Sk::Matrix{Float64} = zeros(para[1])
end

function simubin(para::Union{Tuple,NamedTuple})
    L, β, sector, tth, ts = para
    N = prod(L)
    # prepare ψ0
    ψ0 = zeros(Bool, L)
    for c ∈ sector
        ψ0[c] = true
    end

    ####################[observables]#####################

    # average
    ob = Obs(
        para=para,
        ψ0=ψ0,
    )
    sizehint!(ob.ntrace, tth + ts)

    #################### simulation ###################

    # prepare the string
    Λ::Int = ceil(Int, N * β * 2)
    n::Int = 0
    H::OpString = [Op(i) for i ∈ 1:Λ]
    legs_first = null_legs(ψ0)
    legs_last = null_legs(ψ0)
    for i ∈ 1:10tth
        n = sweep!(H::OpString, legs_first, legs_last, ψ0, Λ, n, 
            annealing(i, β, tth)
        )
        ρ = obs_n̄(ψ0)
        push!(ob.ntrace, ρ)
    end
    for i ∈ 1:ts
        n = sweep!(H, legs_first, legs_last, ψ0, β)

        E = obs_E(n, β, N)
        ρ = obs_n̄(ψ0)
        Sk = obs_Sk(ψ0)
        ninj = obs_ninj(Sk)
        ρflip, ρleaf = obs_leaves(ψ0)
        push!(ob.ntrace, ρ)

        # recording
        ob.Nwl += Nwl
        ob.Nwl2 += Nwl^2
        ob.Nσx += Nσx

        ob.E += E
        ob.ρ += ρ
        ob.ρ2 += ρ^2
        ob.ρflip += ρflip
        ob.ρleaf += ρleaf
        ob.Sk += Sk
        ob.ninj += ninj
    end
    for i ∈ (4:fieldcount(Obs))
        setfield!(ob, i, getfield(ob, i) ./ ts)
    end
    return dumptoNamedTuple(ob)
end
