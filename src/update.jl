function update_ψ0!(ψ0::Matrix{Bool}, legs_last::Matrix{MaybeLeg})
    for i ∈ eachindex(ψ0, legs_last)
        if legs_last[i] isa Leg
            ψ0[i] = legs_last[i].ψ
        end
    end
end

function clear_estimator!(X::Estimator, ψ00::Matrix{Bool})
    X.ψ0 .= ψ00
    X.n = 0
    for h ∈ X.H
        h.flag = 0
    end
    for i ∈ eachindex(X.legs_first, X.legs_last)
        X.legs_first[i] = X.legs_last[i] = nothing
    end
    return nothing
end

###### diag update
# the shift const : μ>0 => 1 ; μ<0 => 1-μ
μshift(μ::Float64)::Float64 = μ < 0.0 ? (1.0 - μ) : 1.0
diag_weight(ψ::Bool, μ::Float64)::Float64 = μshift(μ) + μ * ψ
offdiag_weight(cnt::Int, ξ::Float64)::Float64 = cnt == 1 ? (1.0-ξ) : ξ

function sweep_diag!(X::Estimator)
    X.n, X.Sz1, X.Sz2 = sweep_diag!(X.H,
        X.n, X.β, X.ξ, X.μ,
        X.legs_first, X.legs_last,
        X.ψ0,
    )
end
function sweep_diag!(H::OpString,
    n::Int64, β::Float64, ξ::Float64, μ::Float64,
    legs_first::Matrix{MaybeLeg}, legs_last::Matrix{MaybeLeg},
    ψ0::Matrix{Bool}
    )
    L::Tuple{Int,Int} = size(ψ0)
    N::Int = length(ψ0)
    Λ::Int = length(H)
    
    w::Float64 = 0.
    Sz::Int = sum(ψ0)
    Sz1::Int = 0
    Sz2::Int = 0

    fill!(legs_first, nothing)
    fill!(legs_last, nothing)

    @inbounds for h ∈ H
        # insert diagonal operator
        if h.flag == 0
            ic = rand(1:N)
            w = diag_weight(ψ0[ic], μ)
            if metro((N * w * β) / (Λ - n))
                h.flag = ic
                n += 1
            end
        # remove diagonal operator
        elseif h.flag > 0
            w = diag_weight(ψ0[h.flag], μ)
            if metro((Λ - n + 1) / (N * w * β))
                h.flag = 0
                n -= 1
            end
        elseif h.flag < 0
            ic = -h.flag
            ψ0[ic] ⊻= true
            # also record the total particle number
            Sz += ψ0[ic] ? 1 : -1
        end
        
        # otherwise update config
        # 更新传播过程中的构型(linklist)
        if h.flag ≠ 0
            for (r, j) ∈ enumerate(udlrx(abs(h.flag), L))
                l::Leg = h[r]
                l.j = j
                l.ψ = ψ0[j]
                if legs_first[j] |> isnothing
                    legs_last[j] = l
                    legs_first[j] = l
                    # l.prev = l
                    # l.next = l
                else
                    l.prev = legs_last[j]
                    # l.next = legs_first[j]
                    legs_last[j].next = l
                    # legs_first[j].prev = l
                    legs_last[j] = l
                end
            end
        end
        # ψT[p] = ψ0[iwl] # record only one boson world-line
        Sz1 += Sz
        Sz2 += Sz^2
    end
    for (leg_first,leg_last) ∈ zip(legs_first,legs_last)
        if leg_first isa Leg && leg_last isa Leg
            leg_first.prev = leg_last
            leg_last.next = leg_first
        end
    end
    return n, Sz1/Λ, Sz2/Λ
end

####### off diag update

function wf2wi_tailhead(l::Leg, ξ::Float64, μ::Float64, istail::Bool)::Float64
    ψ0::Bool = l.ψ ⊻ (istail && l.flag < 0)
    wσx::Float64 = offdiag_weight(count(l), ξ)
    wμj::Float64 = diag_weight(ψ0, μ)
    return l.flag < 0 ? (wμj / wσx) : (wσx / wμj)
end

function wf2wi_wormbody(l::Leg, ξ::Float64)::Float64
    if l.flag > 0
        return 1.0
    else
        ci::Int = count(l)
        cf::Int = ci + (l.ψ ? -1 : +1)
        wi::Float64 = offdiag_weight(ci, ξ)
        wf::Float64 = offdiag_weight(cf, ξ)
        return wf/wi
    end
end

function wf2wi_cyclic(l::Leg, μ::Float64)::Float64
    if l.flag < 0
        return 1.0
    else
        ψ = l.ψ
        wi = diag_weight(ψ, μ)
        wf = diag_weight(~ψ, μ)
        return wf/wi
    end
end

function update_ahead!(h0::Op{Leg}, ξ::Float64, μ::Float64)::Bool
    # ensure the entrance is a center
    # I0 cannot flip
    if h0.flag == 0 return false end
    tail::Leg = head::Leg = h0[5]
    wr::Float64 = 1.0
    
    # Find the worm body
    while true
        head = head.next
        if is_center(head) # end this segment anyway
            break
        else
            wr *= wf2wi_wormbody(head, ξ)
            if iszero(wr) return false end
        end
    end

    if head == tail
        if iszero(ξ) && count(head) ≠ 1
            return false
        else
            wr *= wf2wi_cyclic(head, μ)
        end
    else
        wr *= wf2wi_tailhead(tail, ξ, μ, true)
        wr *= wf2wi_tailhead(head, ξ, μ, false)
    end

    # now flip the segment
    if metro(wr)
        ## finally update the configuration
        head.flag *= -1
        tail.flag *= -1
        while true
            tail.ψ ⊻= true
            tail = tail.next
            if tail == head
                return true
            end
        end
    else
        return false
    end
    return false
end

function sweep_off!(H::OpString, ξ::Float64, μ::Float64)
    @inbounds for h ∈ H
        update_ahead!(h, ξ, μ)
    end
end

# function sweep_off!(H::OpString, ξ::Float64, μ::Float64)
#     wl = eachindex(H)
#     l = zip(wl , reverse(wl))
#     @inbounds for (i,j) ∈ l
#         update_ahead!(
#             H[rand((i,j,rand(wl)))],
#             ξ, μ)
#     end
#     return nothing
# end

function sweep_off!(X::Estimator)
    sweep_off!(X.H, X.ξ, X.μ)
    update_ψ0!(X.ψ0, X.legs_last)
    return nothing
end