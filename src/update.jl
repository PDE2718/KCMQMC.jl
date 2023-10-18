function update_ψ0!(ψ0::Matrix{Bool}, legs_last::Matrix{Leg})
    for (i, l) ∈ enumerate(legs_last)
        if l != null_leg
            ψ0[i] = l.ψ
        end
    end
end

function clear_string!(H::OpString)
    for h ∈ H
        h.flag = 0
    end
    return nothing
end

###### diag update
# the shift const : μ>0 => 1 ; μ<0 => 1-μ
μshift(μ::Float64)::Float64 = μ < 0.0 ? (1.0 - μ) : 1.0
diag_weight(ψ::Bool, μ::Float64)::Float64 = μshift(μ) + μ * ψ
offdiag_weight(cnt::Int, ξ::Float64)::Float64 = cnt == 1 ? (1.0-ξ) : ξ

function sweep_diag!(H::OpString, Λ::Int64, n::Int64, β::Float64,
    ξ::Float64, μ::Float64,
    legs_first::Matrix{Leg}, legs_last::Matrix{Leg}, ψ0::Matrix{Bool}, ψT)

    L::Tuple{Int,Int} = size(ψ0)
    N::Int = length(ψ0)
    # iwl::Int = rand(eachindex(ψ0))
    
    # ic::Int = 0 # index of center
    w::Float64 = 0.
    l::Leg = null_leg

    Sz::Int = sum(ψ0)
    Sz1::Int = 0
    Sz2::Int = 0

    for i ∈ eachindex(legs_first, legs_last)
        legs_first[i] = legs_last[i] = null_leg
    end
    for (p,h) ∈ enumerate(H)
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
                l = h[r]
                l.j = j
                l.ψ = ψ0[j]
                if legs_first[j] === null_leg
                    legs_first[j] = l
                    legs_last[j] = l
                    l.prev = l
                    l.next = l
                else
                    l.prev = legs_last[j]
                    l.next = legs_first[j]
                    legs_last[j].next = l
                    legs_first[j].prev = l
                    legs_last[j] = l
                end
            end
        end
        # ψT[p] = ψ0[iwl] # record only one boson world-line
        Sz1 += Sz
        Sz2 += Sz^2
    end
    return n, Sz1/Λ, Sz2/Λ
end

####### off diag update

function wf2wi_tailhead(l::Leg, ξ::Float64, μ::Float64, istail::Bool)::Float64
    ψ0::Bool = l.ψ ⊻ (istail && l.flag < 0)
    wμj::Float64 = diag_weight(ψ0, μ)
    wσx::Float64 = offdiag_weight(count(l), ξ)
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
        # println(confsign(l.op, l.r))
        # println("wf=$(wf), wi=$(wi), flag=$(l.flag)")
        return wf / wi
    end
end

function wf2wi_cyclic(l::Leg, μ::Float64)
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
    if h0.flag == 0
        return false
    end
    tail::Leg = head::Leg = h0[5]
    wr::Float64 = 1.0
    while true
        head = head.next
        if is_center(head)
            if head == tail
                wr *= wf2wi_cyclic(head, μ)
            else
                wr *= wf2wi_tailhead(tail, ξ, μ, true)
                wr *= wf2wi_tailhead(head, ξ, μ, false)
            end
            metro(wr) ? break : return false
        else
            wr *= wf2wi_wormbody(head, ξ)
        end
    end

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
end

function sweep_off!(H::OpString, ξ::Float64, μ::Float64)
    for h ∈ H
        update_ahead!(h, ξ, μ)
        update_ahead!(rand(H), ξ, μ)
    end
    return nothing
end
