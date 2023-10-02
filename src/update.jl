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

function sweep_diag!(H::OpString, Λ::Int64, n::Int64, β::Float64, ξ::Float64,
    legs_first::Matrix{Leg}, legs_last::Matrix{Leg}, ψ0::Matrix{Bool}, ψT)
    L::Tuple{Int,Int} = size(ψ0)
    N::Int = length(ψ0)
    iwl = rand(eachindex(ψ0))

    wtot = 1+ξ
    P1::Float64 = 1/wtot
    Sz::Int64 = sum(ψ0)
    Sz1::Int64 = 0
    Sz2::Int64 = 0

    l::Leg = null_leg

    for i ∈ eachindex(legs_first, legs_last)
        legs_first[i] = legs_last[i] = null_leg
    end
    for (p,h) ∈ enumerate(H)
        if h.flag == 0 && metro((N * wtot * β) / (Λ - n))
            h.flag = rand(1:Nd)
            n += 1
            h.cons = metro(P1)
        elseif h.flag > 0 && metro((Λ - n + 1) / (N * wtot * β))
            h.flag = 0
            n -= 1
        end

        if h.flag ≠ 0
            i_site = abs(h.flag)
            up_range = h.cons ? (1:5) : (5:5)
            nbs = udlrx(i_site, L)
            if h.flag < 0
                ψ0[i_site] ⊻= true
                Sz += ψ0[i_site] ? 1 : -1
            end
            # 更新传播过程中的构型
            for r ∈ up_range
                l = h[r]
                j = nbs[r]

                l.ψ = ψ0[j]
                l.j = j
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
        ψT[p] = ψ0[iwl] # record only one boson world-line
        Sz1 += Sz
        Sz2 += Sz2
    end
    return n, Sz1/Λ, Sz2/Λ
end

####### off diag update

function update_ahead!(h0::Op{Leg})::Bool
    # ensure the entrance is a center
    # I0 cannot flip
    if h0.flag ≠ 0 && flipable(h0) # the entrance is flipable

        tail::Leg = head::Leg = h0[5]
        while true
            head = head.next
            if flipable(head)
                is_center(head) ? break : continue
            else
                return false
            end
        end

        # The update procedure
        if metro(0.99) # return Bool
            head.flag *= -1
            tail.flag *= -1
            while head !== tail
                head = head.prev
                head.ψ ⊻= true
            end
            return true
        else
            return false
        end

    else
        return false
    end
end

function sweep_off!(H::OpString)
    for h ∈ H
        update_ahead!(h)
        update_ahead!(rand(H))
    end
    return nothing
end
