@inline function sortedpair(i::T, j::T) where {T}
    return i < j ? (i, j) : (j, i)
end
pbcid(x::T, L::T) where {T<:Integer} = mod1(x, L)
pbcid(x::NTuple{N,T}, L::NTuple{N,T}) where {N,T<:Integer} = CartesianIndex(mod1.(Tuple(x), L))

@inline function pbcshift(p, s, L)
    return mod1.(Tuple(p) .+ Tuple(s), L) |> CartesianIndex
end

# 2D neighbors search
@inline function udlr(i::Int, r0::Int, c0::Int)::NTuple{4,Int}
    N::Int = r0 * c0
    c::Int, r::Int = divrem(i - 1, r0)
    return (
        r == 0 ? (i + r0 - 1) : (i - 1),
        r == r0 - 1 ? (i - r0 + 1) : (i + 1),
        c == 0 ? (i - r0 + N) : (i - r0),
        c == c0 - 1 ? (i + r0 - N) : (i + r0),
    )
end
@inline udlr(i::Int, L::Tuple{Int,Int})::NTuple{4,Int} = udlr(i, L[1], L[2])
@inline udlr(i::Int, L::Int)::NTuple{4,Int} = udlr(i, L, L)

@inline function udlrx(i::Int, r0::Int, c0::Int)::NTuple{5,Int}
    N::Int = r0 * c0
    c::Int, r::Int = divrem(i - 1, r0)
    return (
        r == 0 ? (i + r0 - 1) : (i - 1),
        r == r0 - 1 ? (i - r0 + 1) : (i + 1),
        c == 0 ? (i - r0 + N) : (i - r0),
        c == c0 - 1 ? (i + r0 - N) : (i + r0),
        i
    )
end
@inline udlrx(i::Int, L::Tuple{Int,Int})::NTuple{5,Int} = udlrx(i, L[1], L[2])
@inline udlrx(i::Int, L::Int)::NTuple{5,Int} = udlrx(i, L, L)

function randsector(L::Tuple{Int64,Int64}, Nseed::Int64)
    N = prod(L)
    lattice = CartesianIndices(L)
    sector = CartesianIndex{2}[]
    nlist = Int64[]
    for i ∈ 1:Nseed
        for (j, q) ∈ enumerate(udlrx(rand(1:N), L))
            if q ∈ nlist
                return randsector(L, Nseed)
            else
                push!(nlist, q)
            end
            if j == 5
                push!(sector, lattice[q])
            end
        end
    end
    return (sector...,)
end
randsector(L0::Int64, Nseed::Int64) = randsector((L0,L0), Nseed)