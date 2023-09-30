# base abstract types
abstract type AbstractOp end
abstract type AbstractLeg end

# Typedef of Op and Leg
mutable struct Op{T<:AbstractLeg} <: AbstractOp
    const leg1::T
    const leg2::T
    const leg3::T
    const leg4::T
    const leg5::T
    const p::Int
    flag::Int
end

mutable struct Leg <: AbstractLeg
    const r::Int # const?
    ψ::Bool
    j::Int
    prev::Leg
    next::Leg
    op::Op{Leg}
    function Leg(r::Integer)
        l = new(r, false, 0)
        l.prev = l
        l.next = l
        return l
    end
end

## some methods for the Op{Leg} and Leg
import Base: getproperty, setproperty!, getindex, count, firstindex, lastindex, eachindex
getproperty(o::Op{Leg}, y::Symbol) = getfield(o::Op{Leg}, y::Symbol)
setproperty!(o::Op{Leg}, y::Symbol, z) = setfield!(o::Op{Leg}, y::Symbol, z)
getindex(o::Op{Leg}, i::Int) = getfield(o, i)

import Base.iterate
iterate(o::Op{Leg}) = iterate(o,1)
iterate(o::Op{Leg}, i::Int) = i == 6 ? nothing : (getfield(o,i), i+1)

firstindex(x::Op{Leg}) = 1
lastindex(x::Op{Leg}) = 5
eachindex(x::Op{Leg}) = 1:5

count(x::Op{Leg})::Int = x[1].ψ + x[2].ψ + x[3].ψ + x[4].ψ
count(l::Leg)::Int = count(x.op)
getproperty(x::Leg, y::Symbol) = hasfield(Leg, y) ? getfield(x, y) : getfield(getfield(x, :op), y)
setproperty!(x::Leg, y::Symbol, z) = hasfield(Leg, y) ? setfield!(x::Leg, y::Symbol, z) : setfield!(getfield(x, :op), y, z)

is_center(l::Leg)::Bool = l.r == 5
is_σx(op::Op{Leg})::Bool = op.flag < 0
is_σx(l::Leg)::Bool = l.flag < 0

flipable(h::Op{Leg})::Bool = is_σx(h) || count(h) == 1
unflipable(h::Op{Leg})::Bool = ~flipable(h)

flipable(l::Leg)::Bool   = is_center(l) ?  flipable(l.op) : ~is_σx(l)
unflipable(l::Leg)::Bool = is_center(l) ? ~flipable(l.op) :  is_σx(l)

# function unflipable(l::Leg)::Bool
#     if is_center(l)
#         return 
#     else
#         return is_σx(l)
#     end
# end

function Op(p::Integer)::Op{Leg}
    h = Op(Leg(1), Leg(2), Leg(3), Leg(4), Leg(5), p, 0)
    for l ∈ h
        l.op = h
    end
    return h
end
global const null_leg::Leg = Leg(0)
global const null_op::Op = Op(fill(null_leg,5)..., 0, 0)
global const OpString::Type = Vector{Op{Leg}}
null_legs(ψ0::AbstractArray{Bool}) = [null_leg for _ ∈ ψ0]

## How to display the OpString

function Base.show(io::IO, l::Leg)
    if !isdefined(l, :op)
        print(io, " leg NULL")
    else
        oprev::Char = l.prev.ψ ? '■' : '□'
        onext::Char = l.ψ ? '■' : '□'
        print(io, " leg$(l.r): ⋯$(l.prev.p)→$(oprev)|$(l.p)|$(onext)→$(l.next.p)⋯ ")
    end
end
function Base.show(io::IO, O::Op{Leg})
    if O.flag ≠ 0
        println(io, "$(O.flag < 0 ? "σx" : "I")($(abs(O.flag))) at p = $(O.p)")
        for l ∈ O
            println(io, l)
        end
    else
        println(io, "I0 at p = $(O.p)")
    end
end
