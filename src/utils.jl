@inline metro(p) = rand() < p

function dumptoNamedTuple(x)
    keys = fieldnames(typeof(x))
    vals = [getfield(x, k) for k ∈ keys]
    return (; zip(keys, vals)...)
end