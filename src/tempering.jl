function annealing(t::Int64, β::Float64, τ::Int64)::Float64
    βmin = 0.5
    r = log(β/βmin)
    τp = 0.15τ
    τc = τ-2τp
    x = (t-τp) / τc
    x = sqrt(clamp(x, 0, 1))
    rt = r * ( x - 0.2sinpi(15x)^2)
    βt = βmin * exp(rt)
    return βt
end

function linear_annealing(t, β, ξ, τ)
    if t > τ
        return (β, ξ)
    else
        x = clamp(t / τ, 0, 1)
        βt = (1-x) * 0.1β + x * β
        ξt = ξ == 0.0 ? 0.0 : ((1-x) * 0.1 + x * ξ)
        return (βt, ξt)
    end
end