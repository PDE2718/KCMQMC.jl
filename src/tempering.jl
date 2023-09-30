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