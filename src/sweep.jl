function sweep!(X::Estimator)
    H = X.H
    ψ0 = X.ψ0
    β, ξ, μ = X.β, X.ξ, X.μ
    X.n, X.Sz1, X.Sz2 = sweep_diag!(H,
        X.Λ, X.n, β, ξ, μ,
        X.legs_first, X.legs_last, ψ0, X.ψT
    )
    sweep_off!(H, ξ, μ)
    update_ψ0!(ψ0, X.legs_last)
    return nothing
end

# function sweep!(H::OpString, legs_first, legs_last, ψ0, Λ, n, β)::Int
#     n = sweep_diag!(H, legs_first, legs_last, ψ0, Λ, n, β)
#     for h ∈ H
#         update_ahead!(h)
#         update_ahead!(rand(H))
#     end
#     update_ψ0!(ψ0, legs_last)
#     return n
# end

# function sweep!(X::Estimator, n_times::Int)
#     for nt ∈ 1:n_times
#         sweep!(X)
#     end
#     return nothing
# end