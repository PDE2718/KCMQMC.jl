function sweep!(X::Estimator)
    sweep_diag!(X)
    sweep_off!(X)
    return nothing
end

function check_increment!(X::Estimator)::Bool
    Λ0::Int = length(X.H)
    Λtol::Int = ceil(Int, X.n * 1.2)
    Λt::Int = ceil(Int, 1.5*Λ0)
    if Λtol > Λ0
        for p ∈ (Λ0+1):Λt
            push!(X.H, Op(p))
        end
        return true
    else
        return false
    end
end