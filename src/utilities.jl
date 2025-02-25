# -*- coding: utf-8 -*-

function cumtrapz(X::T, Y::T) where {T <: AbstractVector}
    @assert length(X) == length(Y)
    
    out = similar(X)
    out[1] = 0
    
    for i in 2:length(X)
        h = X[i] - X[i-1]
        B = Y[i] + Y[i-1]
        out[i] = out[i-1] + h * B / 2
    end
    
    return out
end
