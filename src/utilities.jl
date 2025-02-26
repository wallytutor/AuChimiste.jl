# -*- coding: utf-8 -*-

function heaviside(T, T₀)
    return 1//2 * (sign(T - T₀) + 1)
end

function heaviside(T, T₀, k)
    return 1//2 * (tanh(k * (T - T₀)) + 1)
end

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
