# -*- coding: utf-8 -*-

export LinearProgramTemperature
export set_rate

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

struct LinearProgramTemperature
    T₀::Ref{Float64}
    θ::Ref{Float64}

    function LinearProgramTemperature(T₀::Float64, θ::Float64)
        return new(Ref(T₀), Ref(θ/60.0))
    end
end

function (self::LinearProgramTemperature)(t)
    return self.T₀[] + self.θ[] * t
end

function set_rate(self::LinearProgramTemperature, θ)
    self.θ[] = θ / 60.0
end
