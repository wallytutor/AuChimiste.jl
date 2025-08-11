# -*- coding: utf-8 -*-
export fermilike_viscosity
export FermiLikeViscosity

"""
    fermilike_viscosity(T, μ∞, δ, Θ, Δ)

Evaluation of Fermi-distribution-alike viscosity function.
"""
function fermilike_viscosity(T, μ∞, δ, Θ, Δ)
    return μ∞ + δ / (1 + exp((T - Θ) / Δ))
end

"""
Temperature-dependent viscosity with a Fermi-like distribution dependency.

## Fields

$(TYPEDFIELDS)
"""
struct FermiLikeViscosity <: AbstractViscosity
    "Center temperature of melting range [K]."
    Θ::Float64

    "Spread factor over melting range  [K]."
    Δ::Float64

    "Viscosity change during melting [Pa.s]."
    δ::Float64

    "High temperature viscosity [Pa.s]."
    μ∞::Float64

    "Low temperature viscosity [Pa.s]."
    μ₀::Float64

    "Spread coefficient used to compute Δ."
    κ::Float64

    "Melting start temperature [K]."
    Ts::Float64
    
    "Melting end temperature [K]."
    Te::Float64

    function FermiLikeViscosity(Ts, Te, μ₀, μ∞, κ)
        Θ = (Te + Ts) / 2
        Δ = (Te - Ts) / κ
        δ = μ₀ - μ∞
        return new(Θ, Δ, δ, μ∞, μ₀, κ, Ts, Te)
    end
end

function (obj::FermiLikeViscosity)(T)
    return fermilike_viscosity(T, obj.μ∞, obj.δ, obj.Θ, obj.Δ)
end
