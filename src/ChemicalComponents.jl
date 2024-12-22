# -*- coding: utf-8 -*-
module ChemicalComponents

using ChemicalExceptions
using ChemicalElements
using PhysicalChemistry

#######################################################################
# INTERNALS:
#######################################################################

unpack(c) = vcat(keys(c.data)...), vcat(values(c.data)...)

const ElementAmount = Pair{Symbol, <:Number}

@enum CompositionTypes begin
    Stoichiometry
    MoleProportion
    MassProportion
end

struct Composition{T}
    data::NamedTuple
    scale::Pair{Symbol,<:Number}

    function Composition{T}(;
            scale::Union{Nothing, ElementAmount} = nothing,
            kw...
        ) where T
        # Empty composition keywords is unacceptable:
        isempty(kw) && throw(EmptyCompositionError())

        # Handle absence of scaling for generality:
        scale = something(scale, first(kw))

        # Enforce type of data (no duplicates!)
        data = NamedTuple(kw)

        # Validate existence of scaling element:
        let
            scaler = scale.first
            elements = keys(data)

            scaler in elements || begin
                throw(InvalidScalerError(scaler, elements))
            end
        end

        # Create and return object:
        return new{T}(data, scale)
    end
end

struct ChemicalComponent
    elements::Vector{Symbol}
    coefficients::Vector{Float64}
    mole_fractions::Vector{Float64}
    mass_fractions::Vector{Float64}
    molar_mass::Float64
end

#######################################################################
# EXPORTED:
#######################################################################

function stoichiometry(; kw...)
    return Composition{Stoichiometry}(; kw...)
end

function mole_proportions(; scale = nothing, kw...)
    scale = something(scale, first(kw).first => 1.0)
    return Composition{MoleProportion}(; scale, kw...)
end

function mass_proportions(; scale = nothing, kw...)
    scale = something(scale, first(kw).first => 1.0)
    return Composition{MassProportion}(; scale, kw...)
end

function component(c::T) where {T <: Composition{Stoichiometry}}
    elements, coefs = unpack(c)
    W = atomic_mass.(elements)

    idx = findfirst(x->x==c.scale.first, elements)
    coefs = (c.scale.second / coefs[idx]) .* coefs

    X = coefs ./ sum(coefs)
    Y = get_mass_fractions(X, W)
    M = coefs' * W

    return ChemicalComponent(elements, coefs, X, Y, M)
end

function component(c::T) where {T <: Composition{MoleProportion}}
end

function component(c::T) where {T <: Composition{MassProportion}}
end

function component(spec; kw...)
    c = if spec == :stoichiometry
        stoichiometry(; kw...)
    elseif spec == :mole_proportions
        mole_proportions(; kw...)
    elseif spec == :mass_proportions
        mass_proportions(; kw...)
    else
        error("Invalid composition specification $(spec)")
    end

    return component(c)
end

#######################################################################
# API
#######################################################################

@doc """\

    stoichiometry(; kw...)

Create composition based on elemental stoichiometry.
""" stoichiometry
export stoichiometry

@doc """\

    mole_proportions(; scale = nothing, kw...)

Create composition based on relative molar proportions.
""" mole_proportions
export mole_proportions

@doc """\

    mole_proportions(; scale = nothing, kw...)

Create composition based on relative molar proportions.
""" mass_proportions
export mass_proportions

@doc """\

    component(spec; kw...)

Compile component from given composition specification.
""" component
export component

end # (module ChemicalComponents)