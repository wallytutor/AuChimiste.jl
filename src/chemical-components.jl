# -*- coding: utf-8 -*-

export component
export stoichiometry
export mole_proportions
export mass_proportions
export stoichiometry_map
export mole_fractions_map
export mass_fractions_map

"""
Represents a chemical component.

Fields
======
$(TYPEDFIELDS)

Notes
=====

- This structure is not intended to be called as a constructor,
  safe use of its features require using [`component`](@ref)
  construction in combination with a composition specification.

- The array of elements is unsorted when construction is performed
  through [`component`](@ref) but may get rearranged when composing
  new chemical components through supported algebra.

- Care must be taken when using `molar_mass` because it is given
  for the associated coefficients. That is always the expected
  behavior for molecular components but might not be the case in
  other applications (solids, solutions) when the mean molecular
  mass may be required.
"""
struct ChemicalComponent
    "Array of component symbols."
    elements::Vector{Symbol}

    "Array of stoichiometric coefficients."
    coefficients::Vector{Float64}

    "Array of elemental mole fractions."
    mole_fractions::Vector{Float64}

    "Array of elemental mass fractions."
    mass_fractions::Vector{Float64}

    "Molar mass of corresponding stoichiometry."
    molar_mass::Float64
end

"""
    component(spec; kw...)

Compile component from given composition specification. This function
is a wrapper eliminating the need of calling [`stoichiometry`](@ref),
[`mole_proportions`](@ref) or [`mass_proportions`](@ref) directly. The
value of `spec` must be the symbol representing one of their names.
"""
function component(spec::Symbol; kw...)
    valid = [:stoichiometry, :mole_proportions, :mass_proportions]
    spec in valid || error("Invalid composition specification $(spec)")
    return component(getfield(AuChimiste, spec)(; kw...))
end

"""
    stoichiometry(; kw...)

Create composition based on elemental stoichiometry.
"""
function stoichiometry(; kw...)
    return Composition{Stoichiometry}(; kw...)
end

"""
    mole_proportions(; scale = nothing, kw...)

Create composition based on relative molar proportions. The main
different w.r.t. [`stoichiometry`](@ref) is the presence of a
scaling factor to correct stoichiometry representation of the
given composition.
"""
function mole_proportions(; scale = nothing, kw...)
    scale = something(scale, first(kw).first => 1.0)
    return Composition{MoleProportion}(; scale, kw...)
end

"""
    mass_proportions(; scale = nothing, kw...)

Create composition based on relative molar proportions. This is
essentially the same thing as [`mole_proportions`](@ref) but in
this case the element keywords are interpreted as being the mass
proportions ofa associated elements.
"""
function mass_proportions(; scale = nothing, kw...)
    scale = something(scale, first(kw).first => 1.0)
    return Composition{MassProportion}(; scale, kw...)
end

"""
    stoichiometry_map(c::ChemicalComponent)

Returns component map of elemental stoichiometry.
"""
function stoichiometry_map(c::ChemicalComponent)
   return NamedTuple(zip(c.elements, c.coefficients))
end

"""
    mole_fractions_map(c::ChemicalComponent)

Returns component map of elemental mole fractions.
"""
function mole_fractions_map(c::ChemicalComponent)
   return NamedTuple(zip(c.elements, c.mole_fractions))
end

"""
    mass_fractions_map(c::ChemicalComponent)

Returns component map of elemental mass fractions.
"""
function mass_fractions_map(c::ChemicalComponent)
   return NamedTuple(zip(c.elements, c.mass_fractions))
end

#######################################################################
# INTERNALS
#######################################################################

@enum CompositionTypes begin
    Stoichiometry
    MoleProportion
    MassProportion
end

struct Composition{T}
    data::NamedTuple
    scale::Pair{Symbol,<:Number}

    function Composition{T}(;
            scale::Union{Nothing, Pair{Symbol, <:Number}} = nothing,
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

function component_core(c)
    elems = vcat(keys(c.data)...)
    coefs = vcat(values(c.data)...)
    U     = coefs ./ sum(coefs)
    W     = atomic_mass.(elems)
    idx   = findfirst(x->x==c.scale.first, elems)
    return elems, coefs, U, W, idx
end

function component_close(c, idx, X, W)
    coefs = (c.scale.second / X[idx]) .* X
    return coefs, coefs' * W
end

function component(c::Composition{Stoichiometry})
    elems, coefs, X, W, idx = component_core(c)
    Y = get_mass_fractions(X, W)
    coefs, M = component_close(c, idx, coefs, W)
    return ChemicalComponent(elems, coefs, X, Y, M)
end

function component(c::Composition{MoleProportion})
    elems, coefs, X, W, idx = component_core(c)
    Y = get_mass_fractions(X, W)
    coefs, M = component_close(c, idx, X, W)
    return ChemicalComponent(elems, coefs, X, Y, M)
end

function component(c::Composition{MassProportion})
    elems, coefs, Y, W, idx = component_core(c)
    X = get_mole_fractions(Y, W)
    coefs, M = component_close(c, idx, X, W)
    return ChemicalComponent(elems, coefs, X, Y, M)
end

# Create Quantity
# Operations in quantities act on masses

# MassQuantity(d::Dict) = MassQuantity([n=>v for (n, v) in d])

# function Base.:*(c::Number, s::ChemicalCompound)::MassQuantity
#     return MassQuantity(map(Pair, s.elements, c * s.Y))
# end

# function Base.:*(c::Number, s::MassQuantity)::MassQuantity
#     return MassQuantity(map((e)->e.first => c * e.second, s.amounts))
# end

# function Base.:+(a::MassQuantity, b::MassQuantity)::MassQuantity
#     da, db = Dict(a.amounts), Dict(b.amounts)
#     allkeys = [union(keys(da), keys(db))...]
#     return MassQuantity(map(k->k=>get(da, k, 0)+get(db, k, 0), allkeys))
# end
