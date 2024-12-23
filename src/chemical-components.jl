# -*- coding: utf-8 -*-


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
    c
end

function component(c::T) where {T <: Composition{MassProportion}}
    c
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

# function Base.:*(c::Number, s::Stoichiometry)::Stoichiometry
#     return Stoichiometry(map(x->x[1]=>c*x[2], s.amounts))
# end

# function Base.:*(s::Stoichiometry, c::Number)::Stoichiometry
#     return c * s
# end

# function Base.:+(a::Stoichiometry, b::Stoichiometry)::Stoichiometry
#     da, db = Dict(a.amounts), Dict(b.amounts)
#     allkeys = [union(keys(da), keys(db))...]
#     return Stoichiometry(map(k->k=>get(da, k, 0)+get(db, k, 0), allkeys))
# end

# struct MassQuantity
#     amounts::Vector{ElementalQuantity}
# end

# MassQuantity(d::Dict) = MassQuantity([n=>v for (n, v) in d])

# function Base.:*(c::Number, s::ChemicalCompound)::MassQuantity
#     return MassQuantity(map(Pair, s.elements, c * s.Y))
# end

# function Base.:*(c::Number, s::MassQuantity)::MassQuantity
#     return MassQuantity(map((e)->e.first => c * e.second, s.amounts))
# end

# Base.:*(s::ChemicalCompound, c::Number)::MassQuantity = c * s

# Base.:*(s::MassQuantity, c::Number)::MassQuantity = c * s

# function Base.:+(a::MassQuantity, b::MassQuantity)::MassQuantity
#     da, db = Dict(a.amounts), Dict(b.amounts)
#     allkeys = [union(keys(da), keys(db))...]
#     return MassQuantity(map(k->k=>get(da, k, 0)+get(db, k, 0), allkeys))
# end

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
