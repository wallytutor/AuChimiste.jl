# -*- coding: utf-8 -*-

export mean_molecular_mass
export get_mole_fractions
export get_mass_fractions

"""
    mean_molecular_mass(U, W; basis)

Compute mean molecular mass based on given composition data.
"""
function mean_molecular_mass(U, W; basis)
    basis == :mole && return mean_molecular_mass_x(U, W)
    basis == :mass && return mean_molecular_mass_y(U, W)
    error("Unknown composition basis $(basis).")
end

@doc """
    get_mole_fractions(Y, W)
    get_mole_fractions(Y, W, M)

Get mole fractions from mass fractions.
""" get_mole_fractions

function get_mole_fractions(Y, W)
    return mean_molecular_mass_y(Y, W) * @. Y / W
end

function get_mole_fractions(Y, W, M)
    return M * @. Y / W
end

@doc """
    get_mass_fractions(X, W)
    get_mass_fractions(X, W, M)

Get mass fractions from mole fractions.
""" get_mass_fractions

function get_mass_fractions(X, W)
    return (@. X * W) / mean_molecular_mass_x(X, W)
end

function get_mass_fractions(X, W, M)
    return (@. X * W) / M
end

#######################################################################
# INTERNALS
#######################################################################

"""
    mean_molecular_mass_y(Y, W)

Mean molecular mass computed from mass fractions.
"""
function mean_molecular_mass_y(Y, W)
    return sum(@. Y / W)^(-1.0)
end

"""
    mean_molecular_mass_x(X, W)

Mean molecular mass computed from mole fractions.
"""
function mean_molecular_mass_x(X, W)
    return sum(@. X * W)
end
