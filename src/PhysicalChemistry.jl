# -*- coding: utf-8 -*-
module PhysicalChemistry

mean_molecular_mass_y(Y, W) = sum(@. Y / W)^(-1.0)
mean_molecular_mass_x(X, W) = sum(@. X * W)

function mean_molecular_mass(U, W; basis)
    basis == :mole && return mean_molecular_mass_x(U, W)
    basis == :mass && return mean_molecular_mass_y(U, W)
    error("Unknown composition basis $(basis).")
end

get_mole_fractions(Y, W) = mean_molecular_mass_y(Y, W) * @. Y / W
get_mole_fractions(Y, W, M) = M * @. Y / W

get_mass_fractions(X, W) = (@. X * W) / mean_molecular_mass_x(X, W)
get_mass_fractions(X, W, M) = (@. X * W) / M

#######################################################################
# API
#######################################################################

@doc """\

    mean_molecular_mass(U, W; basis)

Compute mean molecular mass based on given composition data.
""" mean_molecular_mass
export mean_molecular_mass

@doc """\

    get_mole_fractions(Y, W)
    get_mole_fractions(Y, W, M)

Get mole fractions from mass fractions.
""" get_mole_fractions
export get_mole_fractions

@doc """\

    get_mass_fractions(X, W)
    get_mass_fractions(X, W, M)

Get mass fractions from mole fractions.
""" get_mass_fractions
export get_mass_fractions

end # (module PhysicalChemistry)