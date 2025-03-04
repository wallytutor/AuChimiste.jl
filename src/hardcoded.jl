# -*- coding: utf-8 -*-

#######################################################################
# MujumdarFlueProperties
#######################################################################

"""
Provides properties for flue gases as proposed by [Mujumdar2006ii](@cite)
for the simulation of rotary kilns. These properties are provided for
benchmarking against reference model only and are not recommended to be
used in simulations as they are known not to be very accurate and do not
account for composition dependency. This type implements the traits of
[`specific_heat`](@ref), [`thermal_conductivity`](@ref), and
[`viscosity`](@ref).

Thermal conductivity is extracted from [Thunman2002](@cite). The other
properties are tracked back to [Guo2003](@cite). Although the later also
provides a model for thermal conductivity, it is not used by the authors
of [Mujumdar2006ii](@cite), probably due to its linear dependency on
the temperature.
"""
struct MujumdarFlueProperties end

@doc """
    specific_heat(::MujumdarFlueProperties, T)

Provides specific heat capacity for flue gases [Mujumdar2006ii](@cite).
"""
function specific_heat(::MujumdarFlueProperties, T)
    return 0.106T + 1173.0
end

@doc """
    thermal_conductivity(::MujumdarFlueProperties, T)

Provides thermal conductivity for flue gases [Mujumdar2006ii](@cite).
"""
function thermal_conductivity(::MujumdarFlueProperties, T)
    k = 1.581e-17
    k = T * k - 9.463e-14
    k = T * k + 2.202e-10
    k = T * k - 2.377e-07
    k = T * k + 1.709e-04
    k = T * k - 7.494e-03
    return k
end

@doc """
    viscosity(::MujumdarFlueProperties, T)

Provides dynamic viscosity for flue gases [Mujumdar2006ii](@cite).
"""
function viscosity(::MujumdarFlueProperties, T)
    return 1.672e-06sqrt(T) - 1.058e-05
end

#######################################################################
# LawnHfoProperties
#######################################################################

"""
Provides properties for heavy fuel-oil as proposed by [Lawn1987](@cite).
This type implements the traits of [`specific_heat`](@ref).
"""
struct LawnHfoProperties end

@doc """
    specific_heat(::LawnHfoProperties, T, S)

Heavy fuel-oil specific heat estimation in terms of relative density
``S`` as provided by Cragoe (1929). Temperature supplied in kelvin.
"""
function specific_heat(::LawnHfoProperties, T, S)
    return 1000.0 * (1.683 + 0.00339T) / S
end

"""
    enthalpy_net_bs2869(; rho, pct_water, pct_ash, pct_sulphur)
    enthalpy_net_bs2869(rho, x, y, s)

Heavy fuel-oil net energy capacity accordinto to BS2869:1983. Value
is computed in [MJ/kg].

For the key-word interface (recommended) parameters are given as:

- `rho`: HFO density at 15 °C, [kg/m³].
- `pct_water`: Mass percentage of water, [%].
- `pct_ash`: Mass percentage of ashes, [%].
- `pct_sulphur`: Mass percentage of sulphur, [%].

For the positional interface (internal) parameters are given as:

- `rho`: HFO density at 15 °C, [kg/L].
- `x`: Mass fraction of water.
- `y`: Mass fraction of ashes.
- `s`: Mass fraction of sulphur.
"""
function enthalpy_net_bs2869(;
        density,
        pct_water = 0.0,
        pct_ash = 0.0,
        pct_sulphur = 0.0
    )
    # In the reference ρ is given as kg/L.
    rho = density / 1000.0

    # These are mass fractions.
    x = 0.01pct_water
    y = 0.01pct_ash
    s = 0.01pct_sulphur

    return enthalpy_net_bs2869(rho, x, y, s)
end

function enthalpy_net_bs2869(rho, x, y, s)
    # Base enthalpy estimation from empirical data.
    H = 46.423 - rho * (8.792rho - 3.170)

    # Correction factor and additional components.
    return H * (1.0 - x - y - s) + 9.420s - 2.449x
end
