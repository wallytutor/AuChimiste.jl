# -*- coding: utf-8 -*-

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
