# -*- coding: utf-8 -*-

"""
    specific_heat(args...; kwargs...)

Common interface for evaluation of the specific heat of a substance.
Generally this function will take an object with the type for which it
implements the specific heat, a temperature, and maybe an array of
mass fractions. Its return value must be provided in standard SI units
of ``J\\cdotp{}kg^{-1}\\cdotp{}K^{-1}``.
"""
function specific_heat end

"""
    thermal_conductivity(args...; kwargs...)

Common interface for evaluation of the thermal conductivity of a
substance. Generally this function will take an object with the type
for which it implements the thermal conductivity, a temperature, and
maybe an array of mass fractions. Its return value must be provided
in standard SI units of ``W\\cdotp{}m^{-1}\\cdotp{}K^{-1}``.
"""
function thermal_conductivity end

"""
    viscosity(args...; kwargs...)

Common interface for evaluation of the viscosity of a substance.
Generally this function will take an object with the type for which it
implements the viscosity, a temperature, and maybe an array of mass
fractions. Its return value must be provided in standard SI units of
``Pa\\cdotp{}s``.
"""
function viscosity end
