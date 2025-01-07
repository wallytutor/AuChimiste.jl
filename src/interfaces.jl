# -*- coding: utf-8 -*-

"""
    density(args...; kwargs...)

Evaluation of the density of a substance.
Its return value must be in ``kg\\cdotp{}m^{-3}``.
"""
function density end

"""
    specific_heat(args...; kwargs...)

Evaluation of the specific heat of a substance.
Its return value must be in ``J\\cdotp{}kg^{-1}\\cdotp{}K^{-1}``.
"""
function specific_heat end

"""
    enthalpy(args...; kwargs...)

Evaluation of the enthalpy of a substance.
Its return value must be in ``J\\cdotp{}kg^{-1}``.
"""
function enthalpy end

"""
    entropy(args...; kwargs...)

Evaluation of the entropy of a substance.
Its return value must be in ``J\\cdotp{}K^{-1}``.
"""
function entropy end

"""
    thermal_conductivity(args...; kwargs...)

Evaluation of the thermal conductivity of a substance.
Its return value must be in ``W\\cdotp{}m^{-1}\\cdotp{}K^{-1}``.
"""
function thermal_conductivity end

"""
    viscosity(args...; kwargs...)

Evaluation of the viscosity of a substance.
Its return value must be in ``Pa\\cdotp{}s``.
"""
function viscosity end
