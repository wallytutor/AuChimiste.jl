# -*- coding: utf-8 -*-

abstract type ChemicalException       <: Exception end
abstract type ChemicalElementsError   <: ChemicalException end
abstract type ChemicalComponentsError <: ChemicalException end

abstract type AbstractThermodynamicData end
abstract type AbstractThermodynamicsModel end
abstract type AbstractTransportModel end
abstract type AbstractThermalAnalysis end

abstract type AbstractReactorModel end
abstract type AbstractDrumBedModel <: AbstractReactorModel end
