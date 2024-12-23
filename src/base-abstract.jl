# -*- coding: utf-8 -*-

abstract type ChemicalException       <: Exception end
abstract type ChemicalElementsError   <: ChemicalException end
abstract type ChemicalComponentsError <: ChemicalException end
