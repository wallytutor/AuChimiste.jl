# -*- coding: utf-8 -*-

export ELECTRON_MASS
export AVOGADRO
export GAS_CONSTANT
export STEFAN_BOLTZMANN
export JOULE_PER_CALORIE
export P_NORMAL
export T_NORMAL
export C_NORMAL
export T_STANDARD
export JOULE_PER_CALORIE

"Electron Mass ``m_e`` [$(ELECTRON_MASS) kg]"
const ELECTRON_MASS = 9.109_382_915e-31

"Avogadro's Number ``N_\\mathrm{A}`` [$(AVOGADRO) number/kmol]"
const AVOGADRO = 6.022_140_76e26

"Ideal gas constant ``R`` [$(GAS_CONSTANT) J/(mol.K)]."
const GAS_CONSTANT::Float64 = 8.314_462_618_153_24

"Stefan-Boltzmann constant ``\\sigma`` [$(STEFAN_BOLTZMANN) W/(m².K⁴)]."
const STEFAN_BOLTZMANN::Float64 = 5.670_374_419e-08

"Reference atmospheric pressure [$(P_NORMAL) Pa]."
const P_NORMAL::Float64 = 101325.0

"Normal atmospheric temperature [$(T_NORMAL) K]."
const T_NORMAL::Float64 = 273.15

"Normal atmospheric concentration [$(C_NORMAL) mol/m³]. "
const C_NORMAL::Float64 = P_NORMAL / (GAS_CONSTANT * T_NORMAL)

"Standard atmospheric temperature [$(T_STANDARD) K]."
const T_STANDARD::Float64 = 298.15

"Conversion factor from calories to joules [$(JOULE_PER_CALORIE) J/cal]."
const JOULE_PER_CALORIE::Float64 = 4.184
