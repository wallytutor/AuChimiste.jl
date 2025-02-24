# Users

```@meta
CurrentModule = AuChimiste
```

## Elements

```@docs
AuChimiste.has_element
AuChimiste.list_elements
AuChimiste.reset_elements_table
AuChimiste.add_element
AuChimiste.add_isotope
AuChimiste.atomic_mass
AuChimiste.atomic_number
AuChimiste.element_name
AuChimiste.element
AuChimiste.AtomicData
```

## Components

```@docs
AuChimiste.ChemicalComponent
AuChimiste.ComponentQuantity

AuChimiste.component
AuChimiste.stoichiometry
AuChimiste.mole_proportions
AuChimiste.mass_proportions
AuChimiste.stoichiometry_map
AuChimiste.mole_fractions_map
AuChimiste.mass_fractions_map
AuChimiste.quantity
```

The following are not exported but is worth the end-user to known them:

```@docs
AuChimiste.CompositionTypes
AuChimiste.Composition
```

## Interfaces

The following interfaces are provided as a centralization of names for the package. Generally these functions will take an object with the type for which they implement the quantity associated with their name, and other relevant parameters, such as temperature, pressure, and/or an array of mass fractions. In all cases values are expected to be returned in SI units, as documented by each function.

```@docs
AuChimiste.density
AuChimiste.specific_heat
AuChimiste.enthalpy
AuChimiste.entropy
AuChimiste.thermal_conductivity
AuChimiste.viscosity
```

## Kinetics


## Reactors


## Combustion


## Physical Chemistry

```@docs
AuChimiste.mean_molecular_mass_y
AuChimiste.mean_molecular_mass_x
AuChimiste.mean_molecular_mass
AuChimiste.get_mole_fractions
AuChimiste.get_mass_fractions
```

## Thermodynamics

```@docs
AuChimiste.ThermoData
AuChimiste.NASAThermo
AuChimiste.ShomateThermo
```

## Hardcoded

```@docs
AuChimiste.MujumdarFlueProperties
```

## Drummers

Kramers equation [Kramers1952](@cite) describes the slope of bed height ``h`` over rotary kiln axis ``z``, with discharge being given at ``z=0`` where initial condition of granule size is expected to be provided. It accounts for rotation speed ``\omega``, volumetric flow rate ``\Phi_v``, kiln slope ``\alpha``, kiln internal radius ``R``, and solids dynamic repose angle ``\beta``.  In `Auchimiste` its implementation is done by [`DrumMediumKramersChunk`](@ref) and is decomposed as provided in the following equation:

```math
\begin{align*}
\frac{dh}{dz} &= A\left(2r-r^2\right)^{-\frac{3}{2}} - B
\\[6pt]
A &= \frac{3}{4}\frac{\Phi_{v}\tan{\beta}}{\omega\pi{}R^3}
\\[6pt]
B &= \frac{\tan{\alpha}}{\cos{\beta}}
\\[6pt]
r &= \frac{h}{R}
\end{align*}
```

Because thermal effects may impact solids dynamic repose angle ``\beta``, it can be provided as a function of coordinate ``z`` (it is expected the user has solved the thermal model elsewhere and solution can be retrieved in terms of this coordinate); the same applies to volumetric flow ``\Phi_v``. For modeling transitions of radius, ``R`` is also to be provided as a function of ``z``, but that must be done with care to provide a suitable discretization that is compatible with the provided transitions.

```@docs
AuChimiste.DrumMediumKramersSolution
AuChimiste.DrumMediumKramersChunk
AuChimiste.solve_kramers_stack
```

## Database parsing

## Constants

```@docs
AuChimiste.ELECTRON_MASS
AuChimiste.AVOGADRO
AuChimiste.GAS_CONSTANT
AuChimiste.STEFAN_BOLTZMANN
AuChimiste.P_NORMAL
AuChimiste.T_NORMAL
AuChimiste.T_STANDARD
```

## Exception types

Exception types are not exported as they are not intended for other purposes than handling errors inside `AuChimiste`. The documentation below is provided so that the end-user can better understand their occurrence for debugging their own code or reporting bugs:

```@docs
AuChimiste.NoSuchElementError
AuChimiste.NoIsotopeProvidedError
AuChimiste.EmptyCompositionError
AuChimiste.InvalidScalerError
```