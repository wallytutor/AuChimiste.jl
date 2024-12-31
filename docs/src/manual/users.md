# Users

```@meta
CurrentModule = AuChimiste
```

## Chemical Elements

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

## Chemical Components

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

## Chemical Kinetics


## Chemical Reactors


## Combustion Chemistry


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