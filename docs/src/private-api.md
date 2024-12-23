# Private API

```@meta
CurrentModule = AuChimiste
```

This part of the documentation is intended for developers. It might also be useful for standard users trying to understand bugs or propose features. `AuChimiste` aims at having 100% first-level entities documented so that design features can be understood in the future.

## Development rules

- Code written, code documented, code tested.
- Code lines make 72 characters, never more than 79.
- Code is not cluttered and comments are minimal.
- Code abuses of multiple dispatch if needed.
- Code is Julia, nothing else.

## Chemical Elements

```@docs
AuChimiste.ELEMENTS
AuChimiste.USER_ELEMENTS
AuChimiste.handle_element
AuChimiste.find_element
```

## Chemical Components

```@docs
```

## Chemical Kinetics


## Chemical Reactors


## Combustion Chemistry


## Physical Chemistry

```@docs
AuChimiste.mean_molecular_mass_y
AuChimiste.mean_molecular_mass_x
```

## Chemical Thermodynamics

## Exception types

```@docs
AuChimiste.NoSuchElementError
AuChimiste.NoIsotopeProvidedError
AuChimiste.EmptyCompositionError
AuChimiste.InvalidScalerError
```