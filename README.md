# AuChimiste.jl

Please read the [docs](https://wallytutor.github.io/AuChimiste.jl/dev/).

## Project goals and status

*One toolbox, all chemistry.*

- [ ] Provide chemical elements with symbolic support and built-in data; utilities are expected to allow users to define their own elements (*e.g.* isotopes) and retrieve data. This is all implemented in [ChemicalElements.jl](src/ChemicalElements.jl).

- [ ] By making use of `ChemicalElements` we provide [ChemicalComponents.jl](src/ChemicalComponents.jl). This module is called this way because it is intended to include anything from species, compounds, solids, etc., so no other named suited its ends. It is responsible by:

    - [ ] Provide creation of arbitrary compounds from mass or molar composition (try to understand this term in the broader sense) with the other composition being computed, *i.e.* if mass fractions were provided, the compound can access its molar composition, and molecular mass.

    - [ ] Arithmetic of compounds to create complex compositions and manipulation of amounts of matter. This sort of functionality is aimed at computing mixtures for experiments, validation of chemical reactions mass balance, or simply creating new compounds expressed in terms of components, as is often the case in intermetallics or complex oxide systems.

- [ ] Putting `ChemicalComponents` together one can express reactions; with reactions we open the gates to [ChemicalKinetics.jl](src/ChemicalKinetics.jl). This module allows for the expression of symbolic kinetics for ease of integration in reactor models. It provides parsing of [Cantera](https://cantera.org/) mechanism and reusable code generation for simulating mechanisms.

- [ ] Module `ChemicalKinetics` provides the basis for the construction of reactor models in [ChemicalReactors.jl](src/ChemicalReactors.jl). It is built upon [ModelingToolkit](https://docs.sciml.ai/ModelingToolkit/stable/) blocks and allows for chains of reactors, plug-flow reactors,...

-[ ] Supporting the above there is [PhysicalChemistry.jl](src/PhysicalChemistry.jl) which provides the required closure models for the different modules, and [CombustionChemistry.jl](src/CombustionChemistry.jl) a specialized package for the analysis and simulation of combustion systems.
