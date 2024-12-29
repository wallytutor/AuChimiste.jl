# AuChimiste Toolbox

From elements to plain gold (and kinetics), all in plain Julia.

*AuChimiste* is an *alchimiste* wordplay meaning *to the chemist* in French.

Please check the [Getting Started](basics/elements.md) (and the whole suit) in the sidebar if it is your first time here.

If this module was useful in your work, please consider citing us [auchimistejl](@cite).

## Project goals and status

*One toolbox, all chemistry.*

- Provide chemical elements with symbolic support and built-in data; utilities are expected to allow users to define their own elements (*e.g.* isotopes) and retrieve data. This is all implemented in [chemical-elements.jl](https://github.com/wallytutor/AuChimiste.jl/tree/main/src/chemical-elements.jl).

- **WIP:** By making use of chemical elements we provide [chemical-components.jl](https://github.com/wallytutor/AuChimiste.jl/tree/main/src/chemical-cmponents.jl). This file is called this way because it is intended to include anything from species, compounds, solids, etc., so no other name suited its ends. It is responsible by:

    - Provide creation of arbitrary compounds from mass or molar composition (try to understand this term in the broader sense) with the other composition being computed, *i.e.* if mass fractions were provided, the compound can access its molar composition, and molecular mass.

    - **WIP:** Arithmetic of compounds to create complex compositions and manipulation of amounts of matter. This sort of functionality is aimed at computing mixtures for experiments, validation of chemical reactions mass balance, or simply creating new compounds expressed in terms of components, as is often the case in intermetallics or complex oxide systems.

- **WIP:** Putting chemical components together one can express reactions; with reactions we open the gates to [chemical-kinetics.jl](https://github.com/wallytutor/AuChimiste.jl/tree/main/src/chemical-kinetics.jl). This introduces the expression of symbolic kinetics for ease of integration in reactor models. It provides parsing of [Cantera](https://cantera.org/) mechanism and reusable code generation for simulating mechanisms.

- **WIP:** Chemical kinetics provides the basis for the construction of reactor models in [chemical-reactors.jl](https://github.com/wallytutor/AuChimiste.jl/tree/main/src/chemical-reactors.jl). It is built upon [ModelingToolkit](https://docs.sciml.ai/ModelingToolkit/stable/) blocks and allows for chains of reactors, plug-flow reactors,...

- **WIP:** Supporting the above there is [physical-chemistry.jl](https://github.com/wallytutor/AuChimiste.jl/tree/main/src/physical-chemistry.jl) which provides the required closure models for the different models, and [combustion-chemistry.jl](https://github.com/wallytutor/AuChimiste.jl/tree/main/src/combustion-chemistry.jl) a specialized package for the analysis and simulation of combustion systems.

- **WIP:** Going one step further, [chemical-thermodynamics.jl](https://github.com/wallytutor/AuChimiste.jl/tree/main/src/chemical-thermodynamics.jl) makes use of some of the above to provide chemical thermodynamics computations, with focus in phase equilibria and [CALPHAD](https://calphad.org/) approaches.

- **WIP:** For all types for which it applies, provide [recipes](https://korsbo.github.io/Latexify.jl/stable/tutorials/recipes/) for ``\LaTeX`` representation of entities, so that report generation can be automated.

## Related tools

Searching for [chemistry](https://juliapackages.com/packages?search=chemistry), [kinetics](https://juliapackages.com/packages?search=kinetics), or [thermodynamics](https://juliapackages.com/packages?search=thermodynamics) in Julia Packages does not lead to any convincing package or ecosystem in competition with what is aimed here, justifying its existence. 

Some clarifications regarding the design choices of this package:

- It does not intend to replace [Cantera](https://cantera.org/), but to provide similar functionality in a algorithmically-differentiable way for some of its applications. The main difference here is the focus at supporting user-defined models.

- It also does not compete with [pyJac](https://slackha.github.io/pyJac/) as all code generation is aimed to be plain Julia. While `pyJac` uses analytically derived formulas for jacobian evaluation, our intent here is to let the user chose how the AD code will be employed in their simulations.

- Regarding [Catalyst.jl](https://docs.sciml.ai/Catalyst/stable/), our goal is not to analyse kinetics in the same sense, but to use mechanisms (with thermochemistry integrated, what lacks there) in larger simulations.
