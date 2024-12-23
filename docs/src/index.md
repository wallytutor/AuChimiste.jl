# AuChimiste Toolbox

From elements to plain gold (and kinetics), all in plain Julia.

*AuChimiste* is an *alchimiste* wordplay meaning *to the chemist* in French.

Please check the [Getting Started](getting-started.md) guide in the sidebar.

If this module was useful in your work, please consider citing us [auchimistejl](@cite).

## Development rules

- Code written, code documented, code tested.
- Code lines make 72 characters, never more than 79.
- Code is not cluttered and comments are minimal.
- Code abuses of multiple dispatch if needed.
- Code is Julia, nothing else.

## Related tools

Searching for [chemistry](https://juliapackages.com/packages?search=chemistry), [kinetics](https://juliapackages.com/packages?search=kinetics), or [thermodynamics](https://juliapackages.com/packages?search=thermodynamics) in Julia Packages does not lead to any convincing package or ecosystem in competition with what is aimed here, justifying its existence. 

Some clarifications regarding the design choices of this package:

- It does not intend to replace [Cantera](https://cantera.org/), but to provide similar functionality in a algorithmically-differentiable way for some of its applications. The main difference here is the focus at supporting user-defined models.

- It also does not compete with [pyJac](https://slackha.github.io/pyJac/) as all code generation is aimed to be plain Julia. While `pyJac` uses analytically derived formulas for jacobian evaluation, our intent here is to let the user chose how the AD code will be employed in their simulations.

- Regarding [Catalyst.jl](https://docs.sciml.ai/Catalyst/stable/), our goal is not to analyse kinetics in the same sense, but to use mechanisms (with thermochemistry integrated, what lacks there) in larger simulations.