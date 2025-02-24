# Drummers

*Drummers* is how we call rotary drum models in AuChimiste. This more specialized set of functionalities can be used for process estimations and simulations, especially in the field of rotary kilns. In what follows we illustrate how to solve relevant equations and extract useful data with the provided functionalities.

```@setup getting-started-1
using AuChimiste
```

## Solving Kramers equation

The simplest way to get Kramers equation solved is by calling [`solve_kramers_stack`](@ref)

```@example getting-started-1
sol = solve_kramers_stack(;
	grid   = [0, 13.7],
	radius = (_) -> 0.95,
	beta   = (_) -> deg2rad(45.0),
	phiv   = (_) -> 2.88e-03,
	h      = 0.001,
	ω̇      = 0.05,
	α      = 0.0416
)
```

Standardized plotting of [`DrumMediumKramersSolution`](@ref) bed profile is provided bellow. It supports normalization of axes throught keywords `normz` for axial coordinate and `normh` for bed depth.

```@example getting-started-1
fig, ax = AuChimiste.plot(sol)
resize!(fig.scene, 650, 350)
fig
```

```@example getting-started-1
grid = LinRange(0, 2, 50)
grid = vcat(grid, LinRange(1.7, 13.7, 12))

sol = AuChimiste.solve_kramers_stack(;
    grid   = grid,
    radius = (_) -> 0.95,
    beta   = (_) -> deg2rad(45.0),
    phiv   = (_) -> 2.88e-03,
    h      = 0.001,
    ω̇      = 0.05,
    α      = 0.0416
)
```
