# Drummers

*Drummers* is how we call rotary drum models in AuChimiste. This more specialized set of functionalities can be used for process estimations and simulations, especially in the field of rotary kilns. In what follows we illustrate how to solve relevant equations and extract useful data with the provided functionalities.

```@setup getting-started-1
using AuChimiste
using Statistics: mean
```

## Solving Kramers equation

The simplest way to get Kramers equation solved is by calling [`solve_kramers_stack`](@ref) directly; all inputs are provided in international system units as stated in the documentation. Below we provide a partial implementation of the last example provided in Kramers [Kramers1952](@cite) paper.

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

Because the modeling of transitions might require special arrangements of grid, a composite discretization can be provided as follows. It must be emphasized that by no means this is checked agains radius transitions and it is up to the user to provide validation. It can be useful do adopt this approach for later post-processing and assembly in kiln models.

```@example getting-started-1
grid = let
    grid1 = LinRange(0, 2, 50)
    grid2 = LinRange(1.7, 13.7, 12)
    vcat(grid1, grid2)
end

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

In the following dummy example we force a very thick *analytical* bed solution, filling the radius of the rotary drum. Next we confirm the *internal* evaluations of the model match the expected *analytical* values.

```@example getting-started-1
let
    z = collect(0.0:0.1:10.0)
    h = 1.0 * ones(size(z))
    R = 1.0 * ones(size(z))
    β = deg2rad(45.0) * ones(size(z))
    ϕ = 0.01 * ones(size(z))
    
    bed = DrumMediumKramersSolution(z, h, R, β, ϕ)

    A = π * R.^2 / 2;
    V = A * z[end];

    @assert mean(bed.θ) ≈ π
    @assert mean(bed.l) ≈ mean(2R)
    @assert mean(A .- bed.A) ≈ 0
    @assert mean(bed.η) ≈ 0.5
    @assert V[end] / ϕ[end] ≈ bed.τ[end]

    bed
end
```
