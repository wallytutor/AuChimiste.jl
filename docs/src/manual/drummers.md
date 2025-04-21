# Drummers

```@meta
CurrentModule = AuChimiste
```

## Kramers equation

Kramers equation [Kramers1952](@cite) describes the slope of bed height ``h`` over rotary kiln axis ``z``, with discharge being given at ``z=0`` where initial condition of granule size is expected to be provided. It accounts for rotation speed ``\omega``, volumetric flow rate ``\Phi_v``, kiln inclination angle ``\alpha``, kiln internal radius ``R``, and solids dynamic repose angle ``\beta``.  In `Auchimiste` its implementation is done by [`DrumMediumKramersChunk`](@ref) and is decomposed as provided in the following equation:

```math
\begin{align*}
    \frac{dh}{dz} &= A\left(2r-r^2\right)^{-\frac{3}{2}} - B
    \\[6pt]
    %
    A &= \frac{3}{4}\frac{\Phi_{v}\tan{\beta}}{\omega\pi{}R^3}
    \\[6pt]
    %
    B &= \frac{\tan{\alpha}}{\cos{\beta}}
    \\[6pt]
    %
    r &= \frac{h}{R}
\end{align*}
```

Because thermal effects may impact solids dynamic repose angle ``\beta``, it can be provided as a function of coordinate ``z`` (it is expected the user has solved the thermal model elsewhere and solution can be retrieved in terms of this coordinate); the same applies to volumetric flow ``\Phi_v``. For modeling transitions of radius, ``R`` is also to be provided as a function of ``z``, but that must be done with care to provide a suitable discretization that is compatible with the provided transitions.

```@docs
AuChimiste.DrumMediumKramersChunk
```

Raw results of Kramers equation [Kramers1952](@cite) integration are rarely useful without further processing. To this end a post-processing type [`DrumMediumKramersSolution`](@ref) is provided, putting together most quantities of practical use those required to provide the geometry of a bed inside a kiln for process simulation.

The bed view angle ``\theta`` from kiln axis is computed as

```math
\theta = 2\cos^{-1}\left(1 - \dfrac{h}{R}\right)
```

With its value at hand one can compute the local bed chord length ``l_{c}``

```math
l_{c} = 2R\sin\left(\dfrac{\theta}{2}\right)
```

From geometric consideration one shows that bed cross-section area ``A_{b}`` is

```math
A_{s} = \dfrac{\theta R^2 - l_{c}(R - h)}{2}
```

Local fractional loading ``\eta`` can be demonstrated to be

```math
\eta = \dfrac{\theta - \sin(\theta)}{2\pi}
```

Mean loading ``\Eta`` is simply an application of the fundamental theorem of Calculus for ``\eta`` over kiln length ``L``

```math
\Eta = \dfrac{1}{L}\displaystyle\int_{0}^{L}\eta {}dz
```

Since the model does not assume conservation of ``\Phi_{v}``, residence time ``\tau`` needs to be evaluated from local mean drift-velocity ``v_{d}``; it is trivial to show that in interval ``z\in[z_1;z_2]`` the residence is given as

```math
\tau = \displaystyle\int_{z_1}^{z_2}\dfrac{dz}{v_{d}}
\quad\text{where}\quad
v_{d} = \dfrac{A_{b}}{\Phi_{v}}
```

```@docs
AuChimiste.DrumMediumKramersSolution
```

Integration of the model can be done directly with `Common.solve` interface for single chunks (kilns without user-defined discretization) or with [`solve_kramers_stack`](@ref).

```@docs
AuChimiste.solve_kramers_stack
CommonSolve.solve(::DrumMediumKramersChunk)
```