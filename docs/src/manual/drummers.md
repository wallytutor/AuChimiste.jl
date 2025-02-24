# Drummers

```@meta
CurrentModule = AuChimiste
```

## Kramers equation

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
