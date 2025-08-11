# Transport

```@meta
CurrentModule = AuChimiste
```

## Fermi-like viscosity parametrization

In numerical simulation one often faces the task to represent melting of solution solids that soften over a temperature range. An easy way to set this up with a volume-of-fluid (VOF) approach is to have some sort of exponential temperature-dependent viscosity. This structure encapsulates the evaluation of a viscosity function based on Fermi distribution - a sort of sigmoid function - to this end. This can be expressed as (definitions in the structure documentation):

$$
\mu(T) = \mu_{\infty} +
\dfrac{\mu_{0}-\mu_{\infty}}{1+\exp\left(\dfrac{T-\Theta}{\Delta}\right)}
%
\quad\text{where}\quad
%
\begin{cases}
\Theta &= \dfrac{T_{e}+T_{s}}{2}\\[12pt]
\Delta &= \dfrac{T_{e}-T_{s}}{\kappa}
\end{cases}
$$

```@docs
AuChimiste.FermiLikeViscosity
AuChimiste.fermilike_viscosity
```
