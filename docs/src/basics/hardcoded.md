# Hardcoded

```@setup getting-started-1
using AuChimiste
```

## Mujumdar's data

Mujumdar *et al.* [Mujumdar2006ii](@cite) propose temperature dependent properties for flue gases to be used in the simulation of rotary kilns. Although these values are not recommended for practical purposes (especially because their thermal conducivity diverges from reasonable values above 1500 K), they are implemented for benchmarking against the reference literature model. Please, check [`AuChimiste.MujumdarFlueProperties`](@ref) for more details.

```@example getting-started-1
model = AuChimiste.MujumdarFlueProperties()

c = specific_heat(model, T_NORMAL)
k = thermal_conductivity(model, T_NORMAL)
μ = viscosity(model, T_NORMAL)

c, k, μ
```
